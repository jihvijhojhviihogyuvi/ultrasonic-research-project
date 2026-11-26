"""Deterministic experiment runner. Every attempt, including failures, is saved.

Large clouds receive a primary-force screen, not an independence certificate.
Resume skips existing case.json files. Use --rerun to deliberately replace one.
"""
from __future__ import annotations
import argparse,gc,hashlib,json,platform,time,traceback
from pathlib import Path
from dataclasses import asdict,replace
import numpy as np
import torch
from physics import Config,AcousticSystem,geometry,targets,linear_floquet,affine_periodic_orbit
from optimization import Design,optimize,optimize_order,optimize_nonlinear_floquet,equilibrium,sdp_bound

ROOT=Path(__file__).resolve().parent


def clean(value):
    if isinstance(value,dict):return {k:clean(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [clean(v) for v in value]
    if isinstance(value,np.ndarray):return clean(value.tolist())
    if isinstance(value,(np.integer,np.floating)):return clean(value.item())
    if isinstance(value,float) and not np.isfinite(value):return None
    return value


def write(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(clean(value),indent=2,allow_nan=False),encoding='utf-8')


def load_design(folder,config_changes=None):
    meta=json.loads((folder/'case.json').read_text(encoding='utf-8'))
    c=Config(**meta['config'])
    if config_changes:c=replace(c,**config_changes)
    a=np.load(folder/'design.npz')
    s=AcousticSystem(c,a['targets'],a['positions'],a['normals'],a['radii'])
    d=Design(s,meta['architecture'],meta['field_count'],initial=torch.tensor(a['raw'][...,0]+1j*a['raw'][...,1],device=s.device))
    for name in ('raw','timing','gauge','common','residual_scale'):
        getattr(d,name).data.copy_(torch.tensor(a[name],device=s.device))
    return d,meta


def run_case(spec,out,rerun=False):
    folder=out/spec['id'];folder.mkdir(parents=True,exist_ok=True)
    if (folder/'case.json').exists() and not rerun:
        return json.loads((folder/'case.json').read_text(encoding='utf-8'))
    start=time.perf_counter();arch=spec['architecture'];fields=spec.get('fields',1)
    changes=dict(spec.get('config',{}));changes.update(n=spec['n'])
    if arch in ('static','ib','support'):changes['harmonics']=0
    c=Config(**changes);m=spec['m']
    row={'id':spec['id'],'category':spec['category'],'architecture':arch,'field_count':fields,
         'particles':m,'transducers':c.n,'config':asdict(c),'specification':spec,
         'independent_particles':None,'status':'started','seed':c.seed,'equilibrium_error_mm':None,
         'weakest_coupled_stiffness_N_m':None,'floquet_radius':None,'acoustic_power_W':None}
    write(folder/'attempt.json',row)
    print('START',spec['id'],flush=True)
    try:
        pos,norm,rad=geometry(replace(c,seed=spec.get('geometry_seed',c.seed)),spec.get('aperiodic',False))
        q=targets(c,m)*spec.get('cloud_scale',1.)
        s=AcousticSystem(c,q,pos,norm,rad);row['setup_s']=time.perf_counter()-start
        d=Design(s,arch,fields,seed=spec.get('initial_seed',c.seed))
        if spec.get('ib_gravity_compensation'):
            from optimization import ib
            aim=q.copy();corrections=[]
            for _ in range(8):
                qe,fe,ke,ok,cur=equilibrium(s,d.coefficients().detach(),iterations=8)
                error=qe.cpu().numpy()-q
                corrections.append(float(np.linalg.norm(error,axis=1).max())*1000)
                if np.linalg.norm(error,axis=1).max()<c.accuracy*.2:break
                aim-=.7*error
                d.raw.data.copy_(torch.view_as_real(ib(s,aim)[None]))
            row['ib_focus_correction_errors_mm']=corrections
        optimization_s=optimize(d,spec.get('iterations',120),jacobian_weight=spec.get('jacobian_weight',0))
        if spec.get('stage'):
            optimization_s+=optimize(d,spec.get('stage_iterations',8),stage=spec['stage'])
        if spec.get('order_trials'):row['order_trials']=optimize_order(d,spec['order_trials'])
        if spec.get('nonlinear_iterations'):
            optimization_s+=optimize_nonlinear_floquet(d,spec['nonlinear_iterations'],spec.get('nonlinear_steps',16))
        row['optimization_s']=optimization_s
        coef=d.coefficients().detach()
        with torch.no_grad():f,k=s.primary_mean(coef)
        row['weakest_primary_stiffness_N_m']=float(torch.linalg.eigvalsh(k).min())
        row['primary_gravity_residual_fraction']=float(torch.linalg.vector_norm(f,dim=1).max())/(c.mass*9.80665)
        row['acoustic_power_W']=float(s.power(coef))
        phases=torch.arange(1024,device=s.device)*2*np.pi/1024
        command=d.command_coefficients().detach()
        row['maximum_command_envelope']=float((torch.exp(-1j*phases[:,None]*s.ell[None])@command).abs().max())
        row['minimum_separation_mm']=float(np.min(np.linalg.norm(q[:,None]-q[None,:],axis=-1)+np.eye(m)*1e6))*1000
        row['particle_density_per_cm3']=m/(c.volume_side*100)**3
        row['cloud_density_per_cm3']=m/(c.volume_side*spec.get('cloud_scale',1.)*100)**3
        row['status']='primary_screen_rejection'
        qe=s.q.detach();eq_ok=False
        # Full coupled analysis is expensive O(M^3). High-count screens do not
        # claim to measure coupled stability, retention or capacity.
        if m<=16 and row['weakest_primary_stiffness_N_m']>=c.stiffness_required*.65:
            qe,fe,ke,eq_ok,current=equilibrium(s,coef,iterations=8)
            row['equilibrium_converged']=eq_ok
            row['equilibrium_error_mm']=float(torch.linalg.vector_norm(qe-s.q,dim=1).max())*1000
            row['coupled_gravity_residual_fraction']=float(torch.linalg.vector_norm(fe,dim=1).max())/(c.mass*9.80665)
            row['weakest_coupled_stiffness_N_m']=float(torch.linalg.eigvalsh((ke+ke.T)/2).min())
            with torch.no_grad():
                ft,kt,_=current.matrices(coef)
                rho,phi=linear_floquet(kt,c)
                row['floquet_radius']=float(rho)
                try:
                    _,orbit,_=affine_periodic_orbit(ft,kt,c)
                    row['affine_micromotion_mm']=float(torch.linalg.vector_norm(orbit.reshape(-1,m,3),dim=-1).max())*1000
                    np.save(folder/'affine_orbit_m.npy',orbit.cpu().numpy())
                except RuntimeError:row['affine_micromotion_mm']=None
            row['status']='stationary_screen_rejection'
            if (eq_ok and row['equilibrium_error_mm']<=c.accuracy*1000 and
                row['weakest_coupled_stiffness_N_m']>=c.stiffness_required and float(rho)<1):
                row['status']='stationary_candidate_not_motion_validated'
            np.savez_compressed(folder/'coupled.npz',force=ft.cpu().numpy(),stiffness=kt.cpu().numpy(),monodromy=phi.cpu().numpy())
        np.savez_compressed(folder/'design.npz',raw=d.raw.detach().cpu().numpy(),timing=d.timing.detach().cpu().numpy(),
          gauge=d.gauge.detach().cpu().numpy(),common=d.common.detach().cpu().numpy(),
          residual_scale=d.residual_scale.detach().cpu().numpy(),targets=q,positions=pos,normals=norm,radii=rad,
          equilibrium=qe.cpu().numpy(),coefficients=coef.cpu().numpy(),duties=d.duties().detach().cpu().numpy())
        row['duties']=d.duties().detach().cpu().tolist();write(folder/'history.json',d.history)
        if spec.get('sdp'):row['sdp']=sdp_bound(s)
        if spec.get('validate') and row['status']=='stationary_candidate_not_motion_validated':
            from validation import strict_validate
            result,raw=strict_validate(d,qe.cpu().numpy(),degree=8,time_samples=64)
            write(folder/'validation.json',result);write(folder/'frames.json',raw.pop('frames',[]))
            np.savez_compressed(folder/'trajectory_data.npz',**raw)
            row['validation']=result;row['status']=result['status'];row['independent_particles']=result['independent_particles']
    except Exception as error:
        row['status']='execution_error';row['exception']=str(error);row['traceback']=traceback.format_exc()
        if 'd' in locals():
            np.savez_compressed(folder/'failure_state.npz',raw=d.raw.detach().cpu().numpy(),
                timing=d.timing.detach().cpu().numpy(),gauge=d.gauge.detach().cpu().numpy(),
                targets=d.system.q.detach().cpu().numpy())
            write(folder/'history.json',d.history);row['failure_snapshot']='failure_state.npz'
        print(row['traceback'],flush=True)
    row['runtime_s']=time.perf_counter()-start
    row['source_hashes']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.glob('*.py')}
    write(folder/'case.json',row)
    print('END',spec['id'],row['status'],'K',row.get('weakest_coupled_stiffness_N_m') or row.get('weakest_primary_stiffness_N_m'),
          'rho',row.get('floquet_radius'),'seconds',round(row['runtime_s'],2),flush=True)
    gc.collect()
    return row


def plan(group):
    specs=[]
    def add(id,arch,n,m,fields=1,category=None,**kw):
        specs.append(dict(id=id,architecture=arch,n=n,m=m,fields=fields,category=category or group,**kw))
    if group=='comparisons':
        for m in (2,8,16):
            for arch,fields in [('ib',1),('static',1),('multifrequency',2),('singleton',m),('covariance',4),('support',1)]:
                add(f'compare_{m}_{arch}',arch,512,m,fields,iterations=180)
    elif group=='temporal':
        for fields in (2,3,4,8,16):
            add(f'slots_{fields}', 'covariance',512,2,fields,iterations=150)
        add('timing_floquet','covariance',512,2,4,iterations=150,stage='timing',stage_iterations=12,order_trials=8)
        add('joint_floquet','covariance',128,2,4,iterations=150,stage='floquet',stage_iterations=12)
        add('nonlinear_floquet','covariance',64,2,2,iterations=100,nonlinear_iterations=2,nonlinear_steps=16)
        add('permanent_support','common',512,2,4,iterations=180,stage='timing',stage_iterations=6)
        for smooth in (0.,.00004,.00016):
            add(f'smooth_{smooth:.6f}','covariance',512,2,4,iterations=150,config={'smooth_s':smooth})
        for cycle in (.001,.0025,.010):
            add(f'cycle_{cycle:.4f}','covariance',512,2,4,iterations=150,config={'cycle_s':cycle})
        add('control_jacobian','static',128,2,1,iterations=120,jacobian_weight=2.)
        add('sdp32','static',32,2,1,iterations=80,sdp=True)
        for seed in (7,8,9):
            add(f'geometry_{seed}','static',128,2,1,iterations=120,aperiodic=True,config={'seed':seed})
    elif group=='scaling':
        for n in (32,64,128,256,512):
            for m in (2,4,8,16,25,32,50):
                for arch,fields in [('static',1),('covariance',4)]:
                    add(f'scale_N{n}_M{m}_{arch}',arch,n,m,fields,iterations=80 if m<=16 else 50)
        for arch,fields in [('static',1),('covariance',4)]:
            add(f'scale_N512_M100_{arch}',arch,512,100,fields,iterations=40)
    elif group=='power_density':
        for factor in (.25,.5,1.):
            for arch,fields in [('static',1),('covariance',4)]:
                add(f'power_{factor}_{arch}',arch,128,2,fields,iterations=120,config={'power_W':Config().power_W*factor})
        for scale in (.25,.5,1.):
            for arch,fields in [('static',1),('covariance',4)]:
                add(f'density_{scale}_{arch}',arch,128,8,fields,iterations=120,cloud_scale=scale)
    elif group=='refinements':
        add('nonlinear_floquet_regularized','covariance',64,2,2,iterations=100,nonlinear_iterations=2,nonlinear_steps=16)
        for seed in (7,8,9):
            add(f'geometry_matched_{seed}','static',128,2,iterations=120,aperiodic=True,geometry_seed=seed)
        add('gravity_compensated_ib2','ib',512,2,ib_gravity_compensation=True)
        add('sdp32_slow','static',32,2,iterations=80,sdp=True,config={'stiffness_required':1.9378922925187383e-6})
    return specs


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--group',choices=['comparisons','temporal','scaling','power_density','refinements'],required=True)
    parser.add_argument('--only');parser.add_argument('--rerun',action='store_true')
    args=parser.parse_args();out=ROOT/'results'/'experiments'
    specs=plan(args.group)
    if args.only:specs=[s for s in specs if args.only in s['id']]
    write(ROOT/'results'/f'plan_{args.group}.json',specs)
    rows=[run_case(s,out,args.rerun) for s in specs]
    write(ROOT/'results'/f'summary_{args.group}.json',rows)
    print('EXPERIMENT_GROUP_COMPLETE',args.group,len(rows),flush=True)


if __name__=='__main__':main()
