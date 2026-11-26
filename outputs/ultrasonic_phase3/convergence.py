"""Fixed-drive spatial/spectral/time convergence and nonlinear periodic audit."""
import argparse,json,time
from pathlib import Path
from dataclasses import replace
import numpy as np
import torch
from numba import njit
from scipy.optimize import root
from physics import AcousticSystem,source_jets,linear_floquet
from experiments import ROOT,load_design,write
from dynamics import Surrogate,validate_surrogate,radiation_force
from threadpoolctl import threadpool_limits
_limit=threadpool_limits(limits=1)


@njit(cache=False)
def periodic_map(y,coeff,origin,scale,exp,omega,ell,radius,f0,f2,rho,c0,mass,drag,viscosity,cycle,steps,record=False):
    m=len(origin);q=origin+y[:3*m].reshape(m,3)*.0001
    v=y[3*m:].reshape(m,3)*.0001*(2*np.pi*15)
    dt=cycle/steps;history=np.zeros((steps,m,3))
    def acc(qq,vv,ph):
        f=radiation_force(qq,coeff,origin,scale,exp,omega,ell,ph,radius,f0,f2,rho,c0,mass,False)
        for i in range(m):
            reynolds=2*radius*rho*np.sqrt(np.sum(vv[i]**2))/viscosity
            f[i]-=drag*(1+.15*reynolds**.687)*vv[i]
        return f/mass
    for step in range(steps):
        ph=2*np.pi*step/steps
        k1q=v;k1v=acc(q,v,ph)
        k2q=v+.5*dt*k1v;k2v=acc(q+.5*dt*k1q,k2q,ph+np.pi/steps)
        k3q=v+.5*dt*k2v;k3v=acc(q+.5*dt*k2q,k3q,ph+np.pi/steps)
        k4q=v+dt*k3v;k4v=acc(q+dt*k3q,k4q,ph+2*np.pi/steps)
        q=q+dt*(k1q+2*k2q+2*k3q+k4q)/6
        v=v+dt*(k1v+2*k2v+2*k3v+k4v)/6
        history[step]=q
    return np.concatenate(((q-origin).reshape(-1)/.0001,v.reshape(-1)/(.0001*2*np.pi*15))),history


def nonlinear_audit(design,degree=8,steps=128):
    start=time.perf_counter();s=design.system;c=s.c
    surrogate=Surrogate(s,degree=degree)
    check=validate_surrogate(surrogate,design.coefficients().detach())
    coeff=surrogate.polynomial(design.coefficients().detach().cpu().numpy())
    args=(coeff,surrogate.origin,surrogate.scale,surrogate.exp,surrogate.omega,surrogate.ell,
        c.particle_radius,c.f0,c.f2,c.rho,c.sound_speed,c.mass,c.drag,c.viscosity,c.cycle_s,steps)
    def flow(y):return periodic_map(y,*args)[0]
    initial=np.zeros(6*s.m)
    solution=root(lambda y:flow(y)-y,initial,tol=1e-9)
    end,orbit=periodic_map(solution.x,*args)
    columns=[];eps=1e-4
    for j in range(len(initial)):
        delta=np.zeros_like(initial);delta[j]=eps
        columns.append((flow(solution.x+delta)-flow(solution.x-delta))/(2*eps))
    phi=np.stack(columns,axis=1)
    result={'degree':degree,'steps_per_cycle':steps,'shooting_converged':bool(solution.success),
            'closure_scaled_norm':float(np.linalg.norm(end-solution.x)),
            'floquet_radius_nonlinear':float(np.abs(np.linalg.eigvals(phi)).max()),
            'maximum_position_error_mm':float(np.linalg.norm(orbit-surrogate.origin[None],axis=-1).max())*1000,
            'surrogate':check,'runtime_s':time.perf_counter()-start,
            'valid_domain':bool(np.max(np.abs(orbit-surrogate.origin[None]))<=surrogate.scale)}
    return result,{'orbit':orbit,'monodromy':phi,'periodic_initial_scaled':solution.x}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--nonlinear',action='store_true');args=parser.parse_args()
    out=ROOT/'results'/'convergence';out.mkdir(parents=True,exist_ok=True)
    folder=ROOT/'results'/'experiments'/'compare_2_covariance';d,meta=load_design(folder);s=d.system
    rows=[]
    fields=d.fields().detach();duties=d.duties().detach()
    for harmonics in (8,16):
        ss=AcousticSystem(replace(s.c,harmonics=harmonics),s.q.cpu().numpy(),s.positions,s.normals,s.radii)
        coeff=ss.coefficients(fields,duties).detach()
        for samples in (64,128,256):
            start=time.perf_counter()
            with torch.no_grad():f,k,_=ss.matrices(coeff,samples=samples);rho,_=linear_floquet(k,ss.c)
            rows.append({'test':'spectral_time','harmonics':harmonics,'time_samples':samples,
                         'weakest_stiffness_N_m':float(torch.linalg.eigvalsh((k.mean(0)+k.mean(0).T)/2).min()),
                         'floquet_radius':float(rho),'gravity_residual_fraction':float(torch.linalg.vector_norm(f.mean(0),dim=1).max())/(s.c.mass*9.80665),
                         'acoustic_power_W':float(ss.power(coeff)),'runtime_s':time.perf_counter()-start})
            write(out/'results.json',rows)
    coeff=d.coefficients().detach()
    for refinement in (1,2,3):
        start=time.perf_counter();ss=AcousticSystem(s.c,s.q.cpu().numpy(),s.positions,s.normals,s.radii)
        ss.h=[[] for _ in range(4)]
        for f in ss.frequencies:
            jj=source_jets(ss.q,ss.positions,ss.normals,ss.radii,float(f),ss.c,refinement=refinement)
            for i,j in enumerate(jj):ss.h[i].append(j)
        ss.h=[torch.stack(j) for j in ss.h]
        with torch.no_grad():f,k,_=ss.matrices(coeff);rho,_=linear_floquet(k,ss.c)
        rows.append({'test':'disk_quadrature','refinement':refinement,'patches_per_disk':24*refinement**2,
           'weakest_stiffness_N_m':float(torch.linalg.eigvalsh((k.mean(0)+k.mean(0).T)/2).min()),
           'floquet_radius':float(rho),'runtime_s':time.perf_counter()-start})
        write(out/'results.json',rows)
    for method in ('none','born','foldy'):
        s.c.scattering=method
        with torch.no_grad():f,k,_=s.matrices(coeff);rho,_=linear_floquet(k,s.c)
        rows.append({'test':'scattering','method':method,
           'weakest_stiffness_N_m':float(torch.linalg.eigvalsh((k.mean(0)+k.mean(0).T)/2).min()),
           'floquet_radius':float(rho),'gravity_residual_fraction':float(torch.linalg.vector_norm(f.mean(0),dim=1).max())/(s.c.mass*9.80665)})
    s.c.scattering='born'
    for step in (.00005,.000025,.0000125):
        s.c.fd_step=step
        with torch.no_grad():f,k,_=s.matrices(coeff);rho,_=linear_floquet(k,s.c)
        rows.append({'test':'spatial_derivative','step_m':step,
           'weakest_stiffness_N_m':float(torch.linalg.eigvalsh((k.mean(0)+k.mean(0).T)/2).min()),'floquet_radius':float(rho)})
    write(out/'results.json',rows)
    if args.nonlinear:
        for name in ('compare_2_static','compare_2_covariance'):
            design,_=load_design(ROOT/'results'/'experiments'/name)
            for steps in (64,128):
                print('NONLINEAR',name,steps,flush=True)
                result,raw=nonlinear_audit(design,steps=steps)
                write(out/f'{name}_nonlinear_{steps}.json',result)
                np.savez_compressed(out/f'{name}_nonlinear_{steps}.npz',**raw)
                print(json.dumps(result),flush=True)
    print('CONVERGENCE_CALCULATIONS_COMPLETE',len(rows),flush=True)


if __name__=='__main__':main()
