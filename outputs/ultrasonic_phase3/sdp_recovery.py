"""Recover a finite sequence from an executed ideal PSD covariance solution."""
import numpy as np
import torch
from dataclasses import replace
from experiments import ROOT,load_design,write
from physics import AcousticSystem,linear_floquet
from optimization import Design,optimize,equilibrium


def main():
    d,meta=load_design(ROOT/'results'/'experiments'/'sdp32_slow');s=d.system
    qmat=np.array(meta['sdp']['covariance_real'])+1j*np.array(meta['sdp']['covariance_imag'])
    values,vectors=np.linalg.eigh(qmat);order=np.argsort(values)[::-1]
    values=np.maximum(values[order],0);vectors=vectors[:,order]
    ff=torch.zeros((s.m,3));kk=torch.zeros((s.m,3,3))
    for j,value in enumerate(values):
        if value<1e-10:continue
        u=torch.tensor((np.sqrt(value)*vectors[:,j])[None],device=s.device)
        f,k=s.primary_mean(u);f[:,2]+=s.c.mass*9.80665;ff+=f;kk+=k
    ff[:,2]-=s.c.mass*9.80665
    rows=[{'test':'ideal_covariance','positive_eigenvalues':values.tolist(),
       'gravity_residual_fraction':float(torch.linalg.vector_norm(ff,dim=1).max())/(s.c.mass*9.80665),
       'weakest_stiffness_N_m':float(torch.linalg.eigvalsh(kk).min()),
       'required_stiffness_N_m':s.c.stiffness_required,'power_W':meta['sdp']['power_W']}]
    out=ROOT/'results'/'sdp_recovery';out.mkdir(parents=True,exist_ok=True)
    for rank in (1,2,4,8):
        ss=AcousticSystem(replace(s.c,harmonics=8),s.q.cpu().numpy(),s.positions,s.normals,s.radii)
        u=np.sqrt(rank*values[:rank,None])*vectors[:,:rank].T
        dd=Design(ss,'covariance',rank,initial=torch.tensor(u,device=ss.device))
        for stage in ('projected','optimized'):
            if stage=='optimized':optimize(dd,100)
            coef=dd.coefficients().detach()
            with torch.no_grad():f,k,_=ss.matrices(coef);rho,phi=linear_floquet(k,ss.c)
            rows.append({'test':stage,'rank':rank,
               'gravity_residual_fraction':float(torch.linalg.vector_norm(f.mean(0),dim=1).max())/(s.c.mass*9.80665),
               'weakest_stiffness_N_m':float(torch.linalg.eigvalsh((k.mean(0)+k.mean(0).T)/2).min()),
               'floquet_radius':float(rho),'power_W':float(ss.power(coef)),
               'required_stiffness_N_m':s.c.stiffness_required,'duties':dd.duties().detach().cpu().tolist()})
        np.savez_compressed(out/f'rank{rank}.npz',raw=dd.raw.detach().cpu().numpy(),timing=dd.timing.detach().cpu().numpy(),
            gauge=dd.gauge.detach().cpu().numpy(),coefficients=dd.coefficients().detach().cpu().numpy())
        write(out/'results.json',rows)
    print('SDP_RECOVERY_EXECUTED',len(rows),flush=True)


if __name__=='__main__':main()
