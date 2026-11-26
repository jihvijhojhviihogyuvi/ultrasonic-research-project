"""Coupled periodic stability at gravity-balanced density-sweep equilibria."""
import math,time
import numpy as np
import torch
from experiments import ROOT,load_design,write
from optimization import equilibrium
from physics import linear_floquet


def main():
    out=ROOT/'results/density_stability';out.mkdir(parents=True,exist_ok=True);rows=[]
    for scale in (.25,.5,1.):
        for arch in ('static','covariance'):
            name=f'density_{scale}_{arch}';start=time.perf_counter()
            d,meta=load_design(ROOT/'results/experiments'/name);s=d.system;coeff=d.coefficients().detach()
            q,f,k,ok,current=equilibrium(s,coeff,iterations=12)
            with torch.no_grad():ft,kt,_=current.matrices(coeff);rho,phi=linear_floquet(kt,s.c)
            row={'case':name,'architecture':arch,'cloud_density_per_cm3':meta['cloud_density_per_cm3'],
               'equilibrium_converged':ok,'equilibrium_error_mm':float(torch.linalg.vector_norm(q-s.q,dim=1).max())*1000,
               'gravity_residual_fraction':float(torch.linalg.vector_norm(f,dim=1).max())/(s.c.mass*9.80665),
               'weakest_coupled_stiffness_N_m':float(torch.linalg.eigvalsh((k+k.T)/2).min()),
               'floquet_radius':float(rho),'growth_rate_per_s':math.log(float(rho))/s.c.cycle_s,
               'runtime_s':time.perf_counter()-start,'note':'Floquet value is not an equilibrium stability certificate if root did not converge'}
            rows.append(row);write(out/'results.json',rows)
            np.savez_compressed(out/f'{name}.npz',equilibrium=q.cpu().numpy(),force=ft.cpu().numpy(),stiffness=kt.cpu().numpy(),monodromy=phi.cpu().numpy())
            print(row,flush=True)
    print('COUPLED_DENSITY_STABILITY_EXECUTED')


if __name__=='__main__':main()
