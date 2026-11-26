"""Refine the marginal sixteen-particle stationary candidate at fixed drives."""
import time
from dataclasses import replace
import numpy as np
import torch
from experiments import ROOT,load_design,write
from physics import AcousticSystem,source_jets,linear_floquet


def main():
    folder=ROOT/'results'/'experiments'/'compare_16_covariance';d,meta=load_design(folder)
    original=d.system;fields=d.fields().detach();duties=d.duties().detach()
    equilibrium=np.load(folder/'design.npz')['equilibrium'];rows=[]
    out=ROOT/'results'/'convergence'/'sixteen_particle.json'
    for harmonics,samples,refinement,fd in [(8,64,1,.000025),(8,128,1,.000025),
            (16,128,1,.000025),(8,128,2,.000025),(8,128,1,.0000125)]:
        start=time.perf_counter();c=replace(original.c,harmonics=harmonics,time_samples=samples,fd_step=fd)
        s=AcousticSystem(c,equilibrium,original.positions,original.normals,original.radii)
        if refinement>1:
            s.h=[[] for _ in range(4)]
            for f in s.frequencies:
                jets=source_jets(s.q,s.positions,s.normals,s.radii,float(f),c,refinement=refinement)
                for i,j in enumerate(jets):s.h[i].append(j)
            s.h=[torch.stack(v) for v in s.h]
        coef=s.coefficients(fields,duties)
        with torch.no_grad():f,k,_=s.matrices(coef);rho,_=linear_floquet(k,c)
        rows.append({'harmonics':harmonics,'time_samples':samples,'disk_refinement':refinement,'fd_step_m':fd,
           'weakest_stiffness_N_m':float(torch.linalg.eigvalsh((k.mean(0)+k.mean(0).T)/2).min()),
           'floquet_radius':float(rho),'gravity_residual_fraction':float(torch.linalg.vector_norm(f.mean(0),dim=1).max())/(c.mass*9.80665),
           'power_W':float(s.power(coef)),'runtime_s':time.perf_counter()-start})
        write(out,rows);print(rows[-1],flush=True)
    print('SIXTEEN_PARTICLE_REFINEMENT_EXECUTED',flush=True)


if __name__=='__main__':main()
