"""Exactly match realized acoustic power of fixed static/covariance drives."""
import numpy as np
import torch
from experiments import ROOT,load_design,write
from optimization import equilibrium
from physics import linear_floquet


def main():
    rows=[];out=ROOT/'results'/'matched_power';out.mkdir(parents=True,exist_ok=True)
    for m in (2,8,16):
        designs=[load_design(ROOT/'results'/'experiments'/f'compare_{m}_{arch}')[0] for arch in ('static','covariance')]
        power=min(float(d.system.power(d.coefficients().detach())) for d in designs)
        for d in designs:
            s=d.system;coeff=d.coefficients().detach();coeff*=np.sqrt(power/float(s.power(coeff)))
            q,f,k,ok,current=equilibrium(s,coeff,iterations=8)
            with torch.no_grad():ff,kk,_=current.matrices(coeff);rho,_=linear_floquet(kk,s.c)
            rows.append({'particles':m,'architecture':d.architecture,'transducers':s.c.n,
               'matched_realized_power_W':float(s.power(coeff)),
               'weakest_stiffness_N_m':float(torch.linalg.eigvalsh((k+k.T)/2).min()),
               'equilibrium_error_mm':float(torch.linalg.vector_norm(q-s.q,dim=1).max())*1000,
               'equilibrium_converged':ok,'floquet_radius':float(rho),
               'note':'fixed-field attenuation followed by equilibrium search; not a new optimized capacity'})
            np.savez_compressed(out/f'M{m}_{d.architecture}.npz',coefficients=coeff.cpu().numpy(),equilibrium=q.cpu().numpy())
            write(out/'results.json',rows);print(rows[-1],flush=True)
    print('MATCHED_REALIZED_POWER_EXECUTED',flush=True)


if __name__=='__main__':main()
