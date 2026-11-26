"""Refine source-center delay bins on the first saved signed-motion segment."""
import math,json
import numpy as np
from experiments import ROOT,load_design,write
from dynamics import Surrogate,advance_frame
from validation import IncrementalController


def main():
    folder=ROOT/'results'/'motion'/'compare_2_static_cached64'
    data=np.load(folder/'raw.npz');frames=json.loads((folder/'frames.json').read_text())[:24]
    d,_=load_design(ROOT/'results'/'experiments'/'compare_2_static');s=d.system;c=s.c
    surrogate=Surrogate(s,degree=8,bins=16)
    controller=IncrementalController(d,surrogate)
    distance=np.linalg.norm(surrogate.origin[:,None,:]-s.positions[None,:,:],axis=-1)
    low=distance.min();high=distance.max()+1e-10;results=[];all_positions=[]
    for bins,samples in ((8,64),(16,64),(16,128)):
        assignment=np.clip(((distance-low)/(high-low)*bins).astype(int),0,bins-1)
        delays=np.tile((low+(np.arange(bins)+.5)*(high-low)/bins)/c.sound_speed,(s.m,1))
        def binned(command):
            return np.stack([np.einsum('lmnc,ln,mn->lmc',surrogate.transfer,command,assignment==b,optimize=True)
                             for b in range(bins)],axis=2)
        dt=c.cycle_s/samples;steps=round(.05/dt);response=s.response().cpu().numpy()
        command,_=controller.command(np.zeros(6));envelope=binned(command)*response[:,None,None,None]
        history=np.tile(envelope[None],(math.ceil(delays.max()/dt)+4,1,1,1,1));cursor=len(history)-1
        q=np.load(ROOT/'results'/'experiments'/'compare_2_static'/'design.npz')['equilibrium'].copy()
        v=np.zeros_like(q);source=response[:,None]*command;qh=np.tile(q[None],(4,1,1));qc=3
        oldtarget=surrogate.origin.copy();positions=[];maxerror=0.
        for i,row in enumerate(frames):
            controller.z=data['z'][i].copy();command,power=controller.command(controller.z)
            target=np.array(row['targets_m'])
            rr=advance_frame(q,v,envelope,history,cursor,binned(command),source,command,delays/dt,qh,qc,
                steps,dt,i*.05,surrogate.origin,surrogate.scale,surrogate.exp,surrogate.omega,surrogate.ell,
                c.particle_radius,c.f0,c.f2,c.rho,c.sound_speed,c.mass,c.drag,c.viscosity,c.settling_s,c.cycle_s,
                oldtarget,target,np.array([False,True]),surrogate.origin,np.zeros(3),s.power_coeff.max(0).values.cpu().numpy())
            q,v,envelope,history,cursor,source,qh,qc,error,held,peak,failure,energy,maxpower=rr
            assert not failure
            maxerror=max(maxerror,error);positions.append(q.copy());oldtarget=target
        positions=np.array(positions);all_positions.append(positions)
        results.append({'delay_bins':bins,'samples_per_cycle':samples,'duration_s':len(frames)*.05,
            'maximum_error_mm':maxerror*1000,
            'maximum_position_difference_from_first_mm':float(np.linalg.norm(positions-all_positions[0],axis=-1).max())*1000})
    out=ROOT/'results'/'convergence';write(out/'source_delay_bins.json',results)
    np.savez_compressed(out/'source_delay_bins.npz',positions=np.array(all_positions))
    print(json.dumps(results,indent=2),flush=True);print('SOURCE_DELAY_REFINEMENT_EXECUTED')


if __name__=='__main__':main()
