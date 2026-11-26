"""Runnable LOCAL linear dynamics test of independent trajectory tracking.

This uses force/stiffness at each commanded location; it is not a full nonlinear
trajectory simulator. Its validity must be checked with excursions << wavelength.
It includes continuous transducer-envelope settling and command-frame holds.
Spectral beats, camera latency, collisions, and scattering are not modeled.
"""
import argparse
import json
import math
import time
from pathlib import Path
import numpy as np
import torch
from simulator import Config, Design, Field, geometry, targets, optimize, equilibrium


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--architecture",choices=["static_force","covariance","ib","time_single"],default="static_force")
    p.add_argument("--n",type=int,default=128)
    p.add_argument("--particles",type=int,default=2)
    p.add_argument("--motion-hz",type=float,default=2)
    p.add_argument("--frame-hz",type=float,default=50)
    p.add_argument("--seconds",type=float,default=0.5)
    p.add_argument("--amplitude-mm",type=float,default=0.25)
    p.add_argument("--tau-ms",type=float,default=0.16)
    p.add_argument("--steps",type=int,default=120)
    p.add_argument("--out",type=Path,default=Path("tracking_result.json"))
    a=p.parse_args()
    c=Config(n=a.n,tau_s=a.tau_ms/1000,probe_iterations=a.steps)
    if c.n>=512: c.gap=.160;c.volume_side=.060
    pos,_=geometry(c); base=targets(c,a.particles)
    d=Design(a.architecture,base,pos,c); optimize(base,d,300)
    q,f,k,ok=equilibrium(base,d)
    if not ok: raise RuntimeError("initial equilibrium did not converge")
    state=torch.cat((q,torch.zeros_like(q)),dim=1)
    actual_envelope=d.commanded().detach()[0].clone()
    gamma=6*math.pi*c.viscosity*c.particle_radius
    records=[]; solve_times=[]
    dt_frame=1/a.frame_hz
    # Midpoint substeps resolve both switching and mechanical time scales.
    substeps=max(20,int(math.ceil(dt_frame/0.00010)))
    dt=dt_frame/substeps
    for frame in range(math.ceil(a.seconds*a.frame_hz)):
        t=frame*dt_frame
        desired=base.clone()
        desired[0,0]+=a.amplitude_mm/1000*math.sin(2*math.pi*a.motion_hz*t)
        desired[0,2]+=a.amplitude_mm/1000*(math.cos(2*math.pi*a.motion_hz*t)-1)
        start=time.perf_counter(); optimize(desired,d,a.steps)
        solve_times.append(time.perf_counter()-start)
        command=d.commanded().detach()
        field=Field(desired,pos,d.freqs,c)
        for sub in range(substeps):
            now=t+sub*dt
            slot=int((now%c.cycle_s)/c.cycle_s*len(command))%len(command)
            alpha=math.exp(-dt/c.tau_s) if c.tau_s>0 else 0
            half=command[slot]+(actual_envelope-command[slot])*math.sqrt(alpha)
            actual_envelope=command[slot]+(actual_envelope-command[slot])*alpha
            force,stiff=field.force_stiffness([(0,half[None,:],torch.ones(1))])
            force,stiff=force.detach(),stiff.detach()
            # Linearize around the commanded position, then integrate the affine
            # mechanics exactly over dt rather than an unstable Euler step.
            for particle in range(len(base)):
                aa=torch.zeros((7,7));aa[:3,3:6]=torch.eye(3)
                aa[3:6,:3]=-stiff[particle]/c.mass
                aa[3:6,3:6]=-gamma/c.mass*torch.eye(3)
                aa[3:6,6]=(force[particle]+stiff[particle]@desired[particle])/c.mass
                yy=torch.cat((state[particle],torch.ones(1)))
                state[particle]=(torch.matrix_exp(aa*dt)@yy)[:6]
        records.append({"t_s":t+dt_frame,"desired_m":desired.tolist(),"actual_m":state[:,:3].tolist()})
    actual=np.array([r["actual_m"] for r in records]); desired=np.array([r["desired_m"] for r in records])
    error=np.linalg.norm(actual-desired,axis=-1)
    held_drift=np.linalg.norm(actual[:,1:]-base.numpy()[None,1:],axis=-1)
    speed=np.linalg.norm(np.diff(actual,axis=0)/dt_frame,axis=-1)
    result={"model":"Tier 1 local linear transient, no sensing delay/scattering",
            "settings":vars(a)|{"out":str(a.out)},
            "rms_tracking_error_mm":float(np.sqrt(np.mean(error**2))*1000),
            "max_tracking_error_mm":float(error.max()*1000),
            "max_held_particle_drift_mm":float(held_drift.max()*1000) if held_drift.size else 0,
            "max_measured_step_speed_mm_s":float(speed.max()*1000),
            "commanded_path_speed_mm_s":2*math.pi*a.motion_hz*a.amplitude_mm,
            "mean_offline_optimizer_s":float(np.mean(solve_times)),
            "real_time_compute_feasible":max(solve_times)<dt_frame,
            "frames":records}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k!="frames"},indent=2))


if __name__=="__main__":main()
