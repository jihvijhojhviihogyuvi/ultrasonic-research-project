"""Integration controls: force fidelity, controller projection, and contacts."""
import math
import numpy as np
import torch
from physics import Config,AcousticSystem,targets
from optimization import Design
from dynamics import Surrogate,advance_frame,movement_protocol
from validation import IncrementalController


def main():
    c=Config(n=32,harmonics=2)
    s=AcousticSystem(c,targets(c,2));d=Design(s,'covariance',3)
    surrogate=Surrogate(s,degree=6,scale=.001,bins=8)
    controller=IncrementalController(d,surrogate)
    command,power=controller.command(np.zeros(6))
    expected=d.command_coefficients().detach().cpu().numpy()
    assert np.allclose(command,expected,atol=1e-12)
    assert math.isclose(power,float(s.power(d.coefficients())),rel_tol=1e-10)
    assert controller.rank==6
    controller.z[:]=.001
    u,_=controller.command(controller.z)
    fast=controller.polynomial(u);exact=surrogate.polynomial(s.response().cpu().numpy()[:,None]*u)
    assert np.allclose(fast,exact,atol=1e-8,rtol=1e-8)
    protocol=movement_protocol(np.zeros((2,3)))
    labels=set(x[2] for x in protocol)
    assert len(labels)>=25 and len(protocol)>=1200
    for distance in ('0.25','1'):
        for particle in range(2):
            for axis in range(3):
                for sign in (-1,1):assert f'{distance}mm_p{particle}_a{axis}_{sign}' in labels
    # Zero-field falling particle must accelerate down; a contact must fail.
    cc=Config(n=4,total_emitting_area=4*math.pi*.0045**2,harmonics=0)
    ss=AcousticSystem(cc,np.array([[0.,0.,0.],[0.0004,0.,0.]]))
    sur=Surrogate(ss,degree=2,bins=1)
    coeff=np.zeros((1,2,1,len(sur.exp)),complex);history=np.tile(coeff[None],(4,1,1,1,1))
    source=np.zeros((1,4),complex);qh=np.tile(sur.origin[None],(4,1,1))
    rr=advance_frame(sur.origin.copy(),np.zeros((2,3)),coeff.copy(),history,3,coeff,
        source,np.zeros_like(source),np.zeros((2,1)),qh,3,1,1e-5,0.,sur.origin,.001,
        sur.exp,sur.omega,sur.ell,cc.particle_radius,cc.f0,cc.f2,cc.rho,cc.sound_speed,
        cc.mass,cc.drag,cc.viscosity,cc.settling_s,cc.cycle_s,sur.origin,sur.origin,
        np.ones(2,bool),sur.origin,np.zeros(3),np.ones(4))
    assert rr[11]==1 and np.all(rr[1][:,2]<0)
    print('PHASE3_DYNAMICS_CONTROLS_PASSED')


if __name__=='__main__':main()
