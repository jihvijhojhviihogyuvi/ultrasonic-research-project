"""Physics and numerical controls that can fail on sign/factor/derivative errors."""
import math
import numpy as np
import torch
from simulator import (Config, Field, Design, ARCHITECTURES, geometry, targets,
                       propagator, envelope, isolated_piston_power, force_rank)


def check(ok, message):
    if not ok: raise AssertionError(message)


def main():
    c=Config(n=32,gap=0.060,volume_side=0.018)
    pos,_=geometry(c)
    x=targets(c,2)
    # Compare exact quadrature pressure derivatives with a separately computed
    # spatial finite difference. This catches phase or radial-spreading errors.
    h,d=propagator(x,pos,c.frequency,c)
    for axis in range(3):
        shift=torch.zeros(3); shift[axis]=1e-7
        hp,_=propagator(x+shift,pos,c.frequency,c)
        hm,_=propagator(x-shift,pos,c.frequency,c)
        rel=float(torch.linalg.vector_norm((hp-hm)/2e-7-d[:,:,axis])/
                  torch.linalg.vector_norm(d[:,:,axis]))
        check(rel<1e-7,f"propagator derivative {axis}: {rel}")
    many=x.repeat(70,1)
    hb,db=propagator(many,pos,c.frequency,c)
    ha,da=propagator(many[:70],pos,c.frequency,c)
    hc,dc=propagator(many[70:],pos,c.frequency,c)
    check(torch.allclose(hb,torch.cat([ha,hc]),rtol=1e-12,atol=1e-12)
          and torch.allclose(db,torch.cat([da,dc]),rtol=1e-12,atol=1e-12),
          "batched source propagation preserves the physical field")

    # Independently known 1-D standing wave: p=P sin(kz), v=P cos(kz)/(rho c).
    # The Gor'kov potential must have a stable PRESSURE NODE for positive contrast.
    field=Field(torch.tensor([[0.,0.,0.]]),pos,[c.frequency],c)
    k=2*math.pi*c.frequency/c.sound_speed; p0=1000.
    z=field.q[0,:,2]
    field.hs[0]=torch.zeros_like(field.hs[0])
    field.ds[0]=torch.zeros_like(field.ds[0])
    field.hs[0][0,:,0]=p0*torch.sin(k*z).to(torch.complex128)
    field.ds[0][0,:,0,2]=p0*k*torch.cos(k*z).to(torch.complex128)
    u=torch.zeros((1,c.n),dtype=torch.complex128); u[0,0]=1
    f,stiff=field.force_stiffness([(0,u,torch.ones(1))])
    vol=4*math.pi*c.particle_radius**3/3
    f1=1-c.compressibility_ratio
    f2=2*(c.particle_density-c.air_density)/(2*c.particle_density+c.air_density)
    analytical=2*k*k*vol*p0*p0/(c.air_density*c.sound_speed**2)*(f1/4+3*f2/8)
    check(abs(float(stiff[0,2,2])/analytical-1)<0.001,"standing-wave stiffness factor/sign")
    check(abs(float(f[0,2])+c.mass*9.80665)<1e-18,"gravity force direction")
    check(float(stiff[0,0,0])==0 and float(stiff[0,1,1])==0,"plane-wave lateral freedom")
    # Positive trace is insufficient. A known saddle must fail full PSD.
    saddle=torch.diag(torch.tensor([2.,2.,-1.]))
    check(float(torch.trace(saddle))>0 and float(torch.linalg.eigvalsh(saddle).min())<0,
          "negative control: Laplacian cannot certify stability")

    # Incoherent slot potentials add; coherent phasors interfere. Wrongly adding
    # fields before temporal averaging would fail these assertions.
    actual=Field(x,pos,[c.frequency],c)
    one=Design("ib",x,pos,c).commanded()
    incoherent=actual.potentials([(0,torch.cat([one,one]),torch.tensor([.5,.5]))])
    coherent=actual.potentials([(0,2*one,torch.ones(1))])
    single=actual.potentials([(0,one,torch.ones(1))])
    check(torch.allclose(incoherent,single,rtol=1e-12,atol=1e-25),"incoherent averaging")
    check(torch.allclose(coherent,4*single,rtol=1e-12,atol=1e-25),"coherent interference")

    # A drive gradient is independently finite-differenced, including envelope.
    design=Design("covariance",x,pos,c)
    with torch.no_grad(): design.raw.mul_(0.8)  # away from the power projection kink
    f,_=actual.force_stiffness(design.drives())
    val=f[0,0]/(c.mass*9.80665)
    grad=torch.autograd.grad(val,design.raw)[0]
    idx=(0,3,0); orig=float(design.raw[idx].detach()); eps=1e-5
    with torch.no_grad(): design.raw[idx]=orig+eps
    plus=float(actual.force_stiffness(design.drives())[0][0,0].detach())/(c.mass*9.80665)
    with torch.no_grad(): design.raw[idx]=orig-eps
    minus=float(actual.force_stiffness(design.drives())[0][0,0].detach())/(c.mass*9.80665)
    with torch.no_grad(): design.raw[idx]=orig
    check(abs((plus-minus)/(2*eps)-float(grad[idx]))<1e-5,"drive autodiff gradient")

    # Constant mode must remain constant through finite settling. Opposite modes
    # cancel strongly as tau becomes large; blindly averaging commands must fail.
    const=one.repeat(3,1); smooth,w=envelope(const,c)
    check(torch.allclose(smooth,const[0].expand_as(smooth),rtol=1e-12,atol=1e-12),"constant envelope")
    slow=Config(**{**c.__dict__,"tau_s":1.0})
    sm,w=envelope(torch.cat([one,-one]),slow)
    check(float(sm.abs().max().detach())<0.01,"negative control: settling destroys fast phase swaps")

    for arch in ARCHITECTURES:
        d=Design(arch,x,pos,c)
        with torch.no_grad(): d.raw.mul_(10)
        u=d.commanded()
        check(d.effort()<=c.n*c.power_fraction*(1+1e-12),f"power cap {arch}")
        crest=u.abs().sum(0) if arch=="multifrequency" else u.abs().max(0).values
        check(float(crest.max().detach())<=1+1e-12,f"peak/crest cap {arch}")
    check(isolated_piston_power(c)>0,"positive acoustic power estimate")
    rank,cond=force_rank(x,Design("ib",x,pos,c))
    check(0<=rank<=6 and 0<=cond<=1,"force-Jacobian rank bounds")
    print("PHYSICS_CONTROLS_PASSED: derivatives, standing wave, gravity, saddle, interference, envelope, power, rank")


if __name__=="__main__": main()
