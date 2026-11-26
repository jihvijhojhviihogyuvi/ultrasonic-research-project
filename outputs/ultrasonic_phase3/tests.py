"""Analytical and independent numerical controls for the Phase-3 physics."""
import math
from dataclasses import replace
import numpy as np
import torch
from physics import *


def require(condition,message):
    if not condition:raise AssertionError(message)


def test_green():
    q=torch.tensor([[.021,-.017,.063],[.018,.012,.092]])
    k=2*math.pi*40000/343
    jets=green_jets(q,k)
    for axis in range(3):
        offset=torch.zeros(3);offset[axis]=1e-7
        fd=(green_jets(q+offset,k,2)[2]-green_jets(q-offset,k,2)[2])/(2e-7)
        err=float(torch.linalg.vector_norm(fd-jets[3][...,axis])/torch.linalg.vector_norm(fd))
        require(err<1e-7,f"third pressure derivative {axis}: {err}")


def synthetic_system(harmonics=0):
    c=Config(n=4,harmonics=harmonics,scattering="none",time_samples=32,total_emitting_area=4*math.pi*.0045**2)
    s=AcousticSystem(c,np.zeros((1,3)))
    s.h=[torch.zeros_like(x) for x in s.h]
    return c,s


def test_standing_wave():
    c,s=synthetic_system();k=2*math.pi*c.frequency/c.sound_speed;p0=1000.
    s.h[1][0,0,0,2]=p0*k
    s.h[3][0,0,0,2,2,2]=-p0*k**3
    u=torch.zeros((1,c.n),dtype=torch.complex128);u[0,0]=1
    f,kk=s.primary_mean(u)
    expected=2*k*k*c.volume*p0*p0/(c.rho*c.sound_speed**2)*(c.f0/4+3*c.f2/8)
    require(abs(float(kk[0,2,2])/expected-1)<1e-12,"standing-wave stiffness coefficient")
    require(float(f[0,2])==-c.mass*9.80665,"gravity direction")
    ff,dd,_=s.matrices(u)
    require(abs(float(dd[0,2,2])/expected-1)<.0002,"coupled derivative versus analytic stiffness")


def test_beats():
    c,s=synthetic_system(1);p0=700.
    left=0;right=2
    km=2*math.pi*s.frequencies[left]/c.sound_speed;kp=2*math.pi*s.frequencies[right]/c.sound_speed
    for l,sign,k in ((left,-1,km),(right,1,kp)):
        s.h[0][l,0,0]=p0
        s.h[1][l,0,0,0]=sign*1j*k*p0
        s.h[2][l,0,0,0,0]=-k*k*p0
        s.h[3][l,0,0,0,0,0]=-sign*1j*k**3*p0
    coef=torch.zeros((3,c.n),dtype=torch.complex128);coef[[left,right],0]=1
    phases=torch.linspace(0,2*math.pi,33)[:-1]
    f=s.forces(s.q,s.jets(coef),phases)[0,:,0,0]
    expected=-c.volume*p0*p0*(kp+km)/(2*c.rho*c.sound_speed**2)*(c.f0+1.5*c.f2)*torch.sin(2*phases)
    require(torch.allclose(f,expected,rtol=1e-10,atol=1e-18),"cross-frequency beats retained")
    # Retarded modulation phase follows delta omega * travel time.
    r1=torch.tensor([[0.,0.,.050]]);r2=torch.tensor([[0.,0.,.090]])
    phase=(green_jets(r2,kp,0)[0]/green_jets(r2,km,0)[0])/(green_jets(r1,kp,0)[0]/green_jets(r1,km,0)[0])
    predicted=np.exp(1j*(kp-km)*.040)
    require(abs(complex(phase[0])-predicted)<1e-10,"retarded sideband propagation")


def test_power_geometry():
    reference=Config()
    for n in (32,64,128,256,512):
        c=replace(reference,n=n,harmonics=0)
        pos,norm,rad=geometry(c)
        require(abs(np.sum(math.pi*rad**2)-c.total_emitting_area)<1e-12,"equal emitting area")
        require(np.max(np.abs(pos[:,:2])+rad[:,None])<=c.aperture/2+1e-12,"aperture bounds")
    c=Config(n=32,harmonics=2);s=AcousticSystem(c,targets(c,2))
    rng=np.random.default_rng(9)
    raw=torch.tensor(rng.normal(size=(3,c.n))+1j*rng.normal(size=(3,c.n)))
    duties=torch.tensor([.2,.3,.5]);u,a=s.project(raw,duties)
    require(float(u.abs().max())<=1+1e-12,"peak limit")
    require(float(s.power(a))<=c.power_W*(1+1e-12),"radiated-power normalization")
    constant=torch.ones((3,c.n),dtype=torch.complex128)*.2
    aa=s.coefficients(constant,duties)
    require(float(aa[s.ell!=0].abs().max())<1e-14,"constant waveform has no modulation sidebands")
    # Force-drive differentiation versus independent finite difference.
    a=a.detach().clone().requires_grad_(True)
    f=s.primary_mean(a)[0][0,0]
    grad=torch.autograd.grad(f,a)[0]
    plus=a.detach().clone();minus=plus.clone();eps=1e-5
    plus[2,4]+=eps;minus[2,4]-=eps
    fd=float((s.primary_mean(plus)[0][0,0]-s.primary_mean(minus)[0][0,0])/(2*eps))
    require(abs(fd-float(grad[2,4].real))<1e-13,"force autodiff")
    # Born/Foldy comparison at a dilute configuration.
    with torch.no_grad():
        jets=s.jets(a.detach());phases=torch.tensor([.4])
        born=s.forces(s.q,jets,phases)
        s.c.scattering="foldy";foldy=s.forces(s.q,jets,phases)
        require(float(torch.linalg.vector_norm(born-foldy))/(c.mass*9.80665)<.001,"Born/Foldy dilute convergence")


def main():
    for test in (test_green,test_standing_wave,test_beats,test_power_geometry):
        test();print(test.__name__,"passed",flush=True)
    print("PHASE3_PHYSICS_TESTS_PASSED")


if __name__=="__main__":main()
