"""Frequency-resolved finite-aperture Rayleigh acoustic model, SI units.

Peak phasors exp(-i omega t). Finite disk quadrature, filtered command spectra,
frequency-dependent propagation (including envelope delays), coherent beat forces,
Rayleigh particle rescattering, and full coupled force derivatives. Homogeneous
air; streaming, finite-size/viscothermal scattering and walls remain model limits.
"""
from __future__ import annotations
import math
from dataclasses import dataclass,asdict
import numpy as np
import torch
from scipy.special import j1

torch.set_default_dtype(torch.float64)
torch.set_num_threads(1)


@dataclass
class Config:
    n:int=512
    volume_side:float=.060
    aperture:float=.160
    gap:float=.160
    total_emitting_area:float=512*math.pi*.0045**2
    frequency:float=40000.
    rho:float=1.2
    sound_speed:float=343.
    viscosity:float=1.81e-5
    particle_radius:float=.00025
    particle_density:float=30.
    compressibility_ratio:float=.01
    power_W:float=1.492951690156511
    source_gain_reference:float=2.
    settling_s:float=.00016
    cycle_s:float=.005
    harmonics:int=8
    time_samples:int=64
    smooth_s:float=.00008
    stiffness_required:float=1.7441030632668644e-5
    accuracy:float=.0001
    seed:int=7
    device:str="auto"
    memory_mb:int=1024
    scattering:str="born"
    fd_step:float=.000025

    def __post_init__(self):
        for name in ('volume_side','aperture','gap','total_emitting_area','frequency','rho',
                     'sound_speed','viscosity','particle_radius','particle_density','power_W',
                     'settling_s','cycle_s','stiffness_required','accuracy','fd_step'):
            if not np.isfinite(getattr(self,name)) or getattr(self,name)<=0:
                raise ValueError(f'{name} must be finite and positive')
        if self.harmonics<0 or self.time_samples<1 or self.memory_mb<1:
            raise ValueError('invalid spectral, temporal or memory budget')
        if self.scattering not in ('none','born','foldy'):raise ValueError('unknown scattering model')

    @property
    def mass(self):return 4*math.pi*self.particle_radius**3*self.particle_density/3
    @property
    def drag(self):return 6*math.pi*self.viscosity*self.particle_radius
    @property
    def f0(self):return 1-self.compressibility_ratio
    @property
    def f2(self):return 2*(self.particle_density-self.rho)/(2*self.particle_density+self.rho)
    @property
    def volume(self):return 4*math.pi*self.particle_radius**3/3
    @property
    def surface_velocity_reference(self):
        area=math.pi*.0045**2;k=2*math.pi*self.frequency/self.sound_speed
        return self.source_gain_reference*2*math.pi/(self.rho*self.sound_speed*k*area)


def device_of(c):
    return torch.device("cuda" if c.device=="auto" and torch.cuda.is_available() else
                        "cpu" if c.device=="auto" else c.device)


def geometry(c,aperiodic=False):
    """Fixed total emitting area, aperture, and box at every source count."""
    if c.n<4 or c.n%2:raise ValueError("source count must be even and >=4")
    per=c.n//2;side=math.ceil(math.sqrt(per));pitch=c.aperture/side
    xy=(np.arange(side)-(side-1)/2)*pitch
    candidates=np.array([(x,y) for x in xy for y in xy])
    if per<len(candidates):
        # Farthest-point subset retains aperture coverage for 64/256 sources.
        selected=[0]
        while len(selected)<per:
            distance=np.min(np.sum((candidates[:,None]-candidates[selected][None])**2,axis=-1),axis=1)
            distance[selected]=-1;selected.append(int(np.argmax(distance)))
        candidates=candidates[selected]
    pos=np.concatenate([np.column_stack([candidates,np.full(per,z)]) for z in (-c.gap/2,c.gap/2)])
    radius=math.sqrt(c.total_emitting_area/(c.n*math.pi))
    if 2*radius>pitch*(1+1e-8):raise ValueError("emitting disks overlap in proposed geometry")
    if aperiodic:
        rng=np.random.default_rng(c.seed)
        amplitude=min(.0003,(pitch-2*radius)/3)
        pos[:,:2]+=rng.uniform(-amplitude,amplitude,(c.n,2))
        pos[:,:2]=np.clip(pos[:,:2],-c.aperture/2+radius,c.aperture/2-radius)
    normal=np.zeros_like(pos);normal[:,2]=-np.sign(pos[:,2])
    return pos,normal,np.full(c.n,radius)


def targets(c,m,min_spacing=.0045):
    rng=np.random.default_rng(c.seed+m);points=[]
    for _ in range(200000):
        q=rng.uniform(-c.volume_side/2+.0015,c.volume_side/2-.0015,3)
        if all(np.linalg.norm(q-p)>=min_spacing for p in points):
            points.append(q)
            if len(points)==m:return np.array(points)
    raise ValueError("target packing failed at the specified spacing")


def disk_offsets(radius,normal,refinement=1):
    ratio=radius/.0045
    rings=max(3,round(3*ratio))*refinement
    angles=max(8,4*math.ceil(2*ratio))*refinement
    normal=np.asarray(normal)/np.linalg.norm(normal)
    helper=np.array([1.,0.,0.]) if abs(normal[0])<.9 else np.array([0.,1.,0.])
    e1=np.cross(normal,helper);e1/=np.linalg.norm(e1);e2=np.cross(normal,e1)
    return np.array([radius*math.sqrt((r+.5)/rings)*
                     (math.cos(2*math.pi*(j+.5*(r%2))/angles)*e1+
                      math.sin(2*math.pi*(j+.5*(r%2))/angles)*e2)
                     for r in range(rings) for j in range(angles)])


def green_jets(delta,k,order=3):
    """Analytic derivatives through third order of exp(ikr)/r."""
    r=torch.linalg.vector_norm(delta,dim=-1).clamp_min(1e-12)
    nn=delta/r[...,None]
    g=torch.exp(1j*k*r)/r
    result=[g]
    if order>=1:result.append(g[...,None]*(1j*k-1/r)[...,None]*nn)
    if order>=2:
        eye=torch.eye(3,device=delta.device)
        aa=(-k*k-3j*k/r+3/r**2);bb=1j*k/r-1/r**2
        result.append(g[...,None,None]*(aa[...,None,None]*nn[..., :,None]*nn[...,None,:]+bb[...,None,None]*eye))
    if order>=3:
        cc=-1j*k**3+6*k*k/r+15j*k/r**2-15/r**3
        dd=-k*k/r-3j*k/r**2+3/r**3
        cube=nn[..., :,None,None]*nn[...,None,:,None]*nn[...,None,None,:]
        sym=(nn[..., :,None,None]*eye[None,:,:]+nn[...,None,:,None]*eye[:,None,:]
             +nn[...,None,None,:]*eye[:,:,None])
        result.append(g[...,None,None,None]*(cc[...,None,None,None]*cube+dd[...,None,None,None]*sym))
    return result


def source_jets(q,positions,normals,radii,f,c,order=3,refinement=1):
    """Memory-bounded finite disk integration, arbitrary configured normals."""
    dev=q.device;outputs=[[] for _ in range(order+1)]
    k=2*math.pi*f/c.sound_speed
    if np.allclose(radii,radii[0]):
        patches=torch.tensor(np.stack([positions[n]+disk_offsets(radii[n],normals[n],refinement)
                                       for n in range(len(positions))]),device=dev)
        gains=torch.tensor(c.source_gain_reference*(radii/.0045)**2*(f/c.frequency),device=dev)
        for start in range(0,len(q),8):
            pieces=[[] for _ in outputs]
            for source_start in range(0,len(positions),16):
                source=patches[source_start:source_start+16]
                jets=green_jets(q[start:start+8,None,None,:]-source[None],k,order)
                gain=gains[source_start:source_start+16]
                for d,j in enumerate(jets):
                    pieces[d].append(j.mean(2)*gain.reshape((1,-1)+(1,)*d))
            for d in range(order+1):outputs[d].append(torch.cat(pieces[d],dim=1))
        return [torch.cat(parts,dim=0) for parts in outputs]
    for n in range(len(positions)):
        patches=torch.tensor(positions[n]+disk_offsets(radii[n],normals[n],refinement),device=dev)
        pieces=[[] for _ in outputs]
        for start in range(0,len(q),16):
            jets=green_jets(q[start:start+16,None,:]-patches[None,:,:],k,order)
            gain=c.source_gain_reference*(radii[n]/.0045)**2*(f/c.frequency)
            for d,j in enumerate(jets):pieces[d].append(gain*j.mean(1))
        for d in range(order+1):outputs[d].append(torch.cat(pieces[d]))
    return [torch.stack(parts,dim=1) for parts in outputs]


def temporal_weights(harmonics,duties):
    """Exact Fourier coefficients of piecewise constant cyclic field commands."""
    center=torch.cumsum(duties,0)-duties/2
    ell=harmonics[:,None]
    return duties[None,:]*torch.sinc(ell*duties[None,:])*torch.exp(2j*math.pi*ell*center[None,:])


class AcousticSystem:
    def __init__(self,c,q,positions=None,normals=None,radii=None,cycle_s=None,harmonics=None):
        self.c=c;self.device=device_of(c)
        self.q=torch.as_tensor(q,device=self.device,dtype=torch.float64)
        self.m=len(q)
        if positions is None:positions,normals,radii=geometry(c)
        self.positions=np.asarray(positions);self.normals=np.asarray(normals);self.radii=np.asarray(radii)
        self.cycle_s=c.cycle_s if cycle_s is None else cycle_s
        limit=c.harmonics if harmonics is None else harmonics
        limit=min(limit,int(.20*c.frequency*self.cycle_s))
        self.ell=torch.arange(-limit,limit+1,device=self.device)
        self.frequencies=c.frequency+self.ell.detach().cpu().numpy()/self.cycle_s
        if np.any(self.frequencies<=0):raise ValueError("nonpositive acoustic sideband frequency")
        estimate=len(self.ell)*self.m*c.n*40*16/1024**2
        if estimate>c.memory_mb:raise MemoryError(f"source-jet cache requires {estimate:.0f} MiB")
        self.h=[[] for _ in range(4)]
        for f in self.frequencies:
            jets=source_jets(self.q,self.positions,self.normals,self.radii,float(f),c)
            for d,j in enumerate(jets):self.h[d].append(j)
        self.h=[torch.stack(j) for j in self.h]
        self.omega=torch.tensor(2*np.pi*self.frequencies,device=self.device)
        area=math.pi*self.radii**2
        kk=2*np.pi*self.frequencies[:,None]*self.radii[None,:]/c.sound_speed
        resistance=1-2*j1(2*kk)/(2*kk)
        self.power_coeff=torch.tensor(.5*c.rho*c.sound_speed*area[None,:]*
                c.surface_velocity_reference**2*resistance,device=self.device)

    def coefficients(self,fields,duties,smooth_s=None):
        command=self.command_coefficients(fields,duties,smooth_s)
        return self.response()[:,None]*command

    def response(self):
        delta=2*math.pi*self.ell/self.cycle_s
        return 1/(1-1j*delta*self.c.settling_s)

    def command_coefficients(self,fields,duties,smooth_s=None):
        weights=temporal_weights(self.ell,duties)
        delta=2*math.pi*self.ell/self.cycle_s
        smoothing=self.c.smooth_s if smooth_s is None else smooth_s
        smooth=torch.exp(-.5*(delta*smoothing)**2)
        return smooth[:,None]*(weights@fields)

    def power(self,coeff):return (self.power_coeff*coeff.abs().square()).sum()

    def project(self,fields,duties,smooth_s=None):
        fields=fields/fields.abs().clamp_min(1)
        command=self.command_coefficients(fields,duties,smooth_s)
        # Sampled maximum plus a derivative bound controls between-sample peaks
        # in the actually modeled finite-bandwidth command, not just slot values.
        samples=256
        phases=torch.arange(samples,device=self.device)*2*math.pi/samples
        envelope=torch.exp(-1j*phases[:,None]*self.ell[None,:])@command
        derivative=(self.ell.abs()[:,None]*command.abs()).sum(0)
        peak=envelope.abs().max(0).values+math.pi/samples*derivative
        fields=fields/peak.clamp_min(1)[None,:]
        coeff=self.coefficients(fields,duties,smooth_s)
        factor=torch.sqrt(self.c.power_W/self.power(coeff).clamp_min(1e-30)).clamp_max(1)
        return fields*factor,coeff*factor

    def jets(self,coeff):
        return [torch.einsum("lmn,ln->lm",self.h[0],coeff),
                torch.einsum("lmna,ln->lma",self.h[1],coeff),
                torch.einsum("lmnab,ln->lmab",self.h[2],coeff),
                torch.einsum("lmnabc,ln->lmabc",self.h[3],coeff)]

    def external_at(self,positions,jets):
        """Third-order local pressure expansion; exact jets at reference sites."""
        q=positions if positions.ndim==3 else positions[None,:,:]
        d=(q-self.q[None,:,:]).to(torch.complex128)
        p,g,h,t=jets
        pp=p[None]+torch.einsum("lmi,zmi->zlm",g,d)+.5*torch.einsum("lmij,zmi,zmj->zlm",h,d,d)+(
            torch.einsum("lmijk,zmi,zmj,zmk->zlm",t,d,d,d)/6)
        gg=g[None]+torch.einsum("lmij,zmj->zlmi",h,d)+.5*torch.einsum("lmijk,zmj,zmk->zlmi",t,d,d)
        hh=h[None]+torch.einsum("lmijk,zmk->zlmij",t,d)
        return pp,gg,hh

    def scatter(self,q,p,g,h):
        if self.c.scattering=="none" or self.m<2:return p,g,h
        q=q if q.ndim==3 else q[None,:,:]
        mask=torch.eye(self.m,dtype=torch.bool,device=self.device)[None,:,:]
        delta=q[:,:,None,:]-q[:,None,:,:]
        delta=torch.where(mask[...,None],torch.tensor([1.,0.,0.],device=self.device),delta)
        scattered=[]
        for l,f in enumerate(self.frequencies):
            k=2*math.pi*f/self.c.sound_speed
            kernels=green_jets(delta,k,3)
            kernels=[torch.where(mask[(...,)+(None,)*d],torch.zeros_like(j),j) for d,j in enumerate(kernels)]
            g0,g1,g2,g3=kernels
            aa=-self.c.f0*k*k*self.c.particle_radius**3/3
            bb=-self.c.f2*self.c.particle_radius**3/2
            pc=p[:,l];gc=g[:,l];hc=h[:,l]
            if self.c.scattering=="foldy":
                b=len(q);size=4*self.m
                ss=torch.zeros((b,self.m,4,self.m,4),device=self.device,dtype=torch.complex128)
                ss[:,:,0,:,0]=aa*g0
                ss[:,:,0,:,1:]=bb*k*g1
                ss[:,:,1:,:,0]=(aa/k*g1).permute(0,1,3,2)
                ss[:,:,1:,:,1:]=bb*g2.permute(0,1,3,2,4)
                rhs=torch.cat((pc[...,None],gc/k),-1).reshape(b,size,1)
                zz=torch.linalg.solve(torch.eye(size,device=self.device)-ss.reshape(b,size,size),rhs).reshape(b,self.m,4)
                pc=zz[...,0];gc=zz[...,1:]*k
            ps=aa*torch.einsum("zij,zj->zi",g0,pc)+bb*torch.einsum("zija,zja->zi",g1,gc)
            gs=aa*torch.einsum("zija,zj->zia",g1,pc)+bb*torch.einsum("zijab,zjb->zia",g2,gc)
            hs=aa*torch.einsum("zijab,zj->ziab",g2,pc)+bb*torch.einsum("zijabc,zjc->ziab",g3,gc)
            scattered.append((p[:,l]+ps,g[:,l]+gs,hc+hs))
        return tuple(torch.stack([s[d] for s in scattered],dim=1) for d in range(3))

    def forces(self,positions,jets,phases):
        q=positions if positions.ndim==3 else positions[None,:,:]
        p,g,h=self.scatter(q,*self.external_at(q,jets))
        phase=torch.exp(-1j*phases[:,None]*self.ell[None,:])
        pp=torch.einsum("tl,zlm->ztm",phase,p)
        gg=torch.einsum("tl,zlmi->ztmi",phase,g)
        vv=torch.einsum("tl,zlmi->ztmi",phase,g/(1j*self.omega[None,:,None,None]*self.c.rho))
        dv=torch.einsum("tl,zlmij->ztmij",phase,h/(1j*self.omega[None,:,None,None,None]*self.c.rho))
        alpha=self.c.volume*self.c.f0/(4*self.c.rho*self.c.sound_speed**2)
        beta=3*self.c.volume*self.c.f2*self.c.rho/8
        force=-2*alpha*(pp.conj()[...,None]*gg).real+2*beta*(vv.conj()[..., :,None]*dv).real.sum(-2)
        force=force.clone();force[...,2]-=self.c.mass*9.80665
        return force

    def matrices(self,coeff,samples=None,coupled=True):
        samples=1 if len(self.ell)==1 else max(self.c.time_samples if samples is None else samples,2*int(self.ell.abs().max())+1)
        phases=(torch.arange(samples,device=self.device)+.5)*2*math.pi/samples
        jets=self.jets(coeff)
        force=self.forces(self.q,jets,phases)[0]
        columns=[];dim=3*self.m;step=self.c.fd_step
        # Small batches bound memory; all off-diagonal interaction derivatives kept.
        for start in range(0,dim,3):
            indices=list(range(start,min(start+3,dim)))
            dq=torch.eye(dim,device=self.device)[indices].reshape(len(indices),self.m,3)*step
            fp=self.forces(self.q[None]+dq,jets,phases)
            fm=self.forces(self.q[None]-dq,jets,phases)
            columns.append((-(fp-fm)/(2*step)).permute(1,2,3,0).reshape(samples,dim,len(indices)))
        stiffness=torch.cat(columns,dim=-1)
        return force,stiffness,phases

    def primary_mean(self,coeff):
        p,g,h,t=self.jets(coeff);v=g/(1j*self.omega[:,None,None]*self.c.rho)
        dv=h/(1j*self.omega[:,None,None,None]*self.c.rho)
        ddv=t/(1j*self.omega[:,None,None,None,None]*self.c.rho)
        alpha=self.c.volume*self.c.f0/(4*self.c.rho*self.c.sound_speed**2)
        beta=3*self.c.volume*self.c.f2*self.c.rho/8
        force=(-2*alpha*(p.conj()[...,None]*g).real+2*beta*(v.conj()[..., :,None]*dv).real.sum(-2)).sum(0)
        stiffness=(2*alpha*((g.conj()[...,:,None]*g[...,None,:]).real+(p.conj()[...,None,None]*h).real)
            -2*beta*((dv.conj()[..., :, :,None]*dv[..., :,None,:]).real.sum(-3)+
                    (v.conj()[..., :,None,None]*ddv).real.sum(-3))).sum(0)
        force=force.clone();force[:,2]-=self.c.mass*9.80665
        return force,stiffness


def linear_floquet(stiffness,c):
    n=stiffness.shape[-1];samples=len(stiffness);omega_ref=2*math.pi*15
    aa=torch.zeros((samples,2*n,2*n),device=stiffness.device)
    aa[:,:n,n:]=torch.eye(n,device=stiffness.device)*omega_ref
    aa[:,n:,:n]=-stiffness/(c.mass*omega_ref)
    aa[:,n:,n:]=-c.drag/c.mass*torch.eye(n,device=stiffness.device)
    transitions=torch.matrix_exp(aa*c.cycle_s/samples)
    phi=torch.eye(2*n,device=stiffness.device)
    for transition in transitions:phi=transition@phi
    return torch.linalg.eigvals(phi).abs().max(),phi


def affine_periodic_orbit(force,stiffness,c):
    """Full coupled affine periodic orbit, not just the averaged equilibrium."""
    n=stiffness.shape[-1];samples=len(stiffness);ref=2*math.pi*15;scale=.0001
    aa=torch.zeros((samples,2*n+1,2*n+1),device=stiffness.device)
    aa[:,:n,n:2*n]=torch.eye(n,device=stiffness.device)*ref
    aa[:,n:2*n,:n]=-stiffness/(c.mass*ref)
    aa[:,n:2*n,n:2*n]=-c.drag/c.mass*torch.eye(n,device=stiffness.device)
    aa[:,n:2*n,-1]=force.reshape(samples,n)/(c.mass*scale*ref)
    transitions=torch.matrix_exp(aa*c.cycle_s/samples)
    phi=torch.eye(2*n+1,device=stiffness.device)
    for tr in transitions:phi=tr@phi
    rho=torch.linalg.eigvals(phi[:2*n,:2*n]).abs().max()
    fixed=torch.linalg.solve(torch.eye(2*n,device=stiffness.device)-phi[:2*n,:2*n],phi[:2*n,-1])
    y=torch.cat((fixed,torch.ones(1,device=stiffness.device)));orbit=[]
    for tr in transitions:y=tr@y;orbit.append(y[:n]*scale)
    return rho,torch.stack(orbit),phi
