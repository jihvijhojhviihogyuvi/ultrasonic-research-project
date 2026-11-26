"""Compiled nonlinear particle tests using validated spectral pressure surrogates.

Pressure is fit to a high-order polynomial from the finite-disk model; forces
retain all beat phases, gravity, finite-Re drag, and Rayleigh pair scattering.
The source envelope is integrated with its complex periodic state; source-center
travel delays are binned, refined in convergence tests. Periodic aperture delays
are exact in the frequency transfers. Between-frame derivative-of-envelope terms
in the radiation model are a narrowband approximation, explicitly reported.
"""
from __future__ import annotations
import math,time
import numpy as np
import torch
from numba import njit
from physics import source_jets


def exponents(degree):
    return np.array([(a,b,c) for n in range(degree+1) for a in range(n+1)
                    for b in range(n-a+1) for c in [n-a-b]],dtype=np.int64)


@njit(cache=False)
def monomials(q,origin,scale,exp):
    m=len(q);nc=len(exp);degree=int(exp.max())
    mono=np.empty((m,nc));grad=np.zeros((m,nc,3));hess=np.zeros((m,nc,3,3))
    for i in range(m):
        powers=np.ones((3,degree+1))
        for a in range(3):
            x=(q[i,a]-origin[i,a])/scale
            for n in range(1,degree+1):powers[a,n]=powers[a,n-1]*x
        for j in range(nc):
            e=exp[j];mono[i,j]=powers[0,e[0]]*powers[1,e[1]]*powers[2,e[2]]
            for a in range(3):
                if e[a]>0:
                    value=e[a]/scale
                    for b in range(3):value*=powers[b,e[b]-(1 if a==b else 0)]
                    grad[i,j,a]=value
                for b in range(3):
                    if e[a]>0 and e[b]-(1 if a==b else 0)>0:
                        value=e[a]*(e[b]-(1 if a==b else 0))/scale**2
                        for d in range(3):value*=powers[d,e[d]-(1 if a==d else 0)-(1 if b==d else 0)]
                        hess[i,j,a,b]=value
    return mono,grad,hess


@njit(cache=False)
def fields_at(q,coeff,origin,scale,exp):
    mono,dm,hm=monomials(q,origin,scale,exp)
    nl,m,nc=coeff.shape;p=np.zeros((nl,m),np.complex128)
    g=np.zeros((nl,m,3),np.complex128);h=np.zeros((nl,m,3,3),np.complex128)
    for l in range(nl):
        for i in range(m):
            for j in range(nc):
                z=coeff[l,i,j];p[l,i]+=z*mono[i,j]
                for a in range(3):
                    g[l,i,a]+=z*dm[i,j,a]
                    for b in range(3):h[l,i,a,b]+=z*hm[i,j,a,b]
    return p,g,h


@njit(cache=False)
def scattered_fields(q,p,g,h,omega,radius,f0,f2,sound_speed):
    nl,m=p.shape;pp=p.copy();gg=g.copy();hh=h.copy()
    for l in range(nl):
        k=omega[l]/sound_speed;aa=-f0*k*k*radius**3/3;bb=-f2*radius**3/2
        for i in range(m):
            for j in range(m):
                if i==j:continue
                delta=q[i]-q[j];r=np.sqrt(np.sum(delta**2));nn=delta/r
                green=np.exp(1j*k*r)/r
                dg=green*(1j*k-1/r)*nn
                hgreen=np.empty((3,3),np.complex128);tg=np.empty((3,3,3),np.complex128)
                a2=-k*k-3j*k/r+3/r**2;b2=1j*k/r-1/r**2
                a3=-1j*k**3+6*k*k/r+15j*k/r**2-15/r**3;b3=-k*k/r-3j*k/r**2+3/r**3
                for a in range(3):
                    for b in range(3):
                        hgreen[a,b]=green*(a2*nn[a]*nn[b]+b2*(1 if a==b else 0))
                        for d in range(3):
                            tg[a,b,d]=green*(a3*nn[a]*nn[b]*nn[d]+b3*(nn[a]*(1 if b==d else 0)+
                                            nn[b]*(1 if a==d else 0)+nn[d]*(1 if a==b else 0)))
                pp[l,i]+=aa*green*p[l,j]
                for a in range(3):
                    pp[l,i]+=bb*dg[a]*g[l,j,a]
                    gg[l,i,a]+=aa*dg[a]*p[l,j]
                    for b in range(3):
                        gg[l,i,a]+=bb*hgreen[a,b]*g[l,j,b]
                        hh[l,i,a,b]+=aa*hgreen[a,b]*p[l,j]
                        for d in range(3):hh[l,i,a,b]+=bb*tg[a,b,d]*g[l,j,d]
    return pp,gg,hh


@njit(cache=False)
def radiation_force(q,coeff,origin,scale,exp,omega,ell,phase,radius,f0,f2,rho,c0,mass,mean=False):
    p,g,h=fields_at(q,coeff,origin,scale,exp)
    if len(q)>1:p,g,h=scattered_fields(q,p,g,h,omega,radius,f0,f2,c0)
    nl,m=p.shape;force=np.zeros((m,3))
    alpha=(4*np.pi*radius**3/3)*f0/(4*rho*c0*c0);beta=3*(4*np.pi*radius**3/3)*f2*rho/8
    for i in range(m):
        if mean:
            for l in range(nl):
                for a in range(3):
                    force[i,a]-=2*alpha*(np.conj(p[l,i])*g[l,i,a]).real
                    for b in range(3):force[i,a]+=2*beta*(np.conj(g[l,i,b]/(1j*omega[l]*rho))*h[l,i,b,a]/(1j*omega[l]*rho)).real
        else:
            pp=0j;gp=np.zeros(3,np.complex128);vv=np.zeros(3,np.complex128);dv=np.zeros((3,3),np.complex128)
            for l in range(nl):
                ph=np.exp(-1j*ell[l]*phase);pp+=p[l,i]*ph
                for a in range(3):
                    gp[a]+=g[l,i,a]*ph;vv[a]+=g[l,i,a]*ph/(1j*omega[l]*rho)
                    for b in range(3):dv[a,b]+=h[l,i,a,b]*ph/(1j*omega[l]*rho)
            for a in range(3):
                force[i,a]=-2*alpha*(np.conj(pp)*gp[a]).real
                for b in range(3):force[i,a]+=2*beta*(np.conj(vv[b])*dv[b,a]).real
        force[i,2]-=mass*9.80665
    return force


class Surrogate:
    def __init__(self,system,degree=6,scale=.0015,bins=8):
        self.system=system;self.c=system.c;self.origin=system.q.detach().cpu().numpy()
        self.degree=degree;self.scale=scale;self.exp=exponents(degree);self.bins=bins
        nl=len(system.ell);m=system.m;n=system.c.n;nc=len(self.exp)
        estimate=nl*m*n*nc*16/1024**2
        if estimate>system.c.memory_mb*.65:raise MemoryError(f"surrogate transfer requires {estimate:.0f} MiB")
        axis=np.cos(np.pi*(np.arange(degree+1)+.5)/(degree+1))
        grid=np.array([(x,y,z) for x in axis for y in axis for z in axis])
        v=np.prod(grid[:,None,:]**self.exp[None,:,:],axis=-1)
        inverse=np.linalg.pinv(v,rcond=1e-12)
        self.transfer=np.empty((nl,m,n,nc),np.complex128)
        for l,f in enumerate(system.frequencies):
            for i in range(m):
                q=torch.tensor(self.origin[i]+grid*scale,device=system.device)
                h=source_jets(q,system.positions,system.normals,system.radii,float(f),system.c,order=0)[0]
                self.transfer[l,i]=(inverse@h.detach().cpu().numpy()).T
        distances=np.linalg.norm(self.origin[:,None,:]-system.positions[None,:,:],axis=-1)
        low=distances.min();high=distances.max()+1e-10
        self.assignment=np.clip(((distances-low)/(high-low)*bins).astype(int),0,bins-1)
        self.delay_s=np.tile((low+(np.arange(bins)+.5)*(high-low)/bins)/self.c.sound_speed,(m,1))
        self.omega=system.omega.cpu().numpy();self.ell=system.ell.cpu().numpy()
        self.transfer_matrix=np.ascontiguousarray(self.transfer.transpose(0,1,3,2)).reshape(nl,m*nc,n)
        self.bin_matrix=None
        if estimate*(bins+2)<system.c.memory_mb*.60:
            self.bin_matrix=np.stack([self.transfer*(self.assignment==b)[None,:,:,None]
                                       for b in range(bins)],axis=1)
            self.bin_matrix=np.ascontiguousarray(self.bin_matrix.transpose(0,1,2,4,3)).reshape(nl,bins*m*nc,n)

    def polynomial(self,command):
        return np.stack([self.transfer_matrix[l]@command[l] for l in range(len(command))]).reshape(len(command),len(self.origin),len(self.exp))

    def binned(self,command):
        if self.bin_matrix is not None:
            result=np.stack([self.bin_matrix[l]@command[l] for l in range(len(command))])
            return result.reshape(len(command),self.bins,len(self.origin),len(self.exp)).transpose(0,2,1,3)
        return np.stack([np.einsum("lmnc,ln,mn->lmc",self.transfer,command,self.assignment==b,optimize=True)
                         for b in range(self.bins)],axis=2)

    def force(self,q,coeff,phase=0,mean=False):
        c=self.c
        return radiation_force(q,coeff,self.origin,self.scale,self.exp,self.omega,self.ell,phase,
                               c.particle_radius,c.f0,c.f2,c.rho,c.sound_speed,c.mass,mean)

    def stiffness(self,q,coeff):
        m=len(q);cols=[];step=self.c.fd_step
        for j in range(m):
            for a in range(3):
                delta=np.zeros_like(q);delta[j,a]=step
                cols.append((-(self.force(q+delta,coeff,mean=True)-self.force(q-delta,coeff,mean=True))/(2*step)).reshape(-1))
        return np.stack(cols,axis=1)


def exact_force(system,coeff,q,phase):
    """Independent reference at displaced sites using fresh finite-disk derivatives."""
    qq=torch.tensor(q,device=system.device);arrays=[[],[],[]]
    for l,f in enumerate(system.frequencies):
        hs=source_jets(qq,system.positions,system.normals,system.radii,float(f),system.c,order=2)
        uu=coeff[l]
        arrays[0].append(torch.einsum("mn,n->m",hs[0],uu))
        arrays[1].append(torch.einsum("mna,n->ma",hs[1],uu))
        arrays[2].append(torch.einsum("mnab,n->mab",hs[2],uu))
    p,g,h=system.scatter(qq,*[torch.stack(v)[None] for v in arrays])
    ph=torch.exp(-1j*system.ell*phase)
    pp=torch.einsum("l,zlm->zm",ph,p)
    gg=torch.einsum("l,zlmi->zmi",ph,g)
    vv=torch.einsum("l,zlmi->zmi",ph,g/(1j*system.omega[None,:,None,None]*system.c.rho))
    dv=torch.einsum("l,zlmij->zmij",ph,h/(1j*system.omega[None,:,None,None,None]*system.c.rho))
    c=system.c;alpha=c.volume*c.f0/(4*c.rho*c.sound_speed**2);beta=3*c.volume*c.f2*c.rho/8
    f=-2*alpha*(pp.conj()[...,None]*gg).real+2*beta*(vv.conj()[..., :,None]*dv).real.sum(-2)
    f[...,2]-=c.mass*9.80665
    return f[0].detach().cpu().numpy()


def validate_surrogate(surrogate,coeff,trials=6):
    rng=np.random.default_rng(surrogate.c.seed+100);errors=[]
    coeff_np=coeff.detach().cpu().numpy();polynomial=surrogate.polynomial(coeff_np)
    for _ in range(trials):
        q=surrogate.origin+rng.uniform(-.001,.001,surrogate.origin.shape)
        phase=rng.uniform(0,2*np.pi)
        predicted=surrogate.force(q,polynomial,phase)
        exact=exact_force(surrogate.system,coeff,q,phase)
        errors.append(float(np.linalg.norm(predicted-exact,axis=1).max()/(surrogate.c.mass*9.80665)))
    return {"force_error_fractions_of_weight":errors,"maximum_force_error_fraction":max(errors),
            "degree":surrogate.degree,"scale_m":surrogate.scale}


@njit(cache=False)
def delayed_coeff(history,cursor,offset,delay_steps):
    length,nl,m,bins,nc=history.shape;out=np.zeros((nl,m,nc),np.complex128)
    for i in range(m):
        for b in range(bins):
            location=cursor+offset-delay_steps[i,b]
            index=int(np.floor(location));fraction=location-index
            for l in range(nl):
                for j in range(nc):out[l,i,j]+=(1-fraction)*history[index%length,l,i,b,j]+fraction*history[(index+1)%length,l,i,b,j]
    return out


@njit(cache=False)
def advance_frame(q,v,envelope,history,cursor,command,source,source_command,delay_steps,
                  qhistory,qcursor,steps,dt,tstart,origin,scale,exp,omega,ell,
                  radius,f0,f2,rho,c0,mass,drag,viscosity,tau,cycle,target_start,target_end,
                  held,held_origin,wind,power_weights):
    nl=len(ell);lam=1/tau-1j*2*np.pi*ell/cycle
    response=1/(1-1j*2*np.pi*ell/cycle*tau)
    target=command*response[:,None,None,None]
    source_target=source_command*response[:,None]
    source_decay=np.exp(-lam*dt)
    maxerror=0.;maxheld=0.;energy=0.;peak=0.;failure=0;maxpower=0.
    for step in range(steps):
        if step*dt<30*tau:
            for l in range(nl):
                envelope[l]=target[l]+(envelope[l]-target[l])*source_decay[l]
                source[l]=source_target[l]+(source[l]-source_target[l])*source_decay[l]
        else:
            envelope[:]=target;source[:]=source_target
        cursor=(cursor+1)%len(history);history[cursor]=envelope
        a0=delayed_coeff(history,cursor,-1.,delay_steps)
        amid=delayed_coeff(history,cursor,-.5,delay_steps)
        a1=delayed_coeff(history,cursor,0.,delay_steps)
        phase=2*np.pi*(tstart+step*dt)/cycle
        def acceleration(qq,vv,co,ph):
            ff=radiation_force(qq,co,origin,scale,exp,omega,ell,ph,radius,f0,f2,rho,c0,mass,False)
            for i in range(len(qq)):
                rel=vv[i]-wind;re=2*radius*rho*np.sqrt(np.sum(rel**2))/viscosity
                ff[i]-=drag*(1+.15*re**.687)*rel
            return ff/mass
        k1q=v;k1v=acceleration(q,v,a0,phase)
        k2q=v+.5*dt*k1v;k2v=acceleration(q+.5*dt*k1q,k2q,amid,phase+np.pi*dt/cycle)
        k3q=v+.5*dt*k2v;k3v=acceleration(q+.5*dt*k2q,k3q,amid,phase+np.pi*dt/cycle)
        k4q=v+dt*k3v;k4v=acceleration(q+dt*k3q,k4q,a1,phase+2*np.pi*dt/cycle)
        q=q+dt*(k1q+2*k2q+2*k3q+k4q)/6
        v=v+dt*(k1v+2*k2v+2*k3v+k4v)/6
        qcursor=(qcursor+1)%len(qhistory);qhistory[qcursor]=q
        expected=target_start+(target_end-target_start)*(step+1)/steps
        for i in range(len(q)):
            error=np.sqrt(np.sum((q[i]-expected[i])**2));maxerror=max(maxerror,error)
            if held[i]:maxheld=max(maxheld,np.sqrt(np.sum((q[i]-held_origin[i])**2)))
            if np.max(np.abs(q[i]-origin[i]))>scale:failure=2
            for j in range(i):
                if np.sum((q[i]-q[j])**2)<(2*radius)**2:
                    # Contact is a failed independence trial. Inelastic impulse
                    # is still resolved, so trajectories do not pass through.
                    nn=q[i]-q[j];norm=np.sqrt(np.sum(nn**2));nn/=max(norm,1e-12)
                    relative=np.dot(v[i]-v[j],nn)
                    if relative<0:
                        impulse=-(1+.2)*relative/2;v[i]+=impulse*nn;v[j]-=impulse*nn
                    failure=1
        # Normalized source peak; complete periodic power is evaluated from the
        # exact radiated spectral coefficients in the caller.
        power=0.
        for n in range(source.shape[1]):
            amp=0j
            for l in range(nl):amp+=source[l,n]*np.exp(-1j*ell[l]*(phase+2*np.pi*dt/cycle))
            peak=max(peak,abs(amp));power+=power_weights[n]*abs(amp)**2
        maxpower=max(maxpower,power);energy+=power*dt
        if failure:break
    return q,v,envelope,history,cursor,source,qhistory,qcursor,maxerror,maxheld,peak,failure,energy,maxpower


def movement_protocol(origin,frame_hz=20):
    frames=[];m=len(origin);base=origin.copy()
    def segment(destination,seconds,moved,label):
        start=frames[-1][0] if frames else base
        count=max(1,round(seconds*frame_hz))
        for j in range(count):
            # Smooth cubic ramp, avoiding an unphysical instantaneous 1-mm jump.
            x=(j+1)/count;w=x*x*(3-2*x)
            frames.append((start+(destination-start)*w,np.array(moved,dtype=bool),label))
    for distance in (.00025,.001):
        for particle in range(m):
            for axis in range(3):
                for sign in (-1,1):
                    target=base.copy();target[particle,axis]+=sign*distance
                    moved=np.zeros(m,bool);moved[particle]=True
                    ramp=.4 if distance<.001 else 2.
                    label=f"{distance*1000:g}mm_p{particle}_a{axis}_{sign}"
                    segment(target,ramp,moved,label);segment(target,.2,moved,label)
                    segment(base,ramp,moved,label);segment(base,.2,moved,label)
    # Distinct particles move simultaneously in independent directions.
    for sign in (-1,1):
        target=base.copy()
        for i in range(m):target[i,i%3]+=sign*.001
        segment(target,2.,np.ones(m,bool),"simultaneous");segment(target,.2,np.ones(m,bool),"simultaneous")
        segment(base,2.,np.ones(m,bool),"simultaneous");segment(base,.2,np.ones(m,bool),"simultaneous")
    while len(frames)<60*frame_hz:frames.append((base.copy(),np.zeros(m,bool),"retention"))
    return frames
