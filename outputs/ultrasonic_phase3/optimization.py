"""Holographic, force, timing, and direct periodic-stability optimizers."""
from __future__ import annotations
import math,time
from dataclasses import replace
import numpy as np
import torch
from physics import AcousticSystem,Config,linear_floquet,affine_periodic_orbit


def ib(system,points,iterations=40):
    c=system.c;dev=system.device
    # Recompute the pressure transfer at each control point, not particle force.
    from physics import source_jets
    h=source_jets(torch.tensor(points,device=dev),system.positions,system.normals,
                  system.radii,c.frequency,c,order=0)[0]
    u=torch.ones(c.n,device=dev,dtype=torch.complex128)
    for _ in range(iterations):
        p=h@u;back=h.conj().T@(p/p.abs().clamp_min(1e-30))
        u=back/back.abs().clamp_min(1e-30)
    u=torch.where(torch.tensor(system.positions[:,2]>0,device=dev),-u,u)
    return u*.65


class Design:
    def __init__(self,system,architecture="static",fields=1,seed=None,initial=None):
        self.system=system;self.architecture=architecture;self.s=fields
        c=system.c;dev=system.device
        rng=np.random.default_rng(c.seed if seed is None else seed)
        if initial is None:
            modes=[]
            for j in range(fields):
                subset=system.q.detach().cpu().numpy()[j%min(fields,system.m)::min(fields,system.m)]
                if architecture in ("static","ib","support","common"):
                    subset=system.q.detach().cpu().numpy()
                u=ib(system,subset)
                if fields>1:u*=np.exp(1j*rng.uniform(-math.pi,math.pi))
                modes.append(u)
            initial=torch.stack(modes)
        self.raw=torch.nn.Parameter(torch.view_as_real(initial).clone())
        self.timing=torch.nn.Parameter(torch.zeros(fields,device=dev))
        self.gauge=torch.nn.Parameter(torch.zeros(fields,device=dev))
        self.support=torch.where(torch.tensor(system.positions[:,2]<0,device=dev),.4,-.4).to(torch.complex128)
        self.common=torch.nn.Parameter(torch.view_as_real(ib(system,system.q.detach().cpu().numpy())).clone())
        self.residual_scale=torch.nn.Parameter(torch.tensor(-.5,device=dev))
        self.history=[]

    def duties(self):
        return .02/self.s+.98*torch.softmax(self.timing,0)

    def fields(self):
        raw=torch.view_as_complex(self.raw.contiguous())*torch.exp(1j*self.gauge)[:,None]
        if self.architecture in ("ib","singleton"):raw=raw/raw.abs().clamp_min(1e-30)
        elif self.architecture=="support":raw=self.support[None,:]+.45*raw/raw.abs().clamp_min(1)
        elif self.architecture=="common":
            common=torch.view_as_complex(self.common.contiguous())
            raw=common[None,:]+torch.sigmoid(self.residual_scale)*raw
        if self.architecture=="multifrequency":
            raw=raw/raw.abs().sum(0).clamp_min(1)[None,:]
            coeff=self.command_coefficients(raw)*self.system.response()[:,None]
            scale=torch.sqrt(self.system.c.power_W/self.system.power(coeff).clamp_min(1e-30)).clamp_max(1)
            return raw*scale
        return self.system.project(raw,self.duties())[0]

    def command_coefficients(self,fields=None):
        fields=self.fields() if fields is None else fields
        if self.architecture=="multifrequency":
            if self.s!=2:raise ValueError("two-tone control requires two drive rows")
            result=torch.zeros((len(self.system.ell),self.system.c.n),device=self.system.device,dtype=torch.complex128)
            for j,tone in enumerate((-2,2)):
                index=torch.where(self.system.ell==tone)[0]
                if len(index)!=1:raise ValueError("sideband support must contain the two frequencies")
                result[index[0]]=fields[j]
            return result
        return self.system.command_coefficients(fields,self.duties())

    def coefficients(self):return self.command_coefficients()*self.system.response()[:,None]

    def state(self):
        return {"fields_real":self.fields().detach().cpu().real.tolist(),
                "fields_imag":self.fields().detach().cpu().imag.tolist(),
                "duties":self.duties().detach().cpu().tolist(),
                "architecture":self.architecture,"field_count":self.s,"history":self.history}


def mean_loss(system,coeff):
    f,k=system.primary_mean(coeff)
    eig=torch.linalg.eigvalsh(k)
    # Lexicographic-style penalties; no numerical result is assumed feasible.
    force=((f/(system.c.mass*9.80665))**2).mean()
    shortage=torch.nn.functional.softplus((system.c.stiffness_required-eig[:,0])/system.c.stiffness_required)
    return 20*force+5*shortage.square().mean()


def optimize(design,iterations=100,stage="mean",jacobian_weight=0):
    s=design.system;c=s.c
    if design.architecture in ("ib","singleton"):return 0.
    params=[design.raw,design.timing,design.gauge]
    if design.architecture=="common":params += [design.common,design.residual_scale]
    if stage=="timing":params=[design.timing,design.gauge]
    opt=torch.optim.Adam(params,lr=.025 if stage=="mean" else .012)
    best=float("inf");saved=None;start=time.perf_counter()
    for iteration in range(iterations):
        opt.zero_grad();coeff=design.coefficients();loss=mean_loss(s,coeff)
        rho=None
        if stage in ("floquet","timing"):
            force,stiff,_=s.matrices(coeff)
            rho,_=linear_floquet(stiff,c)
            dyn=torch.nn.functional.softplus((rho-.97)*30)/30
            loss=loss+150*dyn.square()
            mean_force=force.mean(0)/(c.mass*9.80665)
            loss+=5*mean_force.square().mean()
        if jacobian_weight and s.m<=4:
            # Normalize force-control sensitivity by weight and source drive units.
            def ff(raw):
                old=design.raw;design.raw=raw
                value=s.primary_mean(design.coefficients())[0].reshape(-1)/(c.mass*9.80665)
                design.raw=old;return value
            b=torch.autograd.functional.jacobian(ff,design.raw,create_graph=True,vectorize=True).reshape(3*s.m,-1)
            gram=b@b.T/(2*c.n*design.s)
            sigma2=torch.linalg.eigvalsh(gram)[0]
            loss+=jacobian_weight*torch.nn.functional.softplus(1-sigma2)
        value=float(loss.detach())
        if np.isfinite(value) and value<best:
            best=value;saved=[p.detach().clone() for p in params]
        if not torch.isfinite(loss):break
        loss.backward()
        norm=torch.nn.utils.clip_grad_norm_(params,10)
        if not torch.isfinite(norm):break
        opt.step()
        if iteration%10==0 or iteration==iterations-1:
            design.history.append({"stage":stage,"iteration":iteration,"loss":value,
                                   "rho":None if rho is None else float(rho.detach())})
    if saved is not None:
        for p,v in zip(params,saved):p.data.copy_(v)
    return time.perf_counter()-start


def nonlinear_cycle(system,coeff,y,steps=64):
    """Differentiable coupled nonlinear particle trajectory over one full cycle.

    The external field uses its local third-order jets; accepted micromotion must
    be small and is checked against the higher-order surrogate in validation.
    """
    c=system.c;ref=2*math.pi*15;scale=.0001;jets=system.jets(coeff)
    dt=c.cycle_s/steps
    def rhs(state,phase):
        q=system.q+state[:3*system.m].reshape(system.m,3)*scale
        v=state[3*system.m:].reshape(system.m,3)*scale*ref
        force=system.forces(q,jets,torch.tensor([phase],device=system.device))[0,0]
        # The second derivative of ||v|| is undefined at v=0. A negligible
        # velocity regularizer keeps the flow-Jacobian gradient finite.
        speed=torch.sqrt(v.square().sum(-1)+1e-16)
        reynolds=2*c.particle_radius*c.rho*speed/c.viscosity
        drag=c.drag*(1+.15*reynolds.clamp_min(1e-12)**.687)
        acc=(force-drag[:,None]*v)/c.mass
        return torch.cat((v.reshape(-1)/scale,acc.reshape(-1)/(scale*ref)))
    state=y
    for j in range(steps):
        phase=2*math.pi*j/steps
        k1=rhs(state,phase);k2=rhs(state+.5*dt*k1,phase+math.pi/steps)
        k3=rhs(state+.5*dt*k2,phase+math.pi/steps);k4=rhs(state+dt*k3,phase+2*math.pi/steps)
        state=state+dt*(k1+2*k2+2*k3+k4)/6
    return state


def optimize_nonlinear_floquet(design,iterations=8,steps=32):
    """Joint shooting closure and nonlinear flow-Jacobian/Floquet optimization."""
    s=design.system
    if s.m>4:raise ValueError("nonlinear gradient shooting budget is restricted to <=4 particles")
    y=torch.nn.Parameter(torch.zeros(6*s.m,device=s.device))
    params=[design.raw,design.timing,design.gauge,y]
    opt=torch.optim.Adam(params,lr=.007)
    start=time.perf_counter();saved=None;best=float("inf")
    for iteration in range(iterations):
        opt.zero_grad();coeff=design.coefficients()
        end=nonlinear_cycle(s,coeff,y,steps)
        phi=torch.autograd.functional.jacobian(lambda yy:nonlinear_cycle(s,coeff,yy,steps),
                    y,create_graph=True,vectorize=True)
        rho=torch.linalg.eigvals(phi).abs().max()
        closure=(end-y).square().mean()
        loss=mean_loss(s,coeff)+50*closure+100*torch.nn.functional.softplus((rho-.97)*20).square()/400
        val=float(loss.detach())
        if val<best:best=val;saved=[p.detach().clone() for p in params]
        loss.backward();torch.nn.utils.clip_grad_norm_(params,5);opt.step()
        design.history.append({"stage":"nonlinear_floquet","iteration":iteration,
                    "loss":val,"rho":float(rho.detach()),"closure":float(closure.detach())})
    if saved:
        for p,v in zip(params,saved):p.data.copy_(v)
    return time.perf_counter()-start


def optimize_order(design,trials=8):
    """Discrete ordering search using actual filtered fields and periodic dynamics."""
    rng=np.random.default_rng(design.system.c.seed)
    original=design.raw.detach().clone();timing=design.timing.detach().clone()
    gauge=design.gauge.detach().clone()
    best=float("inf");best_state=None;rows=[]
    orders=[np.arange(design.s)] + [rng.permutation(design.s) for _ in range(trials-1)]
    for order in orders:
        design.raw.data.copy_(original[order]);design.timing.data.copy_(timing[order])
        design.gauge.data.copy_(gauge[order])
        with torch.no_grad():
            coeff=design.coefficients();force,stiff,_=design.system.matrices(coeff)
            rho,_=linear_floquet(stiff,design.system.c)
            score=float(mean_loss(design.system,coeff))+100*max(0,float(rho)-.97)**2
        rows.append({"order":order.tolist(),"rho":float(rho),"score":score})
        if score<best:best=score;best_state=(design.raw.detach().clone(),design.timing.detach().clone(),design.gauge.detach().clone())
    design.raw.data.copy_(best_state[0]);design.timing.data.copy_(best_state[1])
    design.gauge.data.copy_(best_state[2])
    return rows


def equilibrium(system,coeff,iterations=12):
    """Gravity-balanced coupled equilibrium; rebuild source jets if position shifts."""
    q=system.q.detach().clone();current=system;ok=False
    with torch.no_grad():
        for _ in range(iterations):
            f,d,_=current.matrices(coeff)
            f=f.mean(0);d=d.mean(0)
            if float(torch.linalg.vector_norm(f,dim=1).max())<system.c.mass*9.80665*1e-4:
                ok=True;break
            try:step=torch.linalg.solve(d,f.reshape(-1)).reshape(-1,3)
            except RuntimeError:break
            norm=torch.linalg.vector_norm(step,dim=-1).clamp_min(1e-30)
            q=q+step*(.0002/norm).clamp_max(1)[:,None]
            if float(torch.linalg.vector_norm(q-system.q,dim=1).max())>.003:
                q=current.q.detach().clone();break
            current=AcousticSystem(system.c,q.cpu().numpy(),system.positions,system.normals,system.radii,
                                   system.cycle_s,int(system.ell.abs().max()))
        f,d,_=current.matrices(coeff)
    return q,f.mean(0),d.mean(0),ok,current


def sdp_bound(system):
    """Executed single-carrier PSD bound at tractable N; no dynamic certificate."""
    import cvxpy as cp
    if system.c.n>64:raise ValueError("SDP experiment is limited to 64 sources")
    index=int(torch.where(system.ell==0)[0][0]);h,g,hh,t=[x[index].cpu().numpy() for x in system.h]
    c=system.c;vol=c.volume;alpha=vol*c.f0/(4*c.rho*c.sound_speed**2);beta=3*vol*c.f2*c.rho/8
    v=g/(1j*2*np.pi*c.frequency*c.rho);dv=hh/(1j*2*np.pi*c.frequency*c.rho)
    ddv=t/(1j*2*np.pi*c.frequency*c.rho)
    n=c.n;Q=cp.Variable((n,n),hermitian=True)
    constraints=[Q>>0,cp.real(cp.diag(Q))<=1]
    weights=system.power_coeff[index].cpu().numpy()
    power=cp.real(cp.sum(cp.multiply(weights,cp.diag(Q))))
    constraints.append(power<=c.power_W)
    for m in range(system.m):
        kexpr=[]
        for a in range(3):
            force_matrix=-alpha*(np.outer(h[m].conj(),g[m,:,a])+np.outer(g[m,:,a].conj(),h[m]))
            for l in range(3):
                force_matrix+=beta*(np.outer(v[m,:,l].conj(),dv[m,:,l,a])+np.outer(dv[m,:,l,a].conj(),v[m,:,l]))
            constraints.append(cp.real(cp.trace(force_matrix@Q))/(c.mass*9.80665)==float(a==2))
            row=[]
            for b in range(3):
                mat=alpha*(np.outer(g[m,:,a].conj(),g[m,:,b])+np.outer(g[m,:,b].conj(),g[m,:,a])+
                           np.outer(h[m].conj(),hh[m,:,a,b])+np.outer(hh[m,:,a,b].conj(),h[m]))
                for l in range(3):
                    mat-=beta*(np.outer(dv[m,:,l,a].conj(),dv[m,:,l,b])+np.outer(dv[m,:,l,b].conj(),dv[m,:,l,a])+
                               np.outer(v[m,:,l].conj(),ddv[m,:,l,a,b])+np.outer(ddv[m,:,l,a,b].conj(),v[m,:,l]))
                row.append(cp.real(cp.trace(mat@Q))/c.stiffness_required)
            kexpr.append(row)
        constraints.append(cp.bmat(kexpr)>>np.eye(3))
    problem=cp.Problem(cp.Minimize(power),constraints)
    # Rescale tiny force equalities above through solver tolerance; SCS is an
    # optimistic numerical bound whose residuals are checked explicitly below.
    start=time.perf_counter();problem.solve(solver="CLARABEL",max_iter=200)
    return {"status":problem.status,"power_W":None if Q.value is None else float(power.value),
            "solver_s":time.perf_counter()-start,
            "maximum_constraint_violation":None if Q.value is None else
                max(float(np.max(np.abs(con.violation()))) for con in constraints),
            "covariance_real":None if Q.value is None else Q.value.real.tolist(),
            "covariance_imag":None if Q.value is None else Q.value.imag.tolist(),
            "warning":"fixed-target ideal-field bound; scattering/peak realizability/periodic dynamics excluded"}
