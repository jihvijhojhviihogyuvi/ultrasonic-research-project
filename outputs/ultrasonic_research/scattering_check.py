"""Leading Rayleigh/Born particle-scattering sensitivity and coupled dynamics.

Pressure scattering is converted directly from Silva & Bruus (2014), Eq. 3:
p_sc = -f0*k^2*a^3*p_source*G/3 - f2*a^3*grad(p_source).grad(G)/2.
Only cross terms with the incident field are retained in the interaction
potential, consistent with leading scattering order. Self scattering is excluded.
This is not a full finite-size T-matrix or thermoviscous validation.
"""
import argparse
import json
import math
from pathlib import Path
import torch
from simulator import Config, Field, geometry, envelope


def particle_kernel(delta,k,exclude):
    distance=torch.linalg.vector_norm(delta,dim=-1)
    distance=torch.where(exclude,torch.ones_like(distance),distance)
    direction=delta/distance[...,None]
    g=torch.exp(1j*k*distance)/distance
    dg=g[...,None]*(1j*k-1/distance)[...,None]*direction
    hh=g[...,None,None]*(
        (-k*k-3j*k/distance+3/distance**2)[...,None,None]*direction[..., :,None]*direction[...,None,:]
        +(1j*k/distance-1/distance**2)[...,None,None]*torch.eye(3))
    g=torch.where(exclude,torch.zeros_like(g),g)
    dg=torch.where(exclude[...,None],torch.zeros_like(dg),dg)
    hh=torch.where(exclude[...,None,None],torch.zeros_like(hh),hh)
    return g,dg,hh


class ScatteringField(Field):
    def potentials(self,drives):
        primary=super().potentials(drives)
        c=self.c;m=len(self.x)
        if m<2:return primary
        if float(torch.pdist(self.x).min())<=2*c.particle_radius+2*c.fd_step:
            raise ValueError("particles/contact stencil too close for Rayleigh/Born screening")
        volume=4*math.pi*c.particle_radius**3/3
        f0=1-c.compressibility_ratio
        f2=2*(c.particle_density-c.air_density)/(2*c.particle_density+c.air_density)
        delta=self.q[:,:,None,:]-self.x[None,None,:,:]
        exclude=torch.eye(m,dtype=torch.bool)[:,None,:].expand(m,self.q.shape[1],m)
        correction=torch.zeros_like(primary)
        for fi,u,weights in drives:
            omega=2*math.pi*self.freqs[fi];k=omega/c.sound_speed
            g,dg,hh=particle_kernel(delta,k,exclude)
            p=torch.einsum("mjn,sn->smj",self.hs[fi],u)
            dp=torch.einsum("mjnl,sn->smjl",self.ds[fi],u)
            pc=p[:,:,0];dc=dp[:,:,0,:]
            monopole=-f0*k*k*c.particle_radius**3/3
            dipole=-f2*c.particle_radius**3/2
            ps=monopole*torch.einsum("mjt,st->smj",g,pc)+dipole*torch.einsum("mjtl,stl->smj",dg,dc)
            ds=monopole*torch.einsum("mjtl,st->smjl",dg,pc)+dipole*torch.einsum("mjtlb,stb->smjl",hh,dc)
            up=volume*(f0*(p.conj()*ps).real/(2*c.air_density*c.sound_speed**2)
                -3*f2*c.air_density*(dp.conj()*ds).real.sum(-1)/(4*(omega*c.air_density)**2))
            correction+=torch.einsum("s,smj->mj",weights,up)
        return primary+correction


def coupled_force(q,pos,drives,freqs,c):
    return ScatteringField(q,pos,freqs,c).force_stiffness(drives)[0].detach()


def coupled_stiffness(q,pos,drives,freqs,c):
    # Moving a physical source particle changes both its location and excitation.
    # Rebuild the full force, not just an isolated-particle Hessian block.
    cols=[];h=c.fd_step/2
    for j in range(len(q)):
        for axis in range(3):
            shift=torch.zeros_like(q);shift[j,axis]=h
            plus=coupled_force(q+shift,pos,drives,freqs,c)
            minus=coupled_force(q-shift,pos,drives,freqs,c)
            cols.append((-(plus-minus)/(2*h)).reshape(-1))
    return torch.stack(cols,dim=1)


def selftest():
    c=Config(n=32,gap=.060,particle_radius=.000025,fd_step=.000001)
    pos,_=geometry(c);r=.00030
    q=torch.tensor([[-r/2,0.,0.],[r/2,0.,0.]])
    f=ScatteringField(q,pos,[c.frequency],c)
    k=2*math.pi*c.frequency/c.sound_speed;p0=1000.
    f.hs[0]=torch.zeros_like(f.hs[0]);f.ds[0]=torch.zeros_like(f.ds[0])
    f.hs[0][:,:,0]=p0*torch.sin(k*f.q[:,:,2]).to(torch.complex128)
    f.ds[0][:,:,0,2]=p0*k*torch.cos(k*f.q[:,:,2]).to(torch.complex128)
    u=torch.zeros((1,c.n),dtype=torch.complex128);u[0,0]=1
    force,_=f.force_stiffness([(0,u,torch.ones(1))])
    f2=2*(c.particle_density-c.air_density)/(2*c.particle_density+c.air_density)
    energy=.5*c.air_density*(p0/(c.air_density*c.sound_speed))**2
    expected=3*math.pi*energy*f2*f2*c.particle_radius**6/r**4
    ratio=float(force[0,0])/expected
    if not .95<ratio<1.05 or float(force[1,0])>=0:
        raise AssertionError(f"near-node r^-4 attraction control: ratio {ratio}")
    print(f"SCATTERING_CONTROL_PASSED: near-node attraction / asymptotic = {ratio:.6f}")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results",type=Path,default=Path("results_128"))
    p.add_argument("--architectures",default="static_force,covariance")
    p.add_argument("--out",type=Path,default=Path("scattering_results.json"))
    p.add_argument("--selftest",action="store_true")
    a=p.parse_args()
    if a.selftest:selftest();return
    c=Config(**json.loads((a.results/"metadata.json").read_text())["config"])
    details=json.loads((a.results/"fields.json").read_text())
    rows=[]
    for row in json.loads((a.results/"results.json").read_text()):
        arch=row["architecture"];m=row["requested_particles"]
        if arch not in a.architectures.split(",") or not row["converged"]:continue
        detail=details[f"{arch}_{m}"]
        pos,_=geometry(c,arch=="aperiodic")
        q=torch.tensor(detail["equilibria_m"]);q0=q.clone()
        u=torch.tensor(detail["command_real"])+1j*torch.tensor(detail["command_imag"])
        if arch=="multifrequency":
            freqs=[c.frequency-400,c.frequency+400]
            drives=[(j,u[j:j+1],torch.ones(1)) for j in range(2)]
        else:
            freqs=[c.frequency];uu,w=envelope(u,c);drives=[(0,uu,w)]
        f=coupled_force(q,pos,drives,freqs,c)
        d=coupled_stiffness(q,pos,drives,freqs,c)
        ok=False
        for _ in range(8):
            if float(torch.linalg.vector_norm(f,dim=1).max())<c.mass*9.80665*1e-4:
                ok=True;break
            step=torch.linalg.solve(d,f.reshape(-1)).reshape(-1,3)
            q+=step
            f=coupled_force(q,pos,drives,freqs,c)
            d=coupled_stiffness(q,pos,drives,freqs,c)
        n=3*m;aa=torch.zeros((2*n,2*n))
        aa[:n,n:]=torch.eye(n);aa[n:,:n]=-d/c.mass
        aa[n:,n:]=-6*math.pi*c.viscosity*c.particle_radius/c.mass*torch.eye(n)
        off=d.clone()
        for j in range(m):off[3*j:3*j+3,3*j:3*j+3]=0
        rows.append({"architecture":arch,"particles":m,"converged":ok,
            "max_scattering_equilibrium_shift_um":float(torch.linalg.vector_norm(q-q0,dim=1).max())*1e6,
            "relative_off_diagonal_stiffness_norm":float(torch.linalg.vector_norm(off)/torch.linalg.vector_norm(d)),
            "symmetric_stiffness_min_N_m":float(torch.linalg.eigvalsh((d+d.T)/2).min()),
            "mean_field_dynamic_max_real_eigenvalue_s_inv":float(torch.linalg.eigvals(aa).real.max()),
            "warning":"mean-field sensitivity only; temporal Floquet coupling and finite-size corrections remain"})
        print(json.dumps(rows[-1]),flush=True)
        a.out.write_text(json.dumps(rows,indent=2))


if __name__=="__main__":main()
