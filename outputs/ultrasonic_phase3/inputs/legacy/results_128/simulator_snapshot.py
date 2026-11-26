"""Reproducible Tier-1 acoustic architecture screening; NOT a validated capacity model.

SI units. Peak (not RMS) phasors, exp(-i omega t). Baffled-piston disk
quadrature in homogeneous air; no enclosure, streaming, or particle scattering.
Force and stiffness use central differences of the correct Gor'kov potential.
Optimization is differentiable in complex emitter drives through PyTorch.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import torch

torch.set_default_dtype(torch.float64)
torch.set_num_threads(1)
CDTYPE = torch.complex128


@dataclass
class Config:
    n: int = 128
    pitch: float = 0.010
    gap: float = 0.100
    volume_side: float = 0.030
    frequency: float = 40000.0
    air_density: float = 1.20
    sound_speed: float = 343.0
    viscosity: float = 1.81e-5
    particle_radius: float = 0.00025
    particle_density: float = 30.0
    # Effective compressibility ratio kappa_particle / kappa_air (not c_p guessed).
    compressibility_ratio: float = 0.01
    emitter_radius: float = 0.0045
    source_gain: float = 2.0   # Pa m at f0, an ASSUMED source calibration
    power_fraction: float = 0.50
    mechanical_hz_required: float = 15.0
    fd_step: float = 0.00005
    cycle_s: float = 0.005
    tau_s: float = 0.00016
    envelope_samples: int = 6
    tone_q: float = 20.0
    probe_distance: float = 0.00025
    accuracy: float = 0.00010
    crosstalk_fraction: float = 0.10
    iterations: int = 250
    probe_iterations: int = 100
    seed: int = 7

    @property
    def mass(self):
        return 4 * math.pi * self.particle_radius**3 * self.particle_density / 3

    @property
    def k_required(self):
        return self.mass * (2 * math.pi * self.mechanical_hz_required)**2


ARCHITECTURES = ("ib", "static_force", "multifrequency", "global_local",
                 "time_single", "covariance", "support_residual", "aperiodic")


def geometry(c: Config, aperiodic=False):
    side = int(round(math.sqrt(c.n / 2)))
    if 2 * side * side != c.n:
        raise ValueError("n must be twice a square: 32, 128, 512, 2048 ...")
    xy = (np.arange(side) - (side - 1) / 2) * c.pitch
    points = np.array([(x, y, z) for z in (-c.gap / 2, c.gap / 2)
                       for x in xy for y in xy], dtype=float)
    if aperiodic:
        # Same N, aperture bounds, source area, and enclosing box. This is a
        # fixed jittered candidate, not a geometry optimizer or validated design.
        rng = np.random.default_rng(c.seed)
        jitter = rng.uniform(-0.03, 0.03, (c.n, 2)) * c.pitch
        points[:, :2] = np.clip(points[:, :2] + jitter, xy.min(), xy.max())
    normals = np.zeros_like(points)
    normals[:, 2] = -np.sign(points[:, 2])
    return torch.tensor(points), torch.tensor(normals)


def targets(c: Config, m: int):
    """Same seeded, genuinely volumetric target cloud for all architectures."""
    rng = np.random.default_rng(c.seed + m)
    pts = []
    min_sep = min(0.009, c.volume_side / max(2, m**(1/3)))
    for _ in range(100000):
        q = rng.uniform(-0.5, 0.5, 3) * c.volume_side
        if all(np.linalg.norm(q - x) >= min_sep for x in pts):
            pts.append(q)
            if len(pts) == m:
                return torch.tensor(np.array(pts))
    raise ValueError("cannot pack requested cloud at stated minimum separation")


def disk_sources(pos, radius):
    # Polar midpoint quadrature: 3 equal-area annuli, 8 azimuths = 24 patches.
    # All supplied geometries have normals parallel to +/- z.
    offsets = []
    for ring in range(3):
        r = radius * math.sqrt((ring + 0.5) / 3)
        for j in range(8):
            angle = 2 * math.pi * (j + 0.5 * (ring % 2)) / 8
            offsets.append([r * math.cos(angle), r * math.sin(angle), 0])
    return pos[:, None, :] + torch.tensor(offsets)[None, :, :]


def propagator(x, pos, frequency, c: Config):
    """Disk-integrated e^(ikr)/r and its EXACT spatial derivative for quadrature."""
    src = disk_sources(pos, c.emitter_radius)
    delta = x[:, None, None, :] - src[None, :, :, :]
    distance = torch.linalg.vector_norm(delta, dim=-1)
    if float(distance.min()) < c.emitter_radius:
        raise ValueError("evaluation too near an emitter for this screening model")
    k = 2 * math.pi * frequency / c.sound_speed
    detuning = 2 * c.tone_q * (frequency - c.frequency) / c.frequency
    response = (frequency / c.frequency) / complex(1, detuning)
    hpatch = c.source_gain * response * torch.exp(1j * k * distance) / distance
    dpatch = hpatch[..., None] * (1j*k - 1/distance)[..., None] * delta / distance[..., None]
    return hpatch.mean(dim=2), dpatch.mean(dim=2)


def stencil(h):
    pts = [np.zeros(3)]
    axial = {}
    for axis in range(3):
        axial[axis] = []
        for sign in (1, -1):
            q = np.zeros(3); q[axis] = sign*h
            axial[axis].append(len(pts)); pts.append(q)
    mixed = {}
    for a in range(3):
        for b in range(a + 1, 3):
            mixed[a, b] = []
            for sa, sb in ((1,1),(1,-1),(-1,1),(-1,-1)):
                q = np.zeros(3); q[a] = sa*h; q[b] = sb*h
                mixed[a, b].append(len(pts)); pts.append(q)
    return torch.tensor(np.array(pts)), axial, mixed


class Field:
    def __init__(self, x, pos, freqs, c):
        self.c = c
        self.x = x
        offsets, self.axial, self.mixed = stencil(c.fd_step)
        self.q = x[:, None, :] + offsets[None, :, :]
        flat = self.q.reshape(-1, 3)
        self.hs, self.ds = [], []
        for f in freqs:
            h, d = propagator(flat, pos, f, c)
            self.hs.append(h.reshape(len(x), len(offsets), c.n))
            self.ds.append(d.reshape(len(x), len(offsets), c.n, 3))
        self.freqs = freqs

    def potentials(self, drives):
        """drives = [(frequency index, SxN complex phasors, S weights)]."""
        c = self.c
        vol = 4 * math.pi * c.particle_radius**3 / 3
        f1 = 1 - c.compressibility_ratio
        f2 = 2*(c.particle_density - c.air_density)/(2*c.particle_density + c.air_density)
        result = None
        for fi, u, weights in drives:
            p = torch.einsum("mjn,sn->smj", self.hs[fi], u)
            grad = torch.einsum("mjnl,sn->smjl", self.ds[fi], u)
            velocity2 = grad.abs().square().sum(-1)/(2*math.pi*self.freqs[fi]*c.air_density)**2
            up = vol * (f1*p.abs().square()/(4*c.air_density*c.sound_speed**2)
                        - 3*f2*c.air_density*velocity2/8)
            term = torch.einsum("s,smj->mj", weights, up)
            result = term if result is None else result + term
        return result

    def force_stiffness(self, drives):
        u = self.potentials(drives)
        h = self.c.fd_step
        force = torch.stack([-(u[:, self.axial[a][0]] - u[:, self.axial[a][1]])/(2*h)
                             for a in range(3)], dim=-1)
        entries = {}
        for a in range(3):
            plus, minus = self.axial[a]
            entries[a,a] = (u[:,plus] - 2*u[:,0] + u[:,minus])/h**2
        for (a,b), idx in self.mixed.items():
            entries[a,b] = entries[b,a] = (u[:,idx[0]] - u[:,idx[1]]
                                         - u[:,idx[2]] + u[:,idx[3]])/(4*h*h)
        k = torch.stack([torch.stack([entries[a,b] for b in range(3)], -1)
                         for a in range(3)], -2)
        force = force.clone()
        force[:,2] -= self.c.mass*9.80665
        return force, k


def ib_drives(x, pos, c, iterations=40):
    h, _ = propagator(x, pos, c.frequency, c)
    u = torch.ones(c.n, dtype=CDTYPE)
    for _ in range(iterations):
        p = h @ u
        desired = p / p.abs().clamp_min(1e-30)
        back = h.conj().T @ desired
        u = back / back.abs().clamp_min(1e-30)
    u[pos[:,2] > 0] *= -1
    return u * math.sqrt(c.power_fraction)


def envelope(u, c):
    """Periodic first-order complex envelope; midpoint quadrature of each slot.

    Exposes finite settling during PHASE changes, rather than treating hardware
    updates as instantaneous acoustic updates. All complex coefficients share
    the same explicitly assumed tau; measured impulse responses must replace it.
    """
    s = len(u)
    if s == 1 or c.tau_s <= 0:
        return u, torch.full((s,), 1/s)
    dt = c.cycle_s/s
    a = math.exp(-dt/c.tau_s)
    start = sum((1-a)*a**(s-1-j)*u[j] for j in range(s))/(1-a**s)
    phases = []
    for j in range(s):
        for q in range(c.envelope_samples):
            t = dt*(q+0.5)/c.envelope_samples
            phases.append(u[j] + (start-u[j])*math.exp(-t/c.tau_s))
        start = a*start + (1-a)*u[j]
    return torch.stack(phases), torch.full((s*c.envelope_samples,), 1/(s*c.envelope_samples))


class Design:
    def __init__(self, arch, x, pos, c):
        self.arch, self.c, self.pos = arch, c, pos
        self.freqs = [c.frequency-400, c.frequency+400] if arch == "multifrequency" else [c.frequency]
        m = len(x)
        self.support = torch.zeros(c.n, dtype=CDTYPE)
        sign = torch.where(pos[:,2] < 0, 1., -1.)
        if arch in ("global_local", "support_residual"):
            self.support = sign.to(CDTYPE)*math.sqrt(c.power_fraction)*0.75
        if arch in ("time_single", "global_local"):
            modes = torch.stack([ib_drives(x[j:j+1], pos, c) for j in range(m)])
            if arch == "global_local":
                modes = self.support + 0.45*modes
        elif arch == "covariance":
            groups = min(4,m)
            modes = torch.stack([ib_drives(x[j::groups],pos,c) for j in range(groups)])
        elif arch == "multifrequency":
            modes = torch.stack([ib_drives(x[j::2] if len(x[j::2]) else x, pos, c)/math.sqrt(2)
                                 for j in range(2)])
        else:
            modes = ib_drives(x, pos, c)[None,:]
            if arch == "support_residual":
                modes = (modes-self.support)/0.5
        self.raw = torch.nn.Parameter(torch.view_as_real(modes).clone())
        self.initial = self.raw.detach().clone()

    def commanded(self):
        u = torch.view_as_complex(self.raw.contiguous())
        if self.arch == "support_residual":
            # Bounded residual around a coherent reference; all interference
            # terms are present because phasors are added BEFORE squaring.
            r = u/u.abs().clamp_min(1)
            u = self.support[None,:] + 0.5*r
        if self.arch == "multifrequency":
            # Common peak voltage/crest bound for simultaneous tones.
            crest = u.abs().sum(0).clamp_min(1)
            u = u/crest[None,:]
            power = u.abs().square().sum()
        else:
            u = u/u.abs().clamp_min(1)
            power = u.abs().square().sum()/len(u)
        factor = torch.sqrt((self.c.n*self.c.power_fraction)/power.clamp_min(1e-30)).clamp_max(1)
        return u*factor

    def drives(self):
        u = self.commanded()
        if self.arch == "multifrequency":
            return [(j,u[j:j+1],torch.ones(1)) for j in range(2)]
        effective, weights = envelope(u, self.c)
        return [(0,effective,weights)]

    def effort(self):
        u = self.commanded().detach()
        if self.arch == "multifrequency":
            return float(u.abs().square().sum())
        return float(u.abs().square().sum()/len(u))


def objective(field, design):
    f,k = field.force_stiffness(design.drives())
    weight = field.c.mass*9.80665
    # Smooth PSD surrogate used for learning; eigenspectrum is checked later.
    eig = torch.linalg.eigvalsh(k)
    shortage = torch.nn.functional.softplus((field.c.k_required-eig)/field.c.k_required)
    force_error = (f/weight).square().mean()
    loss = 20*force_error + 2*shortage.square().mean()
    return loss


def optimize(x, design, steps):
    if design.arch in ("ib", "time_single", "global_local"):
        design.raw.data.copy_(Design(design.arch,x,design.pos,design.c).raw.data)
        return 0.0
    field = Field(x, design.pos, design.freqs, design.c)
    opt = torch.optim.Adam([design.raw], lr=0.025)
    best = math.inf
    best_u = design.raw.detach().clone()
    start = time.perf_counter()
    for _ in range(steps):
        opt.zero_grad()
        loss = objective(field, design)
        value = float(loss.detach())
        if value < best:
            best = value; best_u = design.raw.detach().clone()
        loss.backward()
        torch.nn.utils.clip_grad_norm_([design.raw],10)
        opt.step()
    design.raw.data.copy_(best_u)
    return time.perf_counter()-start


def equilibrium(initial, design, max_steps=25):
    q = initial.clone()
    converged = False
    for _ in range(max_steps):
        f,k = Field(q,design.pos,design.freqs,design.c).force_stiffness(design.drives())
        f,k = f.detach(),k.detach()
        if float(torch.linalg.vector_norm(f,dim=1).max()) < design.c.mass*9.80665*1e-4:
            converged = True; break
        try:
            step = torch.linalg.solve(k,f[...,None]).squeeze(-1)
        except RuntimeError:
            break
        norm = torch.linalg.vector_norm(step,dim=1).clamp_min(1e-30)
        step = step*(0.0005/norm).clamp_max(1)[:,None]
        q += step
        if float(torch.linalg.vector_norm(q-initial,dim=1).max()) > 0.004:
            break
    f,k = Field(q,design.pos,design.freqs,design.c).force_stiffness(design.drives())
    return q, f.detach(), k.detach(), converged


def force_rank(x, design):
    field = Field(x,design.pos,design.freqs,design.c)
    def func(raw):
        # Functional substitution preserves the parameter's graph.
        original = design.raw
        design.raw = raw
        force,_ = field.force_stiffness(design.drives())
        design.raw = original
        return force.reshape(-1)/(design.c.mass*9.80665)
    b = torch.autograd.functional.jacobian(func,design.raw,vectorize=True).reshape(3*len(x),-1)
    sv = torch.linalg.svdvals(b)
    cutoff = max(float(sv[0])*1e-6,1e-12)
    return int((sv>cutoff).sum()), float(sv[-1]/sv[0]) if float(sv[0]) else 0


def periodic_stability(q, design):
    """Piecewise-linear Floquet/affine ripple about averaged equilibrium.

    This checks actual sampled envelope dynamics, not just averaged K > 0.
    No interparticle coupling; validity is local. Spectral beats require separate
    two-time-scale radiation-force modeling and are deliberately not fabricated.
    """
    if design.arch == "multifrequency":
        return None, None
    c=design.c
    u,w = envelope(design.commanded().detach(),c)
    field=Field(q,design.pos,design.freqs,c)
    slot=[]
    for j in range(len(u)):
        f,k=field.force_stiffness([(0,u[j:j+1],torch.ones(1))])
        slot.append((f.detach(),k.detach()))
    max_rho=0.0; max_ripple=0.0
    gamma=6*math.pi*c.viscosity*c.particle_radius
    for particle in range(len(q)):
        transforms=[]
        phi=torch.eye(7)
        for j,(f,k) in enumerate(slot):
            a=torch.zeros((7,7))
            a[:3,3:6]=torch.eye(3)
            a[3:6,:3]=-k[particle]/c.mass
            a[3:6,3:6]=-gamma/c.mass*torch.eye(3)
            a[3:6,6]=f[particle]/c.mass
            tr=torch.matrix_exp(a*c.cycle_s*float(w[j]))
            transforms.append(tr); phi=tr@phi
        rho=float(torch.linalg.eigvals(phi[:6,:6]).abs().max())
        max_rho=max(max_rho,rho)
        if rho>=1:
            return max_rho,None
        fixed=torch.linalg.solve(torch.eye(6)-phi[:6,:6],phi[:6,6])
        y=torch.cat((fixed,torch.ones(1)))
        max_ripple=max(max_ripple,float(torch.linalg.vector_norm(y[:3])))
        for tr in transforms:
            y=tr@y
            max_ripple=max(max_ripple,float(torch.linalg.vector_norm(y[:3])))
    return max_rho,max_ripple


def finite_probes(x,q0,design):
    saved=design.raw.detach().clone()
    gains=[]; crosstalk=[]; error=[]; failures=0
    for j in range(len(x)):
        for axis in range(3):
            desired=x.clone(); desired[j,axis]+=design.c.probe_distance
            design.raw.data.copy_(saved)
            optimize(desired,design,design.c.probe_iterations)
            q,f,k,ok=equilibrium(desired,design)
            err=torch.linalg.vector_norm(q-desired,dim=1)
            eig=torch.linalg.eigvalsh(k)
            own=float((q[j,axis]-q0[j,axis])/design.c.probe_distance)
            others=torch.linalg.vector_norm(q-q0,dim=1)
            others[j]=0
            cross=float(others.max()/design.c.probe_distance)
            gains.append(own); crosstalk.append(cross); error.append(float(err.max()))
            if (not ok or float(eig.min())<design.c.k_required
                    or float(err.max())>design.c.accuracy or abs(own-1)>0.2
                    or cross>design.c.crosstalk_fraction):
                failures+=1
    design.raw.data.copy_(saved)
    return {"finite_probes":3*len(x),"failed_probes":failures,
            "gain_min":min(gains),"gain_max":max(gains),
            "max_crosstalk":max(crosstalk),"probe_error_mm":1000*max(error)}


def isolated_piston_power(c):
    """Nominal real acoustic watts at |u|=1, isolated baffled piston.

    Assumed piston velocity inferred from source_gain. Mutual radiation impedance
    and electrical loss are absent, so this is NOT a measured wall-plug budget.
    """
    k=2*math.pi*c.frequency/c.sound_speed
    area=math.pi*c.emitter_radius**2
    velocity=c.source_gain*2*math.pi/(c.air_density*c.sound_speed*k*area)
    z=2*k*c.emitter_radius
    # J1(z), convergent Taylor series at the supplied 40-kHz piston dimensions.
    term=z/2; bessel=term
    for j in range(1,80):
        term*=-(z*z/4)/(j*(j+1)); bessel+=term
        if abs(term)<1e-15: break
    resistance=1-2*bessel/z
    return 0.5*c.air_density*c.sound_speed*area*velocity**2*resistance


def run_case(c,arch,m,probes=True):
    start=time.perf_counter()
    pos,_=geometry(c,arch=="aperiodic")
    x=targets(c,m)
    d=Design(arch,x,pos,c)
    solve_s=optimize(x,d,c.iterations)
    q,f,k,ok=equilibrium(x,d)
    eig=torch.linalg.eigvalsh(k)
    error=torch.linalg.vector_norm(q-x,dim=1)
    rank,condition=force_rank(q,d)
    rho,ripple=periodic_stability(q,d) if ok else (None,None)
    stable=int(((eig[:,0]>=c.k_required)&(error<=c.accuracy)).sum()) if ok else 0
    row={"architecture":arch,"n":c.n,"requested_particles":m,
         "tau_ms":1000*c.tau_s,"converged":ok,"stable_local_traps":stable,
         "max_error_mm":1000*float(error.max()),"min_stiffness_N_m":float(eig.min()),
         "required_stiffness_N_m":c.k_required,"force_rank":rank,
         "force_rank_required":3*m,"force_singular_ratio":condition,
         "drive_effort":d.effort(),"acoustic_power_estimate_W":d.effort()*isolated_piston_power(c),
         "optimization_s":solve_s,"floquet_radius":rho,
         "ripple_mm":None if ripple is None else 1000*ripple,
         "commands_per_cycle":len(d.commanded()),
         "minimum_spacing_mm":1000*float(torch.pdist(x).min()) if m>1 else None,
         "independence_status":"not_tested"}
    if probes:
        p=finite_probes(x,q,d); row.update(p)
        row["independence_status"]=("tier1_local_pass" if stable==m and rank==3*m
             and p["failed_probes"]==0 and rho is not None and rho<1
             and ripple is not None and ripple<c.accuracy else "tier1_fail_or_incomplete")
    row["runtime_s"]=time.perf_counter()-start
    return row,{"targets_m":x.tolist(),"equilibria_m":q.tolist(),
                "command_real":d.commanded().detach().real.tolist(),
                "command_imag":d.commanded().detach().imag.tolist()}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--n",type=int,default=128)
    p.add_argument("--counts",default="2,4")
    p.add_argument("--architectures",default=",".join(ARCHITECTURES))
    p.add_argument("--iterations",type=int,default=250)
    p.add_argument("--probe-iterations",type=int,default=100)
    p.add_argument("--tau-ms",type=float,default=0.16)
    p.add_argument("--cycle-ms",type=float,default=5)
    p.add_argument("--fd-um",type=float,default=50)
    p.add_argument("--source-gain",type=float,default=2)
    p.add_argument("--seed",type=int,default=7)
    p.add_argument("--no-probes",action="store_true")
    p.add_argument("--out",type=Path,default=Path("results"))
    a=p.parse_args()
    c=Config(n=a.n,iterations=a.iterations,probe_iterations=a.probe_iterations,
             tau_s=a.tau_ms/1000,cycle_s=a.cycle_ms/1000,fd_step=a.fd_um/1e6,
             source_gain=a.source_gain,seed=a.seed)
    if c.n>=512:
        c.gap=0.160; c.volume_side=0.060
    if c.n<=32:
        c.gap=0.060; c.volume_side=0.018
    for arch in a.architectures.split(","):
        if arch not in ARCHITECTURES: raise ValueError(arch)
    if c.tau_s<0 or c.cycle_s<=0 or c.fd_step<=0: raise ValueError("invalid timing or stencil")
    a.out.mkdir(parents=True,exist_ok=True)
    meta={"config":asdict(c),"python":platform.python_version(),"torch":torch.__version__,
          "numpy":np.__version__,"platform":platform.platform(),
          "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          "model_level":"Tier 1: dilute, linear, homogeneous, local only",
          "power_note":"nominal isolated-piston acoustic watts; electrical watts unknown"}
    (a.out/"metadata.json").write_text(json.dumps(meta,indent=2))
    rows=[]; fields={}
    for m in map(int,a.counts.split(",")):
        for arch in a.architectures.split(","):
            row,detail=run_case(c,arch,m,not a.no_probes)
            rows.append(row); fields[f"{arch}_{m}"]=detail
            print(json.dumps(row),flush=True)
            (a.out/"results.json").write_text(json.dumps(rows,indent=2,allow_nan=False))
            (a.out/"fields.json").write_text(json.dumps(fields,indent=2))
    columns=list(dict.fromkeys(k for row in rows for k in row))
    with (a.out/"results.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=columns); writer.writeheader(); writer.writerows(rows)


if __name__=="__main__":
    main()
