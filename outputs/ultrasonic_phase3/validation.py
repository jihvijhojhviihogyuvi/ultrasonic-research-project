"""Strict signed/simultaneous motion, retention, feedback, and resource validation."""
import math,time
import numpy as np
import torch
from threadpoolctl import threadpool_limits
_blas_limit=threadpool_limits(limits=1)
from dynamics import Surrogate,validate_surrogate,movement_protocol,advance_frame


class IncrementalController:
    """Cached force Jacobian with constrained nonlinear feedforward correction.

    The raw emitter controls are recomputed through the same phase/amplitude,
    Fourier peak and radiated-power projection as the optimizer. Pressure
    transfer polynomials avoid a fresh large acoustic solve every camera frame.
    """
    def __init__(self,design,surrogate,refresh_jacobian=False):
        self.design=design;self.surrogate=surrogate;self.system=design.system;c=self.system.c
        self.refresh_jacobian=refresh_jacobian
        if design.architecture=='static':
            # Equivalent feasible coordinates remove artificial saturation of
            # the optimizer's unprojected variables. The initial physical field
            # is unchanged; incremental controls use the available headroom.
            effective=design.fields().detach()*torch.exp(-1j*design.gauge.detach())[:,None]
            design.raw.data.copy_(torch.view_as_real(effective))
        self.raw=design.raw.detach().cpu().numpy().copy()
        self.gauge=design.gauge.detach().cpu().numpy();self.duties=design.duties().detach().cpu().numpy()
        self.ell=self.system.ell.cpu().numpy();self.response=self.system.response().cpu().numpy()
        self.weights=(self.duties[None,:]*np.sinc(self.ell[:,None]*self.duties[None,:])*
                np.exp(2j*np.pi*self.ell[:,None]*(np.cumsum(self.duties)-self.duties/2)[None,:]))
        self.weights*=np.exp(-.5*(2*np.pi*self.ell[:,None]/c.cycle_s*c.smooth_s)**2)
        self.power_coeff=self.system.power_coeff.cpu().numpy();self.support=design.support.cpu().numpy()
        self.common=np.asarray(torch.view_as_complex(design.common.detach().contiguous()).cpu())
        self.residual_scale=float(torch.sigmoid(design.residual_scale.detach()))
        self.phases=np.exp(-1j*np.arange(256)[:,None]*2*np.pi/256*self.ell[None,:])
        phases=torch.arange(max(2*int(abs(self.ell).max())+1,32),device=self.system.device)*2*math.pi/max(2*int(abs(self.ell).max())+1,32)
        def forces(raw):
            original=design.raw;design.raw=raw
            coeff=design.coefficients()
            f=self.system.forces(self.system.q,self.system.jets(coeff),phases).mean(1)[0].reshape(-1)/(c.mass*9.80665)
            design.raw=original;return f
        b=torch.autograd.functional.jacobian(forces,design.raw,vectorize=True).reshape(3*self.system.m,-1).detach().cpu().numpy()
        _,dd,_=self.system.matrices(design.coefficients().detach())
        self.dd=dd.mean(0).detach().cpu().numpy()
        singular=np.linalg.svd(b,compute_uv=False)
        self.rank=int(np.sum(singular>max(singular[0]*1e-7,1e-12)))
        self.singular_ratio=float(singular[-1]/singular[0])
        self.scale=.001
        normalized=self.dd*self.scale/(c.mass*9.80665)
        self.basis=np.linalg.lstsq(b,normalized,rcond=1e-7)[0]
        self.inverse=np.linalg.pinv(normalized,rcond=1e-7)
        self.z=np.zeros(3*self.system.m)
        base,_=self.command(self.z);self.command_base=base
        eps=1e-5;mode=[]
        for j in range(len(self.z)):
            delta=np.zeros_like(self.z);delta[j]=eps
            mode.append((self.command(delta)[0]-self.command(-delta)[0])/(2*eps))
        self.command_modes=np.array(mode)
        self.polynomial_base=surrogate.polynomial(self.response[:,None]*base)
        self.polynomial_modes=np.array([surrogate.polynomial(self.response[:,None]*a) for a in mode])
        self.binned_base=surrogate.binned(base)
        self.binned_modes=np.array([surrogate.binned(a) for a in mode])

    def polynomial(self,command):
        predicted=self.command_base+np.einsum('j,jln->ln',self.z,self.command_modes)
        if np.max(np.abs(predicted-command))<1e-10:
            return self.polynomial_base+np.einsum('j,jlmc->lmc',self.z,self.polynomial_modes)
        return self.surrogate.polynomial(self.response[:,None]*command)

    def binned(self,command):
        predicted=self.command_base+np.einsum('j,jln->ln',self.z,self.command_modes)
        if np.max(np.abs(predicted-command))<1e-10:
            return self.binned_base+np.einsum('j,jlmbc->lmbc',self.z,self.binned_modes)
        return self.surrogate.binned(command)

    def command(self,z):
        c=self.system.c
        raw=self.raw+(self.basis@z).reshape(self.raw.shape)
        u=(raw[...,0]+1j*raw[...,1])*np.exp(1j*self.gauge[:,None])
        arch=self.design.architecture
        if arch in ("ib","singleton"):u=u/np.maximum(np.abs(u),1e-30)
        elif arch=="support":u=self.support[None]+.45*u/np.maximum(np.abs(u),1)
        elif arch=="common":u=self.common[None]+self.residual_scale*u
        if arch=="multifrequency":
            u=u/np.maximum(np.sum(np.abs(u),axis=0),1)[None]
            command=np.zeros((len(self.ell),c.n),complex)
            command[np.where(self.ell==-2)[0][0]]=u[0]
            command[np.where(self.ell==2)[0][0]]=u[1]
        else:
            u=u/np.maximum(np.abs(u),1)
            command=self.weights@u
            peak=np.max(np.abs(self.phases@command),axis=0)+np.pi/256*np.sum(np.abs(self.ell[:,None]*command),axis=0)
            u=u/np.maximum(peak,1)[None];command=self.weights@u
        actual=self.response[:,None]*command
        power=float(np.sum(self.power_coeff*np.abs(actual)**2))
        scale=min(1,math.sqrt(c.power_W/max(power,1e-30)))
        return command*scale,power*scale*scale

    def update(self,target,velocity):
        c=self.system.c;start=time.perf_counter()
        target_force=c.drag*velocity
        for _ in range(2):
            command,power=self.command(self.z)
            polynomial=self.polynomial(command)
            force=self.surrogate.force(target,polynomial,mean=True)-target_force
            normalized=force.reshape(-1)/(c.mass*9.80665)
            # Refresh the small force Jacobian inside the cached emitter basis.
            # This permits finite motion instead of relying on its tangent at
            # the original trap. All projections are still applied exactly.
            if self.refresh_jacobian:
                original=self.z.copy();columns=[];eps=1e-4
                for j in range(len(self.z)):
                    self.z=original.copy();self.z[j]+=eps
                    plus=self.surrogate.force(target,self.polynomial(self.command(self.z)[0]),mean=True)
                    self.z=original.copy();self.z[j]-=eps
                    minus=self.surrogate.force(target,self.polynomial(self.command(self.z)[0]),mean=True)
                    columns.append(((plus-minus)/(2*eps*c.mass*9.80665)).reshape(-1))
                self.z=original
                delta=np.linalg.lstsq(np.stack(columns,axis=1),-normalized,rcond=1e-7)[0]
            else:delta=-self.inverse@normalized
            delta=np.clip(delta,-.15,.15)
            self.z=np.clip(self.z+delta,-1.8,1.8)
        command,power=self.command(self.z)
        return command,power,time.perf_counter()-start


def strict_validate(design,initial_q,out_prefix=None,seconds=60,frame_hz=20,degree=6,
                    bins=8,time_samples=64,noise_m=.00001,delay_s=.020,wind=None,max_particles=16):
    s=design.system;c=s.c;m=s.m
    if m>max_particles:
        return {"status":"not_validated_resource_budget","independent_particles":None,
                "reason":f"strict trajectory budget configured for <= {max_particles} particles"},{}
    start=time.perf_counter()
    surrogate=Surrogate(s,degree=degree,bins=bins)
    check=validate_surrogate(surrogate,design.coefficients().detach())
    if check["maximum_force_error_fraction"]>.01:
        return {"status":"surrogate_convergence_failure","independent_particles":None,"surrogate":check},{}
    controller=IncrementalController(design,surrogate)
    if controller.rank<3*m:
        return {"status":"control_rank_failure","independent_particles":0,"control_rank":controller.rank},{}
    origin=s.q.detach().cpu().numpy();q=np.asarray(initial_q).copy();v=np.zeros_like(q)
    protocol=movement_protocol(origin,frame_hz)
    while len(protocol)<seconds*frame_hz:protocol.append((origin.copy(),np.zeros(m,bool),"retention"))
    dt_frame=1/frame_hz;dt=c.cycle_s/time_samples
    steps=round(dt_frame/dt);dt=dt_frame/steps
    delay_steps=surrogate.delay_s/dt
    hist_len=math.ceil(delay_steps.max())+4
    response=s.response().cpu().numpy()
    command,_=controller.command(controller.z)
    envelope=surrogate.binned(command)*response[:,None,None,None]
    history=np.tile(envelope[None],(hist_len,1,1,1,1));cursor=hist_len-1
    source=response[:,None]*command
    qhistory=np.tile(q[None],(math.ceil(delay_s/dt)+4,1,1));qcursor=len(qhistory)-1
    previous_meas=q.copy();estimate=q.copy();velocity_est=np.zeros_like(q);bias=np.zeros_like(q)
    previous_target=origin.copy();held_origin=q.copy();rng=np.random.default_rng(c.seed+900)
    rows=[];zs=[];measurements=[];latencies=[];failure=None;completed_labels=set()
    maximum_error=0.;maximum_held=0.;maximum_peak=0.;minimum_stiffness=float("inf");maximum_power=0.
    transient_drive_energy=0.;maximum_transient_power_bound=0.;maximum_frame_power_bound=0.
    power_weights=s.power_coeff.max(0).values.cpu().numpy()
    wind=np.zeros(3) if wind is None else np.asarray(wind,float)
    # Compile the trajectory kernel before timing the active frame controller.
    warm=advance_frame(q.copy(),v.copy(),envelope.copy(),history.copy(),cursor,surrogate.binned(command),
        source.copy(),command,delay_steps,qhistory.copy(),qcursor,1,dt,0.,origin,surrogate.scale,surrogate.exp,
        surrogate.omega,surrogate.ell,c.particle_radius,c.f0,c.f2,c.rho,c.sound_speed,c.mass,c.drag,c.viscosity,
        c.settling_s,c.cycle_s,origin,origin,np.ones(m,bool),held_origin,wind,power_weights)
    controller.update(origin,np.zeros_like(origin))
    surrogate.stiffness(origin,surrogate.polynomial(response[:,None]*command))
    controller.z[:]=0.
    for frame,(target,moved,label) in enumerate(protocol):
        frame_start=time.perf_counter()
        sample=qhistory[(qcursor-round(delay_s/dt))%len(qhistory)]
        meas=sample+rng.normal(0,noise_m,sample.shape)
        velocity_est=.7*velocity_est+.3*(meas-previous_meas)/dt_frame
        estimate=meas+delay_s*velocity_est;previous_meas=meas
        # Slow camera correction supplements the cached model-based feedforward.
        bias=np.clip(bias+.01*(previous_target-estimate),-.00005,.00005)
        desired_velocity=(target-previous_target)/dt_frame
        command,power,latency=controller.update(target+bias,desired_velocity)
        online_latency=time.perf_counter()-frame_start
        maximum_power=max(maximum_power,power)
        steady=controller.polynomial(command)
        dd=surrogate.stiffness(target,steady)
        stiffness=float(np.linalg.eigvalsh((dd+dd.T)/2).min());minimum_stiffness=min(minimum_stiffness,stiffness)
        if stiffness<c.stiffness_required:
            failure="moving_trap_stiffness_failure";break
        binned=controller.binned(command)
        # Stiffness certification and acoustic-polynomial command transfers are
        # offline simulation oracles. The hardware command is already complete.
        oracle_latency=time.perf_counter()-frame_start
        latencies.append(online_latency)
        result=advance_frame(q,v,envelope,history,cursor,binned,source,command,delay_steps,qhistory,qcursor,
            steps,dt,frame*dt_frame,origin,surrogate.scale,surrogate.exp,surrogate.omega,surrogate.ell,
            c.particle_radius,c.f0,c.f2,c.rho,c.sound_speed,c.mass,c.drag,c.viscosity,c.settling_s,c.cycle_s,
            previous_target,target,~moved,held_origin,wind,power_weights)
        q,v,envelope,history,cursor,source,qhistory,qcursor,error,held,peak,problem,drive_energy,bound=result
        transient_drive_energy+=drive_energy
        # An upper bound in the same isolated-piston resistance model: use the
        # largest retained-band resistance per source. This conservatively
        # handles frequency mixtures during frame-to-frame filter transients.
        maximum_transient_power_bound=max(maximum_transient_power_bound,bound)
        maximum_frame_power_bound=max(maximum_frame_power_bound,drive_energy/dt_frame)
        maximum_error=max(maximum_error,error);maximum_held=max(maximum_held,held);maximum_peak=max(maximum_peak,peak)
        rows.append({"t_s":(frame+1)*dt_frame,"label":label,"stiffness_N_m":stiffness,
                     "position_error_m":error,"held_displacement_m":held,"power_W":power,
                     "controller_s":latencies[-1],"controller_and_oracle_s":oracle_latency,
                     "transient_power_bound_W":bound,
                     "positions_m":q.tolist(),"targets_m":target.tolist()})
        zs.append(controller.z.copy());measurements.append(meas.copy());completed_labels.add(label)
        previous_target=target.copy()
        if problem:failure="collision" if problem==1 else "surrogate_domain_exceeded";break
        if error>c.accuracy:failure="position_error_failure";break
        if held>c.accuracy:failure="held_particle_displacement_failure";break
        if peak>1+1e-6:failure="source_peak_failure";break
        if drive_energy/dt_frame>c.power_W*(1+1e-6):failure="transient_average_power_bound_failure";break
    duration=len(rows)*dt_frame
    physical_pass=failure is None and duration>=60 and len(rows)==len(protocol)
    compute_pass=bool(latencies) and max(latencies)<=dt_frame
    strict=physical_pass and compute_pass
    status="strict_pass" if strict else "compute_rate_failure" if physical_pass else failure or "incomplete"
    result={"status":status,"independent_particles":m if strict else 0,
            "motion_contract_pass":physical_pass,"real_time_compute_feasible":compute_pass,
            "duration_s":duration,"required_protocol_duration_s":len(protocol)*dt_frame,
            "frames_completed":len(rows),"frames_required":len(protocol),"update_hz":frame_hz,
            "maximum_position_error_mm":maximum_error*1000,"maximum_held_displacement_mm":maximum_held*1000,
            "minimum_moving_stiffness_N_m":minimum_stiffness,"maximum_source_peak":maximum_peak,
            "maximum_steady_acoustic_power_W":maximum_power,
            "maximum_transient_power_bound_W":maximum_transient_power_bound,
            "maximum_frame_average_power_bound_W":maximum_frame_power_bound,
            "integrated_acoustic_energy_bound_J":transient_drive_energy,
            "controller_mean_s":float(np.mean(latencies)) if latencies else None,
            "controller_timing_scope":"camera processing, correction and emitted spectral command; offline stiffness/field oracles excluded",
            "controller_max_s":max(latencies) if latencies else None,
            "control_rank":controller.rank,"control_singular_ratio":controller.singular_ratio,
            "surrogate":check,"validation_runtime_s":time.perf_counter()-start,
            "noise_m":noise_m,"camera_delay_s":delay_s,"wind_m_s":wind.tolist(),
            "time_step_s":dt,"delay_bins":bins,"tests_completed_labels":sorted(completed_labels),
            "model_note":"validated pressure polynomial; Born scattering; retarded source-center envelope bins; narrowband feedback transients"}
    raw={"frames":rows,"z":np.array(zs),"measurements_m":np.array(measurements),
         "controller_basis":controller.basis,"raw_base":controller.raw,
         "polynomial_exponents":surrogate.exp,"pressure_transfer_sha":None}
    import hashlib
    raw["pressure_transfer_sha"]=hashlib.sha256(surrogate.transfer.tobytes()).hexdigest()
    return result,raw
