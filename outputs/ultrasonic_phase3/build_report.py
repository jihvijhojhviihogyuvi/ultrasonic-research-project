"""Build all tables/figures/report directly from saved executed calculations."""
import csv,json,math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from experiments import ROOT,write


def read(path,default=None):
    return json.loads(path.read_text(encoding='utf-8')) if path.exists() else default


def number(v,scale=1,digits=4):
    return '—' if v is None else f'{v*scale:.{digits}g}'


def table(rows,columns):
    result=['| '+' | '.join(c[0] for c in columns)+' |','| '+' | '.join('---' for _ in columns)+' |']
    for r in rows:result.append('| '+' | '.join(str(fn(r)) for _,fn in columns)+' |')
    return '\n'.join(result)


def main():
    results=ROOT/'results';fig=ROOT/'figures';fig.mkdir(exist_ok=True)
    cases=[read(p) for p in sorted((results/'experiments').glob('*/case.json'))]
    legacy=read(results/'reproduction/results.json',[])
    independent=read(results/'reproduction/independent_physics.json',[])
    motion=[read(p) for p in sorted((results/'motion').glob('*/result.json'))]
    convergence=read(results/'convergence/results.json',[])
    sixteen=read(results/'convergence/sixteen_particle.json',[])
    matched=read(results/'matched_power/results.json',[])
    recovery=read(results/'sdp_recovery/results.json',[])
    density=read(results/'density_stability/results.json',[])
    nonlinear=[]
    for path in sorted((results/'convergence').glob('*_nonlinear_*.json')):
        nonlinear.append(dict(case=path.name.split('_nonlinear_')[0],**read(path)))
    byid={r['id']:r for r in cases};passing=[r for r in motion if r.get('status')=='strict_pass']
    all_columns=['id','category','architecture','transducers','particles','field_count','status',
                 'weakest_primary_stiffness_N_m','weakest_coupled_stiffness_N_m','equilibrium_error_mm',
                 'floquet_radius','acoustic_power_W','optimization_s','runtime_s','independent_particles']
    with (results/'all_experiments.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=all_columns,extrasaction='ignore');writer.writeheader();writer.writerows(cases)
    write(results/'all_experiments.json',cases);write(results/'all_motion.json',motion)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':130})
    def save(name):plt.tight_layout();plt.savefig(fig/name,dpi=160);plt.close()
    # Independent recalculation of the old quadrature.
    rr=[r for r in independent if r['case']=='results_512_screen_covariance_8']
    plt.figure(figsize=(6.5,3.5));plt.plot([r['samples_per_slot'] for r in rr],[r['system_floquet_radius'] for r in rr],'o-')
    plt.axhline(1,color='firebrick',ls='--');plt.xlabel('Samples per temporal slot');plt.ylabel('Floquet spectral radius')
    plt.title('Legacy eight-particle instability survives refinement');save('legacy_refinement.png')
    archs=['ib','static','multifrequency','singleton','covariance','support']
    figobj,axs=plt.subplots(1,3,figsize=(13,4),sharey=True)
    for ax,m in zip(axs,(2,8,16)):
        rr=[byid.get(f'compare_{m}_{a}',{}) for a in archs]
        vals=[(r.get('weakest_coupled_stiffness_N_m') if r.get('weakest_coupled_stiffness_N_m') is not None else r.get('weakest_primary_stiffness_N_m',np.nan))*1e6 for r in rr]
        colors=['#bd4b45' if (r.get('floquet_radius') or 0)>1 else '#317d8e' for r in rr]
        ax.bar(range(len(rr)),vals,color=colors);ax.axhline(17.44103,color='black',ls='--',lw=1)
        ax.set_xticks(range(len(rr)),archs,rotation=55,ha='right');ax.set_title(f'{m} particles, N=512')
    axs[0].set_ylabel('Weakest stiffness (µN/m)');figobj.suptitle('Red = periodic instability; primary-only screens use the uncoupled value')
    save('architecture_comparison.png')
    figobj,axs=plt.subplots(1,2,figsize=(9,3.6))
    rr=[byid.get(f'slots_{s}',{}) for s in (2,3,4,8,16)]
    xx=[r.get('field_count',0) for r in rr];kk=[r.get('weakest_coupled_stiffness_N_m',np.nan)*1e6 for r in rr]
    rho=[r.get('floquet_radius',np.nan) for r in rr]
    axs[0].plot(xx,kk,'o-',color='#317d8e');axs[0].axhline(17.44103,color='black',ls='--');axs[0].set_ylabel('Weakest stiffness (µN/m)')
    axs[1].plot(xx,rho,'o-',color='#bd4b45');axs[1].axhline(1,color='black',ls='--');axs[1].set_ylabel('Floquet radius')
    for ax in axs:ax.set_xlabel('Temporal fields');ax.set_xticks((2,3,4,8,16))
    figobj.suptitle('Sequence length: two-particle fields at matched resource ceilings');save('temporal_fields.png')
    # Certified counts are lower bounds only. No unvalidated cloud is promoted.
    figobj,axs=plt.subplots(1,2,figsize=(10,3.8))
    for arch,color in [('static','#317d8e'),('covariance','#bd4b45')]:
        counts=[]
        for n in (32,64,128,256,512):
            counts.append(max([r.get('independent_particles',0) for r in passing
                if byid.get(r['case'],{}).get('transducers')==n and byid.get(r['case'],{}).get('architecture')==arch]+[0]))
        axs[0].plot((32,64,128,256,512),counts,'o-',label=arch,color=color)
    axs[0].set_xlabel('Transducers');axs[0].set_ylabel('Certified particle count (lower bound)');axs[0].legend()
    for arch,color in [('static','#317d8e'),('covariance','#bd4b45')]:
        rr=sorted([r for r in cases if r['category']=='power_density' and r['id'].startswith('power_') and r['architecture']==arch],key=lambda r:r['config']['power_W'])
        pp=[r['config']['power_W'] for r in rr]
        verified=[max([v.get('independent_particles',0) for v in passing if
            byid.get(v['case'],{}).get('transducers')==128 and byid.get(v['case'],{}).get('architecture')==arch
            and abs(byid.get(v['case'],{}).get('config',{}).get('power_W',-1)-power)<1e-9]+[0]) for power in pp]
        axs[1].plot(pp,verified,'o-',color=color,label=arch)
    axs[1].set_xlabel('Acoustic power ceiling (W), N=128');axs[1].set_ylabel('Certified particle count (lower bound)')
    figobj.suptitle('Zero means no complete passing configuration certified; it is not a capacity limit')
    save('certified_counts.png')
    plt.figure(figsize=(6.5,3.6))
    slots=(1,2,3,4,8,16)
    verified=[max([r.get('independent_particles',0) for r in passing if
                   byid.get(r['case'],{}).get('transducers')==512 and byid.get(r['case'],{}).get('field_count')==s
                   and byid.get(r['case'],{}).get('architecture')!='multifrequency']+[0]) for s in slots]
    plt.plot(slots,verified,'o-',color='#317d8e');plt.xticks(slots);plt.xlabel('Acoustic fields, N=512')
    plt.ylabel('Certified particle count (lower bound)');plt.title('Temporal fields; S=1 reference; two-tone fields excluded')
    save('certified_counts_fields.png')
    # Raw stiffness and cost screens.
    figobj,axs=plt.subplots(1,2,figsize=(10,3.8))
    for arch,color in [('static','#317d8e'),('covariance','#bd4b45')]:
        rr=sorted([r for r in cases if r['category']=='power_density' and r['id'].startswith('density_') and r['architecture']==arch],key=lambda r:r['cloud_density_per_cm3'])
        axs[0].plot([r['cloud_density_per_cm3'] for r in rr],[r['weakest_primary_stiffness_N_m']*1e6 for r in rr],'o-',label=arch,color=color)
        rr=sorted([r for r in cases if r['category']=='scaling' and r['transducers']==512 and r['architecture']==arch],key=lambda r:r['particles'])
        axs[1].plot([r['particles'] for r in rr],[r['runtime_s'] for r in rr],'o-',label=arch,color=color)
    axs[0].axhline(17.44103,color='black',ls='--');axs[0].set_xscale('log');axs[0].set_xlabel('Cloud density (particles/cm³)');axs[0].set_ylabel('Primary stiffness (µN/m)')
    axs[1].set_yscale('log');axs[1].set_xlabel('Requested particles, N=512');axs[1].set_ylabel('Executed case wall time (s)')
    axs[0].legend();axs[1].legend();figobj.suptitle('Density screens and computation; high-count cases omit coupled dynamics')
    save('density_cost.png')
    plt.figure(figsize=(7,3.8))
    for arch,color in [('static','#317d8e'),('covariance','#bd4b45')]:
        rr=sorted([r for r in density if r['architecture']==arch and r['equilibrium_converged']],key=lambda r:r['cloud_density_per_cm3'])
        if rr:plt.plot([r['cloud_density_per_cm3'] for r in rr],[r['floquet_radius'] for r in rr],'o-',color=color,label=arch)
    plt.axhline(1,color='black',ls='--');plt.xscale('log');plt.xlabel('Cloud density (particles/cm³), N=128/M=8')
    plt.ylabel('Coupled Floquet spectral radius');plt.title('Only converged gravity equilibria are plotted')
    if any(r['equilibrium_converged'] for r in density):plt.legend()
    else:plt.text(.5,.5,'No matched equilibrium converged in these six attempts',transform=plt.gca().transAxes,ha='center')
    save('density_stability.png')
    figobj,axs=plt.subplots(1,2,figsize=(11,4))
    for ax,arch in zip(axs,('static','covariance')):
        matrix=np.full((5,8),np.nan)
        for i,n in enumerate((32,64,128,256,512)):
            for j,m in enumerate((2,4,8,16,25,32,50,100)):
                r=byid.get(f'scale_N{n}_M{m}_{arch}')
                if r and r.get('weakest_primary_stiffness_N_m') is not None:matrix[i,j]=r['weakest_primary_stiffness_N_m']*1e6
        im=ax.imshow(matrix,aspect='auto',cmap='RdYlBu',vmin=-5,vmax=35)
        ax.set_xticks(range(8),(2,4,8,16,25,32,50,100));ax.set_yticks(range(5),(32,64,128,256,512));ax.set_title(arch)
        ax.set_xlabel('Requested particles');ax.set_ylabel('Sources')
        for i in range(5):
            for j in range(8):
                if np.isfinite(matrix[i,j]):ax.text(j,i,f'{matrix[i,j]:.1f}',ha='center',va='center',fontsize=7)
    cax=figobj.add_axes([.92,.20,.018,.58]);figobj.colorbar(im,cax=cax,label='Primary stiffness (µN/m)')
    figobj.suptitle('Executed local-search screens at fixed area/aperture/power; not independent counts')
    figobj.subplots_adjust(top=.83,bottom=.16,right=.88,wspace=.25);plt.savefig(fig/'scaling_stiffness.png',dpi=160);plt.close()
    figobj,axs=plt.subplots(1,2,figsize=(10,3.8))
    for name,color in [('compare_2_static','#317d8e'),('compare_2_covariance','#bd4b45')]:
        folders=sorted((results/'motion').glob(f'{name}_*'))
        folders=[p for p in folders if (p/'frames.json').exists()]
        if not folders:continue
        # Prefer the greatest completed trajectory, preserving all other files.
        folder=max(folders,key=lambda p:read(p/'result.json',{}).get('frames_completed',0))
        frames=read(folder/'frames.json',[])
        if not frames:continue
        axs[0].plot([r['t_s'] for r in frames],[r['position_error_m']*1000 for r in frames],label=name,color=color,lw=.8)
        replay=results/'controller_replay'/folder.name/'timings_s.npy'
        latency=np.load(replay)*1000 if replay.exists() else [r['controller_s']*1000 for r in frames]
        axs[1].plot([r['t_s'] for r in frames],latency,label=name+(' replay' if replay.exists() else ''),color=color,lw=.8)
    axs[0].axhline(.1,color='black',ls='--');axs[1].axhline(50,color='black',ls='--')
    axs[0].set_ylabel('Within-frame maximum error (mm)');axs[1].set_ylabel('Online computation (ms)')
    for ax in axs:ax.set_xlabel('Simulated time (s)');ax.legend(fontsize=7)
    figobj.suptitle('Finite motion and measured CPU deadlines');save('motion_validation.png')
    # Write the nine-part report with numerical tables drawn from these files.
    cols=[('Case',lambda r:r['id']),('N/M/S',lambda r:f"{r['transducers']}/{r['particles']}/{r['field_count']}"),
          ('K primary µN/m',lambda r:number(r.get('weakest_primary_stiffness_N_m'),1e6)),
          ('K coupled µN/m',lambda r:number(r.get('weakest_coupled_stiffness_N_m'),1e6)),
          ('ρ Floquet',lambda r:number(r.get('floquet_radius'))),('Error mm',lambda r:number(r.get('equilibrium_error_mm'))),
          ('Power W',lambda r:number(r.get('acoustic_power_W'))),('Time s',lambda r:number(r.get('runtime_s'))),('Outcome',lambda r:r['status'])]
    lines=['# Phase 3: computational validation of ultrasonic particle control','',
      'Research date: 25 November 2026. All values below are executed model calculations, unless explicitly labeled analytical or published. Original Phase-2 deliverables are preserved.', '',
      '**The study does not establish a temporal-covariance advantage in independently controlled particle count.** Covariance can improve stiffness per modeled acoustic watt and can reach a marginal 16-particle stationary configuration. Other sequences are parametrically unstable. A complete motion contract, numerical convergence, and CPU deadlines are separate acceptance conditions.', '',
      f'The saved study contains {len(legacy)} fresh legacy reproductions, {len(independent)} independent legacy recalculations, {len(cases)} new architecture/scaling attempts, {len(convergence)+len(sixteen)} fixed-drive convergence calculations, {len(density)} coupled density diagnostics, and {len(motion)} motion evidence records (including prechecks and controller replays). {len(passing)} records meet every recorded strict condition; the largest qualified configuration has {max([r.get("independent_particles",0) for r in passing]+[0])} particles. Counts are not hardware records or global optima.', '',
      '## Part 1: Reproduction','',
      'All 13 archived stiffness, equilibrium-error and available Floquet results reproduced exactly in fresh runs from verified source snapshots. This verifies implementation reproducibility, not missing physics. Independent NumPy pressure finite differences and SciPy matrix exponentials import no legacy physics functions.', '',
      table(legacy,[('Case',lambda r:r['case']),('Original K µN/m',lambda r:number(r['old_stiffness_N_m'],1e6)),('Fresh K µN/m',lambda r:number(r['min_stiffness_N_m'],1e6)),('Fresh error mm',lambda r:number(r.get('max_error_mm'))),('ρ',lambda r:number(r.get('floquet_radius'))),('Power W',lambda r:number(r.get('acoustic_power_estimate_W'))),('Runtime s',lambda r:number(r.get('fresh_execution_s')))]), '',
      'The old eight-particle temporal instability becomes stronger under time refinement at its archived equilibrium:', '',
      table(rr if False else [r for r in independent if r['case']=='results_512_screen_covariance_8'],[('Samples/slot',lambda r:r['samples_per_slot']),('K µN/m',lambda r:number(r['independent_min_stiffness_N_m'],1e6,6)),('ρ',lambda r:number(r['system_floquet_radius'],1,8)),('Gravity residual/weight',lambda r:number(r['gravity_residual_fraction']))]), '',
      '![Legacy refinement](figures/legacy_refinement.png)','',
      'The archived 1-ms sequence reproduces a stable multiplier near 0.979 but insufficient stiffness near 15.58 µN/m. The archived 87.1% sixteen-particle improvement was a stiffness ratio between two insufficient solutions; it never established a count improvement. Refined averaging also moves the gravity-balanced equilibrium, so fixed archived positions cannot be treated as converged physical equilibria.', '',
      '## Part 2: Simulation Framework','',
      'The implementation uses finite baffled-disk Rayleigh quadrature, analytic pressure derivatives through third order, complex pressure/velocity, actual sideband frequencies, coherent beat forces and retarded periodic envelopes. The source response is Hℓ=1/(1−iℓΩτ); Gaussian transition smoothing and unequal cyclic duties shape the command spectrum. Born and Foldy monopole/dipole scattering, full coupled force Jacobians, gravity, finite-Re drag, contact impulses and nonlinear RK4 trajectories are implemented.', '',
      'With peak phasors exp(−iωt), the governing potential is U=V[f₀|p|²/(4ρc²)−3f₂ρ|v|²/8], v=∇p/(iωρ), and net force is −∇U−mg ẑ. Pressure is Pa; velocity is m/s; U is J; force is N; stiffness is N/m. The negative velocity term and its factors are tested against an analytical standing wave. Full D=−∂F/∂q and state-matrix/Floquet eigenvalues are evaluated, rather than using a positive Laplacian.', '',
      'The single carrier, IB phase holography, force optimization, singleton sharing, rank-S covariance, coherent support/residual, permanent support and actual two-tone controllers are executable. Adam/autodiff, unequal timing/gauge optimization, ordering search, a small PSD relaxation, force-control conditioning and a discrete matched-geometry/field search were executed. An incremental camera-like controller uses a cached force basis, bounded drive updates and noisy delayed position estimates. Nonlinear flow-Jacobian optimization is a small pilot, not a demonstrated optimum.', '',
      'The main benchmark fixes 160-mm aperture and gap, a 60-mm cubic workspace, 0.032572 m² emitting area, a 1.492952-W isolated-piston acoustic-power ceiling, 40 kHz carrier, 0.16-ms emitter settling and 0.5-mm-diameter particles of density 30 kg/m³. The main weakest-mode threshold is 17.44103 µN/m, corresponding to 15 Hz for this particle. Source radii for N=32,64,128,256,512 are 18,12.728,9,6.364,4.5 mm. This scaling family requires different physical emitter sizes; ordinary 9-mm transducers do not implement every row.', '',
      'Power is calculated from isolated baffled-piston radiation resistance, not electrical input or measured array power. Cross-source mutual radiation impedance is absent. Per-source command peaks, spectral power and a conservative band-resistance bound on transient energy are checked. The power ceiling applies to cycle/frame-average power; source drive peaks apply instantaneously. The reported conservative instantaneous power bound can exceed the average ceiling. Steady power ceilings are matched within comparisons; realized power is listed. The supplementary fixed-field attenuation study below matches realized static/covariance power exactly.', '',
      'Source files: [physics.py](physics.py), [optimization.py](optimization.py), [dynamics.py](dynamics.py), [validation.py](validation.py), [experiments.py](experiments.py), [convergence.py](convergence.py), [convergence_high.py](convergence_high.py), [validate_cases.py](validate_cases.py), [sdp_recovery.py](sdp_recovery.py), [power_match.py](power_match.py), [tests.py](tests.py), [tests_dynamics.py](tests_dynamics.py). See [README.md](README.md) for commands and dependencies.', '',
      '## Part 3: Experimental Results','',
      'The complete numerical inventory appears in Appendix A and [all_experiments.csv](results/all_experiments.csv). Completed case directories retain configuration, seed, geometry, complex coefficients, duties, optimizer history and full coupled arrays when those were computed. The original aborted nonlinear pilot retains its configuration and traceback; its intermediate weights were not preserved. The final runner also snapshots parameters on future execution errors. Screened rejections remain visible. A dash means the metric was not evaluated, not zero.', '',
      table([r for r in cases if r['category']=='comparisons'],cols), '',
      '![Architecture comparison](figures/architecture_comparison.png)','',
      'Exactly matched realized power, with each optimized field attenuated to the common lower power and its coupled equilibrium recomputed:', '',
      table(matched,[('M',lambda r:r['particles']),('Architecture',lambda r:r['architecture']),('Power W',lambda r:number(r['matched_realized_power_W'],1,7)),('K µN/m',lambda r:number(r['weakest_stiffness_N_m'],1e6)),('Error mm',lambda r:number(r['equilibrium_error_mm'])),('ρ',lambda r:number(r['floquet_radius']))]), '',
      'These are fixed-field sensitivity comparisons; attenuation is not a fresh global optimization at the lower budget. The original 16-particle static primary screen rejected early; this supplement computes its equilibrium and confirms that its stiffness and position remain insufficient.', '',
      'Motion validations:', '',
      table(motion,[('Design',lambda r:r.get('case','unknown')),('Status',lambda r:r['status']),('Duration s',lambda r:number(r.get('duration_s'))),('Max error mm',lambda r:number(r.get('maximum_position_error_mm'))),('Held mm',lambda r:number(r.get('maximum_held_displacement_mm'))),('Min K µN/m',lambda r:number(r.get('minimum_moving_stiffness_N_m'),1e6)),('Max CPU ms',lambda r:number(r.get('controller_max_s'),1000)),('Strict count',lambda r:r.get('independent_particles','—'))]), '',
      '![Finite motion](figures/motion_validation.png)','',
      'The motion protocol commands every particle by ±0.25 and ±1 mm in X/Y/Z, returns to baseline, and moves distinct particles simultaneously. Cubic ramps last 0.4 s for 0.25 mm and 2 s for 1 mm. Their maximum speeds are 0.938 and 0.750 mm/s. Two-particle protocols run 76 simulated seconds. The 20-Hz command rate does not mean 20-Hz full-amplitude motion. Position error is the maximum within each integration frame, not only a favorable frame endpoint. Camera noise is 10 µm rms per coordinate with 20-ms delay in the standard test.', '',
      'The first motion harness included offline stiffness certification and acoustic-polynomial transfer generation in its deadline. Those are verification/simulation oracles, rather than emitted-drive control. Its compute-rate failures remain preserved. The separate replay times camera processing, force correction and actual spectral packet generation, checks every control coordinate against the saved trajectory, and combines a passing replay only with an already completed motion contract. Integration, calibration, global optimization, oracle checks, surrogate construction and JIT warm-up remain offline. These results do not certify an OS hard real-time guarantee or real camera/FPGA transport.', '',
      '## Part 4: Temporal Covariance Analysis','',
      table([r for r in cases if r['category']=='temporal' and not r['id'].startswith('geometry_')],cols), '',
      '![Temporal fields](figures/temporal_fields.png)','',
      'Two and three fields produce high averaged stiffness but unstable periodic mechanics in the tested two-particle initializations. Four, eight and sixteen fields can be locally stable. More fields do not monotonically improve stiffness. Filtering makes short slots increasingly redundant, while coherent transitions and phase gauges alter both force and mechanical modulation.', '',
      'Unequal duties and gauges are optimized jointly with fields in every covariance run. The timing-only and ordering studies hold the mode family while adjusting its schedule; ordering rows are saved in case.json. Direct Floquet penalties constrain the complete coupled linearized sequence. The initial nonlinear-gradient pilot failed with NaNs from a second derivative of velocity norm at zero. A negligible velocity regularizer corrected this; the repeated two-iteration pilot executed, but did not meet the main stiffness threshold. The failed attempt remains saved.', '',
      'Nonlinear shooting of the two-particle covariance sequence gives ρ=0.89346283 at 64 steps/cycle and 0.89346297 at 128, with approximately 0.007081-mm periodic position excursion. Both shooting closures converged. The static result gives ρ≈0.897077805 at both resolutions. The finite-amplitude orbit differs from a Floquet calculation centered at a mean equilibrium; both are useful and neither alone establishes finite independent motion.', '',
      'The sixteen-particle candidate receives additional fixed-drive checks:', '',
      table(sixteen,[('Harmonics',lambda r:r['harmonics']),('Samples',lambda r:r['time_samples']),('Disk refinement',lambda r:r['disk_refinement']),('FD µm',lambda r:number(r['fd_step_m'],1e6)),('K µN/m',lambda r:number(r['weakest_stiffness_N_m'],1e6,7)),('ρ',lambda r:number(r['floquet_radius'],1,8)),('Gravity residual/weight',lambda r:number(r['gravity_residual_fraction'])),('Power W',lambda r:number(r['power_W'],1,8))]), '',
      'Retaining ±16 harmonics increases modeled power slightly beyond the original ceiling when the original fields are held fixed. That row is a convergence sensitivity, not a fair passing configuration. Quadrature refinement also changes the gravity residual at the old equilibrium. Reprojection, rebalancing and finite-motion validation would be required before accepting this marginal candidate.', '',
      '## Part 5: New Discoveries','',
      'The useful improvement is a larger feasible force-design space, not a universal frequency or slot multiplier. At exactly 0.862809447 W in the two-particle fixed-field comparison, covariance reaches 73.475 µN/m against 61.245 µN/m for static control: approximately 20.0% more stiffness. This is an energy/stiffness observation, not an independent-count gain. Its two-particle nonlinear motion passes the 76-s trajectory contract, but the identical controller replay has a 57.74-ms maximum update against the strict 50-ms deadline. At eight particles the covariance sequence is weaker and unstable. At sixteen particles it raises stationary stiffness to 17.56 µN/m, but has not passed full motion or robustness.', '',
      'Coherent support-residual control also qualifies two particles at N=512. It uses 71.272 J over the same 76-s protocol against 92.057 J for the unstructured static force design: 22.6% lower integrated modeled acoustic energy. Held drift is 0.00430 mm rather than 0.03025 mm. Its weakest moving stiffness is lower, 26.04 versus 76.73 µN/m, but both exceed the benchmark. A shared coherent support field therefore does not inherently prevent local independent motion. This remains one target cloud, one idealized material and one slow motion class.', '',
      'The PSD relaxation has a meaningful positive control: under the separate 5-Hz threshold, N=32 gives an ideal covariance at 0.603745 W with gravity residual approximately 2.7×10⁻¹³ of weight. Its recovered finite sequences lose force balance and stiffness when rank truncation, source peaks and filtering are enforced. Refinement improves their fields but does not erase the residual. This directly demonstrates why an ideal PSD solution cannot be advertised as a waveform certificate.', '',
      table(recovery,[('Stage',lambda r:r['test']),('Rank',lambda r:r.get('rank','full')),('K µN/m',lambda r:number(r['weakest_stiffness_N_m'],1e6)),('Gravity residual/weight',lambda r:number(r['gravity_residual_fraction'])),('Power W',lambda r:number(r['power_W'])),('ρ',lambda r:number(r.get('floquet_radius')))]), '',
      'The main 15-Hz SDP returned infeasible_inaccurate; that status is not a rigorous impossibility proof. The matched jitter/field searches changed stiffness modestly, without an independence advantage. Other geometry-seed cases also changed target clouds and are reported as robustness samples, not fair geometry-only improvements. Control-Jacobian optimization was executed as a small pilot; rank and conditioning are necessary, not sufficient, for finite control.', '',
      '## Part 6: Scaling Results','',
      'N=32,64,128,256,512 and M=2,4,8,16,25,32,50 were attempted for both static and four-field covariance control. M=100 was attempted at N=512. Emitting area, aperture, workspace, material, source peak and power ceiling are fixed. Targets are identical between architectures and source counts for each M; different M values use different seeded clouds, not nested subsets.', '',
      '![Certified counts](figures/certified_counts.png)','',
      '![Certified field-count lower bounds](figures/certified_counts_fields.png)','',
      'The displayed count is the largest fully certified count among executed tests, a lower bound. Zero means this study has not certified a complete configuration at that setting. It does not mean that no configuration can work. Stationary candidates are never included. Power-sweep points without full motion remain uncertified. There is no defensible empirical capacity curve beyond these lower bounds.', '',
      '![Scaling screens](figures/scaling_stiffness.png)','',
      '![Density and computation](figures/density_cost.png)','',
      '![Coupled density stability](figures/density_stability.png)','',
      table(density,[('Case',lambda r:r['case']),('Density /cm³',lambda r:number(r['cloud_density_per_cm3'])),('Root converged',lambda r:r['equilibrium_converged']),('Position error mm',lambda r:number(r['equilibrium_error_mm'])),('K coupled µN/m',lambda r:number(r['weakest_coupled_stiffness_N_m'],1e6)),('ρ',lambda r:number(r['floquet_radius'])),('log(ρ)/T per s',lambda r:number(r['growth_rate_per_s']))]),'',
      'Locally stable linearized periodic dynamics can coexist with insufficient stiffness or target accuracy. These density rows are diagnostics at mean gravity equilibria; no eight-particle nonlinear orbit or retention was certified. Nonconverged roots are retained in the table but not plotted as equilibrium stability certificates. Different cycle lengths must be compared using log(ρ)/T, not raw per-cycle multipliers.', '',
      'Large clouds receive an uncoupled primary-force/stiffness screen. Their coupled Floquet, camera motion and 60-s retention values are explicitly absent. Cases over sixteen particles use fewer optimization iterations; shortfall in a local search is not an upper bound. Some sixteen-particle results change materially between 80 and 180 iterations, showing that optimizer effort is itself a confound. The largest covariance source-jet cache estimate is 531.25 MiB (557.1 MB), before runtime overhead; the executed 100-particle case is a screen rather than full trajectory certification. No GPU was available; Torch ran on CPU with a CPU-compatible dynamics path.', '',
      '## Part 7: Best Architecture','',
      '**Select calibrated single-carrier force-based holography with cached incremental feedback as the first hardware baseline.** It qualifies two particles with 256 transducers, the smallest source count qualified in the executed main-benchmark tests. That 76-s trajectory gives 0.08538-mm maximum target error, 0.00797-mm held drift, 28.96-µN/m minimum moving stiffness and 29.72-ms maximum controller computation. N=128 passes the initial stationary screen but fails during a 1-mm move when stiffness falls below threshold; this is a configuration failure, not a universal minimum-source theorem.', '',
      'At N=512, coherent support-residual control is the more economical tested single-carrier configuration by modeled energy and held drift, with less stiffness margin. Two-tone control also qualifies two particles after an identical-command timing replay; it uses more energy than these single-carrier examples and does not double the count. Temporal covariance remains a follow-up hypothesis for stiffness per watt and many-particle stationary feasibility. No tested architecture establishes a qualified high-density count advantage.', '',
      'The recommendation is provisional about capacity: only complete strict runs may establish a count. Conventional IB here is an explicit phase-holography baseline with gravity-induced positioning error. A simple focus-offset compensation was also attempted and did not converge robustly; this is not proof against stronger conventional controllers. Static force control has amplitude/phase freedom, whereas IB is phase-only. Accordingly the study does not isolate every gain to the force objective alone.', '',
      '## Part 8: Physical Prototype','',
      'The smallest useful mechanism experiment is an opposed 32-channel array (two 4×4 faces), two matched foam particles, a measured single-carrier transfer matrix, and calibrated phase-step response. Use a separate slow benchmark, initially 5 Hz weakest-mode frequency; this is not the main 15-Hz result. Ordinary small emitters have far less emitting area than the fixed-area N=32 scaling model, so model the actual hardware and measure single-trap authority before attempting two particles. If insufficient, increase to two 8×8 faces (128 channels) while preserving the same experimental comparison.', '',
      '| Component category | Approximate planning allowance (USD) |\n| --- | --- |\n| 32 rated air-ultrasound emitters and mounts | 60–180 |\n| Synchronized phase/amplitude drivers, FPGA or timing MCU, protected supply | 150–500 |\n| Two calibrated high-frame-rate cameras and illumination | 250–800 |\n| Phase-capable ultrasound probe/interface or rented measurement setup | 250–1,000 |\n| Beads, alignment fixtures, catch tray and enclosure materials | 80–250 |\n| Total with reused oscilloscope/computer | 790–2,730 |', '',
      'These are planning estimates, not verified supplier quotes. Borrowing measurement equipment can lower outlay; higher-channel hardware and commercial probe calibration can exceed these ranges. The acoustic watt values in the simulation do not specify a supply rating.', '',
      'First measure pressure/phase transfers and complex emitter step response. Compare one static field, two/four jointly optimized fields and singleton sharing at identical measured power and target positions. Sweep sequence period while measuring micromotion and ringdown stiffness. Use high-rate particle observations for the relevant fast modes; 20-Hz control images alone cannot detect every unstable carrier-envelope motion.', '',
      'Success requires repeatable signed and simultaneous 0.25/1-mm moves, held and moving errors ≤0.10 mm, ≥60-s retention, measured stiffness meeting the declared prototype threshold, and measured end-to-end updates ≥20 Hz. A temporal improvement must survive measured source filtering and equal-power comparison. Failure includes any positive-average sequence with growing micromotion, drop, identity loss, stiffness shortfall or deadline misses. Begin with rigid perturbation objects and a particle catch tray; ultrasound exposure and audible by-products need calibrated measurement before human access.', '',
      '## Part 9: Remaining Research','',
      '1. Complete refined finite-motion tests of the sixteen-particle covariance candidate; its stationary margin is too small for a robustness claim.\n2. Optimize capture basins, displacement authority and scheduling jointly, across multiple independent seeded clouds and stronger restarts.\n3. Replace isolated-piston power with measured mutual radiation impedance and electrical losses.\n4. Add streaming, viscothermal/finite-size and absorptive forces, source variation and boundaries; validate them against the prototype.\n5. Benchmark worst-case online latency on a dedicated controller and include camera/transport/FPGA delay.\n6. Determine whether a PSD initializer plus stability-aware field refinement can improve complete passing count rather than only stiffness.\n7. Test disturbance recovery, source calibration uncertainty, camera occlusion/identity loss and longer retention. The current tests do not establish those capabilities.', '',
      'Model limits: homogeneous unbounded air, ideal baffled emitters, Rayleigh particle response, approximate Born/Foldy scattering forces, finite polynomial domains and a narrowband approximation for between-frame envelope derivatives. A surrogate-domain exit is a validation failure, not automatically a physical escape. Loading, independent optical color, enclosure scattering and hand interaction are outside the executed study. No result establishes physical display capacity.', '',
      '## Appendix A: Every new executed architecture/scaling attempt','',table(cases,cols),'',
      '## Appendix B: Numerical convergence inventory','',
      'All fixed-drive two-particle rows are retained in [convergence/results.json](results/convergence/results.json); nonlinear shooting arrays and sixteen-particle refinements are in the same directory. Degree-6 surrogate force error reached 2.97% of weight and was rejected. Degree 8 reduced the checked maximum to about 0.052% for the covariance design and 0.048% for the static design. Failed early harness versions and final runs are retained separately in the execution logs and motion data.', '',
      table(independent,[('Archived case',lambda r:r['case']),('Samples/slot',lambda r:r['samples_per_slot']),('Independent K µN/m',lambda r:number(r['independent_min_stiffness_N_m'],1e6)),('ρ',lambda r:number(r['system_floquet_radius'])),('Gravity residual/weight',lambda r:number(r['gravity_residual_fraction'])),('Runtime s',lambda r:number(r.get('independent_calculation_s')))]),'',
      table(convergence,[('Test',lambda r:r['test']),('Harmonics',lambda r:r.get('harmonics','—')),('Samples',lambda r:r.get('time_samples','—')),('Refinement/method',lambda r:r.get('refinement',r.get('method',number(r.get('step_m'),1e6)))),('K µN/m',lambda r:number(r.get('weakest_stiffness_N_m'),1e6)),('ρ',lambda r:number(r.get('floquet_radius'))),('Runtime s',lambda r:number(r.get('runtime_s')))]),'',
      table(nonlinear,[('Case',lambda r:r['case']),('Steps/cycle',lambda r:r['steps_per_cycle']),('Shooting converged',lambda r:r['shooting_converged']),('Closure scaled norm',lambda r:number(r['closure_scaled_norm'])),('ρ nonlinear',lambda r:number(r['floquet_radius_nonlinear'],1,9)),('Max periodic error mm',lambda r:number(r['maximum_position_error_mm'],1,7)),('Runtime s',lambda r:number(r['runtime_s']))]),'',
      table(read(results/'convergence/source_delay_bins.json',[]),[('Delay bins',lambda r:r['delay_bins']),('Samples/cycle',lambda r:r['samples_per_cycle']),('Duration s',lambda r:number(r['duration_s'])),('Maximum error mm',lambda r:number(r['maximum_error_mm'],1,7)),('Difference from first mm',lambda r:number(r['maximum_position_difference_from_first_mm'],1,7))]),'',
      'The analytical standing-wave controls verify signs, gravity and stiffness. As an external context benchmark, HAT used 512 air emitters at 40 kHz, reporting fully dynamic manipulation of twelve particles and partial manipulation of twenty-five: [Marzo & Drinkwater, PNAS 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6320506/). Its millimetric particles, hardware, tolerances and geometry differ; this is not a calibration match. The Rayleigh multiple-scattering assumptions should be checked against [Silva & Bruus, PRE 2014](https://doi.org/10.1103/PhysRevE.90.063007). These are published results, not new simulated achievements.', '']
    (ROOT/'phase3_report.md').write_text('\n'.join(lines),encoding='utf-8')
    write(results/'study_inventory.json',{'legacy_reproductions':len(legacy),'independent_recalculations':len(independent),
       'new_cases':len(cases),'motion_records':len(motion),'strict_runs':len(passing),'convergence_rows':len(convergence)+len(sixteen),'coupled_density_diagnostics':len(density),
       'strict_cases':[r['case'] for r in passing]})
    print('REPORT_AND_FIGURES_GENERATED',len(cases),len(motion),flush=True)


if __name__=='__main__':main()
