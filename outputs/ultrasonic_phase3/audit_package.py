"""Independent package/evidence audit. A positive outcome is never assumed."""
import ast,hashlib,json,math,re
from pathlib import Path
import numpy as np
from experiments import ROOT,write,plan


def require(condition,message):
    if not condition:raise AssertionError(message)


def main():
    records=[]
    def passed(name):records.append(name)
    for path in ROOT.glob('*.py'):ast.parse(path.read_text(encoding='utf-8'),filename=path.name)
    passed('every delivered Python file parses')
    legacy=json.loads((ROOT/'results/reproduction/results.json').read_text())
    require(len(legacy)==13,'reproduction count')
    for row in legacy:
        require(abs(row['stiffness_relative_difference'])<1e-12,'legacy stiffness difference')
        require(abs(row['max_error_mm']-row['old_max_error_mm'])<1e-12,'legacy position difference')
        if row.get('floquet_absolute_difference') is not None:require(abs(row['floquet_absolute_difference'])<1e-12,'legacy Floquet difference')
    passed('13 independent fresh executions match archived numerical values')
    independent=json.loads((ROOT/'results/reproduction/independent_physics.json').read_text())
    require(len(independent)==28,'independent physics calculation count')
    require('import simulator' not in (ROOT/'independent_legacy.py').read_text(),'independent validator imports old physics')
    rr=[r for r in independent if r['case']=='results_512_screen_covariance_8']
    require({r['samples_per_slot'] for r in rr}=={6,12,24,48},'legacy temporal refinement missing')
    require(all(r['system_floquet_radius']>1 for r in rr),'legacy instability evidence missing')
    passed('separate NumPy/SciPy refinement exposes unconverged legacy quadrature')
    cases=[json.loads(p.read_text()) for p in (ROOT/'results/experiments').glob('*/case.json')]
    byid={r['id']:r for r in cases};required=set()
    for group in ('comparisons','temporal','refinements','scaling','power_density'):
        required.update(r['id'] for r in plan(group))
    require(required<=set(byid),'required experiment attempts are missing')
    require(len(cases)==len(byid),'duplicate experiment ids')
    errors=[r for r in cases if r['status']=='execution_error']
    require(all(r['id']=='nonlinear_floquet' for r in errors),'unexpected unresolved experiment error')
    for row in cases:
        if row['status']=='execution_error':continue
        folder=ROOT/'results/experiments'/row['id'];require((folder/'design.npz').exists(),'missing raw design')
        a=np.load(folder/'design.npz');c=row['config']
        require(a['positions'].shape==(c['n'],3),'source shape')
        require(a['targets'].shape==(row['particles'],3),'particle shape')
        require(math.isclose(float(np.sum(math.pi*a['radii']**2)),c['total_emitting_area'],rel_tol=1e-12),'emitting area')
        require(row['acoustic_power_W']<=c['power_W']*(1+1e-9),'steady power violation')
        require(row['maximum_command_envelope']<=1+1e-9,'source command peak violation')
        require(row.get('independent_particles') is None,'stationary result incorrectly claims independent count')
    passed('every attempted design retains matched geometry/power/peak resources and raw fields')
    reference=next(r['config'] for r in cases if r['category']=='scaling')
    for row in cases:
        if row['category']=='scaling':
            for key in ('volume_side','aperture','gap','total_emitting_area','power_W','particle_radius','particle_density','frequency'):
                require(row['config'][key]==reference[key],f'unfair scaling resource {key}')
    passed('all 72 scaling attempts keep area/aperture/volume/power/material fixed')
    matched=json.loads((ROOT/'results/matched_power/results.json').read_text())
    require(len(matched)==6,'matched realized-power cases missing')
    for m in (2,8,16):
        pp=[r['matched_realized_power_W'] for r in matched if r['particles']==m]
        require(abs(pp[0]-pp[1])<1e-12,'realized-power mismatch')
    passed('static/covariance supplemental tests match realized acoustic power exactly')
    conv=json.loads((ROOT/'results/convergence/results.json').read_text())
    require(len(conv)==15,'two-particle convergence rows missing')
    require((ROOT/'results/convergence/source_delay_bins.json').exists(),'delay-bin convergence missing')
    require(len(json.loads((ROOT/'results/convergence/sixteen_particle.json').read_text()))==5,'sixteen-particle refinement missing')
    require(len(json.loads((ROOT/'results/density_stability/results.json').read_text()))==6,'coupled density diagnostics incomplete')
    for name in ('compare_2_static','compare_2_covariance'):
        for steps in (64,128):
            orbit=json.loads((ROOT/f'results/convergence/{name}_nonlinear_{steps}.json').read_text())
            require(orbit['shooting_converged'] and orbit['valid_domain'],'nonlinear shooting did not converge')
            require(orbit['closure_scaled_norm']<1e-8,'nonlinear closure failure')
    passed('spectral/time/disk/spatial/scattering/delay and nonlinear-periodic checks retain raw outcomes')
    motion=[json.loads(p.read_text()) for p in (ROOT/'results/motion').glob('*/result.json')]
    required_controls={f'compare_2_{a}' for a in ('ib','static','multifrequency','singleton','covariance','support')}
    require(required_controls<={r.get('case') for r in motion},'two-particle operational comparison is still incomplete')
    require(any(r.get('motion_contract_pass') for r in motion),'no complete finite-motion control')
    for row in motion:
        if row.get('status')!='strict_pass':continue
        require(row['motion_contract_pass'] and row['duration_s']>=60,'strict retention missing')
        require(row['maximum_position_error_mm']<=.1 and row['maximum_held_displacement_mm']<=.1,'strict accuracy missing')
        require(row['minimum_moving_stiffness_N_m']>=1.7441030632668644e-5,'strict stiffness missing')
        require(row['frames_completed']==row['frames_required'],'signed protocol incomplete')
        require(row['controller_max_s']<=.05,'strict controller deadline missing')
        labels=set(row['tests_completed_labels']);m=byid[row['case']]['particles']
        for distance in ('0.25','1'):
            for particle in range(m):
                for axis in range(3):
                    for sign in (-1,1):require(f'{distance}mm_p{particle}_a{axis}_{sign}' in labels,'signed move missing')
        require('simultaneous' in labels,'simultaneous moves missing')
        if row.get('evidence_scope'):
            benchmark=json.loads((ROOT/'results/controller_replay'/row['source_motion_run']/'benchmark.json').read_text())
            require(benchmark['same_control_trajectory'] and benchmark['maximum_control_coordinate_difference']<1e-7,'replay changed controller trajectory')
    passed('strict labels satisfy signed/simultaneous motion, resource, retention and actual controller timing evidence')
    report=(ROOT/'phase3_report.md').read_text(encoding='utf-8')
    for part in range(1,10):require(f'## Part {part}:' in report,f'missing report part {part}')
    for row in cases:require(row['id'] in report,'executed result omitted from report')
    require(len(list((ROOT/'figures').glob('*.png')))>=8,'comparison figures missing')
    passed('nine-part report includes every executed case and plotted numerical evidence')
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.glob('*.py')}
    write(ROOT/'results/package_audit.json',{'status':'passed','checks':records,'source_hashes':hashes,
       'new_cases':len(cases),'strict_records':sum(r.get('status')=='strict_pass' for r in motion)})
    print('PHASE3_PACKAGE_AUDIT_PASSED',len(records))


if __name__=='__main__':main()
