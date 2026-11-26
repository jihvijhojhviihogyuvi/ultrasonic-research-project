"""Measure actual emitted-drive controller on saved camera inputs.

Replay checks every control coordinate against the original trajectory. It does
not rerun or fabricate particle motion. Original oracle-inclusive deadline
failures are preserved and linked as the physical trajectory evidence.
"""
import argparse,json,time
import numpy as np
from experiments import ROOT,load_design,write
from dynamics import Surrogate
from validation import IncrementalController


def main():
    p=argparse.ArgumentParser();p.add_argument('motion_run');a=p.parse_args()
    folder=ROOT/'results'/'motion'/a.motion_run
    result=json.loads((folder/'result.json').read_text());frames=json.loads((folder/'frames.json').read_text())
    data=np.load(folder/'raw.npz')
    changes={'harmonics':result['harmonics_used']} if 'harmonics_used' in result else None
    design,_=load_design(ROOT/'results'/'experiments'/result['case'],changes)
    sur=Surrogate(design.system,degree=result['surrogate']['degree'],bins=result['delay_bins'])
    controller=IncrementalController(design,sur);controller.update(sur.origin,np.zeros_like(sur.origin));controller.z[:]=0.
    assert np.allclose(controller.raw,data['raw_base'],atol=1e-10)
    assert np.allclose(controller.basis,data['controller_basis'],atol=1e-7,rtol=1e-7)
    initial=np.load(ROOT/'results'/'experiments'/result['case']/'design.npz')['equilibrium']
    previous_meas=initial.copy();previous_target=sur.origin.copy();velocity=np.zeros_like(initial);bias=np.zeros_like(initial)
    duration=1/result['update_hz'];delay=result['camera_delay_s'];timings=[];errors=[]
    for i,frame in enumerate(frames):
        start=time.perf_counter();meas=data['measurements_m'][i];target=np.array(frame['targets_m'])
        velocity=.7*velocity+.3*(meas-previous_meas)/duration
        estimate=meas+delay*velocity;previous_meas=meas
        bias=np.clip(bias+.01*(previous_target-estimate),-.00005,.00005)
        desired_velocity=(target-previous_target)/duration
        command,power,_=controller.update(target+bias,desired_velocity)
        packet=command.tobytes()  # materialize the actual drive buffer
        timings.append(time.perf_counter()-start)
        errors.append(float(np.max(np.abs(controller.z-data['z'][i]))))
        previous_target=target.copy()
    same=max(errors,default=0)<1e-7
    meets=same and max(timings,default=float('inf'))<=duration
    out=ROOT/'results'/'controller_replay'/a.motion_run;out.mkdir(parents=True,exist_ok=True)
    benchmark={'case':result['case'],'source_motion_run':a.motion_run,'frames_replayed':len(frames),
       'same_control_trajectory':same,'maximum_control_coordinate_difference':max(errors,default=0),
       'controller_mean_s':float(np.mean(timings)),'controller_max_s':max(timings),'controller_p99_s':float(np.quantile(timings,.99)),
       'compute_deadline_pass':meets,'update_hz':result['update_hz'],
       'scope':'camera processing, correction and spectral packet; offline simulation/verification oracles excluded'}
    write(out/'benchmark.json',benchmark);np.save(out/'timings_s.npy',np.array(timings))
    if result.get('motion_contract_pass'):
        combined=dict(result);combined.update(benchmark)
        combined['status']='strict_pass' if meets else 'compute_rate_failure'
        combined['independent_particles']=design.system.m if meets else 0
        combined['real_time_compute_feasible']=meets
        combined['evidence_scope']='same saved nonlinear motion plus separately measured identical controller replay'
        destination=ROOT/'results'/'motion'/f'{a.motion_run}_controller_replay'
        destination.mkdir(parents=True,exist_ok=True);write(destination/'result.json',combined)
    print(json.dumps(benchmark,indent=2),flush=True)
    print('CONTROLLER_REPLAY_EXECUTED',flush=True)


if __name__=='__main__':main()
