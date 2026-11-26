"""Reproducible sequential pipeline; model rejections are saved observations."""
import argparse,os,subprocess,sys
from pathlib import Path

ROOT=Path(__file__).resolve().parent


def main():
    p=argparse.ArgumentParser();p.add_argument('--skip-reproduction',action='store_true')
    p.add_argument('--skip-motion',action='store_true');a=p.parse_args()
    steps=[['tests.py'],['tests_dynamics.py']]
    if not a.skip_reproduction:steps += [['reproduce.py'],['independent_legacy.py']]
    steps += [['experiments.py','--group',g] for g in ('comparisons','temporal','refinements','scaling','power_density')]
    steps += [['convergence.py','--nonlinear'],['convergence_high.py'],['power_match.py'],['sdp_recovery.py'],['density_stability.py']]
    if not a.skip_motion:
        steps += [['validate_cases.py','compare_2_static','compare_2_covariance','compare_2_support','compare_2_ib','compare_2_singleton','--samples','64','--tag','pipeline'],
                  ['validate_cases.py','compare_2_multifrequency','--harmonics','2','--samples','64','--tag','pipeline'],
                  ['validate_cases.py','scale_N128_M2_static','scale_N256_M2_static','--samples','64','--tag','pipeline'],
                  ['transition_convergence.py']]
    steps += [['build_report.py'],['audit_package.py']]
    env=os.environ.copy();env.setdefault('OPENBLAS_NUM_THREADS','1')
    for step in steps:
        print('Running',' '.join(step),flush=True)
        subprocess.run([sys.executable,'-B',str(ROOT/step[0]),*step[1:]],cwd=ROOT,env=env,check=True)
    print('PHASE3_PIPELINE_COMPLETE')


if __name__=='__main__':main()
