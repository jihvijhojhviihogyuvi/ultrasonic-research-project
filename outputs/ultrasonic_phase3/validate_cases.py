"""Run full signed trajectories on saved designs, retaining partial failures."""
import argparse,time
from pathlib import Path
import numpy as np
from experiments import ROOT,load_design,write
from validation import strict_validate


def main():
    p=argparse.ArgumentParser();p.add_argument('cases',nargs='+')
    p.add_argument('--samples',type=int,default=64);p.add_argument('--bins',type=int,default=8)
    p.add_argument('--degree',type=int,default=8);p.add_argument('--tag',default='final')
    p.add_argument('--noise',type=float,default=1e-5);p.add_argument('--wind',type=float,default=0.)
    p.add_argument('--harmonics',type=int)
    a=p.parse_args()
    for name in a.cases:
        folder=ROOT/'results'/'experiments'/name;out=ROOT/'results'/'motion'/f'{name}_{a.tag}'
        out.mkdir(parents=True,exist_ok=True);print('VALIDATING',name,flush=True)
        try:
            changes=None if a.harmonics is None else {'harmonics':a.harmonics}
            d,meta=load_design(folder,changes);q=np.load(folder/'design.npz')['equilibrium']
            if meta.get('equilibrium_error_mm',0)>d.system.c.accuracy*1000:
                write(out/'result.json',{'case':name,'status':'initial_position_precheck_failure','independent_particles':0,
                    'equilibrium_error_mm':meta['equilibrium_error_mm'],'motion_contract_pass':False})
                continue
            if meta.get('floquet_radius',0)>=1:
                write(out/'result.json',{'case':name,'status':'periodic_instability_precheck_failure','independent_particles':0,
                    'floquet_radius':meta['floquet_radius'],'motion_contract_pass':False})
                continue
            result,raw=strict_validate(d,q,degree=a.degree,time_samples=a.samples,bins=a.bins,
                                      noise_m=a.noise,wind=np.array([a.wind,0.,0.]))
            write(out/'result.json',dict(case=name,harmonics_used=d.system.c.harmonics,**result))
            write(out/'frames.json',raw.pop('frames',[]))
            np.savez_compressed(out/'raw.npz',**raw)
            print(name,result['status'],'duration',result.get('duration_s'),'error',result.get('maximum_position_error_mm'),flush=True)
        except Exception as e:
            import traceback
            write(out/'result.json',{'case':name,'status':'execution_error','exception':str(e),'traceback':traceback.format_exc()})
            traceback.print_exc()
    print('MOTION_VALIDATIONS_EXECUTED',flush=True)


if __name__=='__main__':main()
