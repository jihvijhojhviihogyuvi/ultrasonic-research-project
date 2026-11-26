import sys,json,time,traceback
from pathlib import Path
sys.path.insert(0,str(Path('outputs/ultrasonic_phase3').resolve()))
from physics import *
from optimization import *
from validation import *
out=Path('outputs/ultrasonic_phase3/results/pilot');out.mkdir(parents=True,exist_ok=True)
c=Config(n=512,harmonics=0,time_samples=32)
s=AcousticSystem(c,targets(c,2));d=Design(s,'static',1)
print('OPTIMIZING',flush=True);optimize(d,250)
coef=d.coefficients().detach();q,f,k,ok,cur=equilibrium(s,coef)
print('EQUILIBRIUM',ok,'error_mm',float(torch.linalg.vector_norm(q-s.q,dim=1).max())*1000,
      'stiffness_uN_m',float(torch.linalg.eigvalsh((k+k.T)/2).min())*1e6,flush=True)
torch.save({'config':c,'raw':d.raw.detach(),'gauge':d.gauge.detach(),'timing':d.timing.detach(),
 'q':s.q.cpu(),'equilibrium':q.cpu(),'history':d.history},out/'design.pt')
print('VALIDATING',flush=True)
try:
 result,raw=strict_validate(d,q.cpu().numpy(),degree=6,time_samples=32)
 (out/'validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
 (out/'frames.json').write_text(json.dumps(raw.pop('frames',[])),encoding='utf-8')
 np.savez_compressed(out/'raw.npz',**raw)
 print(json.dumps(result,indent=2),flush=True)
except Exception:
 traceback.print_exc();raise
