import sys,time
from pathlib import Path
sys.path.insert(0,str(Path('outputs/ultrasonic_phase3').resolve()))
from physics import *
from optimization import *
c=Config(n=32,harmonics=0,time_samples=16)
s=AcousticSystem(c,targets(c,2))
d=Design(s,'static',1)
print('optimization_s',optimize(d,80))
coef=d.coefficients();f,k=s.primary_mean(coef)
print('force_ratio',f.detach().numpy()/(c.mass*9.80665),'stiffness',torch.linalg.eigvalsh(k).detach().numpy())
q,f,k,ok,current=equilibrium(s,coef)
print('equilibrium',ok,q.cpu().numpy(),torch.linalg.eigvalsh((k+k.T)/2).cpu().numpy())
print('sdp',sdp_bound(s)['status'])
