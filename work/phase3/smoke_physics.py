import sys
from pathlib import Path
sys.path.insert(0,str(Path('outputs/ultrasonic_phase3').resolve()))
from physics import *
c=Config(n=32,harmonics=2,time_samples=16)
s=AcousticSystem(c,targets(c,2))
u=torch.ones((2,c.n),device=s.device,dtype=torch.complex128)*.5
d=torch.tensor([.4,.6],device=s.device)
u,a=s.project(u,d)
print('power',float(s.power(a)))
f,k=s.primary_mean(a)
print('primary',f.shape,k.shape,torch.linalg.eigvalsh(k))
ff,dd,t=s.matrices(a)
print('coupled',ff.shape,dd.shape,float(linear_floquet(dd,c)[0]))
