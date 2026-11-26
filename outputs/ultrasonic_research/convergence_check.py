"""Reevaluate a fixed archived field with finer stencil/source quadrature."""
import json
import math
from pathlib import Path
import torch
import simulator as sim


def main():
    folder=Path("results_128")
    c=sim.Config(**json.loads((folder/"metadata.json").read_text())["config"])
    detail=json.loads((folder/"fields.json").read_text())["static_force_2"]
    q=torch.tensor(detail["equilibria_m"])
    u=torch.tensor(detail["command_real"])+1j*torch.tensor(detail["command_imag"])
    pos,_=sim.geometry(c)
    drives=[(0,u,torch.ones(1))]
    f0,k0=sim.Field(q,pos,[c.frequency],c).force_stiffness(drives)
    original=sim.disk_sources
    rows=[]
    for label,rings,angles,h in (("half_stencil",3,8,c.fd_step/2),
                                ("96_disk_patches",6,16,c.fd_step/2),
                                ("384_disk_patches",12,32,c.fd_step/2)):
        def refined(pos,radius):
            offsets=[]
            for ring in range(rings):
                rr=radius*math.sqrt((ring+.5)/rings)
                for j in range(angles):
                    angle=2*math.pi*(j+.5*(ring%2))/angles
                    offsets.append([rr*math.cos(angle),rr*math.sin(angle),0])
            return pos[:,None,:]+torch.tensor(offsets)[None,:,:]
        sim.disk_sources=refined
        cf=sim.Config(**{**c.__dict__,"fd_step":h})
        f,k=sim.Field(q,pos,[c.frequency],cf).force_stiffness(drives)
        rows.append({"variant":label,
            "force_change_fraction_of_weight":float(torch.linalg.vector_norm(f-f0,dim=1).max())/(c.mass*9.80665),
            "min_stiffness_relative_change":float(torch.linalg.eigvalsh(k).min()/torch.linalg.eigvalsh(k0).min()-1),
            "linearized_equilibrium_shift_um":float(torch.linalg.vector_norm(torch.linalg.solve(k,f[...,None]).squeeze(-1),dim=1).max())*1e6})
    sim.disk_sources=original
    Path("convergence_results.json").write_text(json.dumps(rows,indent=2))
    print(json.dumps(rows,indent=2))


if __name__=="__main__":main()
