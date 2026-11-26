"""Independent NumPy pressure-FD / SciPy-expm audit of the archived model.

No Phase-2 field, force, Hessian, envelope, or Floquet function is imported.
This independently checks the same simplified physical assumptions; it does not
validate missing scattering, beats, retarded propagation, or real hardware.
"""
import argparse,json,math,time
from pathlib import Path
import numpy as np
from scipy.linalg import expm

ROOT=Path(__file__).resolve().parent


def array_positions(c):
    side=round(math.sqrt(c["n"]/2));xy=(np.arange(side)-(side-1)/2)*c["pitch"]
    pos=np.array([(x,y,z) for z in (-c["gap"]/2,c["gap"]/2) for x in xy for y in xy])
    if c.get("aperiodic",False):
        jitter=np.random.default_rng(c["seed"]).uniform(-.03,.03,(c["n"],2))*c["pitch"]
        pos[:,:2]=np.clip(pos[:,:2]+jitter,xy.min(),xy.max())
    return pos


def transfers(q,c,f):
    pos=array_positions(c)
    disks=np.array([[c["emitter_radius"]*math.sqrt((j+.5)/3)*math.cos(2*math.pi*(a+.5*(j%2))/8),
                     c["emitter_radius"]*math.sqrt((j+.5)/3)*math.sin(2*math.pi*(a+.5*(j%2))/8),0]
                    for j in range(3) for a in range(8)])
    source=pos[:,None,:]+disks[None,:,:]
    k=2*math.pi*f/c["sound_speed"]
    gain=c["source_gain"]*(f/c["frequency"])/(1+1j*2*c["tone_q"]*(f-c["frequency"])/c["frequency"])
    def pressure(x):
        result=[]
        for start in range(0,len(x),32):
            r=np.linalg.norm(x[start:start+32,None,None,:]-source[None,:,:,:],axis=-1)
            result.append(gain*np.mean(np.exp(1j*k*r)/r,axis=2))
        return np.concatenate(result)
    h=pressure(q);derivatives=[];dv=2e-6
    for axis in range(3):
        offset=np.zeros(3);offset[axis]=dv
        derivatives.append((pressure(q+offset)-pressure(q-offset))/(2*dv))
    return h,np.stack(derivatives,axis=-1)/(1j*2*math.pi*f*c["air_density"])


def offsets(h):
    pts=[np.zeros(3)];ax={};mixed={}
    for a in range(3):
        ax[a]=[]
        for s in (1,-1):
            z=np.zeros(3);z[a]=s*h;ax[a].append(len(pts));pts.append(z)
    for a in range(3):
        for b in range(a+1,3):
            mixed[a,b]=[]
            for s,t in ((1,1),(1,-1),(-1,1),(-1,-1)):
                z=np.zeros(3);z[a]=s*h;z[b]=t*h;mixed[a,b].append(len(pts));pts.append(z)
    return np.array(pts),ax,mixed


def field_matrices(q,c,f):
    off,ax,mix=offsets(c["fd_step"])
    points=(q[:,None,:]+off).reshape(-1,3)
    h,v=transfers(points,c,f)
    return h.reshape(len(q),19,-1),v.reshape(len(q),19,-1,3),ax,mix


def moments(mats,drives,c):
    h,v,ax,mix=mats
    pressure=np.einsum("mjn,sn->smj",h,drives)
    velocity=np.einsum("mjnl,sn->smjl",v,drives)
    vol=4*math.pi*c["particle_radius"]**3/3
    f0=1-c["compressibility_ratio"]
    f2=2*(c["particle_density"]-c["air_density"])/(2*c["particle_density"]+c["air_density"])
    u=vol*(f0*np.abs(pressure)**2/(4*c["air_density"]*c["sound_speed"]**2)
           -3*f2*c["air_density"]*np.sum(np.abs(velocity)**2,axis=-1)/8)
    delta=c["fd_step"];force=np.zeros((*u.shape[:2],3));stiff=np.zeros((*u.shape[:2],3,3))
    for a in range(3):
        plus,minus=ax[a]
        force[:,:,a]=-(u[:,:,plus]-u[:,:,minus])/(2*delta)
        stiff[:,:,a,a]=(u[:,:,plus]-2*u[:,:,0]+u[:,:,minus])/delta**2
    for (a,b),idx in mix.items():
        stiff[:,:,a,b]=stiff[:,:,b,a]=(u[:,:,idx[0]]-u[:,:,idx[1]]-u[:,:,idx[2]]+u[:,:,idx[3]])/(4*delta**2)
    return force,stiff


def periodic_envelope(command,c,samples):
    n=len(command)
    if n==1:return command
    dt=c["cycle_s"]/n;a=np.exp(-dt/c["tau_s"])
    # Solve the scalar periodic boundary condition after traversing all slots.
    end=np.zeros(command.shape[1],complex)
    for x in command:end=a*end+(1-a)*x
    start=end/(1-a**n)
    result=[]
    for x in command:
        ts=(np.arange(samples)+.5)*dt/samples
        result.extend(x+(start-x)*np.exp(-ts[:,None]/c["tau_s"]))
        start=a*start+(1-a)*x
    return np.array(result)


def multipliers(stiff,c):
    mass=4*math.pi*c["particle_radius"]**3*c["particle_density"]/3
    drag=6*math.pi*c["viscosity"]*c["particle_radius"]
    radii=[];frequencies=[]
    for i in range(stiff.shape[1]):
        phi=np.eye(6)
        for ks in stiff[:,i]:
            aa=np.block([[np.zeros((3,3)),np.eye(3)],[-ks/mass,-drag/mass*np.eye(3)]])
            phi=expm(aa*c["cycle_s"]/len(stiff))@phi
        radii.append(float(np.max(np.abs(np.linalg.eigvals(phi)))))
        frequencies.append((np.sqrt(np.maximum(np.linalg.eigvalsh(stiff[:,i].mean(0)),0)/mass)/(2*np.pi)).tolist())
    return radii,frequencies


def main():
    p=argparse.ArgumentParser();p.add_argument("--samples",default="6,12,24");a=p.parse_args()
    folder=ROOT/"results/reproduction"
    records=json.loads((folder/"results.json").read_text());details=json.loads((folder/"fields.json").read_text())
    metadata={r["case"]:r for r in json.loads((folder/"provenance.json").read_text())}
    rows=[]
    for row in records:
        key=row["case"];c=metadata[key]["config"];detail=details[key]
        c={**c,"aperiodic":row["architecture"]=="aperiodic"}
        q=np.array(detail["equilibria_m"]);command=np.array(detail["command_real"])+1j*np.array(detail["command_imag"])
        start=time.perf_counter()
        if row["architecture"]=="multifrequency":
            combined=[]
            for j,f in enumerate((c["frequency"]-400,c["frequency"]+400)):
                combined.append(moments(field_matrices(q,c,f),command[j:j+1],c))
            force=sum(x[0][0] for x in combined);stiff=sum(x[1][0] for x in combined)
            radii=None;freq=None;sample_count=0
            attempts=[(0,force,stiff,radii,freq)]
        else:
            mats=field_matrices(q,c,c["frequency"]);attempts=[]
            samples=[1] if len(command)==1 else list(map(int,a.samples.split(",")))
            for s in samples:
                actual=periodic_envelope(command,c,s) if len(command)>1 else command
                ff,kk=moments(mats,actual,c);radii,freq=multipliers(kk,c)
                attempts.append((s,ff.mean(0),kk.mean(0),radii,freq))
        mass=4*math.pi*c["particle_radius"]**3*c["particle_density"]/3
        for sample_count,force,stiff,radii,freq in attempts:
            force=force.copy();force[:,2]-=mass*9.80665
            minimum=float(np.linalg.eigvalsh(stiff).min())
            result={"case":key,"samples_per_slot":sample_count,
                    "independent_min_stiffness_N_m":minimum,
                    "relative_stiffness_difference_from_legacy":minimum/row["min_stiffness_N_m"]-1,
                    "gravity_residual_fraction":float(np.linalg.norm(force,axis=1).max()/(mass*9.80665)),
                    "all_particle_floquet_radii":radii,"system_floquet_radius":max(radii) if radii else None,
                    "mean_stiffness_natural_frequencies_Hz":freq,"independent_calculation_s":time.perf_counter()-start}
            rows.append(result);print(json.dumps(result),flush=True)
        (folder/"independent_physics.json").write_text(json.dumps(rows,indent=2,allow_nan=False))
    print(f"INDEPENDENT_MODEL_AUDIT_EXECUTED: {len(rows)} calculations")


if __name__=="__main__":main()
