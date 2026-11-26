"""Audit result provenance, power/peak constraints, and pass labels from artifacts."""
import ast
import hashlib
import json
import math
import re
from pathlib import Path
import numpy as np


def main():
    root=Path(__file__).resolve().parent
    for source in root.glob("*.py"):
        ast.parse(source.read_text(encoding="utf-8"),filename=str(source))
    report=(root/"research_report.md").read_text(encoding="utf-8")
    sections=re.findall(r"^## (\d+)\. (.+)$",report,re.M)
    assert [int(n) for n,_ in sections]==list(range(1,12)),sections
    count=0;local_pass=0;not_tested=0
    for folder in sorted(root.glob("results_*")):
        meta=json.loads((folder/"metadata.json").read_text())
        rows=json.loads((folder/"results.json").read_text())
        fields=json.loads((folder/"fields.json").read_text())
        assert hashlib.sha256((folder/"simulator_snapshot.py").read_bytes()).hexdigest()==meta["source_sha256"],folder
        config=meta["config"];n=config["n"]
        clouds={}
        keys=set()
        for row in rows:
            arch=row["architecture"];m=row["requested_particles"]
            key=f"{arch}_{m}"
            assert key not in keys;keys.add(key)
            detail=fields[key]
            q=np.array(detail["targets_m"])
            assert q.shape==(m,3) and np.isfinite(q).all()
            if m in clouds:assert np.array_equal(q,clouds[m]),(folder,key,"target mismatch")
            else:clouds[m]=q
            u=np.array(detail["command_real"])+1j*np.array(detail["command_imag"])
            assert u.ndim==2 and u.shape[1]==n and np.isfinite(u).all()
            power=float(np.sum(np.abs(u)**2))
            if arch!="multifrequency":power/=len(u)
            crest=np.sum(np.abs(u),axis=0) if arch=="multifrequency" else np.max(np.abs(u),axis=0)
            assert power<=n*config["power_fraction"]*(1+1e-10),(folder,key,"power")
            assert crest.max()<=1+1e-10,(folder,key,"crest")
            assert abs(power-row["drive_effort"])<1e-8,(folder,key,"effort")
            assert row["force_rank"]<=3*m
            assert 0<=row["stable_local_traps"]<=m
            if row["independence_status"]=="tier1_local_pass":
                local_pass+=1
                assert row["converged"] and row["stable_local_traps"]==m
                assert row["force_rank"]==3*m and row["failed_probes"]==0
                assert row["floquet_radius"]<1 and row["ripple_mm"]<config["accuracy"]*1000
                assert row["min_stiffness_N_m"]>=row["required_stiffness_N_m"]
            if row["independence_status"]=="not_tested":not_tested+=1
            count+=1
        assert keys==set(fields),(folder,"missing result/field")
    motion=json.loads((root/"tracking_static.json").read_text())
    assert motion["real_time_compute_feasible"] is False
    assert len(motion["frames"])==10
    scatter=json.loads((root/"scattering_512.json").read_text())
    assert len(scatter)==4 and all(r["converged"] for r in scatter)
    # Independently recalculate the report's quoted gains from saved executions.
    small=json.loads((root/"results_128/results.json").read_text())
    large=json.loads((root/"results_512_screen/results.json").read_text())
    def stiffness(rows,arch,m):
        return next(r["min_stiffness_N_m"] for r in rows if r["architecture"]==arch and r["requested_particles"]==m)
    gains=[100*(stiffness(small,"covariance",4)/stiffness(small,"static_force",4)-1),
           100*(stiffness(large,"covariance",16)/stiffness(large,"static_force",16)-1)]
    assert abs(gains[0]-43.2)<0.05 and abs(gains[1]-87.1)<0.05
    print(f"PACKAGE_AUDIT_PASSED: {count} cases; {local_pass} local passes; {not_tested} without finite independence probes; 11 sections; gain percentages {gains}")


if __name__=="__main__":main()
