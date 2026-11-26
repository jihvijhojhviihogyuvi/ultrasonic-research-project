"""Fresh execution of archived Phase-2 source, with provenance and comparisons.

This is implementation reproduction, not independent physics validation. The
separate validator will independently recalculate selected forces and dynamics.
"""
import csv
import hashlib
import importlib.util
import json
import platform
import sys
import time
from pathlib import Path

ROOT=Path(__file__).resolve().parent
OLD=ROOT.parent/"ultrasonic_research"
if not OLD.exists():OLD=ROOT/'inputs'/'legacy'


def load(folder):
    path=OLD/folder/"simulator_snapshot.py"
    metadata=json.loads((path.parent/"metadata.json").read_text())
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    if digest!=metadata["source_sha256"]: raise RuntimeError(f"snapshot mismatch: {path}")
    spec=importlib.util.spec_from_file_location("legacy_"+folder,path)
    module=importlib.util.module_from_spec(spec)
    sys.modules[spec.name]=module;spec.loader.exec_module(module)
    return module,metadata,path


def main():
    out=ROOT/"results"/"reproduction";out.mkdir(parents=True,exist_ok=True)
    cases=[("results_128",a,2,True) for a in
           ("ib","static_force","multifrequency","time_single","covariance","support_residual")]
    cases += [("results_512_screen",a,8,False) for a in ("static_force","covariance")]
    cases += [("results_cycle_1ms","covariance",8,False)]
    cases += [("results_512_screen",a,16,False) for a in
              ("static_force","covariance","multifrequency","aperiodic")]
    rows=[];fields={};provenance=[]
    for folder,arch,m,probes in cases:
        module,meta,path=load(folder)
        c=module.Config(**meta["config"])
        old=next(r for r in json.loads((path.parent/"results.json").read_text())
                 if r["architecture"]==arch and r["requested_particles"]==m)
        start=time.perf_counter();fresh,detail=module.run_case(c,arch,m,probes)
        key=f"{folder}_{arch}_{m}"
        row={"case":key,"source_folder":folder,**fresh,
             "old_stiffness_N_m":old["min_stiffness_N_m"],
             "stiffness_relative_difference":fresh["min_stiffness_N_m"]/old["min_stiffness_N_m"]-1,
             "old_max_error_mm":old["max_error_mm"],
             "old_floquet_radius":old["floquet_radius"],
             "fresh_execution_s":time.perf_counter()-start}
        if old["floquet_radius"] is not None and fresh["floquet_radius"] is not None:
            row["floquet_absolute_difference"]=fresh["floquet_radius"]-old["floquet_radius"]
        rows.append(row);fields[key]=detail
        provenance.append({"case":key,"config":meta["config"],
                           "source_sha256":meta["source_sha256"],"python":platform.python_version(),
                           "torch":module.torch.__version__,"numpy":module.np.__version__})
        print(json.dumps(row),flush=True)
        (out/"results.json").write_text(json.dumps(rows,indent=2,allow_nan=False))
        (out/"fields.json").write_text(json.dumps(fields,indent=2,allow_nan=False))
        (out/"provenance.json").write_text(json.dumps(provenance,indent=2))
    columns=list(dict.fromkeys(k for r in rows for k in r))
    with (out/"results.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=columns);w.writeheader();w.writerows(rows)
    print(f"REPRODUCTION_EXECUTED: {len(rows)} fresh cases",flush=True)


if __name__=="__main__":main()
