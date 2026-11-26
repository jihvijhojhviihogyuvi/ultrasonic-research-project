"""Install verified official wheels into the task directory without pip temp files."""
import hashlib,json,re,urllib.request,zipfile
from pathlib import Path,PurePosixPath

root=Path(__file__).resolve().parent
dest=root/"deps";dest.mkdir(exist_ok=True)
wheels=root/"wheels";wheels.mkdir(exist_ok=True)
names=("scipy","matplotlib","cvxpy","numba","llvmlite","clarabel","scs","osqp",
       "contourpy","kiwisolver","pillow","fonttools","cycler","packaging",
       "pyparsing","python-dateutil","six","joblib","jinja2","markupsafe",
       "qdldl","highspy","sparsediffpy")
records=[]
for name in names:
    info=json.load(urllib.request.urlopen(f"https://pypi.org/pypi/{name}/json"))
    candidates=[]
    for f in info["urls"]:
        fn=f["filename"]
        if fn.endswith(".whl"):
            if "cp313-cp313-win_amd64" in fn:score=3
            elif "abi3-win_amd64" in fn and any(int(x)<=313 for x in re.findall(r"cp(\d+)",fn)):score=2
            elif fn.endswith("none-any.whl"):score=1
            else:continue
            candidates.append((score,f))
    if not candidates:raise RuntimeError(f"no compatible wheel for {name}: {info['info']['version']}")
    f=max(candidates,key=lambda x:x[0])[1]
    path=wheels/f["filename"]
    if not path.exists():
        with urllib.request.urlopen(f["url"]) as response,path.open("wb") as out:
            while data:=response.read(1024*1024):out.write(data)
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    if digest!=f["digests"]["sha256"]:raise RuntimeError(f"wheel digest mismatch: {name}")
    with zipfile.ZipFile(path) as z:
        for entry in z.infolist():
            rel=PurePosixPath(entry.filename)
            if rel.is_absolute() or ".." in rel.parts:raise RuntimeError("unsafe wheel path")
            parts=rel.parts
            if parts and parts[0].endswith(".data"):
                if len(parts)>2 and parts[1] in ("purelib","platlib"):parts=parts[2:]
                else:continue
            target=dest.joinpath(*parts)
            if entry.is_dir():target.mkdir(parents=True,exist_ok=True)
            else:target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(z.read(entry))
    records.append({"name":name,"version":info["info"]["version"],"wheel":f["filename"],"sha256":digest})
    (root/"dependency_wheels.json").write_text(json.dumps(records,indent=2))
    print(f"installed {name} {info['info']['version']}",flush=True)
print("TASK_DEPENDENCIES_INSTALLED")
