"""Install SHA-256 verified wheels from the supplied execution manifest.

For the recorded Windows Python 3.13 runtime; use pip/requirements.txt elsewhere.
No pip temporary directories or global installation are needed.
"""
import argparse,hashlib,json,urllib.request,zipfile
from pathlib import Path,PurePosixPath


def main():
    p=argparse.ArgumentParser();p.add_argument('destination',type=Path)
    a=p.parse_args();a.destination.mkdir(parents=True,exist_ok=True)
    records=json.loads((Path(__file__).parent/'dependency_wheels.json').read_text())
    for r in records:
        info=json.load(urllib.request.urlopen(f"https://pypi.org/pypi/{r['name']}/{r['version']}/json"))
        wheel=next(f for f in info['urls'] if f['filename']==r['wheel'])
        data=urllib.request.urlopen(wheel['url']).read()
        assert hashlib.sha256(data).hexdigest()==r['sha256'],f"digest mismatch {r['name']}"
        import io
        with zipfile.ZipFile(io.BytesIO(data)) as z:
            for entry in z.infolist():
                rel=PurePosixPath(entry.filename)
                if rel.is_absolute() or '..' in rel.parts:raise ValueError('unsafe wheel entry')
                parts=rel.parts
                if parts and parts[0].endswith('.data'):
                    if len(parts)>2 and parts[1] in ('purelib','platlib'):parts=parts[2:]
                    else:continue
                target=a.destination.joinpath(*parts)
                if entry.is_dir():target.mkdir(parents=True,exist_ok=True)
                else:target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(z.read(entry))
        print('Installed',r['name'],r['version'],flush=True)


if __name__=='__main__':main()
