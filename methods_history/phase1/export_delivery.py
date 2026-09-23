"""Additive, refuse-overwrite export to the authorized project folder."""
from pathlib import Path
import hashlib,json,shutil
import pandas as pd
R=Path(__file__).resolve().parents[1]
D=Path(r'phase1_runs\20260922_185850')
assert not D.exists(),'Destination exists; refusing to overwrite any delivered run'
assert D.parent.parent.resolve()==Path(r'C:\N01_NK_IFNg').resolve()
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
m=pd.read_csv(R/'artifact_manifest.csv');bad=[]
for _,r in m.iterrows():
 p=R/r['path']
 if p.stat().st_size!=r['bytes'] or sha(p)!=r.sha256:bad.append(r['path'])
assert not bad,bad
D.parent.mkdir(exist_ok=True)
shutil.copytree(R,D)
for _,r in m.iterrows():
 p=D/r['path']
 assert p.stat().st_size==r['bytes'] and sha(p)==r.sha256,r['path']
assert sha(R/'artifact_manifest.csv')==sha(D/'artifact_manifest.csv')
result=dict(destination=str(D),verified_manifest_files=len(m),verified_manifest=True,bytes=int(m['bytes'].sum()),all_previous_runs_untouched=True)
print(json.dumps(result,indent=2))
