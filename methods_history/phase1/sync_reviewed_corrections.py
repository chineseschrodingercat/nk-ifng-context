"""Guarded synchronization of this run's final pre-delivery review edits only."""
from pathlib import Path
import pandas as pd,hashlib,shutil,json
R=Path(__file__).resolve().parents[1];D=Path(r'phase1_runs\20260922_185850')
assert D.parent.parent.resolve()==Path(r'C:\N01_NK_IFNg').resolve()
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
old=pd.read_csv(D/'artifact_manifest.csv').set_index('path');new=pd.read_csv(R/'artifact_manifest.csv').set_index('path')
assert set(old.index)<=set(new.index),'No deletions allowed'
for path,row in old.iterrows():assert sha(D/path)==row.sha256,('Unexpected destination change',path)
changed=[]
for path,row in new.iterrows():
 assert sha(R/path)==row.sha256,('Source mismatch',path)
 if path not in old.index or old.loc[path,'sha256']!=row.sha256:
  target=D/path
  if path not in old.index:assert not target.exists(),('Unexpected file',path)
  target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(R/path,target);changed.append(path)
shutil.copyfile(R/'artifact_manifest.csv',D/'artifact_manifest.csv')
for path,row in new.iterrows():assert sha(D/path)==row.sha256,path
print(json.dumps({'destination':str(D),'verified_files':len(new),'changed_or_added':changed,'bytes':int(new['bytes'].sum())},indent=2))
