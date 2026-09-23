"""Download the four explicitly selected files from the official pinned public archive."""
from pathlib import Path
import requests,json,hashlib,datetime
R=Path(__file__).resolve().parents[1]
p=R/'raw/independent/DepMap24Q4'
a=json.loads((p/'article_manifest.json').read_text(encoding='utf8'))
wanted=['Model.csv','README.txt','OmicsProfiles.csv','OmicsExpressionProteinCodingGenesTPMLogp1.csv']
log=[]
for name in wanted:
 f=next(x for x in a['files'] if x['name']==name)
 out=p/name
 if not out.exists():
  with requests.get(f['download_url'],stream=True,timeout=(30,120)) as response:
   response.raise_for_status()
   with out.with_suffix(out.suffix+'.part').open('wb') as h:
    for chunk in response.iter_content(2**20):h.write(chunk)
  out.with_suffix(out.suffix+'.part').rename(out)
 b=out.read_bytes();md5=hashlib.md5(b).hexdigest();sha=hashlib.sha256(b).hexdigest()
 assert len(b)==f['size'] and md5==f['computed_md5'],name
 row=dict(name=name,url=f['download_url'],bytes=len(b),md5=md5,sha256=sha,md5_verified=True,access_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),release=a['doi'])
 log.append(row);print(json.dumps(row),flush=True)
 (p/'download_log.json').write_text(json.dumps(log,indent=2),encoding='utf8')
