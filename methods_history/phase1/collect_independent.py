from pathlib import Path
import requests,json,hashlib,concurrent.futures
from datetime import datetime,timezone
R=Path(__file__).resolve().parents[1];O=R/'raw/independent';O.mkdir(exist_ok=True)
items={
 'CCLE_RPPA_20181003.csv':'https://data.broadinstitute.org/ccle/CCLE_RPPA_20181003.csv',
 'CCLE_RPPA_Ab_info_20181226.csv':'https://data.broadinstitute.org/ccle/CCLE_RPPA_Ab_info_20181226.csv',
 'CCLE_index.html':'https://data.broadinstitute.org/ccle/',
 'GSM219748_metadata.txt':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM219748&targ=self&form=text&view=quick',
 'PPTP_PMC4209898.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4209898/fullTextXML',
 'CHLA10_datasheet.pdf':'https://www.cccells.org/PDF_Files/Brain/CHLA-10%20Cell%20Line%20Data%20Sheet.pdf'}
def get(item):
 name,url=item;p=O/name
 if p.exists():return dict(file=str(p.relative_to(R)),url=url,status='previous_local',bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
 try:
  r=requests.get(url,timeout=45)
  if r.status_code!=200 and 'json' in r.headers.get('Content-Type',''):p=p.with_suffix('.response.json')
  p.write_bytes(r.content)
  return dict(file=str(p.relative_to(R)),url=url,resolved_url=r.url,status=r.status_code,bytes=len(r.content),content_type=r.headers.get('Content-Type'),sha256=hashlib.sha256(r.content).hexdigest(),utc=datetime.now(timezone.utc).isoformat())
 except Exception as e:return dict(file=name,url=url,error=type(e).__name__+': '+str(e))
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:rows=list(ex.map(get,items.items()))
(O/'acquisition_log.json').write_text(json.dumps(rows,indent=2),encoding='utf-8');print(json.dumps(rows))
