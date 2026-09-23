from pathlib import Path
import requests,json,hashlib,concurrent.futures
from datetime import datetime,timezone
R=Path(__file__).resolve().parents[1];O=R/'raw/independent'
items={
 'CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz':'https://data.broadinstitute.org/ccle/CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz',
 'PPTP_resource.html':'https://gccri.uthscsa.edu/pptp/',
 'PPTP_resource_legacy.html':'https://gccri.uthscsa.edu/PPTP/',
 'CCLE2019_PMC6697103.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6697103/fullTextXML'}
def get(item):
 name,url=item;p=O/name
 if p.exists():return dict(file=str(p.relative_to(R)),url=url,status='previous_local',bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
 try:
  rr=requests.get(url,timeout=(15,55),stream=True);h=hashlib.sha256();size=0
  if rr.status_code!=200 and 'json' in rr.headers.get('Content-Type',''):p=p.with_suffix('.response.json')
  with p.open('wb') as f:
   for b in rr.iter_content(2**20):f.write(b);h.update(b);size+=len(b)
  return dict(file=str(p.relative_to(R)),url=url,resolved_url=rr.url,status=rr.status_code,bytes=size,content_type=rr.headers.get('Content-Type'),sha256=h.hexdigest(),utc=datetime.now(timezone.utc).isoformat())
 except Exception as e:return dict(file=name,url=url,error=type(e).__name__+': '+str(e))
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:rows=list(ex.map(get,items.items()))
(O/'rna_acquisition_log.json').write_text(json.dumps(rows,indent=2),encoding='utf-8');print(json.dumps(rows))
