from pathlib import Path
import concurrent.futures as cf, requests, json, csv, hashlib, datetime, re, subprocess, sys
R=Path(__file__).resolve().parents[1]
for p in ['raw/search','derived','provenance']: (R/p).mkdir(parents=True,exist_ok=True)
OLD=Path('closure_runs/20260923_111617/derived/literature_metadata_final.json')
old=json.loads(OLD.read_text(encoding='utf-8'))
if isinstance(old,dict): old=old.get('records',old.get('results',[]))
old_keys=set()
for d in old:
 for k in ['doi','pmid','pmcid','id']:
  if d.get(k):old_keys.add(str(d[k]).lower())
queries={
 'E01_direct':'(TITLE_ABS:"natural killer" OR TITLE_ABS:"NK cell") AND (TITLE_ABS:"interferon gamma" OR TITLE_ABS:"IFN-gamma" OR TITLE_ABS:"IFNγ" OR TITLE_ABS:"IFN-γ") AND (TITLE_ABS:susceptibility OR TITLE_ABS:resistance OR TITLE_ABS:sensitization OR TITLE_ABS:killing)',
 'E02_death':'(TITLE_ABS:"natural killer" OR TITLE_ABS:"NK cell") AND (TITLE_ABS:Fas OR TITLE_ABS:TRAIL OR TITLE_ABS:"death receptor") AND (TITLE_ABS:interferon OR TITLE_ABS:IFNG)',
 'E03_recent':'(TITLE_ABS:"natural killer" OR TITLE_ABS:"NK cell") AND (TITLE_ABS:"IFN-γ" OR TITLE_ABS:"IFNγ" OR TITLE_ABS:"interferon gamma") AND (TITLE_ABS:tumor OR TITLE_ABS:tumour OR TITLE_ABS:cancer) AND FIRST_PDATE:[2024-01-01 TO 2026-09-23]',
 'E04_data':'(TITLE_ABS:"natural killer" OR TITLE_ABS:"NK cell") AND (TITLE_ABS:screen OR TITLE_ABS:atlas OR TITLE_ABS:dataset OR TITLE_ABS:"single-cell") AND (TITLE_ABS:"IFN-γ" OR TITLE_ABS:"IFNγ" OR TITLE_ABS:"interferon gamma")',
 'E05_context':'(TITLE_ABS:interferon OR TITLE_ABS:IFNG) AND (TITLE_ABS:"NK cytotoxicity" OR TITLE_ABS:"NK-mediated" OR TITLE_ABS:"NK cell-mediated") AND (TITLE_ABS:pretreatment OR TITLE_ABS:pretreated OR TITLE_ABS:conditioning OR TITLE_ABS:blockade)',
 'E06_receptor':'(TITLE_ABS:"NKG2A" OR TITLE_ABS:"HLA-E") AND (TITLE_ABS:"IFN-γ" OR TITLE_ABS:"IFNγ" OR TITLE_ABS:"interferon gamma") AND (TITLE_ABS:tumor OR TITLE_ABS:cancer OR TITLE_ABS:carcinoma)',
 'E07_cancer_death':'(TITLE_ABS:"interferon gamma" OR TITLE_ABS:"IFN-γ" OR TITLE_ABS:"IFNγ") AND (TITLE_ABS:"death competence" OR TITLE_ABS:"death receptor" OR TITLE_ABS:"apoptotic priming") AND (TITLE_ABS:cancer OR TITLE_ABS:tumor)',
 'E08_hist':'(TITLE_ABS:"gamma interferon" OR TITLE_ABS:"IFN-gamma") AND (TITLE_ABS:"natural killer" OR TITLE_ABS:"NK-mediated") AND FIRST_PDATE:[1980-01-01 TO 2006-12-31]'
}
def run(k,q):
 rows=[];logs=[];cursor='*';seen=set()
 for page in range(1,11):
  p=R/'raw/search'/f'{k}_{page:02d}.json'; req={'query':q,'format':'json','resultType':'core','pageSize':500,'cursorMark':cursor}
  log=dict(query_id=k,query=q,page=page,retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
  try:
   if p.exists():b=p.read_bytes();url='cached '+str(p);status=200
   else:
    z=requests.get('https://www.ebi.ac.uk/europepmc/webservices/rest/search',params=req,timeout=45);url=z.url;status=z.status_code;z.raise_for_status();b=z.content;p.write_bytes(b)
   d=json.loads(b);rr=d.get('resultList',{}).get('result',[]);rows.extend(rr)
   log.update(status='success',http_status=status,url=url,records=len(rr),hit_count=d['hitCount'],sha256=hashlib.sha256(b).hexdigest(),file=str(p.relative_to(R)))
   nxt=d.get('nextCursorMark');done=not rr or not nxt or nxt==cursor or len(rows)>=d['hitCount'];log['query_complete']=done;logs.append(log)
   if done:break
   cursor=nxt
  except Exception as e:
   log.update(status='failed',error=str(e),query_complete=False);logs.append(log);break
 return k,rows,logs
allrows={};logs=[]
with cf.ThreadPoolExecutor(max_workers=3) as ex:
 for k,rr,ll in ex.map(lambda kv:run(*kv),queries.items()):
  for d in rr:
   key=(d.get('source'),d.get('id'));d.setdefault('matched_queries',[]).append(k)
   if key in allrows:allrows[key]['matched_queries']+=d['matched_queries']
   else:allrows[key]=d
  logs+=ll;print(json.dumps({'query':k,'records':len(rr),'reported_hits':ll[-1].get('hit_count'),'complete':ll[-1].get('query_complete')}),flush=True)
rows=list(allrows.values())
flat=[]
for d in rows:
 abstract=re.sub('<[^>]+>',' ',d.get('abstractText',''))
 new=not any(str(d.get(k,'')).lower() in old_keys for k in ['doi','pmid','pmcid','id'] if d.get(k))
 score=sum(w in (d.get('title','')+' '+abstract).lower() for w in ['pretreat','sensiti','cytotoxic','resistan','fas','trail','blockade','hla-e','source data'])
 flat.append(dict(source=d.get('source'),id=d.get('id'),pmcid=d.get('pmcid'),doi=d.get('doi'),year=d.get('pubYear'),title=d.get('title'),journal=d.get('journalInfo',{}).get('journal',{}).get('title'),new_relative_to_prior=new,priority_score=score,queries=';'.join(sorted(set(d['matched_queries']))),abstract=abstract))
flat.sort(key=lambda x:(-x['priority_score'],not x['new_relative_to_prior'],str(x['year'])),reverse=False)
(R/'derived/search_records.json').write_text(json.dumps(rows,indent=2),encoding='utf8')
(R/'provenance/search_log.json').write_text(json.dumps(logs,indent=2),encoding='utf8')
with (R/'derived/search_screening.csv').open('w',newline='',encoding='utf-8-sig') as f:
 w=csv.DictWriter(f,fieldnames=list(flat[0]));w.writeheader();w.writerows(flat)
summary=dict(unique_records=len(rows),new_relative_to_prior=sum(x['new_relative_to_prior'] for x in flat),prior_records=len(old),queries=len(queries),completed_queries=sum(any(l['query_id']==k and l.get('query_complete') for l in logs) for k in queries),screening_scope='automated metadata priority, not full-text inclusion')
(R/'derived/search_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf8');print(json.dumps(summary),flush=True)
print(json.dumps([{k:x[k] for k in ['id','pmcid','doi','year','title','priority_score']} for x in flat if x['new_relative_to_prior']][:35],ensure_ascii=False,indent=2),flush=True)
# One complementary PubMed query through the literature skill; preserve its response as an audit.
skill=Path('.codex/plugins/cache/openai-curated-remote/life-sciences-literature/0.1.5/skills/ncbi-entrez-skill/scripts/ncbi_entrez.py')
payload={'endpoint':'esearch','params':{'db':'pubmed','term':'("natural killer"[Title/Abstract] OR "NK cell"[Title/Abstract]) AND ("interferon gamma"[Title/Abstract] OR "IFN-gamma"[Title/Abstract]) AND (tumor[Title/Abstract] OR cancer[Title/Abstract])','retmode':'json','retmax':200},'max_items':10,'save_raw':True,'raw_output_path':str(R/'raw/search/pubmed_ids.json'),'timeout_sec':40}
cp=subprocess.run([sys.executable,str(skill)],input=json.dumps(payload),text=True,capture_output=True,timeout=65)
(R/'raw/search/pubmed_skill_response.json').write_text(cp.stdout or cp.stderr,encoding='utf8')
print('PubMed skill response saved; exit',cp.returncode,flush=True)
