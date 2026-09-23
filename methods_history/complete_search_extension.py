from pathlib import Path
import requests,json,hashlib,datetime,pandas as pd,re,concurrent.futures as cf
R=Path(__file__).resolve().parents[1]
logs=json.loads((R/'provenance/search_log.json').read_text())
failed={x['query_id']:x['query'] for x in logs if x['status']=='failed'}
def repair(k,q):
 cursor='*';rows=[];ll=[]
 for page in range(1,21):
  params=dict(query=q,format='json',resultType='lite',pageSize=200,cursorMark=cursor)
  path=R/'raw/search'/f'{k}_repair_lite{page}.json'
  try:
   if path.exists():b=path.read_bytes();url='cached';status=200
   else:
    z=requests.get('https://www.ebi.ac.uk/europepmc/webservices/rest/search',params=params,timeout=30);z.raise_for_status();b=z.content;url=z.url;status=z.status_code;path.write_bytes(b)
   d=json.loads(b);rr=d.get('resultList',{}).get('result',[]);rows+=rr;nxt=d.get('nextCursorMark');done=not rr or len(rows)>=d['hitCount'] or not nxt or nxt==cursor
   ll.append(dict(query_id=k,query=q,page=page,status='success',result_type='lite',file=str(path.relative_to(R)),sha256=hashlib.sha256(b).hexdigest(),url=url,records=len(rr),hit_count=d['hitCount'],query_complete=done))
   if done:break
   cursor=nxt
  except Exception as e:
   ll.append(dict(query_id=k,query=q,page=page,status='failed',error=str(e),query_complete=False));break
 return k,rows,ll
records=json.loads((R/'derived/search_records.json').read_text());di={(x.get('source'),x['id']):x for x in records};more=[]
with cf.ThreadPoolExecutor(max_workers=2) as ex:
 for k,rr,ll in ex.map(lambda kv:repair(*kv),failed.items()):
  more+=ll
  for x in rr:
   key=(x.get('source'),x['id'])
   if key in di:di[key].setdefault('matched_queries',[]).append(k)
   else:x['matched_queries']=[k];x['metadata_depth']='lite no abstract';di[key]=x
  print(json.dumps({'query':k,'fetched':len(rr),'complete':ll[-1].get('query_complete')}),flush=True)
(R/'provenance/search_repair_log.json').write_text(json.dumps(more,indent=2),encoding='utf8')
(R/'derived/search_records_complete.json').write_text(json.dumps(list(di.values()),indent=2),encoding='utf8')
old=json.loads(Path('closure_runs/20260923_111617/derived/literature_metadata_final.json').read_text(encoding='utf8'))
keys={str(x[k]).lower() for x in old for k in ['doi','id','pmid','pmcid'] if x.get(k)}
new=[x for x in di.values() if not any(str(x[k]).lower() in keys for k in ['doi','id','pmid','pmcid'] if x.get(k))]
summary=dict(unique_records_current=len(di),new_relative_to_prior=len(new),prior_metadata_records=len(old),queries=8,complete_queries=sum(any(y['query_id']==k and y.get('query_complete') for y in logs+more) for k in {y['query_id'] for y in logs}),fulltext_read_claim=False,PubMed_query='complementary 200 IDs retrieved; capped, not exhaustive')
(R/'derived/search_summary_final.json').write_text(json.dumps(summary,indent=2),encoding='utf8');print(json.dumps(summary),flush=True)
flat=[]
for x in di.values():
 text=(x.get('title','')+' '+re.sub('<[^>]+>',' ',x.get('abstractText',''))).lower()
 cancer=any(y in text for y in ['tumor','tumour','cancer','carcinoma','leukemia','leukaemia','myeloma','neuroblastoma','melanoma','lymphoma'])
 treated=any(y in text for y in ['pretreat','sensitiv','susceptib','resistan','sensiti','blockade'])
 infection=any(y in x.get('title','').lower() for y in ['virus','viral','infection','bacteria','salmonella','cmv','hiv','hantaan'])
 flat.append(dict(id=x['id'],source=x.get('source'),pmcid=x.get('pmcid'),doi=x.get('doi'),title=x.get('title'),year=x.get('pubYear'),abstract=re.sub('<[^>]+>',' ',x.get('abstractText','')),metadata_depth=x.get('metadata_depth','core'),cancer_function_priority=cancer and treated and not infection,queries=';'.join(sorted(set(x.get('matched_queries',[]))))))
d=pd.DataFrame(flat);d.to_csv(R/'derived/search_screening_final.csv',index=False)
d[d.cancer_function_priority].to_csv(R/'derived/cancer_function_priority.csv',index=False)
print('cancer function metadata priority',int(d.cancer_function_priority.sum()))
