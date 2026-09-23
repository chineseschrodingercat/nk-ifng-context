from pathlib import Path
from openpyxl import load_workbook
import pandas as pd, numpy as np, json, hashlib
from scipy.stats import spearmanr
R=Path(__file__).resolve().parents[1]
for p in ['analysis','provenance']:(R/p).mkdir(exist_ok=True)
src=R/'raw/Dufva2023/TableS4.xlsx'
sha=hashlib.sha256(src.read_bytes()).hexdigest()
wb=load_workbook(src,read_only=True,data_only=True)
panel=['IFNGR1','IFNGR2','JAK1','JAK2','STAT1','IRF1','B2M','HLA-E','TAP1','ICAM1','CD58','FAS','TNFRSF10A','TNFRSF10B','CASP8','CFLAR','BAX','BAK1','BCL2','BCL2L1','MCL1']
out=[];summ=[]
for ws in wb:
 if ws.title[0] not in 'ABDEFGHIJKLMN':continue
 ws.reset_dimensions();it=ws.iter_rows(values_only=True);header=next(it)
 if ws.title[0] in 'AB':
  df=pd.DataFrame(list(it),columns=header);df.to_csv(R/'provenance'/('Dufva_'+('screen_design' if ws.title[0]=='A' else 'samples')+'.csv'),index=False);continue
 rows=[]
 for i,row in enumerate(it,2):
  d=dict(zip(header,row))
  if not d.get('id'):continue
  def first(keys):
   for k in keys:
    if k in d and d[k] is not None:return d[k],k
   return np.nan,None
  ef,ek=first(['avgfc_ORcorrected','avgfc','neg|lfc','neg.lfc'])
  pf,pk=first(['pos|fdr','pos.fdr']);nf,nk=first(['neg|fdr','neg.fdr'])
  pp,_=first(['pos|p-value','pos.p.value']);npv,_=first(['neg|p-value','neg.p.value'])
  rows.append(dict(screen=ws.title,gene=d['id'],source_row=i,screen_type='GOF' if 'GOF' in ws.title else 'LOF',effector='KHYG1' if 'KHYG1' in ws.title else 'primary_NK',effect=ef,effect_field=ek,pos_p=pp,neg_p=npv,pos_fdr=pf,neg_fdr=nf))
 df=pd.DataFrame(rows)
 df['effect']=pd.to_numeric(df.effect,errors='coerce')
 df['within_screen_percentile']=df.effect.rank(pct=True,method='average')
 df['centered_percentile']=2*df.within_screen_percentile-1
 df['selected_panel']=df.gene.isin(panel)
 out.append(df)
 summ.append(dict(screen=ws.title,rows=len(df),unique_genes=df.gene.nunique(),numeric_effects=int(df.effect.notna().sum()),effect_field=';'.join(df.effect_field.dropna().unique()),duplicate_gene_rows=int(df.gene.duplicated().sum())))
 print(ws.title,len(df),flush=True)
wb.close();allr=pd.concat(out,ignore_index=True)
allr.to_csv(R/'analysis/Dufva_all_gene_results.csv.gz',index=False,compression='gzip')
p=allr[allr.selected_panel].copy();p.to_csv(R/'analysis/Dufva_21gene_panel.csv',index=False)
ambiguous=allr.gene.astype(str).str.match(r'^\d{1,2}-(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)$') | allr.duplicated(['screen','gene'],False)
allr[ambiguous].to_csv(R/'provenance/Dufva_ambiguous_gene_labels.csv',index=False)
resolved=allr[~ambiguous].copy()
resolved[~resolved.screen.str.contains('GOF')].pivot(index='gene',columns='screen',values='effect').to_csv(R/'analysis/Dufva_LOF_effect_matrix.csv.gz',compression='gzip')
pd.DataFrame(summ).to_csv(R/'analysis/Dufva_screen_summary.csv',index=False)
pair=resolved[resolved.screen.isin(['M. KMS11 GOF','N. KMS11 GOF (KHYG1)'])].pivot(index='gene',columns='screen',values='effect').dropna()
rho,pval=spearmanr(pair.iloc[:,0],pair.iloc[:,1])
panelpair=pair.loc[pair.index.intersection(panel)].reset_index();panelpair.to_csv(R/'analysis/Dufva_KMS11_effector_pair.csv',index=False)
summary=dict(source_file=str(src),source_sha256=sha,result_sheets=len(summ),gene_result_rows=len(allr),unique_gene_labels=allr.gene.nunique(),panel_rows=len(p),LOF_screens=int(allr[allr.screen_type=='LOF'].screen.nunique()),GOF_screens=int(allr[allr.screen_type=='GOF'].screen.nunique()),source_design_rows=len(pd.read_csv(R/'provenance/Dufva_screen_design.csv')),source_sample_rows=len(pd.read_csv(R/'provenance/Dufva_samples.csv')),ambiguous_label_rows_excluded_from_gene_joins=int(ambiguous.sum()),ambiguity_rule='duplicate or date-like gene labels retained in raw data; no invented repairs',KMS11_GOF_shared_genes=len(pair),KMS11_GOF_spearman=float(rho),KMS11_GOF_nominal_gene_rank_p=float(pval),inference='descriptive paired gene ranks; genes are dependent and no inferential p value reported; all selection endpoints distinct from acute IFNg effects')
(R/'analysis/Dufva_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf8')
print(json.dumps(summary,indent=2),flush=True)
print(panelpair.to_string(index=False),flush=True)
