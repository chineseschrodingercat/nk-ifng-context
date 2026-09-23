"""Exact-identity independent RNA join and training-only fixed-module validation."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from sklearn.preprocessing import StandardScaler,OneHotEncoder
from sklearn.linear_model import Ridge
from scipy.stats import spearmanr
R=Path(__file__).resolve().parents[1]
P=R/'raw/independent/DepMap24Q4'
modules={
 'mitochondrial':{'pro':['BAX','BAK1','BCL2L11','BID','BBC3','PMAIP1'],'anti':['BCL2','BCL2L1','MCL1','BCL2L2','BCL2A1']},
 'execution':{'pro':['APAF1','CASP3','CASP7','CASP8','CASP9'],'anti':['XIAP','CFLAR']}}
genes=[g for m in modules.values() for v in m.values() for g in v]
norm=lambda s:''.join(c for c in str(s).upper() if c.isalnum())
d=pd.read_csv(R/'analysis/model_feature_target_join.csv')
meta=pd.read_csv(P/'Model.csv');meta['norm_name']=meta.StrippedCellLineName.map(norm)
d['norm_name']=d.model.map(norm)
assert meta[meta.norm_name.isin(d.norm_name)].norm_name.is_unique
relevant_meta=meta[meta.norm_name.isin(d.norm_name)]
ident=d[['model','norm_name','family']].merge(relevant_meta[['norm_name','ModelID','PatientID','CellLineName','CCLEName','RRID']],on='norm_name',how='left',validate='one_to_one')
ident['identity_rule']='case/punctuation normalization of exact official cell-line name; no subline substitutions'
cols=pd.read_csv(P/'OmicsExpressionProteinCodingGenesTPMLogp1.csv',nrows=0).columns
gmap={g:next(c for c in cols if c.startswith(g+' (')) for g in genes}
assert len(gmap)==18
x=pd.read_csv(P/'OmicsExpressionProteinCodingGenesTPMLogp1.csv',usecols=[cols[0]]+list(gmap.values()))
x=x.rename(columns={cols[0]:'ModelID',**{v:k for k,v in gmap.items()}})
assert x.ModelID.is_unique
ident['has_rna']=ident.ModelID.isin(x.ModelID)
ident['status']=np.where(ident.ModelID.isna(),'no exact metadata match',np.where(ident.has_rna,'exact identity and RNA present','metadata matched; RNA absent'))
ident.to_csv(R/'analysis/independent_model_identity.csv',index=False,encoding='utf-8-sig')
rna=d.merge(ident[ident.has_rna][['model','ModelID','PatientID']],on='model',validate='one_to_one').merge(x,on='ModelID',validate='one_to_one')
assert not rna[genes].isna().any().any()
# The source already contains log2(TPM+1); do not transform it a second time.
rna.to_csv(R/'analysis/RNA_matched_feature_target_join.csv',index=False,encoding='utf-8-sig')
matched_patients=rna.groupby('PatientID').model.agg(list)
assert all(set(v)=={'CHLA-9','CHLA-10'} for v in matched_patients if len(v)>1)
spec={'release':'10.25452/figshare.plus.27993248.v1','modules':modules,'gene_columns':gmap,'n_resource_models':len(x),'n_resource_genes':len(cols)-1,'n_matched_models':len(rna),'n_families':rna.family.nunique(),'module_formula':'mean(train-z pro) - mean(train-z anti); module and ligand final scaling also training-only','units':'source log2(TPM+1), already transformed','analysis':'exploratory extension fixed before reading outcome associations','RNA_is_functional_priming':False}
(R/'analysis/molecular_specs.json').write_text(json.dumps(spec,indent=2),encoding='utf8')
rec=['F4B_induced__'+s for s in ['ICAM1','MHC_I','HLA_E']]
names=['mean','lineage','RNA_two_modules','induced_recognition','induced_recognition_plus_RNA',
 'RNA_bounded_z5','induced_recognition_plus_RNA_bounded_z5','RNA_log_balance','induced_recognition_plus_RNA_log_balance']
W=lambda z:(1/z.family.map(z.family.value_counts())).to_numpy()
folds=[];preds=[];metrics=[];coef=[];pre=[];clips=[]
def run_cv(data,name,alpha,split,label):
 data=data.reset_index(drop=True);out=np.empty(len(data));y=data.target_pp.to_numpy(float)
 for held in data[split].unique():
  ti=np.where(data[split].to_numpy()==held)[0];ri=np.where(data[split].to_numpy()!=held)[0]
  a=data.iloc[ri];b=data.iloc[ti];w=W(a);mean=np.average(y[ri],weights=w)
  assert not set(a.family)&set(b.family)
  if name=='mean':pred=np.repeat(mean,len(b))
  else:
   if name=='lineage':
    enc=OneHotEncoder(sparse_output=False,handle_unknown='ignore');X=enc.fit_transform(a[['lineage']]);Z=enc.transform(b[['lineage']]);fn=list(enc.get_feature_names_out())
   elif name=='induced_recognition':X=a[rec].to_numpy();Z=b[rec].to_numpy();fn=rec
   else:
    gs=StandardScaler().fit(a[genes],sample_weight=w)
    az=pd.DataFrame(gs.transform(a[genes]),columns=genes);bz=pd.DataFrame(gs.transform(b[genes]),columns=genes)
    if name.endswith('bounded_z5'):
     for side,zz,models_side in [('train',az,a.model),('test',bz,b.model)]:
      for ii,jj in np.argwhere(abs(zz.to_numpy())>5):clips.append(dict(dataset=label,split=split,alpha=alpha,model_name=name,held_group=held,side=side,tumor_model=models_side.iloc[ii],gene=genes[jj],z_before=zz.iloc[ii,jj],z_after=float(np.clip(zz.iloc[ii,jj],-5,5))))
     az=az.clip(-5,5);bz=bz.clip(-5,5)
    if name.endswith('log_balance'):az=a[genes].reset_index(drop=True);bz=b[genes].reset_index(drop=True)
    X=np.column_stack([az[m['pro']].mean(axis=1)-az[m['anti']].mean(axis=1) for m in modules.values()])
    Z=np.column_stack([bz[m['pro']].mean(axis=1)-bz[m['anti']].mean(axis=1) for m in modules.values()]);fn=list(modules)
    for g,mu,sd in zip(genes,gs.mean_,gs.scale_):pre.append(dict(dataset=label,split=split,alpha=alpha,model_name=name,held_group=held,gene=g,train_mean=mu,train_sd=sd))
    if name.startswith('induced_recognition_plus_RNA'):X=np.column_stack([X,a[rec]]);Z=np.column_stack([Z,b[rec]]);fn+=rec
   sc=StandardScaler().fit(X,sample_weight=w);xs=sc.transform(X);zs=sc.transform(Z)
   fit=Ridge(alpha=alpha).fit(xs,y[ri],sample_weight=w);pred=fit.predict(zs)
   for f,c in zip(fn,fit.coef_):coef.append(dict(dataset=label,split=split,alpha=alpha,model_name=name,held_group=held,feature=f,coefficient=c))
  out[ti]=pred
  folds.append(dict(dataset=label,split=split,alpha=alpha,model_name=name,held_group=held,train_models='|'.join(a.model),test_models='|'.join(b.model),train_patients='|'.join(a.PatientID),test_patients='|'.join(b.PatientID),family_overlap=bool(set(a.family)&set(b.family)),patient_overlap=bool(set(a.PatientID)&set(b.PatientID))))
 return out

def evaluate(data,label,alphas=(10,1,100),splits=('family','lineage')):
 for split in splits:
  results={}
  for alpha in alphas:
   for name in names:
    if name=='mean' and alpha!=10:continue
    results[(alpha,name)]=run_cv(data,name,alpha,split,label)
  y=data.target_pp.to_numpy(float);w=W(data);base=results[(10,'mean')];bmse=np.average((base-y)**2,weights=w)
  for (alpha,name),p in results.items():
   mse=np.average((p-y)**2,weights=w)
   metrics.append(dict(dataset=label,split=split,alpha=alpha,model_name=name,n_models=len(data),n_families=data.family.nunique(),RMSE_family_pp=np.sqrt(mse),MAE_family_pp=np.average(abs(p-y),weights=w),MSE_improvement_vs_fold_mean=1-mse/bmse,Spearman_model=np.nan if name=='mean' else spearmanr(p,y).statistic))
   for (_,row),v in zip(data.iterrows(),p):preds.append(dict(dataset=label,split=split,alpha=alpha,model_name=name,tumor_model=row.model,lineage=row.lineage,family=row.family,patient_id=row.PatientID,observed_pp=row.target_pp,predicted_pp=v,error_pp=v-row.target_pp))
evaluate(rna,'RNA_exact_match')
complete=rna[rna.n_identified_donor_points==4].copy()
evaluate(complete,'RNA_complete_donor',alphas=(10,))
don=pd.read_csv(R/'provenance/donor_values.csv')
don['donor']=pd.to_numeric(don.donor,errors='coerce');don['delta_lysis_pp']=pd.to_numeric(don.delta_lysis_pp,errors='coerce')
for omit in [1,2,3,4]:
 tmp=complete.copy();v=don[(don.donor!=omit)&don.model.isin(tmp.model)].groupby('model').delta_lysis_pp.mean();tmp['target_pp']=tmp.model.map(v)
 assert tmp.target_pp.notna().all()
 evaluate(tmp,f'RNA_omit_donor{omit}',alphas=(10,),splits=('family',))
for filename,rows in [('RNA_prediction_metrics.csv',metrics),('RNA_out_of_fold_predictions.csv',preds),('RNA_fold_membership.csv',folds),('RNA_fold_coefficients.csv',coef),('RNA_fold_gene_preprocessing.csv',pre),('RNA_bounded_z_clipping_events.csv',clips)]:pd.DataFrame(rows).to_csv(R/'analysis'/filename,index=False,encoding='utf-8-sig')

# Protein abundance is retained as an orthogonal descriptive check, not a surrogate NK label.
rp=pd.read_csv(R/'raw/independent/CCLE_RPPA_20181003.csv').rename(columns={'Unnamed: 0':'CCLEName'})
markers=['Bax','Bcl-2','Bcl-xL','Bim(CST2933)','Bim(EP1036)','Bak_Caution','Bid_Caution','Caspase-7_cleavedD198_Caution','Caspase-8_Caution']
prot=ident.merge(rp[['CCLEName']+markers],on='CCLEName',how='left',validate='many_to_one')
prot['has_RPPA']=prot['Bax'].notna();prot.to_csv(R/'analysis/RPPA_identity_and_values.csv',index=False,encoding='utf-8-sig')
ab=pd.read_csv(R/'raw/independent/CCLE_RPPA_Ab_info_20181226.csv');ab[ab.Antibody_Name.isin(markers)].to_csv(R/'analysis/RPPA_antibodies.csv',index=False,encoding='utf-8-sig')
pairs=[]
for mark,g in [('Bax','BAX'),('Bcl-2','BCL2'),('Bcl-xL','BCL2L1'),('Bim(CST2933)','BCL2L11'),('Bim(EP1036)','BCL2L11')]:
 z=prot[['model',mark]].merge(rna[['model',g]],on='model').dropna();rho=spearmanr(z[mark],z[g]).statistic
 pairs.append(dict(protein=mark,gene=g,n_models=len(z),Spearman=rho,interpretation='descriptive RNA-protein agreement; overlapping stocks, not independent functional validation'))
pd.DataFrame(pairs).to_csv(R/'analysis/RNA_RPPA_descriptive_agreement.csv',index=False,encoding='utf-8-sig')
print('RNA',len(rna),'models',rna.family.nunique(),'families; total RNA matrix',x.shape,'RPPA',int(prot.has_RPPA.sum()))
m=pd.DataFrame(metrics);print(m[(m.dataset=='RNA_exact_match')&(m.alpha==10)].to_string(index=False))
