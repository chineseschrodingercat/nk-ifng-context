"""Uniform perturbations inside recorded reading envelopes; not biological CI."""
from pathlib import Path
import json
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[1];rng=np.random.default_rng(20260922)
d=pd.read_csv(R/'analysis/model_feature_target_join.csv')
l=pd.read_csv(R/'analysis/ligand_features_long.csv');l['key']=l.panel+'__'+l.marker
cols=[c for c in d if c.startswith(('F3_','F4B_'))]
lo=l.pivot(index='model',columns='key',values='reading_low').loc[d.model,cols].to_numpy()
hi=l.pivot(index='model',columns='key',values='reading_high').loc[d.model,cols].to_numpy()
models=json.loads((R/'analysis/model_specs.json').read_text())['models'];models={n:v for n,v in models.items() if n!='lineage'}
idx={n:[cols.index(c) for c in cs] for n,cs in models.items()}
w=(1/d.family.map(d.family.value_counts())).to_numpy();splits=[]
for group in d.family.unique():
 te=np.flatnonzero((d.family==group).to_numpy());tr=np.flatnonzero((d.family!=group).to_numpy());splits.append((tr,te))
rows=[]
for repeat in range(200):
 X=rng.uniform(lo,hi);y=rng.uniform(d.mean_lower_reading_pp,d.mean_upper_reading_pp);pr={name:np.empty(len(d)) for name in models}
 for tr,te in splits:
  ww=w[tr];mu=np.average(y[tr],weights=ww);pr['mean'][te]=mu
  for name,ii in idx.items():
   if name=='mean':continue
   a=X[np.ix_(tr,ii)];b=X[np.ix_(te,ii)];means=np.average(a,axis=0,weights=ww);sd=np.sqrt(np.average((a-means)**2,axis=0,weights=ww));sd=np.where(sd==0,1,sd)
   a=(a-means)/sd;b=(b-means)/sd;beta=np.linalg.solve(a.T@(ww[:,None]*a)+10*np.eye(len(ii)),a.T@(ww*(y[tr]-mu)));pr[name][te]=mu+b@beta
 base=np.average((pr['mean']-y)**2,weights=w)
 for name,p in pr.items():
  mse=np.average((p-y)**2,weights=w);rows.append(dict(repeat=repeat,model_name=name,RMSE_family_pp=np.sqrt(mse),MSE_improvement_vs_fold_mean=1-mse/base))
out=pd.DataFrame(rows);out.to_csv(R/'analysis/reading_sensitivity_200.csv',index=False)
summary=out.groupby('model_name').agg(RMSE_min=('RMSE_family_pp','min'),RMSE_max=('RMSE_family_pp','max'),improvement_min=('MSE_improvement_vs_fold_mean','min'),improvement_max=('MSE_improvement_vs_fold_mean','max'),fraction_better_mean=('MSE_improvement_vs_fold_mean',lambda s:float((s>0).mean())))
summary.to_csv(R/'analysis/reading_sensitivity_summary.csv')
(R/'qa/reading_sensitivity_design.json').write_text(json.dumps({'seed':20260922,'repeats':200,'sampling':'independent uniform over recorded image reading envelopes, for numerical sensitivity only','alpha':10,'split':'leave model family out','not_population_CI':True},indent=2),encoding='utf8')
print(summary.to_string())
