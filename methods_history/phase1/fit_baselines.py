"""Prespecified grouped ridge baselines on reconstructed S11 data."""
from pathlib import Path
import json,hashlib,shutil
import numpy as np,pandas as pd
from sklearn.preprocessing import StandardScaler,OneHotEncoder
from sklearn.linear_model import Ridge
from sklearn.metrics import balanced_accuracy_score
from scipy.stats import spearmanr
R=Path(__file__).resolve().parents[1];B=Path(r'phase0b_runs\20260922_170151')
for name,src in [('model_summary.csv',B/'analysis/S11/model_summary.csv'),('donor_values.csv',B/'analysis/S11/reconciled_donor_values.csv'),('ABT263_match.csv',B/'analysis/S12/S11_matched_ABT263.csv')]:
 p=R/'provenance'/name
 if not p.exists():shutil.copyfile(src,p)
features=pd.read_csv(R/'analysis/ligand_features_wide.csv');targets=pd.read_csv(R/'provenance/model_summary.csv')
d=targets.merge(features,on='model',validate='one_to_one');assert len(d)==22
d['family']=np.where(d.model.isin(['CHLA-9','CHLA-10']),'CHLA9_10_same_patient',d.model)
d['target_pp']=d.mean_center_or_interval_midpoint_pp
pre=lambda names:['F3_baseline__'+x for x in names]
post=lambda names:['F4B_induced__'+x for x in names]
recognition=['ICAM1','MHC_I','HLA_E'];death=['FAS_CD95','DR4_CD261','DR5_CD262']
models={'mean':[], 'lineage':['lineage'], 'baseline_recognition':pre(recognition),'induced_recognition':post(recognition),
 'baseline_death_receptors':pre(death),'induced_recognition_plus_baseline_death':post(recognition)+pre(death),
 'induced_recognition_plus_induced_death':post(recognition+death),'all_baseline_ligands':[x for x in d.columns if x.startswith('F3_')],
 'all_induced_ligands':[x for x in d.columns if x.startswith('F4B_')]}
def weights(data):return (1/data.family.map(data.family.value_counts())).to_numpy()
foldlog=[];coeflog=[]
def predict_cv(data,cols,alpha,split='family',run='',capture=True):
 data=data.reset_index(drop=True);y=data.target_pp.to_numpy();out=np.full(len(data),np.nan);grp=data[split].to_numpy()
 for held in pd.unique(grp):
  te=np.where(grp==held)[0];tr=np.where(grp!=held)[0];a=data.iloc[tr];b=data.iloc[te];w=weights(a);mu=np.average(y[tr],weights=w)
  assert not set(a.family)&set(b.family) if split=='family' else not set(a.lineage)&set(b.lineage)
  if not cols:pred=np.full(len(te),mu);co=[];fn=[]
  else:
   if cols==['lineage']:
    enc=OneHotEncoder(sparse_output=False,handle_unknown='ignore');X=enc.fit_transform(a[cols]);Z=enc.transform(b[cols]);fn=list(enc.get_feature_names_out())
   else:X=a[cols].to_numpy(float);Z=b[cols].to_numpy(float);fn=cols
   scaler=StandardScaler().fit(X,sample_weight=w);Xs=scaler.transform(X);Zs=scaler.transform(Z)
   reg=Ridge(alpha=alpha).fit(Xs,y[tr],sample_weight=w);pred=reg.predict(Zs);co=reg.coef_
  out[te]=pred
  if capture:
   foldlog.append(dict(run=run,held_group=held,split=split,n_train_models=len(tr),n_test_models=len(te),train_models='|'.join(a.model),test_models='|'.join(b.model),family_overlap=bool(set(a.family)&set(b.family)),training_mean=mu,alpha=alpha))
   for feature,coef in zip(fn,co):coeflog.append(dict(run=run,held_group=held,feature=feature,standardized_coefficient=coef))
 return out
def metrics(data,pred,meanpred):
 y=data.target_pp.to_numpy();w=weights(data);mse=np.average((y-pred)**2,weights=w);base=np.average((y-meanpred)**2,weights=w)
 return dict(n_models=len(data),n_families=data.family.nunique(),MAE_family_pp=np.average(abs(y-pred),weights=w),RMSE_family_pp=np.sqrt(mse),MAE_model_pp=np.mean(abs(y-pred)),RMSE_model_pp=np.sqrt(np.mean((y-pred)**2)),MSE_improvement_vs_fold_mean=1-mse/base,Spearman_model=float(spearmanr(y,pred).statistic),direction_balanced_accuracy=balanced_accuracy_score(y>0,pred>0,sample_weight=w),positive_models=int(sum(y>0)))
summ=[];predrows=[]
def evaluate(data,label,modeldict,alphas=(10,1,100),splits=('family','lineage')):
 for split in splits:
  base=predict_cv(data,[],10,split,run=label+'_'+split+'_mean',capture=False)
  for alpha in alphas:
   for name,cols in modeldict.items():
    if name=='mean' and alpha!=10:continue
    run=f'{label}__{split}__a{alpha}__{name}';pred=predict_cv(data,cols,alpha,split,run)
    met=metrics(data,pred,base)
    # Leave-one-family-out training means depend algebraically on the held-out y;
    # their OOF rank correlation is not an informative predictive association.
    if name=='mean':met['Spearman_model']=np.nan
    summ.append(dict(dataset=label,split=split,alpha=alpha,model_name=name,features='|'.join(cols),**met))
    for (_,r),pr in zip(data.iterrows(),pred):predrows.append(dict(dataset=label,split=split,alpha=alpha,model_name=name,tumor_model=r.model,lineage=r.lineage,family=r.family,observed_pp=r.target_pp,predicted_pp=pr,error_pp=pr-r.target_pp))
evaluate(d,'all22',models)
full=d[d.n_identified_donor_points==4].copy();assert len(full)==20
evaluate(full,'complete_donor20',models,alphas=(10,))
don=pd.read_csv(R/'provenance/donor_values.csv')
don['donor']=pd.to_numeric(don.donor,errors='coerce');don['delta_lysis_pp']=pd.to_numeric(don.delta_lysis_pp,errors='coerce')
assert don[don.model.isin(full.model)].delta_lysis_pp.notna().all()
for donor in [1,2,3,4]:
 sub=don[(don.donor!=donor)&don.model.isin(full.model)];assert (sub.groupby('model').size()==3).all()
 a=sub.groupby('model').delta_lysis_pp.mean()
 dd=full.copy();dd['target_pp']=dd.model.map(a);evaluate(dd,f'omit_donor{donor}_20',models,alphas=(10,),splits=('family',))
drug=pd.read_csv(R/'provenance/ABT263_match.csv');drug=drug[drug.use_in_correlation.astype(str).str.lower()=='true']
dd=d.merge(drug[['model','log10_relative_sensitivity','censor_status']],on='model',validate='one_to_one');assert len(dd)==20
dd['ABT263_censored']=dd.censor_status.str.contains('ceiling').astype(int)
drugmodels={'mean':[],'lineage':['lineage'],'induced_recognition':post(recognition),'ABT263_proxy':['log10_relative_sensitivity','ABT263_censored'],'induced_recognition_plus_ABT263':post(recognition)+['log10_relative_sensitivity','ABT263_censored']}
evaluate(dd,'drug_matched20',drugmodels)
d.to_csv(R/'analysis/model_feature_target_join.csv',index=False,encoding='utf-8-sig')
dd.to_csv(R/'analysis/drug_feature_target_join.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(summ).to_csv(R/'analysis/prediction_metrics.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(predrows).to_csv(R/'analysis/out_of_fold_predictions.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(foldlog).to_csv(R/'analysis/fold_membership.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(coeflog).to_csv(R/'analysis/fold_coefficients.csv',index=False,encoding='utf-8-sig')
(R/'analysis/model_specs.json').write_text(json.dumps({'alpha_primary':10,'alpha_sensitivity':[1,100],'models':models,'drug_models':drugmodels,'feature_selection':'fixed before modelling; no search','preprocessing':'training-fold weighted standardization','weights':'each family equal; CHLA9/10 split one family weight','target':'IFNg-treated minus untreated source lysis difference, averaged over shared donors','validation':'within-study family holdout; no external validation'},indent=2),encoding='utf-8')
s=pd.DataFrame(summ);print(s[(s.dataset=='all22')&(s.split=='family')&(s.alpha==10)][['model_name','MAE_family_pp','RMSE_family_pp','MSE_improvement_vs_fold_mean','Spearman_model']].to_string(index=False))
print('\nDrug increment');print(s[(s.dataset=='drug_matched20')&(s.split=='family')&(s.alpha==10)][['model_name','RMSE_family_pp','MSE_improvement_vs_fold_mean']].to_string(index=False))
