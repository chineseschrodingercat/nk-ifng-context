"""Independent source-display and numerical checks of the fitted analyses."""
from pathlib import Path
import json
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[1]
f=pd.read_csv(R/'analysis/ligand_features_long.csv')
# Manual source F7 axis and bar-top readings in original embedded-image pixels.
# Each is a second visual display of the same experiment, NOT independent biology.
spec=[('Kasumi-1',376,240,3,313,267),('MOLT-4',377,242,1.5,305,277),('CHLA-136',383,240,1.5,375,263),
 ('BT-12',371,223,4,241,369),('SJ-GBM2',377,235,1,250,315),('NB-1643',378,230,4,253,293),
 ('CHLA-9',816,668,4,694,777),('CHLA-10',813,668,2.5,695,760),('Ramos-RA1',813,670,.8,685,758)]
rows=[]
for model,y0,yt,vm,a,b in spec:
 for marker,yy in [('ICAM1',a),('MHC_I',b)]:
  v=vm*(y0-yy)/(y0-yt);heat=f[(f.model==model)&(f.marker==marker)&(f.panel=='F4B_induced')].value.item()
  rows.append(dict(model=model,marker=marker,axis_zero_y=y0,axis_top_y=yt,axis_max=vm,bar_top_y=yy,F7_bar_reading=v,F4_colour_reading=heat,absolute_difference=abs(v-heat),source_image='derived/S11/p9_Im0.png',reading_note='manual bar/axis reading approx 2 px; paired display concordance only'))
check=pd.DataFrame(rows);check.to_csv(R/'qa/F7_F4_display_concordance.csv',index=False,encoding='utf-8-sig')
assert check.absolute_difference.max()<.1

# Independent weighted ridge normal equations, avoiding the sklearn estimator.
d=pd.read_csv(R/'analysis/model_feature_target_join.csv')
p=pd.read_csv(R/'analysis/out_of_fold_predictions.csv')
rn=pd.read_csv(R/'analysis/RNA_matched_feature_target_join.csv')
rp=pd.read_csv(R/'analysis/RNA_out_of_fold_predictions.csv')
mods=json.loads((R/'analysis/molecular_specs.json').read_text())['modules']
genes=[g for m in mods.values() for gs in m.values() for g in gs]
recognition=['F4B_induced__'+x for x in ['ICAM1','MHC_I','HLA_E']]
def standard(X,Z,w):
 mu=np.average(X,axis=0,weights=w);sd=np.sqrt(np.average((X-mu)**2,axis=0,weights=w));sd=np.where(sd==0,1,sd)
 return (X-mu)/sd,(Z-mu)/sd
numerical=[]
for data,pred,label,name,mode in [(d,p,'all22','induced_recognition','ligand'),(rn,rp,'RNA_exact_match','induced_recognition_plus_RNA','rna'),(rn,rp,'RNA_exact_match','induced_recognition_plus_RNA_bounded_z5','bounded')]:
 for held in data.family.unique():
  a=data[data.family!=held];b=data[data.family==held];w=(1/a.family.map(a.family.value_counts())).to_numpy();y=a.target_pp.to_numpy();mu=np.average(y,weights=w)
  if mode=='ligand':X=a[recognition].to_numpy();Z=b[recognition].to_numpy()
  else:
   ag,bg=standard(a[genes].to_numpy(),b[genes].to_numpy(),w)
   if mode=='bounded':ag=np.clip(ag,-5,5);bg=np.clip(bg,-5,5)
   A=pd.DataFrame(ag,columns=genes);B=pd.DataFrame(bg,columns=genes)
   X=np.column_stack([A[m['pro']].mean(axis=1)-A[m['anti']].mean(axis=1) for m in mods.values()]);Z=np.column_stack([B[m['pro']].mean(axis=1)-B[m['anti']].mean(axis=1) for m in mods.values()])
   X=np.column_stack([X,a[recognition]]);Z=np.column_stack([Z,b[recognition]])
  X,Z=standard(X,Z,w)
  beta=np.linalg.solve(X.T@(w[:,None]*X)+10*np.eye(X.shape[1]),X.T@(w*(y-mu)))
  manual=mu+Z@beta
  saved=pred[(pred.dataset==label)&(pred.split=='family')&(pred.alpha==10)&(pred.model_name==name)].set_index('tumor_model').loc[b.model,'predicted_pp'].to_numpy()
  numerical.append(dict(dataset=label,model_name=name,held_family=held,max_absolute_difference=float(max(abs(manual-saved)))))
num=pd.DataFrame(numerical);num.to_csv(R/'qa/independent_ridge_equation_check.csv',index=False)
assert num.max_absolute_difference.max()<1e-9
fl=pd.read_csv(R/'analysis/fold_membership.csv');rf=pd.read_csv(R/'analysis/RNA_fold_membership.csv')
assert not fl.family_overlap.any() and not rf.family_overlap.any() and not rf.patient_overlap.any()
assert f.record_id.is_unique and len(f)==704
assert p.predicted_pp.notna().all() and rp.predicted_pp.notna().all()
record={'heatmap_records':len(f),'display_checks':len(check),'maximum_display_difference':check.absolute_difference.max(),'mean_display_difference':check.absolute_difference.mean(),'independent_ridge_folds_checked':len(num),'max_ridge_numeric_difference':num.max_absolute_difference.max(),'baseline_fold_records':len(fl),'RNA_fold_records':len(rf),'family_or_patient_leakage':False,'shared_donors_are_external_validation':False,'uncertainty_note':'reading envelopes and omit-donor sensitivity are not population confidence intervals'}
(R/'qa/analysis_validation.json').write_text(json.dumps(record,indent=2),encoding='utf8')
print(json.dumps(record,indent=2))
