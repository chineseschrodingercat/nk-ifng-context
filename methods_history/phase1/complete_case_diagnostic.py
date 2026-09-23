"""Decompose already generated complete-case prediction errors; no refitting."""
from pathlib import Path
import pandas as pd
R=Path(__file__).resolve().parents[1]
p=pd.read_csv(R/'analysis/out_of_fold_predictions.csv')
q=p[(p.dataset=='complete_donor20')&(p.split=='family')&(p.alpha==10)&p.model_name.isin(['mean','baseline_death_receptors'])]
z=q.pivot(index=['tumor_model','family','observed_pp'],columns='model_name',values='predicted_pp').reset_index()
z['squared_error_gain_pp2']=(z['mean']-z.observed_pp)**2-(z.baseline_death_receptors-z.observed_pp)**2
f=z.groupby('family').squared_error_gain_pp2.mean().sort_values(ascending=False)
assert len(z)==20 and len(f)==19
z.to_csv(R/'analysis/complete_case_receptor_error_contributions.csv',index=False)
s=pd.DataFrame({'family':f.index,'squared_error_gain_pp2':f.values})
s['mean_error_gain_after_removing_family_fixed_predictions']=(f.sum()-s.squared_error_gain_pp2)/(len(f)-1)
s.to_csv(R/'analysis/complete_case_receptor_family_diagnostic.csv',index=False)
print('Mean fixed-prediction squared-error advantage:',f.mean())
print('After omitting SJ-GBM2 error contribution:',s.loc[s.family=='SJ-GBM2','mean_error_gain_after_removing_family_fixed_predictions'].item())
