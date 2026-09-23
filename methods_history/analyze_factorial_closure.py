from pathlib import Path
import itertools,json
import numpy as np
import pandas as pd
R=Path(__file__).resolve().parents[1]
SH=pd.read_csv(R/'provenance/Sheard_Fig6_readings.csv')
HT=pd.read_csv(R/'provenance/HT29_Fig3_readings.csv')
records=[]
def combo_bounds(vals,func):
    endpoints=[[r.reading_lower,r.reading_upper] for r in vals]
    out=np.array([func(*v) for v in itertools.product(*endpoints)],float)
    return float(np.nanmin(out)),float(np.nanmax(out))
def summarize(source,model,context,vals):
    # k00: untreated/unblocked; k10: IFNg/unblocked; k01: untreated/blocked; k11: IFNg/blocked.
    k00,k10,k01,k11=[r.value for r in vals]
    f=lambda a,b,c,d:(b-d)-(a-c)
    lo,hi=combo_bounds(vals,f)
    row={'study':source,'model':model,'context':context,
        'baseline_clearance':k00,'IFNg_clearance':k10,'baseline_headroom':100-k00,
        'net_gain':k10-k00,'blockade_sensitive_gain':f(k00,k10,k01,k11),
        'blockade_gain_reading_lower':lo,'blockade_gain_reading_upper':hi,
        'blocked_net_gain':k11-k01,
        'fraction_remaining_eliminated_by_IFNg':(k10-k00)/(100-k00) if k00<100 else np.nan,
        'minimum_alternative_interaction_to_remove_positive_reading_gain':max(0,lo),
        'biological_p_value':np.nan,'uncertainty':'image reading only; no biological replicates reconstructed'}
    # Survival ratio compares blocker effects multiplicatively and reduces dependence on the absolute clearance scale.
    # All four true survival bounds must be positive: otherwise the log calculation is not bounded.
    if all(100-r.reading_upper>0 for r in vals):
        logf=lambda a,b,c,d:np.log((100-d)/(100-b))-np.log((100-c)/(100-a))
        row['log_survival_blockade_interaction']=logf(k00,k10,k01,k11)
        row['log_interaction_lower'],row['log_interaction_upper']=combo_bounds(vals,logf)
    else:
        row.update({'log_survival_blockade_interaction':np.nan,'log_interaction_lower':np.nan,'log_interaction_upper':np.nan})
    records.append(row)
for mod,q in SH[SH.panel=='6C'].groupby('model',sort=False):
    # Survival to clearance, including reversal of uncertainty limits.
    q=q.copy();v=q.value.copy();l=q.reading_lower.copy();h=q.reading_upper.copy()
    q['value']=100-v;q['reading_lower']=100-h;q['reading_upper']=100-l
    for dose in sorted(d for d in q.IFNg_ng_ml.unique() if d>0):
        vals=[q[(q.IFNg_ng_ml==d)&(q.anti_TRAIL.astype(str).str.lower()==str(b).lower())].iloc[0]
              for d,b in [(0,False),(dose,False),(0,True),(dose,True)]]
        summarize('S16',mod,f'IFNg {dose:g} ng/ml',vals)
for (act,et),q in HT.groupby(['effector_activation','ET'],sort=False):
    vals=[q[(q.IFNg==g)&(q.blockade==b)].iloc[0] for g,b in [(False,'none'),(True,'none'),(False,'anti_FAS'),(True,'anti_FAS')]]
    summarize('S17','HT-29',f'{act}; ET {et:g}',vals)
df=pd.DataFrame(records)
df.to_csv(R/'analysis/headroom_and_blockade_sensitivity.csv',index=False,encoding='utf-8-sig')
orth=[]
for (act,et),q in HT.groupby(['effector_activation','ET'],sort=False):
    get=lambda g,b:q[(q.IFNg==g)&(q.blockade==b)].iloc[0]
    none=[get(g,'none') for g in [False,True]];fas=[get(g,'anti_FAS') for g in [False,True]];egta=[get(g,'EGTA_Mg') for g in [False,True]]
    # Operational residual. Additivity of routes is a hypothesis to audit, not an assumption.
    for i,g in enumerate([False,True]):
        v=[none[i],fas[i],egta[i]]
        residual=lambda a,b,c:a-b-c
        lo,hi=combo_bounds(v,residual)
        orth.append({'effector_activation':act,'ET':et,'IFNg':g,
            'unblocked_lysis':none[i].value,'anti_FAS_resistant_lysis':fas[i].value,'EGTA_resistant_lysis':egta[i].value,
            'unblocked_minus_two_resistant_readouts':residual(*[x.value for x in v]),
            'reading_lower':lo,'reading_upper':hi,
            'interpretation':'tests numerical route additivity; EGTA is an operational calcium-dependence perturbation, not a target-specific genetic isolation'})
pd.DataFrame(orth).to_csv(R/'analysis/HT29_orthogonal_perturbation_audit.csv',index=False,encoding='utf-8-sig')
O=pd.DataFrame(orth)
paired=[]
for (act,et),q in O.groupby(['effector_activation','ET'],sort=False):
    q=q.set_index('IFNg');a=q.loc[False];b=q.loc[True]
    paired.append({'effector_activation':act,'ET':et,
        'net_gain':b.unblocked_lysis-a.unblocked_lysis,
        'anti_FAS_sensitive_gain':(b.unblocked_lysis-b.anti_FAS_resistant_lysis)-(a.unblocked_lysis-a.anti_FAS_resistant_lysis),
        'EGTA_resistant_gain':b.EGTA_resistant_lysis-a.EGTA_resistant_lysis,
        'anti_FAS_resistant_gain':b.anti_FAS_resistant_lysis-a.anti_FAS_resistant_lysis})
pd.DataFrame(paired).to_csv(R/'analysis/HT29_orthogonal_paired_gains.csv',index=False,encoding='utf-8-sig')
print(df[['study','model','context','baseline_clearance','net_gain','blockade_sensitive_gain','blockade_gain_reading_lower','log_interaction_lower']].to_string(index=False))
print('ORTHOGONAL',pd.DataFrame(paired).to_string(index=False))
print('ADDIVITY residuals outside image-reading interval:',int(((O.reading_lower>0)|(O.reading_upper<0)).sum()),'of',len(O))
