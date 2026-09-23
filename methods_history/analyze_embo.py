from pathlib import Path
import json, hashlib, itertools
import numpy as np
import pandas as pd
from scipy import stats
R=Path(__file__).resolve().parents[1]
J=json.loads((R/'derived/EMBO_workbooks.json').read_text(encoding='utf-8'))
records=[]; contrasts=[]; summaries=[]; qc=[]
def block(suffix,sheet,start,end,cols,endpoint,unit,cluster=False):
    filename=f'EMBR-23-e54133-{suffix}.xlsx'
    rows={r['row']:r['cells'] for r in J[filename][sheet]}
    data=[]
    for i,r in enumerate(range(start,end+1)):
        v={cond:float(rows[r][str(col)]) for cond,col in cols.items()}
        row={'endpoint':endpoint,'donor_within_panel':i+1,'batch':i//3+1 if cluster else None,**v}
        data.append(row)
        for cond,col in cols.items():
            records.append({'endpoint':endpoint,'donor_within_panel':i+1,'batch':row['batch'],
                'condition':cond,'value':v[cond],'unit':unit,'source_file':f'raw/PMC9346491/{filename}',
                'source_sheet':sheet,'source_row':r,'source_column':col,
                'source_sha256':hashlib.sha256((R/'raw/PMC9346491'/filename).read_bytes()).hexdigest(),
                'pairing':'same published row within panel; donor identities not linked across panels'})
    return pd.DataFrame(data)
def perm_p(d):
    d=np.asarray(d,float)
    signs=np.array(list(itertools.product([-1,1],repeat=len(d))))
    return float(np.mean(np.abs(signs@d)>=abs(d.sum())-1e-10))
def add(df,endpoint,treatment,control,unit,cluster=False):
    d=df[treatment]-df[control]
    for i,x in enumerate(d):
        contrasts.append({'endpoint':endpoint,'donor_within_panel':i+1,'batch':df.iloc[i]['batch'],
            'treatment':treatment,'control':control,'difference':x,'unit':unit})
    primary=pd.DataFrame({'d':d,'batch':df.batch}).groupby('batch').d.mean().to_numpy() if cluster else d.to_numpy()
    n=len(primary);se=np.std(primary,ddof=1)/np.sqrt(n)
    summaries.append({'endpoint':endpoint,'treatment':treatment,'control':control,
        'n_donors':len(d),'n_primary_units':n,'primary_unit':'shared-control batch' if cluster else 'donor',
        'estimate_mean':primary.mean(),'median_donor_difference':d.median(),'sd_primary_difference':np.std(primary,ddof=1),
        't_interval_lower':primary.mean()-stats.t.ppf(.975,n-1)*se,'t_interval_upper':primary.mean()+stats.t.ppf(.975,n-1)*se,
        'p_exact_signflip':perm_p(primary),'p_donor_signflip_sensitivity':perm_p(d),'unit':unit,
        'positive_donors':int((d>0).sum()),'negative_donors':int((d<0).sum()),
        'inference_note':'exploratory; exact sign-flip needs symmetric exchangeable null differences; t interval additionally assumes approximately normal unit means'})
    return d

degran=block('s012','Fig.4',5,14,{'no_antibody':5,'isotype':6,'anti_TRAIL':7},'soluble_block_CD107a','percentage points')
add(degran,'soluble_block_CD107a','anti_TRAIL','isotype','pp')
cyto=block('s003','Fig.6',5,12,{'no_antibody':11,'isotype':12,'anti_TRAIL':13},'soluble_block_IFNg','pg/ml')
add(cyto,'soluble_block_IFNg','anti_TRAIL','isotype','pg/ml')
adh=block('s003','Fig.6',18,25,{'PBS':3,'isotype':4,'anti_TRAIL':5},'soluble_block_conjugates','percentage points')
add(adh,'soluble_block_conjugates','anti_TRAIL','isotype','pp')
plate=block('s003','Fig.6',5,12,{'PBS':2,'isotype':3,'anti_NKp46':4,'anti_TRAIL':5},'plate_bound_IFNg','pg/ml')
add(plate,'plate_bound_IFNg','anti_TRAIL','isotype','pg/ml')
kill=block('s001','Fig.7',5,16,{'isotype':2,'anti_TRAIL':3,'published_delta':5},'soluble_block_remaining_targets','percentage points',True)
d=add(kill,'soluble_block_remaining_targets','anti_TRAIL','isotype','pp',True)
qc.append({'check':'Fig7A numeric published delta equals alphaTRAIL minus isotype (source text header has reversed order)',
           'max_abs_error':float(np.abs(d-kill.published_delta).max()),'pass':bool(np.allclose(d,kill.published_delta,atol=1e-8))})
for start,end,endpoint in [(23,34,'DR4_DR5_KO_relative_clearance'),(41,52,'DR5_overexpression_relative_clearance')]:
    z=block('s001','Fig.7',start,end,{'no_NK_ratio':3,'NK_ratio':4,'published_specific_lysis':6},endpoint,'ratio / percent',True)
    z['derived_specific_lysis']=100*(1-z.no_NK_ratio/z.NK_ratio)
    z['zero']=0
    add(z,endpoint,'derived_specific_lysis','zero','relative clearance (%)',True)
    qc.append({'check':endpoint+' rounded ratios reproduce published normalized percentages',
               'max_abs_error':float(np.abs(z.derived_specific_lysis-z.published_specific_lysis).max()),
               'pass':bool(np.max(np.abs(z.derived_specific_lysis-z.published_specific_lysis))<.2)})
S=pd.DataFrame(summaries)
p=S.p_exact_signflip.to_numpy();order=np.argsort(p);adj=np.minimum(1,np.maximum.accumulate(p[order]*(len(p)-np.arange(len(p)))))
S['p_holm_7_endpoints']=np.nan;S.loc[order,'p_holm_7_endpoints']=adj
pd.DataFrame(records).to_csv(R/'provenance/EMBO_individual_source_values.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(contrasts).to_csv(R/'analysis/EMBO_paired_contrasts.csv',index=False,encoding='utf-8-sig')
S.to_csv(R/'analysis/EMBO_effect_summary.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(contrasts).query('batch == batch').groupby(['endpoint','batch'],as_index=False).difference.mean().to_csv(R/'analysis/EMBO_batch_contrasts.csv',index=False,encoding='utf-8-sig')
(R/'qa/EMBO_numeric_checks.json').write_text(json.dumps(qc,indent=2),encoding='utf-8')
assert all(q['pass'] for q in qc)
print(S[['endpoint','n_donors','n_primary_units','estimate_mean','p_exact_signflip','p_holm_7_endpoints']].to_string(index=False))
