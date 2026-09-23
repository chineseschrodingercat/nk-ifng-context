from pathlib import Path
import sys,json,itertools
import pandas as pd
import numpy as np
from scipy import stats
R=Path(__file__).resolve().parents[1];sys.stdout.reconfigure(encoding='utf-8')
q=pd.read_csv(R/'source_data/EMBO_paired_contrasts.csv');s=pd.read_csv(R/'source_data/EMBO_effect_summary.csv');out=[]
for row in s.itertuples():
    z=q[q.endpoint==row.endpoint]
    d=z.groupby('batch').difference.mean().to_numpy() if row.n_primary_units<row.n_donors else z.difference.to_numpy()
    half=stats.t.ppf(.975,len(d)-1)*d.std(ddof=1)/np.sqrt(len(d));lo=d.mean()-half;hi=d.mean()+half
    assert np.allclose([lo,hi],[row.t_interval_lower,row.t_interval_upper],atol=1e-10)
    out.append({'endpoint':row.endpoint,'n':len(d),'CI_verified':True})
p=s.p_exact_signflip.to_numpy();order=np.argsort(p);adjusted=np.minimum(1,np.maximum.accumulate(p[order]*(len(p)-np.arange(len(p)))))
assert np.allclose(adjusted,s.p_holm_7_endpoints.to_numpy()[order])
(R/'qa/inference_checks.json').write_text(json.dumps({'endpoints':out,'Holm_all_seven_verified':True,'scipy':__import__('scipy').__version__},indent=2))
print('All seven t intervals and Holm-adjusted P values reproduced.')
