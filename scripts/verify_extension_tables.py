"""Recompute retained extension claims from the included source-level tables."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
R=Path(__file__).resolve().parents[1]
(R/'qa').mkdir(exist_ok=True)
checks=[]
def check(name,condition,details=None):
    checks.append(dict(check=name,passed=bool(condition),details=details))
    assert condition,name
for fn in ['Wang2012_Fig6_readings.csv','Wang2012_recognition_factorial_readings.csv','Ovarian2026_Fig5_donor_readings.csv']:
    d=pd.read_csv(R/'analysis'/fn)
    if 'axis_low' in d:
        estimate=d.axis_low+(d.y_low_pixel-d.y_pixel)*(d.axis_high-d.axis_low)/(d.y_low_pixel-d.y_high_pixel)
    else:estimate=d.axis_max*(d.y_zero_pixel-d.y_pixel)/(d.y_zero_pixel-d.y_high_pixel)
    error=float(np.max(abs(estimate-d.value_pp)))
    check('Source-coordinate calibration '+fn,error<1e-10,{'rows':len(d),'max_error_pp':error})
d=pd.read_csv(R/'analysis/Ovarian2026_Fig5_donor_readings.csv')
check('Eight ovarian arms and 40 graphical observations',len(d)==40 and d.groupby(['condition','NK_product']).size().eq(5).all())
g=d.groupby(['NK_product','condition']).value_pp.mean().unstack();effect=g.IFNg-g.control
out=pd.read_csv(R/'analysis/Ovarian2026_NKG2A_interactions.csv').set_index('product')
for product in ['mock','PRAME']:
    check('Ovarian effect modification '+product,np.isclose(effect[product+'_KO']-effect[product+'_WT'],out.loc[product,'interaction_pp'],atol=1e-12))
m=pd.read_csv(R/'analysis/Remsik2025_Fig5i_mouse_values.csv')
check('28 mouse values in four unpaired arms',len(m)==28 and m.groupby(['ifng_overexpression','nk_depletion']).size().eq(7).all())
means=m.groupby(['ifng_overexpression','nk_depletion']).radiance.apply(lambda x:np.log10(x).mean())
host=json.loads((R/'analysis/Remsik2025_Fig5i_summary.json').read_text())
interaction=(means[1,1]-means[0,1])-(means[1,0]-means[0,0])
check('Host descriptive interaction',np.isclose(interaction,host['interaction_depleted_minus_isotype']),{'interaction_log10':interaction})
allr=pd.read_csv(R/'analysis/Dufva_all_gene_results.csv.gz')
check('Full processed genetic-screen table retained',len(allr)==222609 and allr.screen.nunique()==11)
panel=pd.read_csv(R/'analysis/Dufva_21gene_panel.csv')
check('Prespecified gene display complete',len(panel)==231 and panel.groupby('screen').size().eq(21).all())
ambig=allr.gene.astype(str).str.match(r'^\d{1,2}-(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)$')|allr.duplicated(['screen','gene'],False)
check('113 unresolved gene rows retained',int(ambig.sum())==113)
p=allr[~ambig & allr.screen.isin(['M. KMS11 GOF','N. KMS11 GOF (KHYG1)'])].pivot(index='gene',columns='screen',values='effect').dropna()
rho=float(spearmanr(p.iloc[:,0],p.iloc[:,1]).statistic)
check('Descriptive KMS11 correlation',len(p)==18903 and np.isclose(rho,0.062006245341541394),{'shared_labels':len(p),'spearman':rho})
(R/'qa/extension_table_checks.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
print(json.dumps({'checks':len(checks),'passed':True}))
