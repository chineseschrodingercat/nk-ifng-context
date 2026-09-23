from pathlib import Path
from lxml import etree as E
import sys,json,shutil,hashlib,itertools,platform
import numpy as np
import pandas as pd
from openpyxl import load_workbook
sys.stdout.reconfigure(encoding='utf-8')
R=Path(__file__).resolve().parents[1]
P=R; C=R; C0=R; F=R
manifest=[];checks=[]
previous_manifest=R/'provenance/INPUT_MANIFEST.csv'
previous={x['package_file']:x for x in pd.read_csv(previous_manifest).to_dict('records')} if previous_manifest.exists() else {}
manifest=[dict(package_file=str(p.relative_to(R)),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in (R/'source_data').iterdir() if p.is_file()]
don=pd.read_csv(R/'source_data/EMBO_individual_source_values.csv')
books={f:load_workbook(R/'source_data'/f,data_only=True) for f in don.source_file.map(lambda x:Path(x).name).unique()}
delta=[]
for r in don.itertuples():
    wb=books[Path(r.source_file).name]
    v=wb[r.source_sheet].cell(int(r.source_row),int(r.source_column)).value
    delta.append(abs(float(v)-r.value))
checks.append(dict(check='218 values match original published spreadsheet cells',n=len(delta),max_error=max(delta),pass_check=max(delta)<1e-9))
# Repeated no-NK ratios establish the row grouping for four batches, not invented donor identity.
batch_audit=[]
for endpoint in ['DR4_DR5_KO_relative_clearance','DR5_overexpression_relative_clearance']:
    q=don[(don.endpoint==endpoint)&(don.condition=='no_NK_ratio')]
    for b,z in q.groupby('batch'):
        batch_audit.append(dict(endpoint=endpoint,batch=int(b),donors=len(z),n_unique_reference=z.value.nunique(),reference=z.value.iloc[0]))
checks.append(dict(check='Each genetic-control batch shares one published no-NK reference for three donors',pass_check=all(x['donors']==3 and x['n_unique_reference']==1 for x in batch_audit)))
pd.DataFrame(batch_audit).to_csv(R/'analysis/batch_reference_audit.csv',index=False)
q=pd.read_csv(R/'source_data/EMBO_paired_contrasts.csv');s=pd.read_csv(R/'source_data/EMBO_effect_summary.csv')
re=[]
for row in s.itertuples():
    z=q[q.endpoint==row.endpoint]
    d=z.groupby('batch').difference.mean().to_numpy() if row.n_primary_units<row.n_donors else z.difference.to_numpy()
    p=np.mean([abs(np.dot(d,sg))>=abs(d.sum())-1e-10 for sg in itertools.product([-1,1],repeat=len(d))])
    re.append(dict(endpoint=row.endpoint,mean=d.mean(),n=len(d),p=p))
    checks.append(dict(check='Recomputed unit-aware mean and exact sign-flip '+row.endpoint,pass_check=bool(abs(d.mean()-row.estimate_mean)<1e-9 and abs(p-row.p_exact_signflip)<1e-9)))
# Independently reproduce primary leave-family-out ridge with normal equations.
df=pd.read_csv(R/'source_data/model_feature_target_join.csv');spec=json.loads((R/'source_data/model_specs.json').read_text());old=pd.read_csv(R/'source_data/out_of_fold_predictions.csv')
pred=[]
for name,cols in spec['models'].items():
    for fam in df.family.unique():
        a=df[df.family!=fam];b=df[df.family==fam];w=(1/a.family.map(a.family.value_counts())).to_numpy();y=a.target_pp.to_numpy();ym=np.average(y,weights=w)
        if not cols: val=np.full(len(b),ym)
        else:
            if cols==['lineage']:
                levels=sorted(a.lineage.unique());X=np.array([[int(v==k) for k in levels] for v in a.lineage],float);Z=np.array([[int(v==k) for k in levels] for v in b.lineage],float)
            else: X=a[cols].to_numpy(float);Z=b[cols].to_numpy(float)
            mu=np.average(X,axis=0,weights=w);sd=np.sqrt(np.average((X-mu)**2,axis=0,weights=w));sd[sd<1e-12]=1
            X=(X-mu)/sd;Z=(Z-mu)/sd
            beta=np.linalg.solve(X.T@(w[:,None]*X)+10*np.eye(X.shape[1]),X.T@(w*(y-ym)))
            val=ym+Z@beta
        for model,v in zip(b.model,val):pred.append(dict(model_name=name,tumor_model=model,predicted_pp=v))
v=pd.DataFrame(pred).merge(old.query("dataset=='all22' and split=='family' and alpha==10"),on=['model_name','tumor_model'],suffixes=('_new','_old'),validate='one_to_one')
maxerr=float(abs(v.predicted_pp_new-v.predicted_pp_old).max())
checks.append(dict(check='198 primary predictions from independent weighted matrix calculation',max_error=maxerr,pass_check=maxerr<1e-8))
v.to_csv(R/'qa/primary_prediction_recalculation.csv',index=False)
# Recompute the two operational effects directly from each four-arm source.
rows=[]
sh=pd.read_csv(R/'source_data/Sheard_Fig6_readings.csv');ht=pd.read_csv(R/'source_data/HT29_Fig3_readings.csv')
def add(study,model,context,rs,invert=False):
    if invert:
        vals=np.array([100-r.value for r in rs]);bounds=[(100-r.reading_upper,100-r.reading_lower) for r in rs]
    else: vals=np.array([r.value for r in rs]);bounds=[(r.reading_lower,r.reading_upper) for r in rs]
    a,b,c,d=vals;corners=np.array(list(itertools.product(*bounds)))
    net=corners[:,1]-corners[:,0];inter=corners[:,1]-corners[:,0]-corners[:,3]+corners[:,2]
    rows.append(dict(study=study,model=model,context=context,k00=a,k10=b,k01=c,k11=d,net_gain=b-a,net_lower=net.min(),net_upper=net.max(),blockade_sensitive_gain=b-a-d+c,interaction_lower=inter.min(),interaction_upper=inter.max()))
for model,z in sh[sh.panel=='6C'].groupby('model',sort=False):
    for dose in sorted(x for x in z.IFNg_ng_ml.unique() if x>0):
        rs=[z[(z.IFNg_ng_ml==dose0)&(z.anti_TRAIL==block)].iloc[0] for dose0,block in [(0,False),(dose,False),(0,True),(dose,True)]]
        add('Sheard 2013',model,f'{dose:g} ng/ml',rs,True)
for (act,et),z in ht.groupby(['effector_activation','ET'],sort=False):
    rs=[z[(z.IFNg==g)&(z.blockade==block)].iloc[0] for g,block in [(False,'none'),(True,'none'),(False,'anti_FAS'),(True,'anti_FAS')]]
    add('Mori 1997','HT-29',f'{act}; {et}:1',rs)
d=pd.DataFrame(rows);prior=pd.read_csv(R/'source_data/headroom_and_blockade_sensitivity.csv')
checks.append(dict(check='All 16 IFNg-alone four-arm estimates reproduce frozen outputs',pass_check=bool(np.allclose(d.net_gain,prior.net_gain) and np.allclose(d.blockade_sensitive_gain,prior.blockade_sensitive_gain))))
d.to_csv(R/'analysis/factorial_recalculation.csv',index=False)
pd.DataFrame(manifest).to_csv(R/'provenance/INPUT_MANIFEST.csv',index=False)
pd.DataFrame(checks).to_csv(R/'qa/numeric_checks.csv',index=False)
(R/'qa/environment.json').write_text(json.dumps({'python':sys.version,'platform':platform.platform(),'numpy':np.__version__,'pandas':pd.__version__},indent=2))
assert all(x['pass_check'] for x in checks),checks
print(json.dumps({'passed':len(checks),'copied_inputs':len(manifest),'prediction_max_difference':maxerr,'core_condition_means':len(sh)+len(ht)+12,'complete_factorial_contrasts':len(d)}))
