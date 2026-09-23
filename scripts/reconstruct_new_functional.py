from pathlib import Path
import pandas as pd, numpy as np, hashlib, json
from PIL import Image,ImageDraw
R=Path(__file__).resolve().parents[1]
for d in ['analysis','qa','provenance']:(R/d).mkdir(exist_ok=True)
# Each tuple fixes source-specific y calibration, E:T positions, control/IFNg marker centres.
# Image-reading bounds are deterministic measurement bounds, never biological confidence intervals.
specs=[
 ('THP-1',20,70,397,82,[1,4],[153,209.5],[367,356],[280,197]),
 ('HL-60',0,50,397,82,[1,4],[343.5,398.5],[369.5,343],[298,260]),
 ('Kasumi-1',0,50,397,82,[1,4],[529,585.5],[297,263.5],[214,187]),
 ('GDM-1',0,50,397,82,[1,4],[721,775],[380.5,364.5],[340.5,318]),
 ('NCI-H211',0,50,397,82,[1,4],[910,963.5],[380.5,362.5],[347.5,317]),
 ('T47D',0,50,397,82,[1,4],[1098,1152],[379,374],[328.5,322]),
 ('H9',30,80,870,554,[1,3],[157,211],[839.5,651.5],[825,620]),
 ('HeLa',10,70,870,554,[1,4],[346,400],[822,782.5],[808,718]),
 ('K562',25,85,870,560,[1,4],[537,590.5],[852,637],[827,626.5]),
 ('721.221',30,80,870,554,[1,5],[726.5,773.5],[614,593.5],[626,568.5]),
 ('Raji',30,80,870,554,[1,3],[917,970.5],[769,578],[756,588.5]),
 ('Daudi',30,80,870,554,[1,4],[1115,1162],[829,734.5],[799,748]),
 ('293T',0,50,1334,1026,[1,4],[157,211],[1271,1211],[1281.5,1213]),
 ('Jurkat',30,80,1332,1023,[1,5],[350,404],[1123,1092.5],[1218.5,1107.5]),
 ('HH',0,50,1331,1023,[1,4],[540,598],[1185,1106.5],[1250.5,1169]),
 ('MDA-MB-231',0,50,1330,1023,[1,4],[733,793.5],[1320,1320],[1324.5,1324]),
 ('SUP-B15',-10,50,1395,1017,[2,4],[911,963.5],[1356.5,1338],[1331,1337]),
]
f=R/'raw/Wang2012/Wang2012_F6_original1200.png';sha=hashlib.sha256(f.read_bytes()).hexdigest();im=Image.open(f).convert('RGB');draw=ImageDraw.Draw(im)
rows=[]
for model,v0,v1,y0,y1,ets,xs,controls,treated in specs:
 scale=(v1-v0)/(y0-y1)
 for j,et in enumerate(ets):
  for tr,ys in [('control',controls),('IFNg',treated)]:
   y=ys[j];bound=4.0 if model in ['K562','721.221','Raji','293T','MDA-MB-231','SUP-B15'] else 3.0
   value=v0+(y0-y)*scale
   rows.append(dict(study='Wang2012',panel='Fig6',model=model,effector='NK-92MI',ET=et,condition=tr,value_pp=value,lower_reading_pp=value-bound*scale,upper_reading_pp=value+bound*scale,x_pixel=xs[j],y_pixel=y,pixel_bound=bound,y_low_pixel=y0,y_high_pixel=y1,axis_low=v0,axis_high=v1,source_file=str(f.relative_to(R)),source_sha256=sha,unit='displayed mean of technical triplicates; representative experiment',biological_n_recovered=0))
   draw.ellipse((xs[j]-4,y-4,xs[j]+4,y+4),outline='red' if tr=='control' else 'blue',width=2)
im.save(R/'qa/Wang2012_digitization_overlay.png')
df=pd.DataFrame(rows);df.to_csv(R/'analysis/Wang2012_Fig6_readings.csv',index=False)
contr=[]
for (m,et),d in df.groupby(['model','ET'],sort=False):
 c=d[d.condition=='control'].iloc[0];t=d[d.condition=='IFNg'].iloc[0]
 contr.append(dict(model=m,ET=et,control_pp=c.value_pp,IFNg_pp=t.value_pp,delta_pp=t.value_pp-c.value_pp,delta_lower_reading_pp=t.lower_reading_pp-c.upper_reading_pp,delta_upper_reading_pp=t.upper_reading_pp-c.lower_reading_pp,comparison_set='common_1to1' if et==1 else 'SUP_primary_2to1' if m=='SUP-B15' and et==2 else 'higher_ET_sensitivity'))
pd.DataFrame(contr).to_csv(R/'analysis/Wang2012_Fig6_contrasts.csv',index=False)
# Two source-complete IFNg factorials; other displayed cytokine arms retained but not mixed.
fact=[]
for panel,fn,y0,y1,top,xs,ys,modnames,other in [
 ('Fig2B','Wang2012_F2_original.png',970,772,75,[149,186,222,260,297,334],[898,868,857,798,789,797],['IgG','GL183'],'IFNg_TNFa'),
 ('Fig4E','Wang2012_F4_original.png',219,32,40,[790,818,846,875,903,931],[183,127,89,106,104,80],['Empty_vector','ICAM1_overexpression'],'TNFa')]:
 src=R/'raw/Wang2012'/fn;sh=hashlib.sha256(src.read_bytes()).hexdigest()
 for k,(x,y) in enumerate(zip(xs,ys)):
  val=(y0-y)*top/(y0-y1);err=2*top/(y0-y1)
  fact.append(dict(study='Wang2012',panel=panel,model='THP-1',modifier=modnames[k//3],condition=['control','IFNg',other][k%3],value_pp=val,lower_reading_pp=val-err,upper_reading_pp=val+err,x_pixel=x,y_pixel=y,y_zero_pixel=y0,y_high_pixel=y1,axis_max=top,source_file=str(src.relative_to(R)),source_sha256=sh,unit='representative published mean; biological paired values unavailable'))
pd.DataFrame(fact).to_csv(R/'analysis/Wang2012_recognition_factorial_readings.csv',index=False)
fcontr=[]
for panel in ['Fig2B','Fig4E']:
 d=pd.DataFrame(fact);d=d[(d.panel==panel)&d.condition.isin(['control','IFNg'])]
 for mod,g in d.groupby('modifier',sort=False):
  c=g[g.condition=='control'].iloc[0];t=g[g.condition=='IFNg'].iloc[0]
  fcontr.append(dict(panel=panel,modifier=mod,control_pp=c.value_pp,IFNg_pp=t.value_pp,delta_pp=t.value_pp-c.value_pp,delta_lower_reading_pp=t.lower_reading_pp-c.upper_reading_pp,delta_upper_reading_pp=t.upper_reading_pp-c.lower_reading_pp))
pd.DataFrame(fcontr).to_csv(R/'analysis/Wang2012_recognition_factorial_contrasts.csv',index=False)
# Ovarian2026: all 40 displayed donor markers in panels A/B, no inferred cross-panel donor identity.
# Anonymous shape records describe graphical observations only. Analyses use unpaired group means.
ospec={
 ('control','mock_WT'):[('v',117,787),('D',127,724),('o',147,724),('s',137,779),('^',158,786)],
 ('control','mock_KO'):[('v',203,720),('D',203,674),('o',203,783),('s',203,767),('^',203,693)],
 ('control','PRAME_WT'):[('v',268,751),('D',268,700),('o',268,733),('s',289,754),('^',247,756)],
 ('control','PRAME_KO'):[('v',347,705),('D',333,594),('o',347,745),('s',320,714),('^',320,751)],
 ('IFNg','mock_WT'):[('v',503,827),('D',489,796),('o',503,765),('s',517,789),('^',503,875)],
 ('IFNg','mock_KO'):[('v',587,640),('D',558,637),('o',573,810),('s',587,756),('^',558,743)],
 ('IFNg','PRAME_WT'):[('v',620,785),('D',620,773),('o',635,771),('s',662,774),('^',648,782)],
 ('IFNg','PRAME_KO'):[('v',696,641),('D',710,498),('o',710,693),('s',724,630),('^',710,719)]}
src=R/'raw/vanHees2026/Ovarian2026_Fig5_full.png';sha=hashlib.sha256(src.read_bytes()).hexdigest();im=Image.open(src).convert('RGB');draw=ImageDraw.Draw(im);ov=[]
for (cond,product),pts in ospec.items():
 y0,y1=(886,413) if cond=='control' else (883,411)
 for shape,x,y in pts:
  val=(y0-y)*80/(y0-y1);err=2*80/(y0-y1)
  ov.append(dict(study='vanHees2026',panel='Fig5A' if cond=='control' else 'Fig5B',model='SK-OV-3',condition=cond,NK_product=product,shape_id=shape,value_pp=val,lower_reading_pp=val-err,upper_reading_pp=val+err,x_pixel=x,y_pixel=y,y_zero_pixel=y0,y_high_pixel=y1,axis_max=80,pixel_bound=2,source_file=str(src.relative_to(R)),source_sha256=sha,unit='graphical donor observation',donor_id='not supplied',cross_panel_pairing='not assumed'))
  draw.ellipse((x-4,y-4,x+4,y+4),outline='blue',width=2)
im.save(R/'qa/Ovarian2026_digitization_overlay.png');od=pd.DataFrame(ov);od.to_csv(R/'analysis/Ovarian2026_Fig5_donor_readings.csv',index=False)
og=od.groupby(['NK_product','condition'],sort=False).agg(mean_pp=('value_pp','mean'),sd_pp=('value_pp','std'),n=('value_pp','size'),mean_lower_reading_pp=('lower_reading_pp','mean'),mean_upper_reading_pp=('upper_reading_pp','mean')).reset_index();og.to_csv(R/'analysis/Ovarian2026_Fig5_group_means.csv',index=False)
oc=[]
for product,g in og.groupby('NK_product',sort=False):
 c=g[g.condition=='control'].iloc[0];t=g[g.condition=='IFNg'].iloc[0]
 oc.append(dict(NK_product=product,control_pp=c.mean_pp,IFNg_pp=t.mean_pp,delta_pp=t.mean_pp-c.mean_pp,delta_lower_reading_pp=t.mean_lower_reading_pp-c.mean_upper_reading_pp,delta_upper_reading_pp=t.mean_upper_reading_pp-c.mean_lower_reading_pp))
oc=pd.DataFrame(oc);oc.to_csv(R/'analysis/Ovarian2026_Fig5_treatment_contrasts.csv',index=False)
ints=[]
for prod in ['mock','PRAME']:
 w=oc[oc.NK_product==prod+'_WT'].iloc[0];k=oc[oc.NK_product==prod+'_KO'].iloc[0]
 ints.append(dict(product=prod,contrast='IFNg effect in KO minus IFNg effect in WT',interaction_pp=k.delta_pp-w.delta_pp,lower_reading_pp=k.delta_lower_reading_pp-w.delta_upper_reading_pp,upper_reading_pp=k.delta_upper_reading_pp-w.delta_lower_reading_pp,inferential_p='not calculated; cross-panel donor identity unavailable'))
pd.DataFrame(ints).to_csv(R/'analysis/Ovarian2026_NKG2A_interactions.csv',index=False)
# Exact model identity only; cross-study context comparison, not paired samples or effector causality.
prior=pd.read_csv(R/'source_data/model_feature_target_join.csv')
wc=pd.DataFrame(contr);over=wc[wc.ET==1].merge(prior[['model','target_pp','mean_lower_reading_pp','mean_upper_reading_pp']],on='model',how='inner')
over.to_csv(R/'analysis/exact_model_cross_study_comparison.csv',index=False)
primary=wc[wc.ET==1]
summary=dict(Wang_means=len(df),Wang_models=len(specs),Wang_common_ET_models=len(primary),Wang_positive_reading_bounds=int((primary.delta_lower_reading_pp>0).sum()),Wang_negative_reading_bounds=int((primary.delta_upper_reading_pp<0).sum()),Wang_crossing_zero=int(((primary.delta_lower_reading_pp<=0)&(primary.delta_upper_reading_pp>=0)).sum()),Wang_primary_ET_rule='16 models at1:1; SUP-B15 separate2:1; retain model-specific higher ET sensitivity',Ovarian_donor_points=len(od),Ovarian_arms=len(og),Ovarian_points_per_arm=5,exact_model_overlap=list(over.model),Wang_extra_factorial_means=len(fact))
(R/'analysis/new_functional_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf8')
print(json.dumps(summary,indent=2));print(oc.to_string(index=False));print(pd.DataFrame(ints).to_string(index=False));print(over[['model','delta_pp','target_pp']].to_string(index=False))
