"""Decode the complete source heatmaps using their actual colour scales."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from PIL import Image,ImageDraw
R=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
markers=['ICAM1','FAS_CD95','DR4_CD261','DR5_CD262','HVEM_CD270','PDL1','MHC_I','HLA_E','HLA_DR','CD244_ligands','NKG2D_ligands','DNAM1_ligands','NKp30_ligands','NKp44_ligands','NKp46_ligands','NKp80_ligands']
models=['RD','Rh41','Rh30','BT-12','SJ-GBM2','CHLA-266','CHLA-9','CHLA-10','CHLA-258','TC-71','NB-1643','NB-EBc1','CHLA-90','CHLA-136','NALM-6','RS4;11','COG-LL-317','MOLT-4','CCRF-CEM','Kasumi-1','Karpas-299','Ramos-RA1']
indmarkers=['PDL1','ICAM1','HLA_DR','MHC_I','FAS_CD95','HVEM_CD270','DR4_CD261','NKp80_ligands','HLA_E','NKG2D_ligands','NKp44_ligands','CD244_ligands','NKp46_ligands','NKp30_ligands','DR5_CD262','DNAM1_ligands']
indmodels=['CHLA-9','MOLT-4','CHLA-10','NB-EBc1','Kasumi-1','Ramos-RA1','CHLA-136','Rh30','COG-LL-317','CHLA-258','Rh41','TC-71','CHLA-90','CCRF-CEM','RD','CHLA-266','NALM-6','RS4;11','Karpas-299','BT-12','SJ-GBM2','NB-1643']
specs=[dict(panel='F3_baseline',file='derived/S11/p6_Im0.jpg',x0=335,x1=1158,y0=214,y1=813,nx=22,ny=16,cols=models,rows=markers,palette_x0=665,palette_x1=813,palette_y=842,v0=0,v1=5.49,orientation='marker_rows'),
 dict(panel='F4B_induced',file='derived/S11/p7_Im0.jpg',x0=1315,x1=1959,y0=166,y1=1012,nx=16,ny=22,cols=indmarkers,rows=indmodels,palette_x0=1539,palette_x1=1743,palette_y=1069,v0=-3.72,v1=3.72,orientation='model_rows')]
out=[];diagnostics=[]
for s in specs:
 im=Image.open(R/s['file']).convert('RGB');a=np.array(im,dtype=float);overlay=im.copy();draw=ImageDraw.Draw(overlay)
 px=np.arange(s['palette_x0'],s['palette_x1']+1);pal=np.median(a[s['palette_y']-3:s['palette_y']+4,px],axis=0)
 vals=np.linspace(s['v0'],s['v1'],len(px)); np.savetxt(R/'derived/S11'/f"{s['panel']}_palette.csv",np.c_[px,vals,pal],delimiter=',',header='x,source_scale,red,green,blue',comments='')
 for ri,row in enumerate(s['rows']):
  for ci,col in enumerate(s['cols']):
   x=s['x0']+(ci+.5)*(s['x1']-s['x0'])/s['nx'];y=s['y0']+(ri+.5)*(s['y1']-s['y0'])/s['ny'];xx,yy=round(x),round(y)
   color=np.median(a[yy-4:yy+5,xx-4:xx+5].reshape(-1,3),axis=0)
   dist=np.sqrt(np.mean((pal-color)**2,axis=1));idx=int(dist.argmin());good=np.where(dist<=dist[idx]+5)[0]
   lo,hi=vals[good.min()],vals[good.max()]
   model,marker=(col,row) if s['orientation']=='marker_rows' else (row,col)
   out.append(dict(record_id=f"S11_{s['panel']}_{model}_{marker}",panel=s['panel'],model=model,marker=marker,value=vals[idx],reading_low=lo,reading_high=hi,
    rgb_r=color[0],rgb_g=color[1],rgb_b=color[2],nearest_palette_RMSE=dist[idx],x_px=x,y_px=y,source_image=s['file'],source_sha256=sha(R/s['file']),
    unit='author arcsinh heatmap scale; not absolute MFI or untransformed fold change',independent_n='not supplied for CyTOF; one source heatmap cell',raw_missing=False))
   draw.rectangle((x-5,y-5,x+5,y+5),outline=(255,60,60),width=2)
   if dist[idx]>12:diagnostics.append(dict(model=model,marker=marker,panel=s['panel'],RMSE=float(dist[idx])))
 overlay.save(R/'qa'/f"{s['panel']}_reading_overlay.png")
df=pd.DataFrame(out);assert len(df)==704 and not df.record_id.duplicated().any()
df.to_csv(R/'analysis/ligand_features_long.csv',index=False,encoding='utf-8-sig')
w=df.pivot(index='model',columns=['panel','marker'],values='value');w.columns=['__'.join(c) for c in w.columns];w.to_csv(R/'analysis/ligand_features_wide.csv',encoding='utf-8-sig')
(R/'analysis/heatmap_calibration.json').write_text(json.dumps(specs,indent=2),encoding='utf-8')
(R/'qa/heatmap_extraction_qc.json').write_text(json.dumps(dict(records=len(df),models=df.model.nunique(),markers=df.marker.nunique(),RMSE_median=float(df.nearest_palette_RMSE.median()),RMSE_max=float(df.nearest_palette_RMSE.max()),high_colour_mismatch=diagnostics,uncertainty='colour and image reading envelope; not biological CI'),indent=2),encoding='utf-8')
print(json.dumps({'records':len(df),'RMSE_median':float(df.nearest_palette_RMSE.median()),'RMSE_max':float(df.nearest_palette_RMSE.max()),'high_mismatch':diagnostics}))
