"""Traceable figures of actual held-out predictions and molecular diagnostics."""
from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
R=Path(__file__).resolve().parents[1];O=R/'figures';O.mkdir(exist_ok=True)
font=font_manager.findfont(font_manager.FontProperties(family='Arial'),fallback_to_default=False)
plt.rcParams.update({'font.family':'Arial','font.size':12,'font.weight':'bold','axes.labelweight':'bold','axes.titleweight':'bold','axes.titlesize':15,'axes.linewidth':1.4,'axes.grid':False,'xtick.major.width':1.4,'ytick.major.width':1.4,'xtick.major.size':5,'ytick.major.size':5,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
palette={'Ewing sarcoma':('#0072B2','o'),'leukemia':('#D55E00','s'),'lymphoma':('#CC79A7','D'),'neuroblastoma':('#009E73','^'),'rhabdomyosarcoma':('#E69F00','v'),'brain tumor':('#6C42B4','P')}
def style(ax):
 ax.grid(False,which='both')
 for s in ax.spines.values():s.set_visible(True);s.set_linewidth(1.4);s.set_color('#171717')
 ax.tick_params(labelsize=11)
def letter(fig,x,y,t):fig.text(x,y,t,fontsize=25,fontweight='bold',va='top')
records=[]
def save(fig,name,axes):
 fig.canvas.draw()
 for ext in ['pdf','svg','png']:fig.savefig(O/f'{name}.{ext}',dpi=300)
 records.append({'figure':name,'width_inches':fig.get_size_inches()[0],'height_inches':fig.get_size_inches()[1],'axes_inches':[{'width':a.get_position().width*fig.get_size_inches()[0],'height':a.get_position().height*fig.get_size_inches()[1]} for a in axes],'font':font,'grid':False,'all_spines':True})
 plt.close(fig)

pred=pd.read_csv(R/'analysis/out_of_fold_predictions.csv')
p=pred[(pred.dataset=='all22')&(pred.split=='family')&(pred.alpha==10)&(pred.model_name=='induced_recognition')]
met=pd.read_csv(R/'analysis/prediction_metrics.csv');m=met[(met.dataset=='all22')&(met.split=='family')&(met.alpha==10)].copy()
p.to_csv(O/'Figure1A_source.csv',index=False);m.to_csv(O/'Figure1B_source.csv',index=False)
fig=plt.figure(figsize=(14.5,7.8));a=fig.add_axes([.075,.30,.31,.565]);b=fig.add_axes([.715,.23,.255,.635])
for lineage,(color,marker) in palette.items():
 z=p[p.lineage==lineage];a.scatter(z.observed_pp,z.predicted_pp,s=86,c=color,marker=marker,edgecolors='#171717',linewidths=.9,zorder=3)
a.plot([-35,20],[-35,20],color='#707070',linestyle='--',linewidth=1.5,zorder=1)
a.set(xlim=(-35,20),ylim=(-35,20),xlabel='Observed IFNγ effect (pp)',ylabel='Held-out predicted effect (pp)',title='Induced recognition: 22 models')
handles=[Line2D([0],[0],marker=mark,color='none',markerfacecolor=c,markeredgecolor='#171717',markersize=8,label=k.capitalize()) for k,(c,mark) in palette.items()]
fig.legend(handles=handles,loc='lower left',bbox_to_anchor=(.035,.065),ncol=2,frameon=False,fontsize=11,columnspacing=1.4,handletextpad=.5)
labels=['Training mean','Lineage','Baseline recognition','Induced recognition','Baseline death receptors','Induced recognition\n+ baseline death receptors','Induced recognition\n+ induced death receptors','All 16 baseline ligands','All 16 induced ligands']
ys=np.arange(len(m));v=m.RMSE_family_pp.to_numpy();colors=['#555555']+['#0072B2']*8
b.scatter(v,ys,s=90,c=colors,edgecolor='#171717',linewidth=.8,zorder=3)
for yy,vv in zip(ys,v):b.annotate(f'{vv:.2f}',(vv,yy),xytext=(9,0),textcoords='offset points',va='center',fontsize=11)
b.axvline(v[0],color='#555555',linestyle='--',linewidth=1.5)
b.set(xlim=(0,17),ylim=(8.6,-.6),yticks=ys,yticklabels=labels,xlabel='Held-out RMSE (pp)',title='Family-held-out error')
b.tick_params(axis='y',labelsize=10.5,pad=8)
for ax in [a,b]:style(ax)
letter(fig,.018,.965,'A');letter(fig,.430,.965,'B')
save(fig,'Figure1_grouped_prediction',[a,b])

rna=pd.read_csv(R/'analysis/RNA_matched_feature_target_join.csv')
rm=pd.read_csv(R/'analysis/RNA_prediction_metrics.csv');rm=rm[(rm.dataset=='RNA_exact_match')&(rm.alpha==10)].copy()
rna[['model','lineage','BCL2A1']].to_csv(O/'Figure2A_source.csv',index=False);rm.to_csv(O/'Figure2B_source.csv',index=False)
fig=plt.figure(figsize=(14.5,8.4));a=fig.add_axes([.12,.17,.255,.67]);b=fig.add_axes([.715,.17,.255,.67])
for y,(_,r) in enumerate(rna.iterrows()):
 c,mk=palette[r.lineage];a.scatter(r.BCL2A1,y,s=85,c=c,marker=mk,edgecolor='#171717',linewidth=.9,zorder=3)
a.set(yticks=np.arange(len(rna)),yticklabels=rna.model.tolist(),ylim=(len(rna)-.4,-.6),xlabel='BCL2A1 expression\n[log2(TPM + 1)]',title='RNA coverage: 15 models')
a.set_xlim(-.2,max(rna.BCL2A1)*1.16)
order=['mean','lineage','RNA_two_modules','induced_recognition','induced_recognition_plus_RNA','RNA_bounded_z5','induced_recognition_plus_RNA_bounded_z5','RNA_log_balance','induced_recognition_plus_RNA_log_balance']
labels=['Training mean','Lineage','RNA modules: original z','Induced recognition','Recognition + RNA:\noriginal z','RNA modules: bounded z','Recognition + RNA:\nbounded z','RNA modules: log balance','Recognition + RNA:\nlog balance']
for split,col,marker,dy in [('family','#0072B2','o',-.14),('lineage','#D55E00','D',.14)]:
 z=rm[rm.split==split].set_index('model_name').loc[order]
 b.scatter(z.RMSE_family_pp,np.arange(9)+dy,s=77,c=col,marker=marker,edgecolor='#171717',linewidth=.8,zorder=3)
 if split=='family':
  for yy,vv in zip(np.arange(9),z.RMSE_family_pp):b.annotate(f'{vv:.2f}',(vv,yy-.14),xytext=(9,0),textcoords='offset points',va='center',fontsize=10.5)
b.set(yticks=np.arange(9),yticklabels=labels,ylim=(8.6,-.6),xlim=(0,54),xlabel='Held-out RMSE (pp)',title='Matched-subset prediction')
b.tick_params(axis='y',pad=8)
fig.legend(handles=[Line2D([0],[0],marker='o',color='none',markerfacecolor='#0072B2',markeredgecolor='#171717',markersize=8,label='Family holdout'),Line2D([0],[0],marker='D',color='none',markerfacecolor='#D55E00',markeredgecolor='#171717',markersize=8,label='Lineage holdout')],loc='lower center',bbox_to_anchor=(.70,.035),ncol=2,frameon=False,fontsize=11)
for ax in [a,b]:style(ax)
letter(fig,.018,.945,'A');letter(fig,.430,.945,'B')
save(fig,'Figure2_independent_RNA',[a,b])
(R/'qa/figure_style_manifest.json').write_text(json.dumps(records,indent=2),encoding='utf8')
print(json.dumps(records,indent=2))
