from pathlib import Path
import json,sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.legend_handler import HandlerTuple
from figure_quality import inspect_figure
R=Path(__file__).resolve().parents[1];D=R/'source_data';F=R/'figures'
font=font_manager.findfont('Arial',fallback_to_default=False)
plt.rcParams.update({'font.family':'Arial','font.size':9,'font.weight':'bold','axes.labelsize':10,'axes.labelweight':'bold','axes.titlesize':10,'axes.titleweight':'bold','axes.linewidth':1.3,'xtick.major.width':1.2,'ytick.major.width':1.2,'axes.grid':False,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
BLUE='#1764D6';ORANGE='#ED7420';GREEN='#009E73';PURPLE='#A338B7';GRAY='#77828F';BLACK='#202329'
manifest=[]
def style(ax,letter,title=None):
    ax.grid(False,which='both');ax.tick_params(direction='out',length=3.5)
    for spine in ax.spines.values():spine.set_visible(True);spine.set_linewidth(1.3);spine.set_color(BLACK)
    ax.text(-.13,1.04,letter,transform=ax.transAxes,fontsize=18,fontweight='bold',va='bottom')
    if title:ax.set_title(title,pad=11)
def legend(ax,handles,loc='upper right',**kwargs):
    ax.legend(handles=handles,loc=loc,frameon=True,facecolor='white',edgecolor='none',framealpha=1,fontsize=8,handlelength=kwargs.pop('handlelength',1.4),**kwargs)
def save(fig,name,tables):
    fig.canvas.draw();geo=[];quality=inspect_figure(fig)
    for ax in fig.axes:
        b=ax.get_position();lg=ax.get_legend();bb=lg.get_window_extent(fig.canvas.get_renderer()) if lg else None;ab=ax.get_window_extent()
        geo.append({'axes_width_mm':b.width*fig.get_figwidth()*25.4,'axes_height_mm':b.height*fig.get_figheight()*25.4,'four_spines':all(s.get_visible() for s in ax.spines.values()),'grid_visible':any(x.get_visible() for x in ax.get_xgridlines()+ax.get_ygridlines()),'own_legend':lg is not None,'legend_inside':bool(bb is not None and ab.contains(bb.x0,bb.y0) and ab.contains(bb.x1,bb.y1))})
    for ext in ['pdf','svg','png']:fig.savefig(F/f'{name}.{ext}',dpi=300)
    fig.savefig(R/'qa'/f'{name}_preview.png',dpi=145)
    manifest.append({'figure':name,'width_inches':fig.get_figwidth(),'height_inches':fig.get_figheight(),'no_overall_title':fig._suptitle is None,'font':font,'panels':geo,'quality_checks':quality,'source_tables':tables})
    plt.close(fig)
def errors(z,v,lo,hi):return np.vstack([z[v]-z[lo],z[hi]-z[v]])

# Figure 1: all models stay visible, including two interval-constrained means.
m=pd.read_csv(D/'model_feature_target_join.csv').sort_values('target_pp');don=pd.read_csv(D/'S11_donor_values.csv')
don['delta_lysis_pp']=pd.to_numeric(don.delta_lysis_pp,errors='coerce')
met=pd.read_csv(D/'prediction_metrics.csv')
fig=plt.figure(figsize=(7.1,6.9));gs=fig.add_gridspec(1,2,width_ratios=[1,1.03]);axs=[fig.add_subplot(gs[0,i]) for i in range(2)]
fig.subplots_adjust(left=.17,right=.97,bottom=.12,top=.91,wspace=.66)
ax=axs[0]
for i,row in enumerate(m.itertuples()):
    z=don[don.model==row.model]
    for j,r in enumerate(z[z.delta_lysis_pp.notna()].itertuples()):
        ax.plot(r.delta_lysis_pp,i+(j-1.5)*.12,'o',color=BLUE,ms=3.8,mec=BLACK,mew=.4)
    if row.n_identified_donor_points<4:
        ax.errorbar(row.target_pp,i,xerr=[[row.target_pp-row.mean_lower_reading_pp],[row.mean_upper_reading_pp-row.target_pp]],fmt='s',ms=5,color=ORANGE,mfc='white',mew=1.1,capsize=2,lw=1.3)
    else:ax.plot(row.target_pp,i,'D',ms=4.5,color=BLACK)
ax.set_yticks(range(len(m)),m.model,fontsize=8);ax.set_ylim(-1,len(m)+3.7);ax.set_xlim(-43,31);ax.axvline(0,color=GRAY,ls='--',lw=1)
ax.set_xlabel('IFNγ-associated lysis change\n(percentage points)');style(ax,'A')
legend(ax,[Line2D([],[],marker='o',ls='',color=BLUE,label='Donor value'),Line2D([],[],marker='D',ls='',color=BLACK,label='Complete mean'),Line2D([],[],marker='s',ls='',color=ORANGE,mfc='white',label='Interval mean')],loc='upper left')
ax=axs[1];names=['mean','lineage','baseline_recognition','induced_recognition','baseline_death_receptors','induced_recognition_plus_baseline_death','induced_recognition_plus_induced_death','all_baseline_ligands','all_induced_ligands']
labels=['Fold mean','Lineage','Baseline\nrecognition','Induced\nrecognition','Baseline\ndeath receptors','Induced recognition\n+ baseline death','Induced recognition\n+ induced death','All baseline\nligands','All induced\nligands']
z=met.query("dataset=='all22' and split=='family' and alpha==10").set_index('model_name').loc[names]
y=np.arange(9);ax.scatter(z.RMSE_family_pp,y,c=[BLACK]+[PURPLE]*8,s=45,edgecolor=BLACK,lw=.6)
for i,v in enumerate(z.RMSE_family_pp):ax.text(v+.12,i,f'{v:.2f}',va='center',fontsize=8)
ax.axvline(z.iloc[0].RMSE_family_pp,color=GRAY,ls='--',lw=1)
ax.set_yticks(y,labels,fontsize=8);ax.set_ylim(8.7,-2.2);ax.set_xlim(10.3,14.7);ax.set_xlabel('Family-held-out RMSE\n(percentage points)');style(ax,'B')
legend(ax,[Line2D([],[],marker='o',ls='',color=PURPLE,label='Fixed model'),Line2D([],[],ls='--',color=GRAY,marker='o',mfc=BLACK,mec=BLACK,label='Mean baseline')],loc='upper left')
save(fig,'Figure4_Response_and_prediction',['model_feature_target_join.csv','S11_donor_values.csv','prediction_metrics.csv'])

# Figure 2: operational net and interaction effects on the same source scale.
d=pd.read_csv(R/'analysis/factorial_recalculation.csv')
fig,axs=plt.subplots(2,1,figsize=(7.1,7.0),gridspec_kw={'height_ratios':[1,1.55]});fig.subplots_adjust(left=.23,right=.96,bottom=.105,top=.94,hspace=.43)
for ax,z,letter,title in [(axs[0],d[(d.study=='Sheard 2013')&(d.context=='10 ng/ml')],'A','Neuroblastoma  |  anti-TRAIL'),(axs[1],d[d.study=='Mori 1997'],'B','HT-29  |  anti-Fas')]:
    z=z.reset_index(drop=True);y=np.arange(len(z))
    for col,lo,hi,off,c,mark in [('net_gain','net_lower','net_upper',-.15,GREEN,'o'),('blockade_sensitive_gain','interaction_lower','interaction_upper',.15,PURPLE,'s')]:
        ax.errorbar(z[col],y+off,xerr=errors(z,col,lo,hi),fmt=mark,color=c,ms=5,mec=BLACK,mew=.5,capsize=2,elinewidth=1.3)
    labs=z.model if letter=='A' else z.context.str.replace('IL2_IL12','IL-2 + IL-12',regex=False).str.replace('IL2','IL-2',regex=False).str.replace('IFNa','IFNα',regex=False).str.replace('; ','; E:T ',regex=False)
    ax.set_yticks(y,labs,fontsize=9);ax.set_ylim(len(z)-.4,-1.8);ax.set_xlim(-13,31);ax.axvline(0,color=GRAY,ls='--',lw=1);style(ax,letter,title)
    legend(ax,[Line2D([],[],marker='o',color=GREEN,label='Net treatment effect'),Line2D([],[],marker='s',color=PURPLE,label='Blockade-sensitive effect')],loc='upper left',ncol=2)
    ax.set_xlabel('Change (percentage points)')
save(fig,'Figure2_Death_ligand_effects',['factorial_recalculation.csv'])

# Figure 4: direct published donor rows, not values inferred from error bars.
e=pd.read_csv(D/'EMBO_individual_source_values.csv');ec=pd.read_csv(D/'EMBO_paired_contrasts.csv');es=pd.read_csv(D/'EMBO_effect_summary.csv')
fig,axs=plt.subplots(2,2,figsize=(7.1,6.3));fig.subplots_adjust(left=.12,right=.97,bottom=.12,top=.92,wspace=.43,hspace=.53)
for ax,endpoint,letter,ylabel,title in [(axs[0,0],'soluble_block_CD107a','A','CD107a+ NK cells (%)','Soluble antibody'),(axs[0,1],'soluble_block_IFNg','B','IFNγ (pg/ml)','Soluble antibody'),(axs[1,0],'plate_bound_IFNg','C','IFNγ (pg/ml)','Immobilized antibody')]:
    z=e[e.endpoint==endpoint].pivot(index='donor_within_panel',columns='condition',values='value')
    for i,row in z.iterrows():ax.plot([0,1],[row['isotype'],row['anti_TRAIL']],color=GRAY,lw=1,alpha=.7);ax.scatter([0,1],[row['isotype'],row['anti_TRAIL']],c=[BLUE,ORANGE],s=28,edgecolor=BLACK,lw=.5,zorder=3)
    top=z[['isotype','anti_TRAIL']].max().max();ax.set_ylim(0,top*1.36);ax.set_xlim(-.35,1.35);ax.set_xticks([0,1],['Isotype','Anti-TRAIL']);ax.set_ylabel(ylabel);style(ax,letter,title)
    legend(ax,[Line2D([],[],color=GRAY,lw=1,label=f'Paired donors (n = {len(z)})')],loc='upper left')
ax=axs[1,1];z=ec[ec.endpoint=='soluble_block_remaining_targets'];cs=[BLUE,ORANGE,GREEN,PURPLE]
for b0,q in z.groupby('batch'):
    x=float(b0);ax.scatter(x+np.linspace(-.13,.13,len(q)),q.difference,color=cs[int(b0)-1],s=28,edgecolor=BLACK,lw=.5)
    ax.plot(x,q.difference.mean(),'D',color=BLACK,ms=6)
ax.axhline(0,color=GRAY,ls='--',lw=1);ax.set_xlim(.5,4.5);ax.set_ylim(-3,20);ax.set_xticks([1,2,3,4]);ax.set_xlabel('Shared-control batch');ax.set_ylabel('Remaining target difference (pp)');style(ax,'D','Soluble antibody')
donor_key=tuple(Line2D([],[],marker='o',ls='',color=c,ms=4,mec=BLACK,mew=.4) for c in cs)
legend(ax,[donor_key,Line2D([],[],marker='D',ls='',color=BLACK)],labels=['Donor','Batch mean'],handler_map={tuple:HandlerTuple(ndivide=None,pad=.2)},handlelength=3,loc='upper left',ncol=2)
save(fig,'Figure3_Effector_readouts',['EMBO_individual_source_values.csv','EMBO_paired_contrasts.csv','EMBO_effect_summary.csv'])


(F/'RETAINED_FIGURE_MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf8')
assert all(p['four_spines'] and not p['grid_visible'] and p['own_legend'] and p['legend_inside'] for f in manifest for p in f['panels'])
print('Three retained figures rebuilt.')
