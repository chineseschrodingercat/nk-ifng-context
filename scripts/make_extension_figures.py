from pathlib import Path
import json,shutil
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from figure_quality import inspect_figure, panel_note
R=Path(__file__).resolve().parents[1];A=R/'analysis';D=R/'source_data';F=R/'figures'
font=font_manager.findfont('Arial',fallback_to_default=False)
plt.rcParams.update({'font.family':'Arial','font.size':9,'font.weight':'bold','axes.labelsize':10,'axes.labelweight':'bold','axes.titlesize':10,'axes.titleweight':'bold','axes.linewidth':1.3,'xtick.major.width':1.2,'ytick.major.width':1.2,'axes.grid':False,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
BLUE='#1764D6';ORANGE='#ED7420';GREEN='#009E73';PURPLE='#A338B7';GRAY='#77828F';BLACK='#202329'
manifest=[]
def style(ax,letter,title=None):
 ax.grid(False,which='both');ax.tick_params(direction='out',length=3.5)
 for s in ax.spines.values():s.set_visible(True);s.set_linewidth(1.3);s.set_color(BLACK)
 if letter:ax.text(-.14,1.04,letter,transform=ax.transAxes,fontsize=18,fontweight='bold',va='bottom')
 if title:ax.set_title(title,pad=11)
def key(ax,handles,loc='upper left',**kwargs):ax.legend(handles=handles,loc=loc,frameon=True,facecolor='white',edgecolor='none',framealpha=1,fontsize=8,handlelength=1.3,**kwargs)
def save(fig,name,tables):
 fig.canvas.draw();geo=[];quality=inspect_figure(fig)
 for ax in fig.axes:
  lg=ax.get_legend();bb=lg.get_window_extent(fig.canvas.get_renderer()) if lg else None;ab=ax.get_window_extent();p=ax.get_position()
  geo.append(dict(width_mm=p.width*fig.get_figwidth()*25.4,height_mm=p.height*fig.get_figheight()*25.4,grid=any(x.get_visible() for x in ax.get_xgridlines()+ax.get_ygridlines()),own_legend=lg is not None,legend_inside=bool(bb is not None and ab.contains(bb.x0,bb.y0) and ab.contains(bb.x1,bb.y1))))
 for ext in ['pdf','svg','png']:fig.savefig(F/f'{name}.{ext}',dpi=300)
 fig.savefig(R/'qa'/f'{name}_preview.png',dpi=145)
 manifest.append(dict(figure=name,font=font,width_inches=fig.get_figwidth(),height_inches=fig.get_figheight(),no_overall_title=fig._suptitle is None,panels=geo,quality_checks=quality,source_tables=tables));plt.close(fig)
def err(z,v,lo,hi):return np.array([z[v]-z[lo],z[hi]-z[v]])
colorkey=[Line2D([],[],marker='o',ls='',color=BLUE,label='Control'),Line2D([],[],marker='s',ls='',color=ORANGE,label='IFNγ')]
# F1 contains all eight ovarian arms and the two retained recognition factorials.
o=pd.read_csv(A/'Ovarian2026_Fig5_donor_readings.csv');g=pd.read_csv(A/'Ovarian2026_Fig5_group_means.csv');c=pd.read_csv(A/'Ovarian2026_Fig5_treatment_contrasts.csv')
products=['mock_WT','mock_KO','PRAME_WT','PRAME_KO'];labels=['Mock\nWT','Mock\nKO','PRAME\nWT','PRAME\nKO']
fig,axs=plt.subplots(2,2,figsize=(7.1,6.8));fig.subplots_adjust(left=.115,right=.975,bottom=.11,top=.91,wspace=.38,hspace=.60)
ax=axs[0,0]
for j,p in enumerate(products):
 for k,(cond,color,off) in enumerate([('control',BLUE,-.15),('IFNg',ORANGE,.15)]):
  z=o[(o.NK_product==p)&(o.condition==cond)];y=g[(g.NK_product==p)&(g.condition==cond)].iloc[0].mean_pp
  jitter=np.linspace(-.055,.055,len(z));ax.scatter(j+off+jitter,z.value_pp,s=17,c=color,marker='o' if cond=='control' else 's',edgecolor=BLACK,lw=.4,zorder=3)
  ax.plot([j+off-.10,j+off+.10],[y,y],lw=2.2,color=BLACK,zorder=4)
ax.set_xticks(range(4),labels,fontsize=8);ax.set_ylim(0,89);ax.set_ylabel('Cytotoxicity (%)');style(ax,'A','SK-OV-3  |  five NK donors');key(ax,colorkey,ncol=2)
ax=axs[0,1];z=c.set_index('NK_product').loc[products]
ax.bar(range(4),z.delta_pp,color=[GRAY,GREEN,PURPLE,ORANGE],edgecolor=BLACK,lw=.8,yerr=err(z,'delta_pp','delta_lower_reading_pp','delta_upper_reading_pp'),capsize=3)
ax.axhline(0,color=BLACK,lw=1);ax.set_xticks(range(4),labels,fontsize=8);ax.set_ylim(-13,20);ax.set_ylabel('IFNγ effect (pp)');style(ax,'B','Net response by NK product')
panel_note(ax,'Whiskers: reading bounds')
for ax,panel,letter,title,mods,labs in [(axs[1,0],'Fig2B','C','THP-1  |  inhibitory KIR',['IgG','GL183'],['Isotype','GL183']), (axs[1,1],'Fig4E','D','THP-1  |  ICAM-1',['Empty_vector','ICAM1_overexpression'],['Vector','ICAM-1↑'])]:
 z=pd.read_csv(A/'Wang2012_recognition_factorial_readings.csv');z=z[(z.panel==panel)&z.condition.isin(['control','IFNg'])]
 for cond,color,off,marker in [('control',BLUE,-.08,'o'),('IFNg',ORANGE,.08,'s')]:
  q=z[z.condition==cond].set_index('modifier').loc[mods]
  ax.errorbar(np.arange(2)+off,q.value_pp,yerr=err(q,'value_pp','lower_reading_pp','upper_reading_pp'),fmt=marker+'-',color=color,lw=1.7,ms=5,mec=BLACK,mew=.5,capsize=2)
 ax.set_xticks(range(2),labs);ax.set_ylim(0,96 if panel=='Fig2B' else 40);ax.set_xlim(-.35,1.35);ax.set_ylabel('Specific lysis (%)');style(ax,letter,title);key(ax,colorkey,ncol=2)
save(fig,'Figure1_Recognition_context',['Ovarian2026_Fig5_donor_readings.csv','Ovarian2026_Fig5_treatment_contrasts.csv','Wang2012_recognition_factorial_readings.csv'])
# F5 all target models remain visible; higher E:T ratios are labelled separately.
w=pd.read_csv(A/'Wang2012_Fig6_contrasts.csv');lo=w[w.comparison_set.isin(['common_1to1','SUP_primary_2to1'])].sort_values('delta_pp');hi=w[w.comparison_set=='higher_ET_sensitivity'].set_index('model')
fig=plt.figure(figsize=(7.1,6.6));gs=fig.add_gridspec(1,2,width_ratios=[1.45,1]);axs=[fig.add_subplot(gs[0,i]) for i in range(2)];fig.subplots_adjust(left=.18,right=.97,bottom=.13,top=.92,wspace=.57)
ax=axs[0]
for i,row in enumerate(lo.itertuples()):
 other=hi.loc[row.model]
 ax.errorbar(row.delta_pp,i-.11,xerr=[[row.delta_pp-row.delta_lower_reading_pp],[row.delta_upper_reading_pp-row.delta_pp]],fmt='o',color=BLUE,ms=4.5,mec=BLACK,mew=.4,capsize=2)
 ax.plot(other.delta_pp,i+.13,'s',color=ORANGE,ms=4.4,mec=BLACK,mew=.4)
ax.set_yticks(range(len(lo)),[f'{r.model} ({int(hi.loc[r.model,"ET"])}:1)' + ('*' if r.model=='SUP-B15' else '') for r in lo.itertuples()],fontsize=8)
ax.set_ylim(-.8,len(lo)+2.7);ax.set_xlim(-22,30);ax.axvline(0,color=GRAY,ls='--',lw=1);ax.set_xlabel('IFNγ-associated lysis change\n(percentage points)');style(ax,'A','Wang target-cell panel')
key(ax,[Line2D([],[],marker='o',ls='',color=BLUE,label='1:1 (*2:1)'),Line2D([],[],marker='s',ls='',color=ORANGE,label='Higher E:T in row label')])
ax=axs[1];p=pd.read_csv(D/'model_feature_target_join.csv');v=p[p.model=='Kasumi-1'].iloc[0];z=w[w.model=='Kasumi-1'].sort_values('ET')
ys=[z.iloc[0].delta_pp,z.iloc[1].delta_pp,v.target_pp];x=[0,1,2]
ax.scatter(x,ys,c=[BLUE,BLUE,ORANGE],s=65,edgecolor=BLACK,lw=.6,zorder=3)
for i,val in enumerate(ys):ax.text(i,val+1.5,f'{val:+.1f}',ha='center',fontsize=10)
ax.axhline(0,color=GRAY,lw=1,ls='--');ax.set_xticks(x,['Wang\n1:1','Wang\n4:1','Aquino\nSix E:T'],fontsize=8);ax.set_ylim(-24,26);ax.set_xlim(-.55,2.55);ax.set_ylabel('IFNγ-associated lysis change (pp)');style(ax,'B','Kasumi-1 across studies')
key(ax,[Line2D([],[],marker='o',ls='',color=BLUE,label='Wang'),Line2D([],[],marker='o',ls='',color=ORANGE,label='Aquino')])
save(fig,'Figure5_Functional_context_transfer',['Wang2012_Fig6_contrasts.csv','exact_model_cross_study_comparison.csv'])
# F6 ranks describe within-screen relative selection, not acute killing.
s=pd.read_csv(A/'Dufva_21gene_panel.csv');genes=['IFNGR1','IFNGR2','JAK1','JAK2','STAT1','IRF1','B2M','HLA-E','TAP1','ICAM1','CD58','FAS','TNFRSF10A','TNFRSF10B','CASP8','CFLAR','BAX','BAK1','BCL2','BCL2L1','MCL1'];screens=list(s.screen.unique());mat=s.pivot(index='gene',columns='screen',values='centered_percentile').loc[genes,screens]
labs=['K562','MOLM14','SUDHL4','NALM6','MM1S','LP1','KMS11\nKHYG1','MM1S','LP1\nKHYG1','KMS11','KMS11\nKHYG1']
assert all(' LOF' in name for name in screens[:7]) and all(' GOF' in name for name in screens[7:])
fig,ax=plt.subplots(figsize=(7.1,7.0));fig.subplots_adjust(left=.20,right=.80,bottom=.17,top=.90)
hm=ax.imshow(mat.values,aspect='auto',cmap='RdBu_r',vmin=-1,vmax=1,interpolation='none');ax.set_yticks(range(21),genes,fontsize=9);ax.set_xticks(range(11),labs,rotation=55,ha='right',fontsize=8)
ax.set_ylim(20.5,-4.2);ax.set_xlim(-.5,10.5);style(ax,None)
for left,right,label in [(-.5,6.5,'Loss of function'),(6.5,10.5,'Gain of function')]:
 x0,x1=left+.12,right-.12
 ax.plot([x0,x0,x1,x1],[-.75,-1.2,-1.2,-.75],color=BLACK,lw=1.3,clip_on=False)
 ax.text((left+right)/2,-1.7,label,ha='center',va='bottom',fontsize=9)
ax.plot([6.5,6.5],[-.5,20.5],color='white',lw=3,zorder=3)
ax.plot([6.5,6.5],[-.5,20.5],color=BLACK,lw=1.1,zorder=4)
key(ax,[Patch(facecolor='#67a9cf',label='Lower selection rank'),Patch(facecolor='#ef8a62',label='Higher selection rank')],ncol=2)
cbax=fig.add_axes([.83,.33,.024,.40]);cb=fig.colorbar(hm,cax=cbax)
cb.set_ticks([-1,0,1],labels=['−1','0 (median)','+1']);cb.ax.tick_params(labelsize=8)
cb.ax.yaxis.set_label_position('right');cb.set_label('Within-screen centered percentile',fontsize=9,labelpad=8)
cb.ax.grid(False,which='both');cb.outline.set_linewidth(1.3);cb.outline.set_edgecolor(BLACK)
save(fig,'Figure6_Genetic_context_map',['Dufva_21gene_panel.csv','Dufva_all_gene_results.csv.gz'])
# F7 separates no-NK background from a whole-host factorial.
b=pd.read_csv(D/'Breast_Fig2D_readings.csv');bc=pd.read_csv(D/'Breast_factorial_contrasts.csv');bn=pd.read_csv(D/'Breast_normalized_group_means.csv')
mouse=pd.read_csv(D/'Remsik2025_Fig5_mouse_source_cells.csv');mouse=mouse[mouse.panel=='Fig5i'];mouse.to_csv(A/'Remsik2025_Fig5i_mouse_values.csv',index=False)
fig,axs=plt.subplots(2,2,figsize=(7.1,6.6));fig.subplots_adjust(left=.12,right=.97,bottom=.12,top=.92,wspace=.44,hspace=.60)
ax=axs[0,0]
for cy,color,off in [(False,BLUE,-.17),(True,ORANGE,.17)]:
 z=b[(b.cytokines==cy)&b.arm.isin(['no_NK','NK'])].set_index('arm').loc[['no_NK','NK']]
 ax.bar(np.arange(2)+off,z.value,.30,color=color,edgecolor=BLACK,lw=.7,yerr=err(z,'value','reading_lower','reading_upper'),capsize=2)
ax.set_xticks([0,1],['No NK','With NK']);ax.set_ylim(0,109);ax.set_ylabel('AnnV+ or 7AAD+ targets (%)');style(ax,'A','Mixed cytokine background');key(ax,[Patch(facecolor=BLUE,label='Control'),Patch(facecolor=ORANGE,label='IFNγ + TNFα')])
ax=axs[0,1];z=pd.concat([bc[(bc.arm=='NK')&(bc.contrast=='cytokine_total_nonviability_change')],bc[(bc.arm=='NK')&(bc.contrast=='no_NK_adjusted_additive_change')]])
n=bn[bn.arm=='NK'].sort_values('cytokines');vv=n.iloc[1].normalized_NK_associated_nonviability-n.iloc[0].normalized_NK_associated_nonviability
ax.bar([0,1,2],[*z.estimate_pp,vv],color=[ORANGE,GREEN,PURPLE],edgecolor=BLACK,lw=.7);ax.set_xticks([0,1,2],['Total\nchange','Additive\nchange','Normalized\nchange'],fontsize=8);ax.set_ylim(0,32);ax.set_ylabel('Cytokine effect (pp)');style(ax,'B','Response definition');panel_note(ax,'Distinct response measures')
ax=axs[1,0];groups=[(0,0),(1,0),(0,1),(1,1)]
for j,(ifn,depl) in enumerate(groups):
 z=mouse[(mouse.ifng_overexpression==ifn)&(mouse.nk_depletion==depl)];vals=z.log10_radiance.to_numpy();color=ORANGE if ifn else BLUE
 ax.scatter(j+np.linspace(-.12,.12,len(z)),vals,c=color,s=26,edgecolor=BLACK,lw=.5);ax.plot([j-.20,j+.20],[vals.mean()]*2,color=BLACK,lw=2)
ax.set_xticks(range(4),['GFP\nIso','IFNγ\nIso','GFP\nDepl','IFNγ\nDepl'],fontsize=8);ax.set_ylim(2.7,7.8);ax.set_ylabel('Cranial radiance (log10)');style(ax,'C','Leptomeningeal mouse model');key(ax,[Line2D([],[],marker='o',ls='',color=BLUE,label='GFP'),Line2D([],[],marker='o',ls='',color=ORANGE,label='IFNγ')],ncol=2)
ax=axs[1,1];v=[]
for depl in [0,1]:
 z=mouse[mouse.nk_depletion==depl];v.append(z[z.ifng_overexpression==1].log10_radiance.mean()-z[z.ifng_overexpression==0].log10_radiance.mean())
ax.bar([0,1],v,color=[GREEN,PURPLE],edgecolor=BLACK,lw=.8);ax.axhline(0,color=BLACK,lw=1);ax.set_xticks([0,1],['Isotype','Depletion']);ax.set_ylim(-1.7,.48);ax.set_ylabel('IFNγ/GFP effect (log10 ratio)');style(ax,'D','Host perturbation changes effect');panel_note(ax,'Descriptive contrasts')
for j,val in enumerate(v):ax.text(j,val-.1,f'{10**val:.3f}×',ha='center',va='top')
save(fig,'Figure7_Background_and_host_context',['Breast_Fig2D_readings.csv','Breast_factorial_contrasts.csv','Breast_normalized_group_means.csv','Remsik2025_Fig5i_mouse_values.csv'])
(F/'EXTENSION_FIGURE_MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf8')
print(json.dumps({'new_figures':len(manifest),'panels':sum(len(x['panels']) for x in manifest),'font':font}))
