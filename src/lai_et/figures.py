"""Publication-friendly black-label plots from reproducible result tables."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from .io import save

def _setup():
    plt.rcParams.update({'font.family':'DejaVu Serif', 'text.color':'black',
                         'axes.labelcolor':'black','xtick.color':'black','ytick.color':'black',
                         'axes.titlecolor':'black','savefig.dpi':300})

def gradient_figure(gradient,file):
    _setup()
    g=gradient.copy()
    if g.empty:return
    ais=[c for c in ['AI_A','AI_B','AI_C'] if c in set(g.ai)]
    if not ais:return
    metrics=[(k,f'{source}_{k}_local_ICE') for k in ['ET','T','Ei'] for source in []]
    fig,axs=plt.subplots(1,3,figsize=(11,3.5),sharex=True,constrained_layout=True)
    for j,key in enumerate(['ET','T','Ei']):
        ax=axs[j]
        for source,clr in [('PML','#176a85'),('GLEAM','#d35f39')]:
            v=g[g.metric==f'{source}_{key}_local_ICE'].set_index('ai').reindex(ais)
            if v.empty or v.rho.isna().all():continue
            yy=v.rho.to_numpy(float)
            bot=np.maximum(0,yy-v.lo95.to_numpy(float))
            top=np.maximum(0,v.hi95.to_numpy(float)-yy)
            ax.errorbar(np.arange(len(ais)),yy,yerr=[bot,top],color=clr,marker='o',
                        capsize=3,label=source)
        ax.axhline(0,color='0.6',ls='--',lw=.8)
        ax.set_xticks(range(len(ais)),ais)
        ax.set_title(['(a) Total ET','(b) Transpiration','(c) Interception'][j],loc='left')
        ax.spines[['top','right']].set_visible(False)
    axs[0].set_ylabel('10° block-median Spearman $\\rho$')
    h,l=axs[-1].get_legend_handles_labels()
    if h:fig.legend(h,l,ncol=2,loc='lower center',bbox_to_anchor=(.5,-.09),frameon=False)
    file=Path(file);file.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(file,bbox_inches='tight');plt.close(fig)

def class_heatmap(pft,ai='AI_A',endpoint='ET',file='outputs/pft_heatmap.png'):
    _setup()
    metric=f'GLEAM_{endpoint}_local_ICE'
    d=pft.dropna(subset=[ai,metric,'PFT']).copy()
    if d.empty:return
    d['zone']=pd.cut(d[ai],[-np.inf,.2,.5,.65,np.inf],
                     labels=['Arid','Semi-arid','Dry sub-humid','Humid'],right=False)
    order=['Forest','Shrubland','Savanna/Grassland','Cropland']
    cols=['Arid','Semi-arid','Dry sub-humid','Humid']
    m=d.groupby(['PFT','zone'],observed=True)[metric].mean().unstack().reindex(index=order,columns=cols)
    fig,ax=plt.subplots(figsize=(7,4.1),constrained_layout=True)
    im=ax.imshow(m.to_numpy(),cmap='YlGnBu',aspect='auto')
    ax.set_xticks(range(4),cols)
    ax.set_yticks(range(4),['Forest','Shrubland','Grass/savanna','Cropland'])
    for (i,j),value in np.ndenumerate(m.to_numpy()):
        if np.isfinite(value):ax.text(j,i,f'{value:.0f}',ha='center',va='center',color='black',fontsize=9)
    fig.colorbar(im,ax=ax,label='Model-derived LAI–ET sensitivity (mm yr$^{-1}$ LAI$^{-1}$)')
    file=Path(file);file.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(file,bbox_inches='tight');plt.close(fig)
