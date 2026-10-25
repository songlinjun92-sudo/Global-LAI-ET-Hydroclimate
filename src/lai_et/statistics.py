"""Block-bootstrap hydroclimate gradients, PFT contrasts and symmetric decomposition."""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from .io import blocks,classes,ZONES

MAIN_PFT=['Shrubland','Savanna/Grassland','Cropland']

def block_medians(d,metric,ai,deg):
    x=d[['lon','lat',ai,metric]].dropna()
    x=x[x[ai]>=0]
    if len(x)<5:return pd.DataFrame(columns=['block','ai','value'])
    x=x.assign(block=blocks(x.lat,x.lon,deg))
    g=x.groupby('block',as_index=False).agg(ai=(ai,'median'),value=(metric,'median'))
    return g

def rho_ci(d,metric,ai,degrees=10,reps=2000,seed=1234):
    g=block_medians(d,metric,ai,degrees)
    if len(g)<5 or g.ai.nunique()<3 or g.value.nunique()<3:
        return {'rho':np.nan,'lo95':np.nan,'hi95':np.nan,'n_blocks':len(g),'n_pixels':int(d[[metric,ai]].dropna().shape[0])}
    x=g.ai.to_numpy();y=g.value.to_numpy()
    point=float(spearmanr(x,y).statistic)
    rng=np.random.default_rng(seed)
    vals=[]
    for _ in range(reps):
        ids=rng.integers(0,len(x),len(x))
        if len(np.unique(x[ids]))>2 and len(np.unique(y[ids]))>2:
            val=float(spearmanr(x[ids],y[ids]).statistic)
            if np.isfinite(val):vals.append(val)
    lo,hi=np.quantile(vals,[.025,.975]) if len(vals)>=reps*0.75 else (np.nan,np.nan)
    return {'rho':point,'lo95':float(lo),'hi95':float(hi),
            'n_blocks':len(g),'n_pixels':int(d[[metric,ai]].dropna().shape[0])}

def paired_rho_delta(d, col1, col2, ai,degrees=10,reps=2000,seed=1234):
    z=d.dropna(subset=[ai,col1,col2]).copy()
    z['block']=blocks(z.lat,z.lon,degrees)
    g=z.groupby('block').agg(ai=(ai,'median'),m1=(col1,'median'),m2=(col2,'median')).dropna()
    if len(g)<5:return {'difference':np.nan,'lo95':np.nan,'hi95':np.nan,'n_blocks':len(g)}
    ar=g[['ai','m1','m2']].to_numpy()
    def delta(a):
        if any(len(np.unique(a[:,k]))<3 for k in range(3)):return np.nan
        return float(spearmanr(a[:,0],a[:,1]).statistic-spearmanr(a[:,0],a[:,2]).statistic)
    point=delta(ar);rng=np.random.default_rng(seed)
    vals=[delta(ar[rng.integers(0,len(ar),len(ar))]) for _ in range(reps)]
    vals=np.asarray(vals,dtype=float);vals=vals[np.isfinite(vals)]
    lo,hi=np.quantile(vals,[.025,.975]) if len(vals)>=0.75*reps else (np.nan,np.nan)
    return {'difference':point,'lo95':float(lo),'hi95':float(hi),'n_blocks':len(g)}

def kitagawa(df,col,ai,pfts=MAIN_PFT,min_pixels=10,min_blocks=3):
    """Symmetric Kitagawa decomposition: humid minus semi-arid sensitivity.

    Returns NaN if any of the three PFTs is unsupported in either climate class.
    """
    d=df[df.PFT.isin(pfts)].dropna(subset=[ai,col,'lat','lon']).copy()
    d['zone']=classes(d[ai])
    d=d[d.zone.isin(['Semi-arid','Humid'])].copy()
    if d.empty:return {'total':np.nan,'within':np.nan,'composition':np.nan,'supported':False}
    d['weight']=np.cos(np.deg2rad(d.lat)).clip(0)
    d['block']=blocks(d.lat,d.lon,10)
    v={}
    for zone in ['Semi-arid','Humid']:
        z=d[d.zone==zone]
        total_weight=z.weight.sum()
        for pft in pfts:
            q=z[z.PFT==pft]
            if len(q)<min_pixels or q.block.nunique()<min_blocks or q.weight.sum()==0:
                return {'total':np.nan,'within':np.nan,'composition':np.nan,
                        'supported':False,'reason':f'{pft}/{zone} has inadequate spatial support'}
            v[(zone,pft)]=(float(q.weight.sum()/total_weight),float(np.average(q[col],weights=q.weight)))
    within=composition=0.0
    for pft in pfts:
        w0,m0=v[('Semi-arid',pft)];w1,m1=v[('Humid',pft)]
        within+=(w0+w1)/2*(m1-m0)
        composition+=(m0+m1)/2*(w1-w0)
    return {'total':within+composition,'within':within,'composition':composition,'supported':True}

def kitagawa_bootstrap(d,col,ai,cfg):
    """Efficient, geographically paired block bootstrap of Kitagawa terms.

    Resamples whole blocks by multiplicity using block-level sufficient statistics;
    identical block draws can be reproduced for different products using the same seed.
    """
    s=cfg['statistics'];degree=int(s['block_degrees']);reps=int(s['bootstrap_replicates'])
    sub=d[d.PFT.isin(MAIN_PFT)].dropna(subset=[col,ai,'lon','lat']).copy()
    sub['zone']=classes(sub[ai])
    sub=sub[sub.zone.isin(['Semi-arid','Humid'])].copy()
    if sub.empty:return []
    sub['block']=blocks(sub.lat,sub.lon,degree)
    sub['w']=np.cos(np.deg2rad(sub.lat)).clip(0)
    sub['wy']=sub['w']*sub[col]
    point=kitagawa(sub,col,ai,min_pixels=int(s['min_pft_pixels']),
                   min_blocks=int(s['min_pft_blocks']))
    unique=np.sort(sub.block.unique()); nb=len(unique)
    aggregates=sub.groupby(['block','zone','PFT'],observed=True)[['w','wy']].sum()
    ws=np.zeros((nb,2,3));ys=np.zeros((nb,2,3))
    index={b:i for i,b in enumerate(unique)}
    for (block,zone,pft),row in aggregates.iterrows():
        ii=index[block];jj=['Semi-arid','Humid'].index(str(zone));kk=MAIN_PFT.index(pft)
        ws[ii,jj,kk]=row.w;ys[ii,jj,kk]=row.wy
    rng=np.random.default_rng(int(s['seed']))
    vals=np.empty((reps,3),dtype=float);vals.fill(np.nan)
    for t in range(reps):
        draw=rng.multinomial(nb,np.ones(nb)/nb)
        w=np.einsum('i,ijk->jk',draw,ws)
        wy=np.einsum('i,ijk->jk',draw,ys)
        if np.any(w<=0):continue
        norms=w/w.sum(axis=1,keepdims=True)
        means=wy/w
        within=np.sum((norms[0]+norms[1])*.5*(means[1]-means[0]))
        composition=np.sum((means[0]+means[1])*.5*(norms[1]-norms[0]))
        vals[t]=[within+composition,within,composition]
    rows=[]
    for j,term in enumerate(['total','within','composition']):
        v=vals[:,j];v=v[np.isfinite(v)]
        lo,hi=np.quantile(v,[.025,.975]) if point.get('supported') and len(v)>=.75*reps else (np.nan,np.nan)
        rows.append({'term':term,'estimate':point.get(term,np.nan),'lo95':float(lo),'hi95':float(hi),
                     'n_blocks':nb,'n_effective_resamples':len(v),
                     'supported':bool(point.get('supported',False))})
    return rows

def summarize(ice_tables,cfg):
    """Compute continuous gradients on the intersection of local-ICE supports."""
    if not ice_tables:raise ValueError('Missing ICE tables.')
    base=ice_tables[0]
    for q in ice_tables[1:]:
        add=[x for x in q if x.endswith('_local_ICE') and x not in base]
        base=base.merge(q[['lon','lat']+add],on=['lon','lat'],validate='one_to_one',how='inner')
    s=cfg['statistics'];degree=int(s['block_degrees']);reps=int(s['bootstrap_replicates']);seed=int(s['seed'])
    aicols=[c for c in ['AI_A','AI_B','AI_C'] if c in base]
    metrics=[c for c in base if c.endswith('_local_ICE')]
    gradient=[];paired=[];pft=[];decomp=[]
    for ai in aicols:
        for j,metric in enumerate(metrics):
            d=base.dropna(subset=[ai,metric])
            z=rho_ci(d,metric,ai,degree,reps,seed+j)
            gradient.append({'ai':ai,'metric':metric,**z})
            if 'PFT' in d:
                for pf in ['Forest']+MAIN_PFT:
                    v=d[d.PFT==pf]
                    zz=rho_ci(v,metric,ai,degree,reps,seed+j)
                    pft.append({'ai':ai,'metric':metric,'PFT':pf,**zz})
        for endpoint in cfg['targets']:
            a=f'PML_{endpoint}_local_ICE';b=f'GLEAM_{endpoint}_local_ICE'
            if a in base and b in base:
                result=paired_rho_delta(base,a,b,ai,degree,reps,seed)
                paired.append({'ai':ai,'endpoint':endpoint,**result})
        if 'PFT' in base:
            for endpoint in ['ET']:
                for source in ['GLEAM','PML']:
                    col=f'{source}_{endpoint}_local_ICE'
                    if col in base:
                        for r in kitagawa_bootstrap(base,col,ai,cfg):
                            decomp.append({'ai':ai,'source':source,'endpoint':endpoint,**r})
    return base,pd.DataFrame(gradient),pd.DataFrame(paired),pd.DataFrame(pft),pd.DataFrame(decomp)
