"""Local historical LAI ICE; yearly finite-difference diagnostics for fixed models."""
from __future__ import annotations
import numpy as np
import pandas as pd
from .io import blocks


def local_ice(model,panel,cfg):
    """OLS slope of 20 RF predictions over each pixel's historical q05-q95 LAI range.

    At each pixel hold the *multiannual-mean* climate covariates fixed.
    Pixels with insufficient LAI range produce NaN instead of extrapolated slopes.
    """
    predictors=cfg['predictors'];rules=cfg['ice']
    cols=['lon','lat',*predictors]
    grp=panel.groupby(['lon','lat'],sort=False)
    mean=grp[predictors].mean()
    qq=grp['LAI_gs'].quantile([float(rules['q_low']),float(rules['q_high'])]).unstack()
    qq.columns=['q_low','q_high']
    p=mean.join(qq).reset_index()
    extras=[x for x in ['PFT','AI_A','AI_B','AI_C'] if x in panel]
    if extras:
        at=grp[extras].first().reset_index()
        p=p.merge(at,on=['lon','lat'],validate='one_to_one')
    p['area_weight_coslat']=np.cos(np.deg2rad(p.lat)).clip(0)
    p['block10']=blocks(p.lat,p.lon,10)
    p['sensitivity']=np.nan
    good=(p[predictors+['q_low','q_high']].notna().all(axis=1)
          &(p.q_high-p.q_low>=float(rules['minimum_range'])))
    ids=np.flatnonzero(good.to_numpy())
    frac=np.linspace(0,1,int(rules['points']),dtype=np.float32)
    if len(frac)<3:raise ValueError('ICE requires at least three points.')
    bs=int(rules['batch_pixels'])
    if bs<1:raise ValueError('batch_pixels must be positive')
    lai_index=predictors.index('LAI_gs')
    for lo in range(0,len(ids),bs):
        ii=ids[lo:lo+bs]
        x=p.iloc[ii][predictors].to_numpy(np.float32)
        left=p.iloc[ii].q_low.to_numpy(np.float32);right=p.iloc[ii].q_high.to_numpy(np.float32)
        grid=left[:,None]+(right-left)[:,None]*frac[None,:]
        xx=np.repeat(x,len(frac),axis=0)
        xx[:,lai_index]=grid.reshape(-1)
        yy=model.predict(xx).reshape(len(ii),len(frac))
        xx0=grid-grid.mean(axis=1,keepdims=True)
        yy0=yy-yy.mean(axis=1,keepdims=True)
        p.loc[ii,'sensitivity']=(xx0*yy0).sum(axis=1)/(xx0**2).sum(axis=1)
    return p

def annual_fixed_model_sensitivity(model,panel,cfg):
    """Finite difference at each pixel-year; separate estimand from local ICE slopes."""
    predictors=cfg['predictors'];step=float(cfg['temporal']['step'])
    if step<=0:raise ValueError('temporal.step must be positive.')
    data=panel.dropna(subset=predictors).copy()
    x=data[predictors].to_numpy(np.float32)
    col=predictors.index('LAI_gs');n=len(x)
    derivatives=np.empty(n,dtype=np.float64)
    batch=8192
    for start in range(0,n,batch):
        row=x[start:start+batch].copy()
        high=row.copy();low=row.copy()
        high[:,col]+=step;low[:,col]=np.maximum(0,low[:,col]-step)
        delta=(high[:,col]-low[:,col]).astype(float)
        derivatives[start:start+batch]=(model.predict(high)-model.predict(low))/delta
    data['finite_sensitivity']=derivatives
    data['weight']=np.cos(np.deg2rad(data.lat)).clip(0)
    out=[]
    for year,g in data.groupby('year'):
        out.append({'year':int(year), 'mean_finite_sensitivity':float(np.average(g.finite_sensitivity,weights=g.weight)),
                    'n_grid_years':len(g)})
    return pd.DataFrame(out)

def yearly_trend(series):
    """HAC time-series trend, NOT a geographical block-bootstrap statement."""
    import statsmodels.api as sm
    x=series.year.to_numpy(float);y=series.mean_finite_sensitivity.to_numpy(float)
    if len(series)<6:return {'slope':np.nan,'hac_p':np.nan}
    result=sm.OLS(y,sm.add_constant(x-x.mean())).fit(cov_type='HAC',cov_kwds={'maxlags':3})
    return {'slope':float(result.params[1]),'hac_p':float(result.pvalues[1])}
