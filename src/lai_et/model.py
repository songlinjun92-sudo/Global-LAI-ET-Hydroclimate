"""Matched RF modeling, geography-block CV, and independent-product comparisons."""
from __future__ import annotations
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
from sklearn.model_selection import GroupKFold
import joblib
from .io import save,write_json,blocks


def _pair_panel(panel, target_cols, cfg):
    p=panel.dropna(subset=cfg['predictors']+list(target_cols)).copy()
    if cfg['period']['require_balanced']:
        all_years=int(cfg['period']['end'])-int(cfg['period']['start'])+1
        count=p.groupby(['lon','lat'])['year'].nunique()
        eligible=count[count==all_years].index
        if len(eligible)==0:raise ValueError('No pixels have balanced observations. Consider period.require_balanced=false for a documented sensitivity experiment.')
        p=p.set_index(['lon','lat']).loc[eligible].reset_index()
    if len(p)<10:raise ValueError('Too few valid observations.')
    return p.reset_index(drop=True)

def _sample(ix,cap,seed):
    ix=np.asarray(ix,dtype=int)
    if cap and len(ix)>int(cap):
        return np.sort(np.random.default_rng(seed).choice(ix,size=int(cap),replace=False))
    return ix

def _forest(cfg,trees,seed):
    m=cfg['model']
    return RandomForestRegressor(n_estimators=int(trees),
        min_samples_leaf=int(m['min_samples_leaf']),n_jobs=int(m['n_jobs']),
        random_state=int(seed),max_features=1.0)

def _metrics(y,pred):
    return {'r2':float(r2_score(y,pred)), 'rmse':float(np.sqrt(mean_squared_error(y,pred))),
            'mae':float(mean_absolute_error(y,pred))}

def _cv(p,columns,target,cfg,seed):
    """Spatial CV on held-out 10-degree cells; forward temporal CV on held-out future years."""
    m=cfg['model'];x=p[cfg['predictors']].to_numpy(np.float32)
    y=p[target].to_numpy(np.float32)
    yy=p['year'].to_numpy(int)
    groups=blocks(p.lat,p.lon,int(cfg['statistics']['block_degrees']))
    uniq=np.unique(groups)
    if len(uniq)<2:raise ValueError('Spatial CV needs at least two distinct blocks.')
    folds=min(int(m['cv_folds']),len(uniq))
    rows=[]
    if m['skip_cv']:return pd.DataFrame(rows)
    splits=[(f'space_{k}',tr,te,'spatial') for k,(tr,te) in enumerate(GroupKFold(n_splits=folds).split(x,y,groups),1)]
    for label,upto,lo in [('time_2001_2010',2000,(2001,2010)),('time_2011_2020',2010,(2011,2020))]:
        tr=np.flatnonzero(yy<=upto)
        te=np.flatnonzero((yy>=lo[0])&(yy<=lo[1]))
        if len(tr)>10 and len(te)>10:splits.append((label,tr,te,'temporal'))
    for k,(label,tr,te,kind) in enumerate(splits,1):
        tr=_sample(tr,int(m['train_cap']),seed+k)
        for label_model,featurecols in [('climate_only',cfg['predictors'][1:]),('climate_lai',cfg['predictors'])]:
            idx=[cfg['predictors'].index(j) for j in featurecols]
            mod=_forest(cfg,m['cv_trees'],seed+k)
            mod.fit(x[tr][:,idx],y[tr])
            pred=mod.predict(x[te][:,idx]); stat=_metrics(y[te],pred)
            stat.update(source=columns[0],target=columns[1],split=label,split_type=kind,
                        model=label_model,n_train=len(tr),n_test=len(te))
            if kind=='temporal':
                df=pd.DataFrame({'lon':p.lon.iloc[te].values,'lat':p.lat.iloc[te].values,
                                 'obs':y[te],'pred':pred})
                # Calendar-held-out anomalies are not forecast anomalies; descriptive only.
                mu=df.groupby(['lon','lat'])[['obs','pred']].transform('mean')
                obs=(df.obs-mu.obs).to_numpy();pr=(df.pred-mu.pred).to_numpy()
                stat['within_test_period_anomaly_r2']=float(1-np.sum((obs-pr)**2)/np.sum(obs**2)) if np.sum(obs**2)>0 else np.nan
            rows.append(stat)
    return pd.DataFrame(rows)

def train_models(panel,cfg,out):
    from .sensitivity import local_ice
    out=Path(out);(out/'models').mkdir(parents=True,exist_ok=True)
    m=cfg['model'];predictors=cfg['predictors'];pairs=cfg['pml_pairs'];targets=cfg['targets']
    all_ice=[];all_cv=[];all_models={}
    for key,gleam in targets.items():
        target_set=[('GLEAM',gleam)]
        if pairs.get(key) in panel:target_set.append(('PML',pairs[key]))
        # To compare GLEAM with PML, all fit rows for the endpoint are matched.
        p=_pair_panel(panel,[x[1] for x in target_set],cfg)
        p['block10']=blocks(p.lat,p.lon,10)
        ix=_sample(np.arange(len(p)),int(m['train_cap']),int(m['seed'])+109)
        per_target=[]
        for i,(source,col) in enumerate(target_set):
            seed=int(m['seed'])+i
            cv=_cv(p,(source,key),col,cfg,seed)
            if len(cv):all_cv.append(cv)
            model=_forest(cfg,m['trees'],int(m['seed'])+501)
            model.fit(p.iloc[ix][predictors].to_numpy(np.float32),p.iloc[ix][col].to_numpy(np.float32))
            model_file=out/'models'/f'rf_{source.lower()}_{key.lower()}.joblib'
            joblib.dump({'model':model,'predictors':predictors,'source':source,'target':key,
                         'train_rows':len(ix),'matched_rows':len(p),'seed':int(m['seed']),
                         'n_estimators':int(m['trees'])},model_file,compress=3)
            all_models[(source,key)]=model
            print(f'RF {key}/{source}: {len(ix)} training records, {len(p)} matched records',flush=True)
            res=local_ice(model,p,cfg)
            res=res.rename(columns={'sensitivity':f'{source}_{key}_local_ICE'})
            per_target.append(res)
            if m['do_shap']:
                explain_shap(model,p,source,key,cfg,out)
        r=per_target[0]
        for s in per_target[1:]:
            r=r.merge(s[['lon','lat',f'{target_set[1][0]}_{key}_local_ICE']],on=['lon','lat'],validate='one_to_one')
        all_ice.append(r)
        save(r,out/f'ice_{key.lower()}.csv.gz')
        # Warn about target-specific support differences; final combined intersection is explicit.
    if all_cv:save(pd.concat(all_cv,ignore_index=True),out/'rf_cross_validation.csv')
    return all_ice,all_models

def explain_shap(model, panel,source,key,cfg,out):
    import shap
    n=min(len(panel),int(cfg['model']['shap_sample']))
    rng=np.random.default_rng(int(cfg['model']['seed']))
    loc=rng.choice(len(panel),n,replace=False)
    sample=panel.iloc[loc][cfg['predictors']]
    explainer=shap.TreeExplainer(model)
    values=explainer.shap_values(sample)
    df=pd.DataFrame({'feature':cfg['predictors'],'mean_abs_shap':np.mean(np.abs(values),axis=0)})
    save(df,Path(out)/f'shap_{source.lower()}_{key.lower()}.csv')

def flux_agreement(panel,cfg):
    out=[]
    for key,gleam in cfg['targets'].items():
        pml=cfg['pml_pairs'].get(key)
        if pml not in panel:continue
        d=panel.dropna(subset=[gleam,pml])
        if len(d)<3:continue
        a=d[pml].to_numpy(float);b=d[gleam].to_numpy(float)
        da=d[pml]-d.groupby(['lon','lat'])[pml].transform('mean')
        db=d[gleam]-d.groupby(['lon','lat'])[gleam].transform('mean')
        out.append({'endpoint':key,'matched_grid_years':len(d),
           'annual_correlation':float(np.corrcoef(a,b)[0,1]),
           'annual_anomaly_correlation':float(da.corr(db)),
           'bias_pml_minus_gleam':float((a-b).mean()),
           'rmse':float(np.sqrt(np.mean((a-b)**2)))})
    return pd.DataFrame(out)
