from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .config import path

FEATS=('LAI_gs','Ta_ann','P_ann','Rs_ann','VPD_ann')
PFT_MAP={1:'Forest',2:'Shrubland',3:'Savanna/Grassland',4:'Cropland',
         5:'Wetland',6:'Barren/Sparse'}
ZONES=['Arid','Semi-arid','Dry sub-humid','Humid']

def csv_read(pathobj):
    return pd.read_csv(pathobj,low_memory=False)

def save(df,file):
    f=Path(file);f.parent.mkdir(parents=True,exist_ok=True)
    df.to_csv(f,index=False,compression='gzip' if f.suffix=='.gz' else None)
    return f

def write_json(obj,file):
    f=Path(file);f.parent.mkdir(parents=True,exist_ok=True)
    def default(x):
        if isinstance(x,(np.integer,)):return int(x)
        if isinstance(x,(np.floating,)):return float(x)
        if isinstance(x,(np.bool_,)):return bool(x)
        if isinstance(x,Path):return str(x)
        raise TypeError(type(x).__name__)
    f.write_text(json.dumps(obj,indent=2,ensure_ascii=False,default=default),encoding='utf-8')

def blocks(lat,lon,degrees=10):
    lat=np.asarray(lat,float);lon=np.asarray(lon,float)
    return (np.floor((lat+90)/degrees).astype('int64')*(int(np.ceil(360/degrees))+1)
           +np.floor((lon+180)/degrees).astype('int64'))

def classes(ai):
    arr=np.asarray(ai,float)
    return pd.Categorical(pd.cut(arr,[-np.inf,0.2,0.5,0.65,np.inf],
                            labels=ZONES,right=False), categories=ZONES,ordered=True)

def validate_panel(df,cfg):
    """Validate a unique cell-year table without silently inventing missing variables."""
    required={'year','lon','lat',*cfg['predictors'],*cfg['targets'].values()}
    missing=required-set(df.columns)
    if missing:raise ValueError(f'Annual panel missing required columns: {sorted(missing)}')
    d=df.copy()
    for c in required:d[c]=pd.to_numeric(d[c],errors='coerce')
    if d[list(required)].isna().all(axis=1).any():
        raise ValueError('Panel contains entirely empty records.')
    if d.duplicated(['year','lon','lat']).any():raise ValueError('Duplicate year/lon/lat records.')
    start=int(cfg['period']['start']);end=int(cfg['period']['end'])
    d=d[d.year.between(start,end)].copy()
    if d.empty:raise ValueError('No records within configured period.')
    if not np.allclose(d.year,np.round(d.year)):raise ValueError('Years must be integers.')
    if not d.lon.between(-180,180).all() or not d.lat.between(-90,90).all():
        raise ValueError('Invalid longitude/latitude.')
    d['year']=d.year.astype(int)
    d['block10']=blocks(d.lat,d.lon,10)
    return d

def load_inputs(cfg):
    panel_file=path(cfg,cfg['inputs']['annual_panel'])
    if not panel_file.exists():raise FileNotFoundError(f'Annual panel not found: {panel_file}')
    panel=validate_panel(csv_read(panel_file),cfg)
    attrs_file=path(cfg,cfg['inputs']['pixel_attributes'])
    if attrs_file and attrs_file.exists():
        attrs=csv_read(attrs_file)
        if not {'lon','lat','PFT'}.issubset(attrs.columns):
            raise ValueError('Pixel attributes must contain lon, lat, PFT.')
        if attrs.duplicated(['lon','lat']).any():raise ValueError('Duplicate pixel attributes.')
        attrs=add_ai(attrs)
        overlap=[x for x in attrs if x in panel and x not in ('lon','lat')]
        attrs=attrs.drop(columns=overlap)
        panel=panel.merge(attrs,on=['lon','lat'],how='left',validate='many_to_one')
    elif 'PFT' not in panel:
        raise FileNotFoundError('A pixel_attributes CSV with PFT and AI is required.')
    if 'AI_A' not in panel:
        panel=add_ai(panel)
    pf=path(cfg,cfg['inputs']['pml_panel'])
    if pf is not None:
        pml=csv_read(pf)
        if pml.duplicated(['year','lon','lat']).any():raise ValueError('Duplicate PML grid-year records.')
        for c in ['year','lon','lat']:
            if c not in pml:raise ValueError(f'PML panel missing {c}')
        keep=['year','lon','lat']+[x for x in cfg['pml_pairs'].values() if x in pml]
        pml=pml[keep]
        panel=panel.merge(pml,on=['year','lon','lat'],how='left',validate='one_to_one',suffixes=('','_external'))
    if 'PFT' in panel:
        panel['PFT']=panel['PFT'].map(PFT_MAP).fillna(panel['PFT'].astype(str))
    return panel

def add_ai(d):
    d=d.copy()
    # Three physically distinct AI definitions; do not force same numerical climate-bin composition.
    definitions={'AI_A':('P_ERA','Ep_GLEAM'),'AI_B':('P_ERA','PET_Terra'),
                 'AI_C':('P_Terra','PET_Terra')}
    for ai,(numerator,denom) in definitions.items():
        if ai not in d and numerator in d and denom in d:
            p=pd.to_numeric(d[numerator],errors='coerce')
            pet=pd.to_numeric(d[denom],errors='coerce')
            d[ai]=np.where(pet>0,p/pet,np.nan)
    if 'AI_A' not in d and 'AI' in d:d['AI_A']=d['AI']
    for ai in [x for x in definitions if x in d]:
        d[ai]=pd.to_numeric(d[ai],errors='coerce')
        d.loc[(d[ai]<0)|(~np.isfinite(d[ai])),ai]=np.nan
    return d
