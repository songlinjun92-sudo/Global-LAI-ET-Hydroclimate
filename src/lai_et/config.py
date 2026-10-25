from __future__ import annotations
from pathlib import Path
import copy
import yaml

DEFAULTS = {
    'period': {'start': 1982, 'end': 2020, 'require_balanced': True},
    'inputs': {'annual_panel': 'data/processed/annual_panel.csv.gz',
               'pixel_attributes': 'data/processed/pixel_attributes.csv',
               'pml_panel': None},
    'outputs': {'directory': 'outputs/publication_run'},
    'predictors': ['LAI_gs', 'Ta_ann', 'P_ann', 'Rs_ann', 'VPD_ann'],
    'targets': {'ET':'ET_ann','T':'T_ann','Ei':'Ei_ann','Eb':'Eb_ann'},
    'pml_pairs': {'ET':'PML_ET','T':'PML_Ec','Ei':'PML_Ei','Eb':'PML_Es'},
    'model': {'trees':250, 'cv_trees':100, 'n_jobs':2, 'train_cap':150000,
              'cv_folds':5, 'seed':202610, 'min_samples_leaf':2,
              'skip_cv':False, 'do_shap':False, 'shap_sample':200},
    'ice': {'q_low':0.05,'q_high':0.95,'points':20,'minimum_range':0.10,'batch_pixels':128},
    'statistics': {'block_degrees':10,'bootstrap_replicates':2000,'min_pft_pixels':30,
                   'min_pft_blocks':8,'seed':1234},
    'temporal': {'enabled':False,'step':0.05},
    'raster_inputs': {'annual': {}, 'pixel': {}},
    'pml_netcdf': {'root':None,'start_year':1982,'end_year':2020},
}

def recursive_merge(a,b):
    x=copy.deepcopy(a)
    for k,v in b.items():
        if isinstance(v,dict) and isinstance(x.get(k),dict):
            x[k]=recursive_merge(x[k],v)
        else:x[k]=v
    return x

def load_config(file):
    f=Path(file).expanduser().resolve()
    raw=yaml.safe_load(f.read_text(encoding='utf-8')) or {}
    if not isinstance(raw,dict):raise ValueError('Config YAML must be an object.')
    cfg=recursive_merge(DEFAULTS,raw)
    cfg['_config_path']=f
    cfg['_base_dir']=f.parent.parent if f.parent.name=='configs' else f.parent
    return cfg

def path(cfg,value):
    if value is None:return None
    p=Path(value).expanduser()
    return p if p.is_absolute() else (cfg['_base_dir']/p).resolve()

def output_dir(cfg):
    return path(cfg,cfg['outputs']['directory'])
