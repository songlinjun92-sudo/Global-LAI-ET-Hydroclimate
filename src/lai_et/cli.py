from __future__ import annotations
import argparse
from pathlib import Path
import sys
import numpy as np
from .config import load_config,path,output_dir
from .io import load_inputs,save,write_json
from .panel import raster_panel,pml_panel
from .model import train_models,flux_agreement
from .sensitivity import annual_fixed_model_sensitivity,yearly_trend
from .statistics import summarize
from .figures import gradient_figure,class_heatmap


def run(cfg):
    out=output_dir(cfg);out.mkdir(parents=True,exist_ok=True)
    data=load_inputs(cfg)
    print(f'Annual rows: {len(data):,}; pixels: {len(data.drop_duplicates(["lon","lat"])):,}',flush=True)
    ice_tables,models=train_models(data,cfg,out)
    common,grad,paired,pft,decomp=summarize(ice_tables,cfg)
    save(common,out/'pixel_ice_common.csv.gz')
    save(grad,out/'continuous_ai_gradients.csv')
    save(paired,out/'paired_product_differences.csv')
    save(pft,out/'pft_ai_gradients.csv')
    save(decomp,out/'kitagawa_decomposition.csv')
    save(flux_agreement(data,cfg),out/'pml_gleam_flux_agreement.csv')
    gradient_figure(grad,out/'figures/independent_et_gradients.png')
    if 'PFT' in common:class_heatmap(common,file=out/'figures/pft_ai_heatmap.png')
    if cfg['temporal']['enabled']:
        rows=[]
        for (source,key),model in models.items():
            target=cfg['targets'][key] if source=='GLEAM' else cfg['pml_pairs'][key]
            dat=data.dropna(subset=cfg['predictors']+[target])
            ts=annual_fixed_model_sensitivity(model,dat,cfg)
            ts.insert(0,'source',source);ts.insert(1,'endpoint',key)
            save(ts,out/f'temporal_{source.lower()}_{key.lower()}.csv')
            rows.append({'source':source,'endpoint':key,**yearly_trend(ts)})
        save(__import__('pandas').DataFrame(rows),out/'temporal_trend_hac.csv')
    meta={'period':cfg['period'],'predictors':cfg['predictors'],'targets':cfg['targets'],
          'model':cfg['model'],'ice':cfg['ice'],'statistics':cfg['statistics'],
          'data_were_not_included':True,'sensitivity_is_model_association_not_causality':True,
          'independent_product_validation_shares_some_LAI_input':True,
          'common_pixel_count':len(common),'ai_definitions':{
             'AI_A':'ERA5 precipitation / GLEAM potential evaporation',
             'AI_B':'ERA5 precipitation / TerraClimate PET',
             'AI_C':'TerraClimate precipitation / TerraClimate PET'},
          'kitagawa_contrast':'Semi-arid (0.20–0.50) vs Humid (>=0.65), supported nonforest PFTs'}
    write_json(meta,out/'method_manifest.json')
    print('ANALYSIS COMPLETE. Results:',out,flush=True)
    return out

def main(argv=None):
    parser=argparse.ArgumentParser(prog='lai-et',description='Global LAI-ET hydroclimate reproducibility toolkit')
    sub=parser.add_subparsers(dest='command',required=True)
    for name in ['run','make-panel','prepare-pml']:
        sp=sub.add_parser(name);sp.add_argument('--config',default='configs/config.yaml')
    demo=sub.add_parser('demo');demo.add_argument('--root',default='examples/demo_run')
    args=parser.parse_args(argv)
    if args.command=='demo':
        from .demo import generate
        from .config import DEFAULTS
        import copy
        root=Path(args.root).expanduser().resolve();generate(root)
        cfg=copy.deepcopy(DEFAULTS)
        cfg['_base_dir']=root;cfg['_config_path']=root/'demo.yaml'
        cfg['period']={'start':1982,'end':1993,'require_balanced':True}
        cfg['inputs']['pml_panel']='data/processed/pml_panel.csv.gz'
        cfg['outputs']['directory']='outputs/demo'
        cfg['model'].update({'trees':12,'cv_trees':8,'train_cap':4000,'cv_folds':3,'n_jobs':2,
                             'min_samples_leaf':2,'skip_cv':False})
        cfg['statistics'].update({'bootstrap_replicates':50,'min_pft_pixels':3,'min_pft_blocks':1})
        cfg['temporal']['enabled']=False
        out=run(cfg)
        print('DEMO / SMOKE TEST ONLY; NEVER USE SYNTHETIC OUTPUTS IN A PAPER.')
        return
    cfg=load_config(args.config)
    if args.command=='make-panel':raster_panel(cfg);print('Annual panel assembled successfully.')
    elif args.command=='prepare-pml':pml_panel(cfg);print('Annual PML NetCDF files harmonized successfully.')
    elif args.command=='run':run(cfg)

if __name__=='__main__':main()
