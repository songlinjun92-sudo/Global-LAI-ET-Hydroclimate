"""Small deterministic synthetic example. These values are NOT study results."""
from pathlib import Path
import copy
import numpy as np
import pandas as pd
from .config import DEFAULTS,output_dir
from .io import save

def generate(root,seed=42):
    root=Path(root);root.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(seed)
    lon_arr=np.arange(-137.5,112.5,10.0)
    lat_arr=np.arange(-47.5,52.5,10.0)
    cells=[]
    labels=['Shrubland','Savanna/Grassland','Cropland','Forest']
    for i,lat in enumerate(lat_arr):
        for j,lon in enumerate(lon_arr):
            t=(i*len(lon_arr)+j)
            ai=[.12,.35,.56,1.0,1.65][(3*j+i)%5]
            pft=labels[(j+i)%len(labels)]
            cells.append((lon,lat,pft,ai,ai*1.2,ai*.83))
    at=pd.DataFrame(cells,columns=['lon','lat','PFT','AI_A','AI_B','AI_C'])
    annual=[];pml=[]
    for year in range(1982,1994):
        n=len(at)
        ai=at.AI_A.to_numpy()
        gs_lai=1.0+0.22*(year-1982)/11.0 + .2*np.log1p(ai)+rng.normal(0,.19,n)
        ta=20-np.abs(at.lat.to_numpy())*.12+rng.normal(0,.3,n)
        prec=300+350*ai+rng.normal(0,25,n)
        rs=180+15*np.cos(np.deg2rad(at.lat.to_numpy()))+rng.normal(0,8,n)
        vpd=1.8/(.25+ai)+rng.normal(0,.05,n)
        et=180+(.5*prec)+gs_lai*(115/(.3+ai))+rng.normal(0,18,n)
        trans=.55*et+gs_lai*12/(.4+ai)+rng.normal(0,10,n)
        ei=7+6*ai+gs_lai*(4+8*ai)+rng.normal(0,4,n)
        eb=et-trans-ei
        annual.append(at[['lon','lat']].assign(year=year,LAI_gs=gs_lai,
            Ta_ann=ta,P_ann=prec,Rs_ann=rs,VPD_ann=vpd,ET_ann=et,T_ann=trans,Ei_ann=ei,Eb_ann=eb))
        pml.append(at[['lon','lat']].assign(year=year,PML_ET=et+rng.normal(0,20,n),
            PML_Ec=trans*.84+rng.normal(0,10,n),PML_Ei=ei*.78+rng.normal(0,3,n),
            PML_Es=eb+15+rng.normal(0,10,n)))
    data=root/'data/processed';data.mkdir(parents=True,exist_ok=True)
    save(at,data/'pixel_attributes.csv')
    save(pd.concat(annual,ignore_index=True),data/'annual_panel.csv.gz')
    save(pd.concat(pml,ignore_index=True),data/'pml_panel.csv.gz')
    return len(at)
