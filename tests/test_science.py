import numpy as np
import pandas as pd
import pytest
from lai_et.io import classes,blocks,add_ai,validate_panel
from lai_et.config import DEFAULTS
from lai_et.statistics import kitagawa,rho_ci,paired_rho_delta
from lai_et.sensitivity import local_ice
from sklearn.linear_model import LinearRegression


def test_climate_boundaries():
    x=classes([0.0,0.1999,0.2,0.4999,0.5,0.649,0.65,3.55])
    assert list(map(str,x)) == ['Arid','Arid','Semi-arid','Semi-arid','Dry sub-humid','Dry sub-humid','Humid','Humid']


def test_aridity_formulas():
    d=pd.DataFrame({'P_ERA':[500.0],'Ep_GLEAM':[1000.0],
                    'PET_Terra':[250.0],'P_Terra':[200.0]})
    a=add_ai(d)
    assert a.AI_A.iloc[0] == .5
    assert a.AI_B.iloc[0] == 2.0
    assert a.AI_C.iloc[0] == .8


def test_block_key_is_spatial():
    assert blocks([12,18],[-90,-89],10).tolist() == [blocks([12],[-90],10)[0]]*2
    assert blocks([12,23],[-90,-90],10)[0] != blocks([12,23],[-90,-90],10)[1]


def test_balanced_panel_duplicate_rejected():
    data={key:[1.,2.] for key in ['LAI_gs','Ta_ann','P_ann','Rs_ann','VPD_ann','ET_ann','T_ann','Ei_ann','Eb_ann']}
    d=pd.DataFrame({'year':[1982,1982],'lat':[0.,0.],'lon':[0.,0.],**data})
    with pytest.raises(ValueError,match='Duplicate'):
        validate_panel(d,DEFAULTS)


def test_kitagawa_closure_and_sign():
    d=[]
    pfts=['Shrubland','Savanna/Grassland','Cropland']
    for j,pft in enumerate(pfts):
        for zone,ai in [('dry',.3),('wet',1.)]:
            for k in range(20):
                d.append({'lon':(-160+j*30)+(k%4)*2,'lat':-50+10*(k//4),
                          'PFT':pft,'AI_A':ai,'slope':100+10*j+(0 if zone=='dry' else -20)})
    df=pd.DataFrame(d)
    r=kitagawa(df,'slope','AI_A',min_pixels=3,min_blocks=1)
    assert r['supported']
    assert abs(r['total']-r['within']-r['composition'])<1e-9
    assert r['total']<0 and r['within']<0
    assert abs(r['composition'])<1e-6


def test_block_bootstrap_and_paired_difference():
    rng=np.random.default_rng(34)
    n=120
    ai=np.linspace(.05,1.8,n)
    d=pd.DataFrame({'lon':np.linspace(-177,177,n),'lat':np.linspace(-48,48,n),
                    'AI_A':ai,'pml_slope':250-110*ai+rng.normal(0,4,n),
                    'gleam_slope':180-80*ai+rng.normal(0,4,n)})
    r=rho_ci(d,'pml_slope','AI_A',10,60,12)
    p=paired_rho_delta(d,'pml_slope','gleam_slope','AI_A',10,60,12)
    assert r['rho']<-.7
    assert np.isfinite(r['lo95'])
    assert np.isfinite(p['difference'])


def test_ice_approximates_known_linear_slopes():
    rng=np.random.default_rng(7)
    rows=[]
    for k in range(6):
        for year in range(1982,1995):
            lai=.8+k*.05+.04*(year-1982)
            row={'lon':float(k),'lat':float(k),'LAI_gs':lai,'Ta_ann':18.,'P_ann':440.,'Rs_ann':180.,'VPD_ann':1.,'PFT':'Cropland','AI_A':.45}
            rows.append(row)
    df=pd.DataFrame(rows)
    X=df[list(DEFAULTS['predictors'])].to_numpy()
    y=3.2*X[:,0]+.1*X[:,1]
    model=LinearRegression().fit(X,y)
    result=local_ice(model,df,DEFAULTS)
    assert len(result)==6
    assert np.allclose(result.sensitivity.dropna(),3.2,atol=1e-4)
