from pathlib import Path
import numpy as np
import pytest
from lai_et.panel import _annual_nc_to_grid,_read_global_one_degree


def test_annual_pml_nc_aggregation(tmp_path):
    import xarray as xr
    lat=np.array([89.75,89.25],dtype=float)
    lon=np.array([-179.75,-179.25],dtype=float)
    da=xr.DataArray(np.array([[1.,2.],[3.,4.]],dtype='float32'),
                    dims=('lat','lon'),coords={'lat':lat,'lon':lon},attrs={'units':'mm/year'})
    file=tmp_path/'PML_ET_1982.nc'
    xr.Dataset({'ET':da}).to_netcdf(file,engine='scipy')
    ar=_annual_nc_to_grid(file,'ET')
    assert ar.shape==(180,360)
    assert np.isclose(ar[0,0],2.5,atol=.001)
    assert np.isnan(ar[0,1])


def test_bad_pml_units_fail(tmp_path):
    import xarray as xr
    lat=np.array([.5,-.5]);lon=np.array([.5,1.5])
    da=xr.DataArray(np.ones((2,2)),dims=('lat','lon'),coords={'lat':lat,'lon':lon},attrs={'units':'mm/day'})
    file=tmp_path/'PML_Ei_1982.nc';xr.Dataset({'Ei':da}).to_netcdf(file,engine='scipy')
    with pytest.raises(ValueError,match='units'):
        _annual_nc_to_grid(file,'Ei')


def test_gee_extra_edge_crop_without_resampling(tmp_path):
    import rasterio
    from rasterio.transform import from_origin
    array=np.zeros((181,361),dtype='float32')
    array[:180,:360]=7.0
    array[-1,:]=300.0
    array[:,-1]=300.0
    file=tmp_path/'terra_export_1degree.tif'
    with rasterio.open(file,'w',driver='GTiff',height=181,width=361,count=1,
                       dtype='float32',crs='EPSG:4326',transform=from_origin(-180,90,1,1)) as ds:
        ds.write(array,1)
    out=_read_global_one_degree(file)
    assert out.shape==(180,360)
    assert np.all(out==7.0)


def test_integer_pft_tif_mask_nodata(tmp_path):
    import rasterio
    from rasterio.transform import from_origin
    ar=np.full((180,360),3,dtype=np.int16)
    ar[0,0]=-9999
    file=tmp_path/'pft_1degree.tif'
    with rasterio.open(file,'w',driver='GTiff',height=180,width=360,count=1,
                       dtype='int16',nodata=-9999,crs='EPSG:4326',transform=from_origin(-180,90,1,1)) as ds:
        ds.write(ar,1)
    out=_read_global_one_degree(file)
    assert np.isnan(out[0,0])
    assert np.isclose(out[0,1],3.0)
