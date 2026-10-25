"""Optional harmonization entry points for 1-degree annual GeoTIFFs and annual PML NetCDF."""
from __future__ import annotations
from pathlib import Path
import re
import numpy as np
import pandas as pd
from .config import path
from .io import save,add_ai

NROWS,NCOLS=180,360

def _read_global_one_degree(file):
    import rasterio
    from rasterio.transform import from_origin
    f=Path(file)
    if not f.is_file():raise FileNotFoundError(f)
    with rasterio.open(f) as ds:
        expected=from_origin(-180,90,1,1)
        if ds.crs is None or ds.crs.to_epsg()!=4326 or not ds.transform.almost_equals(expected):
            raise ValueError(f'{f.name}: expected EPSG:4326 1-degree grid starting at upper-left (-180,90).')
        if ds.shape not in ((NROWS,NCOLS),(NROWS+1,NCOLS+1)):
            raise ValueError(f'{f.name}: expected 180x360 or 181x361 grid; got {ds.shape}.')
        # Some 1-degree GEE exports contain a redundant southeast edge row/column.
        # Crop just that extra row/column WITHOUT shifting pixels or interpolating values.
        arr=ds.read(1,window=rasterio.windows.Window(0,0,NCOLS,NROWS),masked=True).astype('float64').filled(np.nan)
    return arr

def raster_panel(cfg):
    """Build annual 1-degree panel from one TIFF per variable and year."""
    rasters=cfg.get('raster_inputs',{}).get('annual',{})
    req=set(cfg['predictors']) | set(cfg['targets'].values())
    if req-set(rasters):raise ValueError(f'raster_inputs.annual missing: {sorted(req-set(rasters))}')
    coords=np.indices((NROWS,NCOLS))
    lat=(89.5-coords[0]).reshape(-1)
    lon=(-179.5+coords[1]).reshape(-1)
    frames=[]
    for year in range(int(cfg['period']['start']),int(cfg['period']['end'])+1):
        dat={'year':np.full(NROWS*NCOLS,year,dtype=int),'lon':lon,'lat':lat}
        for field,pattern in rasters.items():
            filename=path(cfg,str(pattern).replace('{year}',str(year)))
            dat[field]=_read_global_one_degree(filename).reshape(-1)
        d=pd.DataFrame(dat)
        # Only valid predictor/response cells are included; missing values are not imputed.
        frames.append(d.dropna(subset=sorted(req)))
    out=pd.concat(frames,ignore_index=True)
    dest=path(cfg,cfg['inputs']['annual_panel'])
    save(out,dest)
    pixels=cfg.get('raster_inputs',{}).get('pixel',{})
    if pixels:
        attr={'lon':lon,'lat':lat}
        for name,file in pixels.items():attr[name]=_read_global_one_degree(path(cfg,file)).reshape(-1)
        attributes=add_ai(pd.DataFrame(attr))
        keep_cols=['lon','lat','PFT']+[c for c in ['AI_A','AI_B','AI_C','P_ERA','Ep_GLEAM','PET_Terra','P_Terra'] if c in attributes]
        if 'PFT' not in attributes:raise ValueError('raster_inputs.pixel requires a PFT grid.')
        save(attributes[keep_cols].dropna(subset=['PFT']),path(cfg,cfg['inputs']['pixel_attributes']))
    return out

def _annual_nc_to_grid(file,field):
    """Average valid 0.1-degree PML annual fluxes onto EPSG:4326 1-degree grid.

    Input units must already be annual millimetres (mm/year). No days multiplier.
    """
    import rasterio
    from rasterio.enums import Resampling
    from rasterio.warp import reproject
    from rasterio.transform import from_origin
    import xarray as xr
    with xr.open_dataset(file,mask_and_scale=True,decode_cf=True) as ds:
        names=[v for v in ds.data_vars if v.lower()==field.lower()]
        if len(names)!=1:raise ValueError(f'{file}: expected NetCDF variable {field}, found {list(ds.data_vars)}')
        da=ds[names[0]]
        lat_dim=next((v for v in ['lat','latitude','y'] if v in da.dims),None)
        lon_dim=next((v for v in ['lon','longitude','x'] if v in da.dims),None)
        if lat_dim is None or lon_dim is None:raise ValueError(f'{file}: 1D lat/lon coordinates required.')
        for other in set(da.dims)-{lat_dim,lon_dim}:
            if da.sizes[other]!=1:raise ValueError(f'{file}: multiple slices in {other}; expected annual data.')
            da=da.isel({other:0})
        da=da.transpose(lat_dim,lon_dim)
        lat=np.asarray(da[lat_dim].values,float)
        lon=np.asarray(da[lon_dim].values,float)
        if not(np.all(np.diff(lat)>0) or np.all(np.diff(lat)<0)):raise ValueError('Latitude not monotonic.')
        if not(np.all(np.diff(lon)>0) or np.all(np.diff(lon)<0)):raise ValueError('Longitude not monotonic.')
        arr=np.asarray(da.values,dtype=np.float32)
        if lat[0]<lat[-1]:arr=arr[::-1,:];lat=lat[::-1]
        if lon[0]>lon[-1]:arr=arr[:,::-1];lon=lon[::-1]
        dx=float(np.median(np.diff(lon)))
        dy=float(abs(np.median(np.diff(lat))))
        if not 0<dx<=1.01 or not 0<dy<=1.01:raise ValueError('Invalid source grid resolution.')
        units=str(da.attrs.get('units','')).lower().replace(' ','')
        if not units or not any(u in units for u in ('mm/year','mmyr','mmyear','mm/a','mm/yr')):
            raise ValueError(f'{file}: units={units!r}; expected annual millimetres. Check source first.')
        # Force masked source pixels to a nodata sentinel, so reproject average ignores them.
        src=np.where(np.isfinite(arr),arr,-9999.0).astype('float32')
        out=np.full((NROWS,NCOLS),-9999.,dtype='float32')
        reproject(source=src,destination=out,
                  src_transform=from_origin(float(lon[0]-dx/2),float(lat[0]+dy/2),dx,dy),
                  src_crs='EPSG:4326',dst_transform=from_origin(-180,90,1,1),
                  dst_crs='EPSG:4326',src_nodata=-9999.,dst_nodata=-9999.,
                  resampling=Resampling.average,num_threads=1)
        out[out==-9999]=np.nan
        return out

def pml_panel(cfg):
    root=path(cfg,cfg['pml_netcdf']['root'])
    if root is None or not root.exists():raise ValueError('Configure pml_netcdf.root with ET/Ec/Ei/Es yearly subfolders.')
    s=int(cfg['pml_netcdf']['start_year']);e=int(cfg['pml_netcdf']['end_year'])
    coords=np.indices((NROWS,NCOLS))
    lat=(89.5-coords[0]).reshape(-1)
    lon=(-179.5+coords[1]).reshape(-1)
    frames=[]
    for year in range(s,e+1):
        frame={'year':year,'lon':lon,'lat':lat}
        for variable in ('ET','Ec','Ei','Es'):
            directory=root/variable
            fs=sorted(x for x in directory.glob('*.nc*') if re.search(rf'(?<!\d){year}(?!\d)',x.stem))
            if len(fs)!=1:raise ValueError(f'Expected exactly one {variable} NetCDF for year {year} in {directory}; found {len(fs)}.')
            frame['PML_'+variable]=_annual_nc_to_grid(fs[0],variable).reshape(-1)
        frames.append(pd.DataFrame(frame).dropna(subset=['PML_ET']))
        print(f'PML {year} harmonized',flush=True)
    combined=pd.concat(frames,ignore_index=True)
    dest=path(cfg,cfg['inputs']['pml_panel'])
    if dest is None:raise ValueError('Set inputs.pml_panel output filename in the config.')
    save(combined,dest)
    return combined
