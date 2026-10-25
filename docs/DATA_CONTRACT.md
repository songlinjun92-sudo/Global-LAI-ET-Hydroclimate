# Data contracts and scientific conventions

## Input granularity

Research period 1982–2020 inclusive; one row per year and 1° grid center. No hidden files or user-specific computer paths are required. Every annual value must be checked for its native units and aggregation definition *before* insertion into the panel. The independent PML annual NetCDF parser expects provider-supplied `mm/year` values and performs spatial averaging only.

## Harmonized response columns

- `ET_ann`: annual GLEAM total evapotranspiration, mm/year
- `T_ann`: annual GLEAM transpiration, mm/year
- `Ei_ann`: annual GLEAM canopy interception evaporation, mm/year
- `Eb_ann`: GLEAM bare-soil evaporation, mm/year
- Optional `Es_ann`: GLEAM sublimation (not PML Es; avoid confusing names)
- Optional `Ew_ann`: GLEAM open-water evaporation
- `PML_ET`, `PML_Ec`, `PML_Ei`, `PML_Es`: PML-V2.2b independent annual ET, vegetation transpiration, interception and soil evaporation

`PML_Es` and `GLEAM_Eb` are only the closest available target analogues; they are not identical definitions. Annual fluxes should not be scaled by the number of days again. Different products may have different spatial support and component closure. No residual-balance identity is imposed on PML component fields.

## GIMMS LAI / explanatory variable definition

The predictor `LAI_gs` is the already harmonized growing-season GIMMS LAI4g measure. Do not substitute annual mean LAI or a generic fixed northern-hemisphere growing season without updating the scientific method. The scripts accept this annual predictor as provided and do not infer per-pixel phenology from raw monthly LAI.

## Aridity and vegetation

`AI_A`: annual-mean ERA5-Land P divided by annual-mean GLEAM Ep (original AI definition).
`AI_B`: ERA5-Land P divided by TerraClimate PET.
`AI_C`: TerraClimate precipitation divided by TerraClimate PET.

All three numerator/denominator rasters should be period-matched, annualized and on the same 1° grid. Missing or nonpositive denominators are invalid (NaN). Never label `AI_A`, `AI_B`, `AI_C` as interchangeable quantities; they are alternate definitions of moisture background.

Climate class thresholds:

| Class | Inclusive lower | Exclusive upper |
| --- | --- | --- |
| Arid | 0.00 | 0.20 |
| Semi-arid | 0.20 | 0.50 |
| Dry sub-humid | 0.50 | 0.65 |
| Humid | 0.65 | unbounded |

PFT codes: 1 Forest; 2 Shrubland; 3 Savanna/Grassland; 4 Cropland; 5 Wetland; 6 Barren/Sparse. Do not treat an empty or sparsely sampled PFT × climate cell as evidence of a zero response.

## Statistical definitions

- `local ICE`: OLS slope over repeated predictions with LAI varying from pixel historical q05 to q95, all other features held at pixel multiannual mean.
- `continuous_ai_gradients`: block-median Spearman rho; geographical blocks are **10° cells** by default.
- `paired_product_differences`: rho(PML) − rho(GLEAM), with *identical geographical blocks resampled together*.
- `kitagawa_decomposition`: symmetric decomposition of **humid minus semi-arid** ET sensitivity into within-PFT response and relative-composition components among the **three supported nonforest PFTs**.
- `annual fixed-model derivative`: centered finite-difference on each annual record; not the same estimand as historical ICE. HAC time p values do not replace spatial-block uncertainty assessment.

## What is not included

Automatic download/decompression and licensing of each third-party archive; undocumented conversion details for raw satellite and reanalysis products; ready-to-redistribute third-party observations; manuscript-specific figure design or final PDF typesetting. See the main README for input-generation entry points.
