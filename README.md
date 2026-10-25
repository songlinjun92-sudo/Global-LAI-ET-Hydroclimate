Global LAI–ET Hydroclimate Analysis
Reproducible Python workflow to examine how hydroclimate and plant functional types (PFTs) organize model-derived growing-season LAI sensitivities of annual evapotranspiration (ET). The scientific analysis period is 1982–2020, with harmonized global 1° × 1° annual records. GLEAM is the primary ET product; PML-V2.2b is an independent ET-product sensitivity check.
Release scope: This repository is an independently structured, executable research-code implementation of the documented methods. Its synthetic demo and automated tests run without private files. Exact numerical reproduction of a given manuscript version additionally requires the identical upstream harmonized 1° input arrays, masks, training settings and product releases. This repository does not ship the observational datasets. The original conversion of upstream source products (e.g., growing-season LAI construction, ERA5-Land derivations, GLEAM v3.8a selection, remote-download steps) must follow the data processing methods and metadata in the paper; generic annual input contracts are provided here, not an undocumented reconstruction of every source-specific preprocessing choice.

1. Research workflow
1. Prepare growing-season LAI and annual ET/T/Ei/Eb and climate predictors as a 1° matched annual panel, or construct that panel from previously processed, exactly aligned annual 1° GeoTIFFs.
2. Build static pixel metadata: PFT and three aridity indices AI = P / PET:
   - AI_A = ERA5 P / GLEAM Ep;
   - AI_B = ERA5 P / TerraClimate PET;
   - AI_C = TerraClimate P / TerraClimate PET.
3. For each endpoint, train Random Forest (RF) models using predictors LAI_gs, Ta_ann, P_ann, Rs_ann, VPD_ann. Check spatially grouped CV and forward temporal holdouts, with and without LAI.
4. Compute pixel-local ICE over each pixel's historical q05–q95 LAI range, holding its multiannual average climate conditions fixed. Estimate the RF response slope using 20 ICE predictions. Exclude pixels whose historical LAI range is narrower than 0.10.
5. Analyze continuous AI gradients with Spearman correlations of 10° geographical-block medians, block bootstrap (2,000 replicates by default), PFT-specific results, and the semi-arid–humid symmetric Kitagawa decomposition of within-PFT versus composition contributions.
6. Independently harmonize PML-V2.2b annual NetCDF fluxes (ET, Ec, Ei, Es) from 1982–2020 to the same grid. Match PML ET ↔ GLEAM ET, PML Ec ↔ GLEAM T, PML Ei ↔ GLEAM Ei, and PML Es ↔ GLEAM Eb (last pair is approximate, not physically identical). Fit both products to the same grid-years and predictors for each target; compare AI trends and paired block resampling.
7. Produce reproducibility manifests, summary tables, and reproducible diagnostic comparison figures. The diagnostic figures are not exact reproductions of the manuscript layout.
Crucial scientific interpretation: RF/ICE quantities are conditional fitted-model sensitivities, not causal LAI interventions. The total ET and transpiration gradient shows cross-product agreement in the examined analysis; the interception gradient is product-dependent. High spatial CV skill does not establish strong year-to-year anomaly skill. PML and GLEAM are distinct ET products but may share LAI-related inputs; this is not entirely independent LAI-based observational validation.
2. Install
Conda (recommended, Windows/macOS/Linux)
conda env create -f environment.yml
conda activate lai-et-reproducible
Or create an existing Python 3.10+ environment and run:
python -m pip install -e ".[test]"
For optional SHAP reports:
python -m pip install -e ".[shap,test]"
3. Verify the code without research data
python -m pytest -q
python -m lai_et demo --root examples/demo_run
A deterministic, synthetic (not publishable) dataset is generated, and all main modeling/ICE/hydroclimate modules are exercised. The outputs are saved under examples/demo_run/outputs/demo/. Never cite the demo's synthetic numbers as study results.
4. Supply your research inputs
Create your real config:
cp configs/config.example.yaml configs/config.yaml
On Windows, you can manually duplicate config.example.yaml; no Unix utilities are necessary.

Edit configs/config.yaml and supply files in data/processed/. All relative paths resolve from the repository root; absolute machine-specific paths do not appear in the repository. The default values require 1982–2020 coverage and fully balanced per-target pixels; pixels without complete annual records for a modeled target are excluded rather than silently filled. For a GLEAM-only run without a PML panel, set inputs.pml_panel: null in configs/config.yaml; the example config otherwise points to a PML file and the loader expects that file to exist.

Annual panel: data/processed/annual_panel.csv.gz containing one record per (year, lon, lat). Required columns:
Column(s)	Definition
year, lon, lat	Calendar year and WGS84 pixel center coordinates (degrees)
LAI_gs	Growing-season LAI, dimensionless, from the specified GIMMS LAI4g processing
ET_ann, T_ann, Ei_ann, Eb_ann	GLEAM annual flux components, mm yr⁻¹
Ta_ann	ERA5-Land annual air temperature; consistent units throughout
P_ann	ERA5-Land annual precipitation, mm yr⁻¹
Rs_ann	Annual shortwave radiation, consistent units throughout
VPD_ann	Annual vapor pressure deficit, consistent units throughout


Optional annual fields (for expanded postprocessing) include Es_ann, Ew_ann, area weights and original grid identifiers. They are not required for the four main endpoints. See docs/DATA_CONTRACT.md for exclusions, units and data licensing. Keep original NaN/mask patterns; no implicit imputation is performed.

Pixel attributes: data/processed/pixel_attributes.csv with lon,lat,PFT,AI_A,AI_B,AI_C (or constituent annual mean P_ERA,Ep_GLEAM,PET_Terra,P_Terra, from which the code derives the three AIs). PFT accepts the official project codes 1 Forest, 2 Shrubland, 3 Savanna/Grassland, 4 Cropland, 5 Wetland, 6 Barren/Sparse or the equivalent category names.

Independent ET annual panel: data/processed/pml_panel.csv.gz with year,lon,lat,PML_ET,PML_Ec,PML_Ei,PML_Es, annual mm yr⁻¹, not daily rates. See step 5 for direct NetCDF support.
Optional: construct the annual panel from 1° GeoTIFFs
Uncomment/edit the variable patterns in raster_inputs.annual and raster_inputs.pixel in configs/config.yaml, then run:
python -m lai_et make-panel --config configs/config.yaml
Raster contract: EPSG:4326, 360 columns × 180 rows (or a 361 × 181 GEE export with only one redundant southeastern edge row and column), upper-left corner (-180°,90°), 1° pixels, north-up. Values are already annual or growing-season according to column names. Only the verified extra bottom row and rightmost column are cropped when present; no interpolation occurs. A different mismatch raises an error—this command does not silently resample LAI, ET or climate products. Constructing growing-season LAI from raw monthly source data is intentionally an upstream documented step, not a generic monthly mean assumed here.
Optional: directly aggregate official annual PML-V2.2b NetCDFs
Place provider files under:
data/raw/pml_annual/
  ET/   ...1982.nc ...2020.nc
  Ec/   ...1982.nc ...2020.nc
  Ei/   ...1982.nc ...2020.nc
  Es/   ...1982.nc ...2020.nc

Make sure each annual filename includes its correct four-digit year and each NetCDF contains its matching 2D annual variable and lon/lat coordinates. Set pml_netcdf.root and run:
python -m lai_et prepare-pml --config configs/config.yaml
Output is the configurable inputs.pml_panel CSV. The code uses valid-cell spatial averaging on the WGS84 1° grid and does not multiply annual flux values by 365 days. Any non-mm/year or ambiguous unit needs explicit inspection rather than silent conversion. Confirm component units and availability against the official PML-V2.2b documentation.
5. Run the harmonized-input analysis
python -m lai_et run --config configs/config.yaml
Key outputs (under outputs/publication_run/, excluded from Git):
models/rf_gleam_et.joblib                # plus other target/product models
rf_cross_validation.csv                # spatial/temporal, climate-only vs full
ice_et.csv.gz                          # plus T, Ei, Eb
pixel_ice_common.csv.gz                # merged pixel table; unsupported ICE slopes remain NaN
continuous_ai_gradients.csv            # 10° block rho & bootstrap CI
paired_product_differences.csv         # paired PML minus GLEAM rho CIs
pft_ai_gradients.csv                   # within-PFT gradients
kitagawa_decomposition.csv             # within-type vs composition
pml_gleam_flux_agreement.csv
method_manifest.json
figures/independent_et_gradients.png
figures/pft_ai_heatmap.png

The pft_ai_heatmap.png generator computes unweighted within-class pixel means; it is a diagnostic plot and must not be presented as the manuscript’s area-weighted PFT–AI heatmap without a matching weighted calculation. The merged pixel_ice_common.csv.gz retains NaN slopes for unsupported pixels; filter to nonmissing slopes for each particular endpoint before using it as a strict local-ICE support table.

model.train_cap: 0 fits on all matched rows (slower). Different training samples, PFT masks, or model settings can change the results. The reproducible defaults are an open-source workflow configuration, not a claim that all original manuscript model fits used exactly these settings. A new method_manifest.json records the actual settings of every run.

The optional temporal.enabled: true performs fixed-model finite difference per year, distinct from historical pixel-local ICE. HAC trends on year-level means are reported, but must not be presented as spatial-block-significant without a separate block-bootstrap temporal test.
6. Important scope and support rules
- Climate classes are Arid (AI < 0.20), Semi-arid (0.20 ≤ AI < 0.50), Dry sub-humid (0.50 ≤ AI < 0.65), Humid (AI ≥ 0.65). The continuous AI range extends beyond 1 where P/PET exceeds 1.
- The primary Kitagawa decomposition uses Semi-arid versus Humid, not Arid-to-Humid, because the very dry PFT groups have insufficient independent geographical support. It is restricted to Shrubland, Savanna/Grassland and Cropland. Unsupported groups are reported as NaN, not silently extrapolated.
- Different AI formulas can classify different cells into the same labeled category; sensitivity contrasts across AI schemes therefore do not automatically share identical samples.
- Missing years or masked cells are not filled. require_balanced: true imposes the configured year window for per-target modeling. Pixel support can differ between targets; pixel_ice_common is the explicit intersection.
- For matched PML–GLEAM comparisons, both models are fitted on identical matched records. A new RF model is fitted for GLEAM as well as PML, rather than comparing an archived GLEAM model with a differently supported PML model.
- SHAP is optional (model.do_shap: true) and is not itself a causal attribution method.
7. Public release checklist
- [ ] Replace the generic author entry in CITATION.cff, choose the definitive license with all contributors, and add the article DOI when available.
- [ ] Check configs/config.yaml and all text files for personally identifiable paths, credentials, emails and institution-specific file names.
- [ ] Confirm GLEAM, GIMMS LAI4g, PML-V2.2b, TerraClimate and ERA5-Land licensing and required citations. Do not upload their raw datasets without permission.
- [ ] Re-run python -m pytest -q and python -m lai_et demo on a fresh environment.
- [ ] Re-run the real analysis using the precisely documented source-product versions and attach an input-manifest / checksums to your release.
- [ ] Record which outputs correspond to manuscript tables/figures and inspect whether all expected numbers match the archived analysis.
- [ ] Add GitHub URL to your manuscript only after the repository is public; create a versioned GitHub Release and Zenodo DOI if possible.
8. Suggested paper text (after release)
Code availability. Python scripts for annual data harmonization, random forest modeling, local ICE sensitivity estimation, hydroclimatic and PFT stratification, spatial-block bootstrap inference, and independent ET-product comparison are available at [repository URL]. The repository provides a software environment specification, documented input contracts, a synthetic reproducibility test, and diagnostic figure-generation scripts. Original third-party datasets must be obtained from their respective providers subject to their access conditions.

9. Data and method references
Cite the exact dataset versions used, not only this repository. In particular, follow the PML-V2.2b provider's official data citation instructions and include the relevant methodology papers (e.g., the PML-V2 framework paper by Zhang et al., 2019) as applicable. For GLEAM, GIMMS LAI4g, ERA5-Land and TerraClimate, cite the official dataset descriptions and primary method references from the final manuscript's bibliography. This repository does not redistribute third-party data.
10. License
The source code in this repository is released under the MIT License. See the LICENSE file for details.

Third-party datasets, including GIMMS LAI4g, GLEAM, ERA5-Land, TerraClimate, and PML-V2.2b, are subject to their respective data access and licensing terms. These datasets are not redistributed under the MIT License.
