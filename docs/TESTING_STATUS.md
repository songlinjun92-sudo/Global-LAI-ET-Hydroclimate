# Validation and known reproducibility boundaries

## Implemented and exercised

- Automated tests cover exact AI thresholds, three AI formulas, spatial block identifiers, duplicate panel rejection, symmetric Kitagawa closure, paired geographical bootstrap, pixel-local ICE on a known linear model, annual NetCDF to 1° averaging, and unit-check rejection.
- Eleven automated tests passed in the build environment.
- End-to-end synthetic demo exercised all four targets for both GLEAM and PML (8 RF fits) and wrote ICE, correlations, paired differences, PFT tables, decomposition, flux agreement and plots.
- Compatibility smoke run used the **column structure** of the pre-existing 1982–2020 annual panel on 300 sampled pixels / 11,700 grid-years; all GLEAM endpoints completed.

## Important pending validation by the corresponding author

- A full 1982–2020 run on the original large-area harmonized inputs with the definitive masking decisions and model settings.
- Confirm official-source units, masks and regridding for every upstream satellite and climate product; checksums for local raster downloads.
- Full-year PML-V2.2b NetCDF aggregation (the aggregation method and a small synthetic NetCDF fixture are tested; the complete provider archive is not bundled).
- Numerical comparison with every manuscript table, figure, historical archived model, and the original interpretation of yearly sensitivity.
- Release-level GitHub packaging, licensing and authorship consent.

The tests guarantee that the supplied code paths *run for tested inputs* and satisfy selected scientific invariants; they do not prove scientific equivalence with a prior unpublished local script set or guarantee that an unrelated version of every third-party dataset will be supported unchanged.
