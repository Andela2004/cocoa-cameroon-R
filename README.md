# Statistical analysis of cocoa production in Cameroon (2000-2024)

Exploratory, national-level analysis of annual cocoa production, harvested area,
mean temperature and precipitation, using first-differenced variables and GLS
models with AR(1)/AR(2) errors (R, package `nlme`).

## Data
- Production and harvested area: FAOSTAT (FAO), annual country data.
- Temperature and precipitation: World Bank Climate Change Knowledge Portal (ERA5).
- Data file: `data/cocoa_cameroon_2000_2024.csv`
  with columns: `year, harvested_area_ha, cocoa_production_t, temperature_c, precipitation_mm`.

## Model
dlog(Production_t) = b0 + b1 dTemperature_t + b2 dlog(Precipitation_t) + b3 dlog(HarvestedArea_t) + e_t
with e_t following AR(1) or AR(2). Models are fitted by maximum likelihood and compared with AIC.

## Run
```
Rscript scripts/analysis.R
```
Requires R with packages `nlme`, `lmtest`, `sandwich`.
Figures are written to `figures/`, numerical results to `outputs/model_results.txt`.

## Limitations
25 annual national observations; no regional or farm-level data; ERA5 is a gridded
reanalysis, not farm measurements; FAOSTAT area and production may be partly
estimated, so their relationship may be partly mechanical; results are associations,
not causal effects.
