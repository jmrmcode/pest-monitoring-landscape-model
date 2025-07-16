# Landscape-scale Simulation and Evaluation of Insect Pest Abundance Sampling Strategies

This repository contains R code and spatial data associated with the manuscript on simulating and evaluating insect pest sampling strategies in greenhouse-dominated landscapes. The study supports the development of spatially informed monitoring systems using simulated abundance data and Bayesian spatial models.

## Contents

- `Requena_Mullor_et_al_R_code_clusterSampling.R`: R script for simulating and analyzing cluster sampling under comparable greenhouse landscape scenarios.
- `Requena_Mullor_et_al_R_code_randomSampling.R`: R script for simulating and analyzing random sampling under comparable greenhouse landscape scenarios.
- `shapefiles.zip`: Compressed folder containing all required shapefiles to reproduce the simulation and modeling results across different greenhouse configurations (dense, moderate, and scarce coverage).

## Requirements

Make sure you have the following installed:

- **R** version 4.0.0 or later
- R packages:
  - `INLA`
  - `INLAspacetime`
  - `raster`
  - `sf`

> **Note**: The `INLA` package is not available on CRAN. It must be installed from its official repository. You can install it using:
> ```r
> install.packages("INLA", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)
> ```

## Instructions

1. Clone or download this repository.
2. Extract the contents of `shapefiles.zip` into a directory (e.g., `sampling_simulation/`).
3. Open the R script corresponding to your sampling strategy:
   - `Requena_Mullor_et_al_R_code_clusterSampling.R`
   - `Requena_Mullor_et_al_R_code_randomSampling.R`
4. Edit the paths to the shapefiles in the script and output directory to match your local directory structure.
5. Edit the *range.fraction* attribute in *inla.barrier.pcmatern()* to select the desired level of greenhouse resistance to pest dispersal, i.e., 0.2 (strong) or 0.8 (weak).
6. Run the script in R. Model fitting and evaluation will be performed, and outputs will be saved as `.csv` files in the specified output folder.
7. To analyze different landscapes (e.g., scarce, dense, and moderate greenhouse density), modify the shapefile names and mesh parameters accordingly within the script.

## Output

Each script will generate performance metrics across sample sizes:
- **MAE** (Mean Absolute Error)
- **Posterior Predictive p-values**
These metrics are saved in `.csv` files named according to the landscape and barrier range used in the simulation.

## Citation

If you use this code or reproduce our results, please cite:

> Requena-Mullor, J. M. (2025). *Landscape-scale simulation and evaluation of insect pest sampling strategies around greenhouse environments*. DOI: [DOI placeholder]

## References

> Bakka, H., Vanhatalo, J., Illian, J.B., Simpson, D., Rue, H., 2019. Non-stationary Gaussian models with physical barriers. Spatial Statistics 29, 268–288. https://doi.org/10.1016/j.spasta.2019.01.002
> 
> Rue, H., Martino, S., Chopin, N., 2009. Approximate Bayesian inference for latent Gaussian models by using integrated nested Laplace approximations. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 71, 319–392. https://doi.org/10.1111/j.1467-9868.2008.00700.x

## License

This repository is licensed under the [MIT License](LICENSE).

---

Developed and maintained by Juan M. Requena-Mullor.

