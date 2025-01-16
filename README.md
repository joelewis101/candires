# Candires

This is a repo to share my analysis of the candires data. For data protection 
reasons, no data are stored here, nor any output of the analysis. The analysis 
code is R code which contained in quarto documents.

The code assumes that the data is in a folder `data-raw/` - the script
`load_mic_data.R` loads and cleans the latest raw data file from this folder -
the filename of the data it loads is in that script, and
`modelling_helper_functions.R` contains a bunch of functions to mung the data
into the formats needed for the models.

There are two models coded in Stan:

- `hmm_bin_2s.stan` is a two state hidden Markov model in
Stan, with a matrix  implementation of the forward algorithm following Zucchini
2016 (Hidden Markov Models for Time Series 2nd ed.) and one covariate. The quarto file
`fit_stan_markov_mods.qmd` will fit this model to the data, save the fitted
models to the `data/` folder and generate a html with some model diagnostics.
The file `sim_from_fitted_models.R` will load the fitted models, run some sims
from the posterior, and save the results in `data/`

- `mic_model_v2_baseline.stan` is a regression of MIC data with a time-varying effect of
  antifungal exposure. See `candires_draft_final_analysis.qmd` for details. As
above, the quarto file `fit_stan_mic_mods_baseline.qmd` will fit the models, save the
output to `data/` and generate a html with some model diagnostics.

The file `candires_draft_final.qmd` loads the saved models,
simulations, and generates analysis plots and tables as a pdf.

The `archive/` folder contains previous iterations of models, scripts. These are
by no means guaranteed to run.

There are also:

  * A 3 state hmm model with support for censored states. This wasn't
included in the final analysis but the model file is `himm3_v2.stan` and can be
fit with `fit_stan_3s_markov_mods.qmd`

  * A MIC regression model that uses a more complex within-particiapnt
  correlation `gp_model.stan` that is fit with `fit_stan_gp_mic_mods.qmd` and is
  described in `candires_gp_model.qmd`

 
