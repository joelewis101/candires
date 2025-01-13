# Aim to fit my stan gaussian process model to the candires data


library(tidyverse)
library(msm)
library(here)
library(ivs)
library(bayesplot)
library(patchwork)

source(here("load_mic_data.R"))
source(here("modelling_helper_functions.R"))

library(cmdstanr)


species = c("albicans", "glabrata", "dubliniensis", "parapsilosis")
drug = c("Azole", "Echinocandin")

current_species <- species[[1]]
current_drug <- drug[[1]]

if (current_drug == "Azole") {
  current_mic_drug <- "fluconazole"
} else {
  current_mic_drug <- "anidulafungin"
}

df <-
  generate_mic_model_data(
    df_mics,
    df_drug,
    spec = current_species,
    exposure_drug_class = current_drug,
    mic_drug = current_mic_drug,
    code_antifungal_exposure_as = "zero_for_onngoing_exposure"
  ) |>
  group_by(psn) |>
  mutate(t = as.numeric(
    difftime(swb_date, min(swb_date), unit = "days")
  ))

list(
  n = nrow(df),
  t = scale(df$t)[,1],
  n_participants = df |>
    pull(psn) |>
    unique() |>
    length(),
  N = df |>
    group_by(psn) |>
    summarise(n = n()) |>
    pull(n),
  n_covariates = 1,
  t_e = df |>
    mutate(days_since_drug_exposure = days_since_drug_exposure/sd(df$t)) |>
    pull(days_since_drug_exposure) |>
    matrix(ncol = 1),
  y = df |>
    pull(mic) |>
    log()
) -> gp_mod_list

gp_mod_list$y <- scale(gp_mod_list$y)[,1] 


mic_mod <- cmdstan_model(here("gp_model.stan"))

gp_mod_list <-
  make_mic_gp_model_list_for_stan(df)

model_fit <- mic_mod$sample(
  data = gp_mod_list,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  refresh = 500,
  # adapt_delta = adapt_delta,
  # step_size = step_size,
  # max_treedepth = 15,
  parallel_chains = 4
  # output_dir = here("data")
)

d <- model_fit$draws()

mcmc_trace(d, regex_pars = "sigma|beta|tau|length|alpha")

mcmc_dens(d, regex_pars = "sigma|beta|tau|length|alpha")

# plot correlation

mod_sum_df <- model_fit$draws(format = "df")

cross_join(
  tibble(x = seq(0, 10, by = 0.1)),
  mod_sum_df |>
    select(c(length_scale, alpha, .draw))
) |>
  mutate(f = alpha^2 * exp(-x^2 / (2 * length_scale^2))) |>
  ggplot(aes(x, f, group = .draw)) +
  geom_line(alpha = 0.1)

  # group_by(x) |>
  # summarise(ll = quantile(f,0.025),
  #   l = quantile(f,0.25),
  #   m =quantile(f,0.5),
  #   h = quantile(f,0.75),
  #   hh = quantile(f,0.975))

