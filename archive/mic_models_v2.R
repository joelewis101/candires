


library(tidyverse)
library(msm)
library(here)
library(ivs)
library(bayesplot)
library(cmdstanr)
source(here("load_mic_data.R"))
source(here("modelling_helper_functions.R"))



generate_mic_model_data(
  df_mics,
  df_drug,
  spec = c("albicans", "glabrata", "dubliniensis", "parapsilosis"),
  exposure_drug_class = "Azole",
  mic_drug = "fluconazole",
  code_antifungal_exposure_as = "zero_at_first_exposure"
) |>
    group_by(psn) |>
  ggplot(aes(days_since_first_drug_exposure,
    log(mic),
    group = psn,
    # group = any_drug_exposure,
    color = any_drug_exposure
  )) +
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.5) +
  # geom_density2d(position = position_jitter(width = 0.1)) +
  # geom_line() +
  # geom_smooth(method = "lm", se = FALSE, size = 0.5, alpha = 0.3) +
  facet_grid(any_drug_exposure ~ species, scales = "free") +
  theme_bw()


generate_mic_model_data(
  df_mics,
  df_drug,
  spec = c("albicans", "glabrata", "dubliniensis", "parapsilosis"),
  exposure_drug_class = "Echinocandin",
  mic_drug = "anidulafungin",
  code_antifungal_exposure_as = "zero_at_first_exposure"
) |>
  ggplot(aes(days_since_first_drug_exposure,
    log(mic),
    group = psn,
    # group = any_drug_exposure,
    color = any_drug_exposure
  )) +
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.5) +
  # geom_density2d(position = position_jitter(width = 0.1)) +
  # geom_line() +
  # geom_smooth(method = "lm", se = FALSE, size = 0.5, alpha = 0.3) +
  facet_grid(any_drug_exposure ~ species, scales = "free") +
  theme_bw()


generate_mic_model_data(
  df_mics,
  df_drug,
  spec = c("albicans", "glabrata", "dubliniensis", "parapsilosis"),
  exposure_drug_class = "Azole",
  mic_drug = "fluconazole",
  code_antifungal_exposure_as = "zero_at_first_exposure"
) |>
  mutate(mic = case_when(
    species == "glabrata" & mic > 16 ~ 1,
    species == "glabrata" & mic <= 16 ~ 0,
    mic > 4 ~ 1,
    mic <= 4 ~ 0,
    TRUE ~ NA_real_
  )) |>
  ggplot(aes(days_since_first_drug_exposure,
    mic,
    group = any_drug_exposure,
    color = any_drug_exposure
  )) +
  geom_point(alpha = 0.5) +
  geom_smooth() +
  facet_grid(any_drug_exposure ~ species) 


generate_mic_model_data(
  df_mics,
  df_drug,
  spec = c("albicans", "glabrata", "dubliniensis", "parapsilosis"),
  exposure_drug_class = "Echinocandin",
  mic_drug = "anidulafungin",
  code_antifungal_exposure_as = "zero_at_first_exposure"
) |>
  mutate(mic = case_when(
    species == "albicans" & mic > 0.03 ~ 1,
    species == "albicans" & mic <= 0.03 ~ 0,
    species %in% c("glabrata", "tropicalis") & mic > 0.06 ~ 1,
    species %in% c("glabrata", "tropicalis") & mic <= 0.06 ~ 0,
    species == "parapsilosis" & mic > 4 ~ 1,
    species == "parapsilosis" & mic <= 4 ~ 0,
    TRUE ~ NA_real_
  )) |>
  ggplot(aes(days_since_first_drug_exposure,
    mic,
    group = any_drug_exposure,
    color = any_drug_exposure
  )) +
  geom_point(alpha = 0.5) +
  geom_smooth() +
  facet_grid(any_drug_exposure ~ species)


# ok model this mutha

generate_mic_model_data(
  df_mics,
  df_drug,
  spec = "albicans",
  exposure_drug_class = "Echinocandin",
  mic_drug = "anidulafungin",
  code_antifungal_exposure_as = "zero_for_onngoing_exposure"
) |>
  make_mic_model_list_for_stan() -> stan_data2


mic_mod <- cmdstan_model(here("mic_model_v2.stan"))



fit_anidula_glabrata2 <-
  mic_mod$sample(
    data = stan_data2,
    chains = 4,
    parallel_chains = 4
  )
