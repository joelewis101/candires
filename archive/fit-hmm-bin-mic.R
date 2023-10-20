
library(tidyverse)
library(msm)
library(here)
library(ivs)

source(here("load_mic_data.R"))
source(here("modelling_helper_functions.R"))


ast_lookup <- tibble(
  species = c(
    "glabrata", "glabrata", "glabrata",
    "albicans", "albicans", "albicans"
  ),
  drug = c(
    "micafungin", "fluconazole", "anadulafungin",
    "micafungin", "fluconazole", "anadulafungin"
  ),
  r = c(0.03, 16, 0.06, 0.016, 4, 0.03)
)

generate_model_df(df_raw,
  df_mics,
  df_drugs,
  spec = "glabrata",
  mic_drug = "fluconazole",
  exposure_drug_class = "Azole"
) |>
  mutate(mic = case_when(
    is.na(mic) ~ NA_integer_,
    mic > 16 ~ 1,
    mic <= 16 ~ 0
  )) |>
  filter(sp_present == 1) |>
  remove_NA_initial_state(state_var = mic) |>
  # a couple of repeated mics - removed
  select(pid, t, mic, covariate) |>
  unique() |>
  group_by(pid,t) |>
  arrange(desc(mic)) |>
  slice(1) |>
  check_long_modelling_df(state_var = mic) |>
  ungroup() |> 
  make_start_stop_df_from_long(state_var = mic) |>
  check_start_stop_modelling_df() |>
  rename(state_start = "mic") ->
df_stan_bin

make_stan_data_list_from_df(df_stan_bin) -> data_list


mod_bin <- cmdstan_model(here("hmm_bin_2s.stan"))

fit_bin <- mod_bin$sample(
  data = data_list,
  chains = 1,
  iter_warmup = 1000,
  iter_sampling = 1000,
  refresh = 100,
  adapt_delta = 0.999,
  step_size = 0.001,
  max_treedepth = 15,
  # parallel_chains = 4
)

