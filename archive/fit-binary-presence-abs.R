library(tidyverse)
library(msm)
library(here)
library(ivs)

source(here("load_mic_data.R"))
source(here("modelling_helper_functions.R"))


ast_lookup <- tibble(
  species = c("glabrata", "glabrata", "glabrata"),
  drug = c("micafungin", "fluconazole", "anadulafungin"),
  r = c(0.03, 16, 0.06)
)

  
generate_model_df(df_raw,
  df_mics,
  df_drugs,
  spec = "glabrata",
  mic_drug = "fluconazole",
  exposure_drug_class = "Azole"
) |>
  remove_NA_initial_state(state_var = sp_present) |>
  # a couple of repeated mics - removed
  select(pid, t, sp_present, covariate) |>
  unique() |>
  check_long_modelling_df() |>
  make_start_stop_df_from_long() |>
  check_start_stop_modelling_df() |>
  rename(state_start = "sp_present") ->
df_stan_bin

make_stan_data_list_from_df(df_stan_bin) -> data_list


library(cmdstanr)
 
mod_bin <- cmdstan_model(here("hmm_bin_2s.stan"))

fit_bin <- mod_bin$sample(
  data = data_list,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  refresh = 100,
  # adapt_delta = 0.999,
  # step_size = 0.001,
  # max_treedepth = 15,
  parallel_chains = 4
)

# okkkk - simulate some data using msms

df_mod2 |>
  filter(!is.na(t)) |>
  transmute(
    subject = pid,
    time = t,
    censor = if_else(mic == 99999, -1, 0),
    cov = covariate) -> df_sim

start <- rbinom(n = length(unique(df_sim$pid)), size = 1, p = 0.2) + 1



simmulti.msm(
  df_sim,
  qmatrix = rbind( c(0,0.5), c(0.3,0)),
  covariates= list(cov = c(0,-2)),
  hmodel = list(hmmLNorm(mean = -4, sd = 1), hmmLNorm(mean = 1, sd = 1)), 
  start = start
) |>
  transmute(
    pid = pid,
    t = time,
    mic = obs, #if_else(censor == -1, -1, obs),
    covariate = cov) -> dfsim2



df_stan_sim <-
  dfsim2 |>
  group_by(pid) |>
  mutate(n_mic_vals = sum(!is.na(mic))) |>
  filter(n_mic_vals > 1) |>
  group_by(pid,t) |>
  arrange(pid,t,mic) |>
  slice(1) |>
  group_by(pid) |>
  mutate(t = t - min(t)) |>
  arrange(pid, t) |>
  mutate(
    t_end = lead(t),
    mic_end = lead(mic)
  ) |>
  filter(!is.na(t_end)) |>
  unique()

# check for dups

df_stan_sim |> group_by(pid,t,t_end) |> count() |> filter(n > 1)

df_stan_sim

n_participants <- length(unique(df_stan_sim$pid))

n_obs_pairs_per_participant <- df_stan_sim |>
  group_by(pid) |>
  summarise(n = n()) |>
  pull(n)

t_start <- df_stan_sim$t
t_end <- df_stan_sim$t_end
mic_start <- df_stan_sim$mic
mic_end <- df_stan_sim$mic_end
cov1 <- df_stan_sim$covariate

data_list_sim <- list(
  N = n_participants,
  n = nrow(df_stan_sim),
  n_obs = n_obs_pairs_per_participant,
  t_start = t_start,
  t_end = t_end,
  mic_start = ifelse(is.na(mic_start), -1, mic_start),
  mic_end = ifelse(is.na(mic_end), -1, mic_end),
  covs = matrix(cov1),
  debug = 0
)


sim_fit <- mod$sample(
  data = data_list_sim,
  chains = 1
  # iter_warmup = 1000,
  # iter_sampling = 1000,
  # adapt_delta = 0.999,
  # step_size = 0.01
  # max_treedepth = 10,
  # parallel_chains = 4
)

# try censoring 5%

cens <- rbinom(n = nrow(dfsim2), 1, 0.05)

dfsim2_cens <- dfsim2

dfsim2_cens$mic[cens == 1] <- -1

df_stan_sim_cens <-
  dfsim2_cens |>
  group_by(pid) |>
  mutate(n_mic_vals = sum(!is.na(mic))) |>
  filter(n_mic_vals > 1) |>
  group_by(pid,t) |>
  arrange(pid,t,mic) |>
  slice(1) |>
  group_by(pid) |>
  mutate(t = t - min(t)) |>
  arrange(pid, t) |>
  mutate(
    t_end = lead(t),
    mic_end = lead(mic)
  ) |>
  filter(!is.na(t_end)) |>
  unique()

# check for dups

df_stan_sim_cens |> group_by(pid,t,t_end) |> count() |> filter(n > 1)

df_stan_sim_cens

n_participants <- length(unique(df_stan_sim_cens$pid))

n_obs_pairs_per_participant <- df_stan_sim_cens |>
  group_by(pid) |>
  summarise(n = n()) |>
  pull(n)

t_start <- df_stan_sim_cens$t
t_end <- df_stan_sim_cens$t_end
mic_start <- df_stan_sim_cens$mic
mic_end <- df_stan_sim_cens$mic_end
cov1 <- df_stan_sim_cens$covariate

data_list_sim_cens <- list(
  N = n_participants,
  n = nrow(df_stan_sim_cens),
  n_obs = n_obs_pairs_per_participant,
  t_start = t_start,
  t_end = t_end,
  mic_start = ifelse(is.na(mic_start), -1, mic_start),
  mic_end = ifelse(is.na(mic_end), -1, mic_end),
  covs = matrix(cov1),
  debug = 0
)



sim_fit_cens <- mod$sample(
  data = data_list_sim_cens,
  chains = 1,
  # iter_warmup = 1000,
  # iter_sampling = 1000,
  adapt_delta = 0.999,
  step_size = 0.01
  # max_treedepth = 10,
  # parallel_chains = 4
)
