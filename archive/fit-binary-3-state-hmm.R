library(tidyverse)
library(msm)
library(here)
library(ivs)

source(here("load_mic_data.R"))


df_mics

df_drug

ast_lookup <- tibble(
  species = c("glabrata", "glabrata", "glabrata"),
  drug = c("micafungin", "fluconazole", "anadulafungin"),
  r = c(0.03, 16, 0.06)
)

# take a df of test/state results
# and a df of exposure results
# and generate a data frame with one row per pair of test results
# with start and stop dtaes of covariates added as censored states 

df_raw |>
  select(psn, date_enrol, swb_date, isol_species) |>
  unique() |>
  group_by(psn, swb_date, date_enrol) |>
  summarise(
    sp_present = if_else("glabrata" %in% isol_species, 1, 0)
  ) |>
  left_join(
    df_mics |>
      filter(
        species == "albicans",
        drug == "fluconazole"
      ) |>
      select(
        psn, swb_date, mic
      )
  ) |>
  ungroup() |>
  transmute(
    pid = psn,
    t = as.numeric(difftime(swb_date, date_enrol, units = "days")),
    sp_present = sp_present,
    mic = mic
  )

df_mod <- 
  full_join(
    # filter mic data and switch to relative time
df_raw |>
  select(psn, date_enrol, swb_date, isol_species) |>
  unique() |>
  group_by(psn, swb_date, date_enrol) |>
  summarise(
    sp_present = if_else("albicans" %in% isol_species, 1, 0)
  ) |>
  left_join(
    df_mics |>
      filter(
        species == "albicans",
        drug == "fluconaozle"
      ) |>
      select(
        psn, swb_date, mic
      )
  ) |>
  ungroup() |>
  transmute(
    pid = psn,
    t = as.numeric(difftime(swb_date, date_enrol, units = "days")),
    sp_present = sp_present,
    mic = mic
  ),
    # filter drug data
    df_drug |>
      filter(
        antifung_class == "Azole",
        !is.na(date_enrol)
      ) |>
      # add one to end dates - we are considering the first day of no exposure
      mutate(
        date_antifung_end = date_antifung_end + days(1), date_antifung_end
      ) |>
      # merge overlapping intervals with the ivs package
      mutate(range = iv(date_antifung_start, date_antifung_end)) |>
      reframe(
        range = iv_groups(range, abutting = TRUE),
        .by = c(psn, date_enrol, antifung_class)
      ) |>
      # make long
      transmute(
        pid = psn,
        t_start = 
          as.numeric(difftime(iv_start(range), date_enrol, units = "days")),
        t_end = 
          as.numeric(difftime(iv_end(range), date_enrol, units = "days"))
      ) |>
      pivot_longer(
        -c(pid),
        names_to = "covariate",
        values_to = "t"
      ) |>
      mutate(covariate = if_else(
        covariate == "t_start", 1, 0
      )),
    by = c("pid", "t")
  ) |>
  # if first covariate is unknown set to 0
  # and fill down other NAs
  arrange(pid, t) |>
  group_by(pid) |>
  mutate(pos = seq_len(n())) |>
  mutate(covariate = case_when(
    pos == 1 & is.na(covariate) ~ 0,
    TRUE ~ covariate
  )) |>
  select(-pos) |>
  fill(covariate) |>
  # zero t on lowest value
  mutate(t = t - min(t))




#

#
# df_mod |>
#   group_by(pid) |>
#   arrange(pid,t) |>
#   mutate(index = seq_len(n())) |>
#   filter((index == 1 & is.na(mic)))
#
# df_mod |>
#   group_by(pid) |>
#   filter(sum(!is.na(mic)) > 1) |>
#   arrange(pid, t) |>
#   mutate(index = seq_len(n())) |>
#   filter(!(index == 1 & is.na(mic))) -> df_stan
# # #
# # #
# df_stan |>
#   group_by(pid) |>
#   arrange(pid, t) |>
#   mutate(index = seq_len(n())) |>
#   filter(!(index == 1 & is.na(mic))) |>
#   select(-index) -> df_stan
# #
#
# df_stan |>
#   group_by(pid) |>
#   arrange(pid, t) |>
#   mutate(index = seq_len(n())) |>
#   filter((index == 1 & is.na(mic)))

breakpoint <- 4

df_mod |>
  transmute(
    pid = pid,
    t = t,
    covariate = covariate,
    state = case_when(
      sp_present == 0 ~ 0,
      sp_present & mic <= breakpoint ~ 1,
      sp_present & mic > breakpoint ~ 2,
      sp_present & is.na(mic) ~ 3,
      is.na(sp_present) ~ 4
    )
  ) |>
  unique() -> df_3s

# any dups at a timepoint?

df_3s |>
  group_by(pid) |>
  filter(n() > 1) |>
  group_by(pid, t) |>
  filter(n() > 1)




df_stan_3s <-
  df_3s |>
  group_by(pid) |>
  filter(n() > 1)|>
  arrange(pid, t) |>
  mutate(
    t_end = lead(t),
    state_end = lead(state)
  ) |>
  filter(!is.na(t_end)) |>
  unique() |>
  ungroup()

# check for dups

df_stan_3s |> group_by(pid,t,t_end) |> count() |> filter(n > 1)

# check for messed up segments

df_stan_3s |> group_by(pid) |> filter(t_end > lead(t))
df_stan_3s |> group_by(pid) |> filter(t_end < t)
df_stan_3s |> group_by(pid) |> filter(t_end != lead(t))

df_stan_3s


n_participants <- length(unique(df_stan_3s$pid))

n_obs_pairs_per_participant <- df_stan_3s |>
  group_by(pid) |>
  summarise(n = n()) |>
  pull(n)

t_start <- df_stan_3s$t
t_end <- df_stan_3s$t_end
state_start <- df_stan_3s$state
state_end <- df_stan_3s$state_end
cov1 <- df_stan_3s$covariate

data_list <- list(
  N = n_participants,
  n = nrow(df_stan_3s),
  n_obs = n_obs_pairs_per_participant,
  t_start = t_start,
  t_end = t_end,
  state_start = state_start,
  state_end = state_end,
  covs = matrix(cov1),
  debug = 0)

library(cmdstanr)
 
mod3s <- cmdstan_model(here("hmm3.stan"))

fit_3s <- mod3s$sample(
  data = data_list,
  chains = 1,
  iter_warmup = 500,
  iter_sampling = 500,
  refresh = 10
  # adapt_delta = 0.999,
  # step_size = 0.001,
  # max_treedepth = 15,
  # parallel_chains = 4
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
