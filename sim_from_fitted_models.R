library(tidyverse)
library(here)
library(patchwork)
library(viridis)
library(cmdstanr)
library(bayesplot)
library(kableExtra)
library(expm)

# simulate candid prevalence from fitted models and save sims

fitted_mods <-
  readRDS(here("data/fitted_stan_models.rda"))

# make transition matrix, q, from 
# state transition and antibiotic effects
q <- function(q01, q10, b01, b10) {
  q <- matrix(
    data = c(-q01 * exp(b01), q01 * exp(b01), q10 * exp(b10), -q10 * exp(b10)),
    byrow = TRUE,
    nrow = 2,
    ncol = 2
  )
  return(q)
}

#given an inital probability of being in state 1 initial_p1,
#q matrix q and time t, output the probabiity of being in state 1 at time t
pr_state1 <- function(initial_p1, q, t) {
  f <- matrix(
    data = c(1 - initial_p1, initial_p1),
    ncol = 2, nrow = 1, byrow = TRUE
  )
  return((f %*% expm(q * t))[1, 2])
}


# calculate stationary values - start at this probability

initial_p1s <- 
bind_rows(
  lapply(
    fitted_mods,
    function(x) {
      bind_cols(
        species = x$species,
        drug = x$drug,
        bayesplot::mcmc_intervals_data(x$model_fit$draws(),
          prob_outer = 0.95
        )
      )
    }
  )
) |>
  filter(grepl("q", parameter)) |>
  select(species, drug, parameter, m) |>
  pivot_wider(
    id_cols = c(species, drug),
    names_from = parameter,
    values_from = m
  ) |>
  janitor::clean_names() |>
  rowwise() |>
  mutate(p1_initial = pr_state1(
    0.5,
    q(q_1, q_2, 0, 0),
    100
  )) 

# iterate over fitted models


t_steps <- 0:50
abx_days <- 10

sims_list_full_posterior <- list()
sims_list_summary <- list()
print("Simulating from posterior ...")

for (i in seq_len(length(fitted_mods))) {
  print(
    paste0("Model ", i, " of ", length(fitted_mods))
  )


  p1_initial <-
    initial_p1s |>
    filter(
      species == fitted_mods[[i]]$species,
      drug == fitted_mods[[i]]$drug
    ) |>
    pull(p1_initial)

  if (length(p1_initial) > 1) {
    stop("More than one initial p1?")
  }


  sim_df <-
    expand_grid(
      tibble(
        t = t_steps,
        abx = c(rep(1, abx_days), rep(0, length(t) - abx_days))
      ),
      fitted_mods[[i]]$model_fit$draws(
        format = "draws_matrix"
      ) |>
        as.data.frame()
    ) |>
    mutate(draw = rep(seq_len(nrow(d)), length(t_steps))) |>
    janitor::clean_names() |>
    rowwise() |>
    mutate(pr_1 = pr_state1(
      initial_p1 = p1_initial,
      q = q(q_1, q_2, beta_01_1 * abx, beta_10_1 * abx),
      t = t
    )) |>
    ungroup() |>
    mutate(
      species = fitted_mods[[i]]$species,
      drug = fitted_mods[[i]]$drug
    )

  sim_df_summary <-
    sim_df |>
    ungroup() |>
    group_by(species, drug, t) |>
    summarise(
      pr_1_med = median(pr_1),
      lci = quantile(pr_1, 0.025),
      uci = quantile(pr_1, 0.975),
      .groups = "keep"
    )

  sims_list_full_posterior[[i]] <- sim_df
  sims_list_summary[[i]] <- sim_df_summary
}

sims_list_summary_df <-
  bind_rows(sims_list_summary)

write_rds(
  sims_list_full_posterior,
  here("data/sims_list_full_posterior.rda"),
  compress = "gz"
)


write_rds(
  sims_list_summary_df,
  here("data/sims_list_summary_df.rda"),
)

