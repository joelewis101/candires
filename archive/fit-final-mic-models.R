library(tidyverse)
library(cmdstanr)
library(bayesplot)
library(here)
library(patchwork)

source(here("load_mic_data.R"))


# helper function to get the data in the right format

prepare_data_in_stan_model_format <-
  function(df_mics = df_mics,
           df_drug = df_drug,
           exposure_drug_class, outcome_drug,
           output = "list",
           scale_outcome = FALSE) { 
  
  df_out <-
    df_mics |>
      left_join(
        df_drug |>
          select(-date_enrol) |>
          filter(antifung_class == {{ exposure_drug_class }}),
        join_by(psn, closest(swb_date >= date_antifung_start))
      ) |>
      mutate(
        days_since_drug_exposure =
          case_when(
            is.na(date_antifung_start) ~ NA_integer_,
            date_antifung_end >= swb_date ~ 0,
            date_antifung_end < swb_date ~ as.numeric(
              difftime(swb_date, date_antifung_end, units = "days")
            ),
            TRUE ~ -999
          )
      ) |>       
      filter(drug == {{ outcome_drug }}) 

  species_lookup <-
    tibble(
      species_id = seq_len(length(unique(df_out$species))),
      species = unique(fct_infreq(df_out$species))
    ) 

  data_list_out <- list(
    n = nrow(df_out),
    n_participants = length(unique(df_out$psn)),
    n_species = length(unique(df_out$species)),
    pid =
      left_join(
        df_out |>
          select(psn),
        tibble(
          psn = unique(df_out$psn),
          id = seq_along(psn)
        )
      ) |>
        pull(id),
    species_id =
      left_join(
        df_out |>
          select(species),
        species_lookup
      ) |>
        pull(species_id),
    y = log(df_out$mic),
    t_e = if_else(is.na(df_out$days_since_drug_exposure), -1,
      df_out$days_since_drug_exposure
    )
  )

  if (scale_outcome == TRUE) {
    data_list_out$y <- scale(data_list_out$y)[,1]
  }

  if (output == "list") {
      return(data_list_out)
    } else if (output == "df") {
    return(df_out) 
  } else if (output == "species lookup") {
      return(species_lookup)
  } else {
    stop("unknown output type in prepare_data_in_stan_model_format()")
  }

}

# helper function to plot the outputs

plot_model_output <- function(mcmc_draws, exposure_drug_class, outcome_drug) {
    
  species_lookup <-
        prepare_data_in_stan_model_format(
          df_mics = df_mics,
          df_drug = df_drug,
          exposure_drug_class = exposure_drug_class,
          outcome_drug = outcome_drug,
          output = "species lookup"
        ) 

  params <-
    mcmc_intervals_data(mcmc_draws,
      prob_outer = 0.95
    )

  plot_a <-
    params |>
    filter(grepl("beta.*,2]", parameter)) |>
    left_join(
      species_lookup |>
        mutate(parameter = paste0("beta[", species_id, ",2]")) |>
        select(species, parameter)
    ) |>
    ggplot(aes(fct_reorder(species, m, .desc = TRUE), m, ymin = ll, ymax = hh)) +
    geom_point() +
    geom_errorbar(width = 0) +
    coord_flip() +
    labs(
      title = "Drug effect",
      y = "Change in mean log(MIC) during drug exposure", x = ""
    )

  plot_b <-
    params |>
    filter(grepl("beta.*,1]", parameter)) |>
    left_join(
      species_lookup |>
        mutate(parameter = paste0("beta[", species_id, ",1]")) |>
        select(species, parameter)
    ) |>
    ggplot(aes(fct_reorder(species, m, .desc = TRUE), m, ymin = ll, ymax = hh)) +
    geom_point() +
    geom_errorbar(width = 0) +
    coord_flip() +
    labs(
      title = "Species effect",
      y = "Baseline (intercept) mean log(MIC) for species", x = ""
    )

  plot_c <-
    params |>
    filter(grepl("gamma\\[", parameter)) |>
    ggplot(aes(fct_reorder(parameter, m, .desc = TRUE), m, ymin = ll, ymax = hh)) +
    geom_point() +
    geom_errorbar(width = 0) +
    coord_flip() +
    labs(title = "Individual effect", y = "Individual effect on mean log(MIC)", x = "pid")

  plot_d <-
    params |>
    filter(grepl("sigma|Omega\\[2,1", parameter)) |>
    mutate(
      parameter =
        case_when(
          parameter == "sigma" ~ "SD of log(MIC)\ndistribution",
          parameter == "gamma_sigma" ~ "SD of within-individual\neffect",
          parameter == "beta_sigma[1]" ~ "SD of species\neffect",
          parameter == "beta_sigma[2]" ~ "SD of drug\neffect",
          parameter == "Omega[2,1]" ~ "Correlation coefficient of\nspecies and drug effects",
          TRUE ~ parameter
        )
    ) |>
    ggplot(aes(fct_reorder(parameter, m, .desc = TRUE), m, ymin = ll, ymax = hh)) +
    geom_point() +
    geom_errorbar(width = 0) +
    coord_flip() +
    labs(title = "Other model parameters", y = "Value", x = "")

  plot_e <-
    params |>
    filter(grepl("tau", parameter)) |>
    mutate(across(where(is.numeric), \(x) x * log(2))) |>
    mutate(
      parameter =
        case_when(
          parameter == "tau" ~ "Half life (days) of\n antibiotic effect",
          TRUE ~ parameter
        )
    ) |>
    ggplot(aes(fct_reorder(parameter, m), m, ymin = ll, ymax = hh)) +
    geom_point() +
    geom_errorbar(width = 0) +
    coord_flip() +
    labs(title = "Other model parameters", y = "Value", x = "")

  (plot_a + plot_b) / (plot_c + (plot_d / plot_e)) + 
    plot_annotation(title = paste0(exposure_drug_class,
      " exposure and effect on ",
      outcome_drug,
      " log(MIC)"))
}

model <- cmdstan_model(here("mic_model_final.stan"))

# azole- fluconazole

fit_azole_fluconazole <-
  model$sample(
    data =
      prepare_data_in_stan_model_format(
        df_mics = df_mics,
        df_drug = df_drug,
        exposure_drug_class = "Azole",
        outcome_drug = "fluconazole",
        output = "list",
    scale_outcome = TRUE
      ),
    chains = 4,
    parallel_chains = 4
  )

draws_azole_fluconazole <- fit_azole_fluconazole$draws()

plot_model_output(draws_azole_fluconazole, "Azole","fluconazole")


fit_azole_voriconazole <-
  model$sample(
    data =
      prepare_data_in_stan_model_format(
        df_mics = df_mics,
        df_drug = df_drug,
        exposure_drug_class = "Azole",
        outcome_drug = "voriconazole",
        output = "list",
    scale_outcome = TRUE
      ),
    chains = 4,
    parallel_chains = 4
  )

draws_azole_voriconazole <- fit_azole_voriconazole$draws()

plot_model_output(draws_azole_voriconazole, "Azole","voriconazole")



fit_echinocandin_micafungin <-
  model$sample(
    data =
      prepare_data_in_stan_model_format(
        df_mics = df_mics,
        df_drug = df_drug,
        exposure_drug_class = "Echinocandin",
        outcome_drug = "micafungin",
        output = "list",
        scale_outcome = FALSE
      ),
    chains = 4,
    parallel_chains = 4,
    adapt_delta = 0.9999,
    step_size = 0.001,
    max_treedepth = 20
  )

draws_echinocandin_micafungin <- fit_echinocandin_micafungin$draws()

plot_model_output(draws_echinocandin_micafungin, "Echinocandin",
  "micafungin")



fit_azole_micafungin <-
  model$sample(
    data =
      prepare_data_in_stan_model_format(
        df_mics = df_mics,
        df_drug = df_drug,
        exposure_drug_class = "Azole",
        outcome_drug = "micafungin",
        output = "list",
        scale_outcome = FALSE
      ),
    chains = 4,
    parallel_chains = 4,
    adapt_delta = 0.9999,
    step_size = 0.001,
    max_treedepth = 20
  )

draws_azole_micafungin <- fit_azole_micafungin$draws()

plot_model_output(draws_azole_micafungin, "Azole",
  "micafungin")



fit_echinocandin_fluconazole <-
  model$sample(
    data =
      prepare_data_in_stan_model_format(
        df_mics = df_mics,
        df_drug = df_drug,
        exposure_drug_class = "Echinocandin",
        outcome_drug = "fluconazole",
        output = "list",
        scale_outcome = FALSE
      ),
    chains = 4,
    parallel_chains = 4,
    adapt_delta = 0.9999,
    step_size = 0.001,
    max_treedepth = 20
  )

draws_echinocandin_fluconazole <- fit_echinocandin_fluconazole$draws()

plot_model_output(draws_echinocandin_fluconazole, "Echinocandin",
  "fluconazole")



fit_nystatin_fluconazole <-
  model$sample(
    data =
      prepare_data_in_stan_model_format(
        df_mics = df_mics,
        df_drug = df_drug,
        exposure_drug_class = "Nystatin",
        outcome_drug = "fluconazole",
        output = "list",
        scale_outcome = FALSE
      ),
    chains = 4,
    parallel_chains = 4,
    adapt_delta = 0.9999,
    step_size = 0.001,
    max_treedepth = 20
  )

draws_nystatin_fluconazole <- fit_nystatin_fluconazole$draws()

plot_model_output(draws_nystatin_fluconazole, "nystatin",
  "fluconazole")


fit_ampho_b_ampho_b <-
  model$sample(
    data =
      prepare_data_in_stan_model_format(
        df_mics = df_mics,
        df_drug = df_drug,
        exposure_drug_class = "Amphotericin B",
        outcome_drug = "amphotericin_b",
        output = "list",
        scale_outcome = FALSE
      ),
    chains = 4,
    parallel_chains = 4,
    adapt_delta = 0.9999,
    step_size = 0.001,
    max_treedepth = 20
  )

draws_ampho_b_ampho_b <- fit_ampho_b_ampho_b$draws()

plot_model_output(draws_ampho_b_ampho_b, "Amphotericin B",
  "amphotericin_b")


fit_nystatin_ampho_b <-
  model$sample(
    data =
      prepare_data_in_stan_model_format(
        df_mics = df_mics,
        df_drug = df_drug,
        exposure_drug_class = "Nystatin",
        outcome_drug = "amphotericin_b",
        output = "list",
        scale_outcome = FALSE
      ),
    chains = 4,
    parallel_chains = 4,
    adapt_delta = 0.9999,
    step_size = 0.001,
    max_treedepth = 20
  )

draws_nystatin_ampho_b <- fit_nystatin_ampho_b$draws()

plot_model_output(draws_nystatin_ampho_b, "Nystatin",
  "amphotericin_b")


