require(tidyverse)

# take a df of test/state results
# and a df of exposure results
# and generate a data frame with one row per pair of test results
# with start and stop dtaes of covariates added as censored states 
# for multistate modelling
generate_model_df <- function(
    df_raw,
    df_mics,
    df_drugs,
    spec = "albicans",
    exposure_drug_class = "Azole",
    mic_drug = "fluconazole") {


  df_mod <- 
    full_join(
      # filter mic data and switch to relative time
  df_raw |>
    select(psn, date_enrol, swb_date, isol_species) |>
    unique() |>
    group_by(psn, swb_date, date_enrol) |>
    summarise(
      sp_present = if_else(spec %in% isol_species, 1, 0),
        .groups = "keep"
    ) |>
    left_join(
      df_mics |>
        filter(
          species == spec &
          drug == mic_drug,
        ) |>
        select(
          psn, swb_date, mic
        ),
        by = join_by(psn, swb_date)
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
          antifung_class == exposure_drug_class,
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

  return(df_mod)

}


# strip out 1st states where the first state_var is NA
# ffor a given participant
remove_NA_initial_state <- function(df, state_var = sp_present) {

n_rows_to_rm <- 
df |>
  group_by(pid) |>
  arrange(pid, t) |>
  mutate(index = seq_len(n())) |>
  filter((index == 1 & is.na({{state_var}}))) |>
  select(-index) |>
  nrow()

while(n_rows_to_rm > 0) {
  df |>
    group_by(pid) |>
    arrange(pid, t) |>
    mutate(index = seq_len(n())) |>
    filter(!(index == 1 & is.na({{state_var}}))) |>
    select(-index) -> df
 
  n_rows_to_rm <- 
  df |>
    group_by(pid) |>
    arrange(pid, t) |>
    mutate(index = seq_len(n())) |>
    filter((index == 1 & is.na({{state_var}}))) |>
    select(-index) |>
    nrow()

  }

return(df)
}

# check the long state modelling df (the one with one
# row per obs) for NA state var first row
# and duplicated pid-timepoints
check_long_modelling_df <- function(df, state_var = sp_present) {
  cat("checking for NA first row ...\n")

  pids_na_first_state <-
    df |>
    group_by(pid) |>
    arrange(pid, t) |>
    mutate(index = seq_len(n())) |>
    filter((index == 1 & is.na({{state_var}}))) |>
    pull(pid)

  if (length(pids_na_first_state > 0)) {
    stop("NA in first states, pid: ", paste(pids_na_first_state, collapse = " "))
  }
  cat("OK! \n")

  cat("checking for duplicated pid-timepoints...\n")

  pids_dup_t <-
  df |>
    group_by(pid) |>
    filter(n() > 1) |>
    group_by(pid, t) |>
    filter(n() > 1) |>
    pull(pid)

  if (length(pids_dup_t > 0)) {
    stop("Duplicated pid-timepoints, pid: ", paste(pids_dup_t, collapse = " "))
  }
  cat("OK! \n")
  return(df)
}


make_start_stop_df_from_long <- function(df, state_var = sp_present) {

df <-  
  df |>
  group_by(pid) |>
  filter(n() > 1)|>
  arrange(pid, t) |>
  mutate(
    t_end = lead(t),
    state_end = lead({{state_var}})
  ) |>
  filter(!is.na(t_end)) |>
  unique() |>
  ungroup()
   
  return(df)

}


  # check hmm start-stop modeling frame for
  # impossible segments that overlap
  # it doesn't worry about state vars
check_start_stop_modelling_df <- function(df) {

  cat("checking for duplicated segments ...\n")

  pids_dup_seg <-
  df |> 
    group_by(pid,t,t_end) |> filter(n() > 1) |>
    pull(pid)

  if (length(pids_dup_seg > 0)) {
    stop("Duplicated segments", paste(pids_na_first_state, collapse = " "))
  }
  
  cat("OK! \n")

  cat("checking for overlapping segments...\n")

  pids_overlap <-
  c(
  df |> 
    group_by(pid) |> 
    filter(t_end > lead(t)) |>
      pull(pid),
  df |> 
    group_by(pid) |>
    filter(t_end < t) |>
    pull(pid),
  df |> 
    group_by(pid) |> 
    filter(t_end != lead(t)) |>
    pull(pid)
    )


  if (length(pids_overlap > 0)) {
    stop("Messed up segments", paste(pids_overlap, collapse = " "))
  }
  cat("OK! \n")
  return(df)
}


# make stan data list from data frame for hmm
make_stan_data_list_from_df <- function(df, 
                                        cens_val = -1,
                                        debug = 0,
                                        covariates = c("covariate")) {
  
  data_list <- list(
    N = length(unique(df$pid)),
    n = nrow(df),
    n_obs = df |>
      group_by(pid) |>
      summarise(n = n(), .groups = "keep") |>
      pull(n),
    t_start = df$t,
    t_end = df$t_end,
    state_start = if_else(is.na(df$state_start), cens_val, df$state_start),
    state_end = if_else(is.na(df$state_end), cens_val, df$state_end),
    covs = as.matrix(df[covariates]),
    debug = debug
  )
  return(data_list)
}

####################################
## THESE ARE FOR MIC MODELS .... ###
####################################

generate_mic_model_data <- function(
    df_mics,
    df_drug,
    spec = c(
      "albicans",
      "glabrata", "dubliniensis",
      "parapsilosis"
    ),
    exposure_drug_class,
    mic_drug,
    code_antifungal_exposure_as) {
  if (!code_antifungal_exposure_as %in%
    c(
      "zero_at_first_exposure",
      "zero_for_onngoing_exposure"
    )) {
    stop("This function can code antifungal exposure as
    either start the clock at first exposure or
    start the clock when exposure finishes - if the latter
    then the days_since_last_drug_exposure will be 0 whilst
    exposure is ongoing, and -1 if no exposure.
    If the former, then days_since_first_drug_exposure
    will just be set to days since enroll_date to allow
    plotting.")
  }

  df_mics |>
    filter(
      species %in% {{ spec }},
      drug == {{ mic_drug }}
    ) |>
    left_join(
      df_drug |>
        select(-date_enrol) |>
        filter(antifung_class == {{ exposure_drug_class }}), !is.na(date_enrol) |>
        mutate(range = iv(date_antifung_start, date_antifung_end)) |>
        reframe(
          range = iv_groups(range, abutting = TRUE),
          .by = c(psn, date_enrol, antifung_class)
        ),
      by = join_by(psn, closest(swb_date >= date_antifung_start))
    ) ->
  df_out

  if (code_antifungal_exposure_as == "zero_at_first_exposure") {
    df_out |>
      group_by(psn) |>
      mutate(
        any_drug_exposure = any(!is.na(date_antifung_start)),
        zero_date_for_unexposed = date_enrol,
        first_antifung_exp = case_when(
          any_drug_exposure ~ min(date_antifung_start, na.rm = TRUE),
          TRUE ~ NA
        )
        # date_enrol +
        #   (difftime(max(swb_date), date_enrol, units = "days") / 2)
      ) |>
      ungroup() |>
      mutate(
        days_since_first_drug_exposure =
          case_when(
            !any_drug_exposure ~ as.numeric(
              difftime(swb_date, zero_date_for_unexposed, units = "days")
            ),
            any_drug_exposure ~ as.numeric(
              difftime(swb_date, first_antifung_exp, units = "days")
            )
          )
      ) -> df_out
  }

  if (code_antifungal_exposure_as == "zero_for_onngoing_exposure") {
    df_out |>
      ungroup() |>
      mutate(
        days_since_drug_exposure =
          case_when(
            is.na(date_antifung_start) ~ -1,
            date_antifung_end >= swb_date & date_antifung_start <= swb_date ~ 0,
            date_antifung_end < swb_date ~ as.numeric(
              difftime(swb_date, date_antifung_end, units = "days")
            )
          )
      ) |>
      arrange(psn, swb_date) -> df_out
  }

  return(df_out)
}

make_mic_model_list_for_stan <- function(df, scale = FALSE) {
    
  list(
    n = nrow(df),
    n_participants = df |>
      pull(psn) |>
      unique() |>
      length(),
    pid = df |>
      left_join(
        df |>
        select(psn) |>
        unique() |>
        ungroup() |>
        mutate(pid = seq_along(1:n()))
      ) |>
      pull(pid),
    t_e = df |>
      pull(days_since_drug_exposure),
    y = df |>
      pull(mic) |> 
      log()
  ) -> list_out

if (scale) {
    list_out$y <- scale(list_out$y)[,1]
  }

  return(list_out)

}

make_mic_gp_model_list_for_stan <- function(df, scale = TRUE) {
  if (!scale) {
    stop("not scaling is not yet implemented!")
  }

  df <-
    df |>
    group_by(psn) |>
    mutate(t = as.numeric(
      difftime(swb_date, min(swb_date), unit = "days")
    ))

  list_out <-
    list(
      n = nrow(df),
      t = scale(df$t)[, 1],
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
        mutate(days_since_drug_exposure = days_since_drug_exposure / sd(df$t)) |>
        pull(days_since_drug_exposure) |>
        matrix(ncol = 1),
      y = df |>
        pull(mic) |>
        log()
    )

  list_out$y <- scale(list_out$y)[, 1]

  return(list_out)
}

cat("\n")
cat("---------------------------\n")
cat("> loaded helper functions from < \n")
cat("> modelling_helper_functions.R < \n")
cat("---------------------------\n")
