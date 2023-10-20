
require(tidyverse)
require(here)
# split off micr measurements and drugs in to seperate dataframes 

path <- here("data-raw/long_df_mic_061023.csv")

df_raw <- read_csv(path)

df_mics <-
  df_raw |>
  select(c(psn, date_enrol, swb_date, isol_species, starts_with("mic"))) |>
  select(-mic_tested) |>
  unique() |>
  pivot_longer(
    cols = starts_with("mic"),
    names_to = "drug",
    names_pattern = "mic_(.*)",
    values_to = "mic",
    values_drop_na = TRUE
  ) |>
  rename(species = isol_species) 

df_drug <-
  df_raw |>
  select(
    psn,
    date_enrol,
    antifung_class,
    date_antifung_start,
    date_antifung_end
  ) |>
  filter(!is.na(antifung_class)) |>
  unique() |>
  mutate(date_antifung_start = date(date_antifung_start))

cat("\n")
cat("---------------------------\n")
cat("> loaded data from < \n")
cat(paste0("> ", path, " < \n"))
cat("---------------------------\n")
cat("> raw dataframe now in df_raw <\n")
cat("> mic data now in df_mics <\n")
cat("> drug data now in df_drug <\n")
cat("---------------------------\n")
