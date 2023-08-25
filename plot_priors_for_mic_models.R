library(tidyverse)
library(ggdist)
library(here)

source(here("load_mic_data.R"))

# assume that mics are lognormally distibuted 
# 

df_mics |>
  filter(species %in% unique(fct_infreq(species))[1:5]) |>
  ggplot(aes(log(mic))) +
  geom_histogram() +
  facet_grid(species ~ drug, scales = "free_y")

# all log(mic) values lie within -5 to +table5

# students t with 3 df and variance 5 looks like
tibble(
  x = seq(from = -8, to = 8, by = 0.1),
  dens = dstudent_t(x, mu = 0, sigma = 5, df = 3)
) |>
  ggplot(aes(x, dens)) +
  geom_line()


tibble(
  x = seq(from = -8, to = 8, by = 0.1),
  dens = pstudent_t(x, mu = 0, sigma = 5, df = 3)
) |>
  ggplot(aes(x, dens)) +
  geom_line()
