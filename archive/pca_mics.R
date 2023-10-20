library(tidyverse)
library(here)
library(patchwork)
library(FactoMineR)
library(pheatmap)
library(factoextra)
library(ggpubr)


source(here("load_mic_data.R"))

df_mics |>
  select(-date_enrol) |>
  group_by(psn,swb_date, species, drug) |>
  slice(1) |>
  ungroup() |>
  mutate(mic = log(mic)) |>
  pivot_wider(id_cols = c(psn, swb_date, species),
    names_from = drug,
    values_from = mic) -> df

pheatmap(
  cor(select(df, -c(psn, swb_date, species)), use = "pair")
  ) -> p_c

PCA(select(df, -c(psn, swb_date)), quali.sup = 1, graph = FALSE) -> p



fviz_pca_var(p, repel = TRUE,axes = c(1,2)) -> p_a

fviz_pca_ind(p, label = "none", habillage = 1, addEllipses = TRUE) -> p_b

p_c

p_a + p_b + plot_layout(ncol = 2)
