library(magrittr)
library(ggplot2)

load("vignettes/bias_res.RData")

t <- bias_res %>% dplyr::filter(!is.na(true_ftmic) & !is.na(model_ftmic)) %>%
  dplyr::mutate(bias = model_ftmic - true_ftmic)


quantile(t$bias, probs = c(.025,.5,.95,.975), na.rm=T)
