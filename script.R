## Load packages ---------------------------------------------------------------

library(tidyverse)
devtools::load_all()

## Simulated Example -----------------------------------------------------------

set.seed(103)

d_sim <- ltm_sim(
  ns = 500, nk = 2,
  vmu = matrix(c(.5,.5), nrow = 2),
  mPhi = diag(2) * c(.99, .99),
  mSigs = c(.1,.1),
  dsig = .15,
  vd = matrix(c(.4,.4), nrow = 2)
)

# adding zeroed beta
d_sim$mx <- cbind(d_sim$mx, runif(500)-.5)
d_sim$mb <- cbind(d_sim$mb, 0)

p_sim <- d_sim$mb %>%
  as.data.frame() %>%
  as_tibble() %>%
  set_names(paste0("beta", 1:(ncol(.)))) %>%
  rowid_to_column() %>%
  gather(beta, val, -rowid) %>%
  ggplot(aes(x = rowid, y = val)) +
  geom_line() +
  facet_wrap(~beta, ncol = 2) +
  geom_hline(yintercept = c(-1,1) * .4, linetype = 2) +
  ggtitle("Simulated")

p_sim

result <- ltm_mcmc(d_sim$mx, d_sim$vy, burnin = 2000, iter = 8000, K = 3)
# readr::write_rds(result, "data-raw/result.rds", compress = "xz")

## Results ---------------------------------------------------------------------

# Results after 2000 burnin and 8000 iterations.

result <- read_rds("data-raw/result.rds")

### Summary statistics ---------------------------------------------------------

# (like Table 1 in the paper)

summary_fun <- function(.x) {
  v1 <- map_dbl(.x, first)
  q <- quantile(v1, c(.05, .95))
  tibble(mean = mean(v1), sd = sd(v1), q05 = q[1], q95 = q[2])
}

summary_table <- map_dfr(result[-1], summary_fun, .id = "parm") %>%
  mutate(true = c(.5, .99, .4, .15, .1)) %>%
  select(parm, true, everything()) %>%
  mutate_if(is.numeric, round, 4)

knitr::kable(summary_table)

## MCMC Chains -----------------------------------------------------------------
map_dfc(result[-1], ~map_dbl(.x, first)) %>%
  rowid_to_column() %>%
  gather(parm, val, -rowid) %>%
  ggplot(aes(x = rowid, y = val)) +
  geom_line(alpha = .9) +
  geom_hline(aes(yintercept = mean), data = summary_table, colour = "red") +
  facet_wrap(~parm, scales = "free_y") +
  labs(x = "Iteration", y = "Value")

## Estimated betas -------------------------------------------------------------

mediana <- apply(simplify2array(result$beta), c(1,2), median, simplify = FALSE)
qt1 <- apply(simplify2array(result$beta), c(1,2), quantile, .05, simplify = FALSE)
qt2 <- apply(simplify2array(result$beta), c(1,2), quantile, .95, simplify = FALSE)
beta_stats <- imap_dfr(list(mediana = mediana, qt1 = qt1, qt2 = qt2), ~{
  .x %>%
    as.data.frame() %>%
    as_tibble() %>%
    rowid_to_column() %>%
    gather(beta, valor, -rowid) %>%
    mutate(tipo = .y)
})

p2 <- beta_stats %>%
  spread(tipo, valor) %>%
  ggplot(aes(x = rowid, y = mediana)) +
  geom_line() +
  facet_wrap(~beta, ncol = 2) +
  geom_ribbon(aes(ymin = qt1, ymax = qt2), alpha = .2) +
  ggtitle("Estimated")

p2


## TODO ------------------------------------------------------------------------

# - [ ] Write pure R script (today)
# - [ ] Write for $b_t$ with $y_{it}$ ()
# - [ ] Test on real data
# - [ ] Hyperparameter control
# - [ ] Better output to use bayesplot
# - [ ] Better Documentation

