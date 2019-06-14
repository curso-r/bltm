## Load packages ---------------------------------------------------------------

library(tidyverse)
devtools::load_all()

## Simulated Example -----------------------------------------------------------

set.seed(103)

d_sim <- ltm_sim(
  ns = 500, nk = 2, ni = 2,
  vmu = matrix(c(.5,.5), nrow = 2),
  mPhi = diag(2) * c(.99, .99),
  mSigs = c(.1,.1),
  dsig = .15,
  vd = matrix(c(.4,.4), nrow = 2)
)

# adding zeroed beta
# d_sim$mx <- cbind(d_sim$mx, runif(500)-.5)
# d_sim$mb <- cbind(d_sim$mb, 0)

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

pega_chain <- function(.z, i, nk) {
  dx <- dim(.z[[1]])
  if (is.null(dx)) {
    if (length(.z[[1]]) == 1) {
      chain <- map_dbl(.z, ~.x[1])
    } else {
      chain <- map_dbl(.z, ~.x[i])
    }
  } else if (all(dx == c(nk, 1))) {
    chain <- map_dbl(.z, ~.x[i,1])
  } else if (all(dx == c(nk, nk))) {
    chain <- map_dbl(.z, ~.x[i,i])
  }
  chain
}

summary_fun <- function(.z, i, nk) {
  v1 <- pega_chain(.z, i, nk)
  q <- quantile(v1, c(.05, .95))
  tibble(mean = mean(v1), sd = sd(v1), q05 = q[1], q95 = q[2])
}

summary_table <- function(i, variaveis, nk) {
  map_dfr(result[-1], summary_fun, i = i, nk = nk, .id = "parm") %>%
    mutate(true = c(.5, .99, .4, .15, .1)) %>%
    select(parm, true, everything()) %>%
    mutate_if(is.numeric, round, 4) %>%
    filter(parm %in% variaveis)
}

mcmc_trace <- function(.z, nk) {
  chains <- purrr::map(seq_len(nk), ~pega_chain(.z, .x, nk))
  d_medianas <- tibble::tibble(
    mediana = purrr::map_dbl(chains, median),
    key = paste0("V", seq_along(mediana))
  )
  chains %>%
    bind_cols() %>%
    rowid_to_column() %>%
    gather(key, value, -rowid) %>%
    inner_join(d_medianas, "key") %>%
    ggplot(aes(x = rowid, y = value)) +
    geom_line() +
    facet_wrap(~key, scales = "free_y") +
    geom_hline(aes(yintercept = mediana), colour = "red")
}

montar_tabelas <- function(result, variaveis = c("mu", "sig_eta", "phi")) {
  nk <- nrow(result$mu[[1]])
  purrr::map(seq_len(nk), summary_table, variaveis = variaveis, nk = nk)
}
montar_graficos <- function(result, variaveis = c("mu", "sig_eta", "phi")) {
  nk <- nrow(result$mu[[1]])
  purrr::map(result[variaveis], mcmc_trace, nk = nk)
}

montar_graficos_beta <- function(result) {
  threshold <- montar_tabelas(result, "d") %>%
    bind_rows(.id = "key") %>%
    transmute(key = paste0("V", key), d = mean)
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
  }) %>%
    inner_join(threshold, c("beta" = "key"))
  beta_stats %>%
    spread(tipo, valor) %>%
    ggplot(aes(x = rowid, y = mediana)) +
    geom_line() +
    facet_wrap(~beta, ncol = 2) +
    geom_ribbon(aes(ymin = qt1, ymax = qt2), alpha = .2) +
    geom_hline(aes(yintercept = d), linetype = 2) +
    geom_hline(aes(yintercept = -d), linetype = 2)

}

tabelas <- montar_tabelas(result)
graficos_mcmc <- montar_graficos(result)
graficos_beta <- montar_graficos_beta(result)

# TABELAS RESUMO
# no seu caso iria de 1 a 9
knitr::kable(tabelas[[1]])

# MCMC CHAINS
# mu, sig_eta, phi
graficos$mu

# BETAS
# betas
graficos_beta


## TODO ------------------------------------------------------------------------

# - [x] Write pure R script (today)
# - [x] Write for $b_t$ with $y_{it}$
# - [ ] Test on real data
# - [ ] Hyperparameter control
# - [ ] Better output to use bayesplot
# - [ ] Better Documentation

