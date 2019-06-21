## Load packages ---------------------------------------------------------------
library(tidyverse)
# install.packages(mcmcr)
if (!require(ltm)) devtools::install()
devtools::load_all()

## Simulated Example -----------------------------------------------------------
set.seed(1)

simular <- function() {
  d_sim <- ltm_sim(
    ns = 500, nk = 2, ni = 10,
    vmu = matrix(c(.5,.5), nrow = 2),
    mPhi = diag(2) * c(.99, .99),
    mSigs = c(.1,.1),
    dsig = .15,
    vd = matrix(c(.4,.4), nrow = 2),
    alpha = 0
  )
  binder <- array(runif(500)-.5, c(10, 500, 1))
  d_sim$mx <- abind::abind(d_sim$mx, binder, along = 3)
  d_sim$mb <- cbind(d_sim$mb, 0)
  d_sim
}

d_sim <- simular()

# faz 4 chains do mesmo modelo em paralelo
library(future)
plan(multiprocess, workers = 4)
res2 <- furrr::future_imap(1:4, ~{
  ltm_mcmc(d_sim$mx, d_sim$vy, burnin = 100, iter = 100, K = 3)
})

a_array <- res2 %>%
  abind::abind(along = 3) %>%
  aperm(c(3, 1, 2))

nm_all <- dimnames(a_array)[[3]]
nm <- nm_all[str_detect(nm_all, "d\\[phi|sig_eta\\[")]
rh <- a_array[,,nm] %>%
  `class<-`("mcmcarray") %>%
  mcmcr::rhat("term") %>%
  setNames(nm)

# grafico
bayesplot::mcmc_rhat(rh) +
  bayesplot::yaxis_text(hjust = 1) +
  bayesplot::theme_default()
