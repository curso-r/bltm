## Load packages ---------------------------------------------------------------
library(tidyverse)
if (!require(ltm)) devtools::install()
devtools::load_all()

## Simulated Example -----------------------------------------------------------
set.seed(1)

simular <- function() {
  d_sim <- ltm_sim(
    ns = 500, nk = 2, ni = 5,
    vmu = matrix(c(.5,.5), nrow = 2),
    mPhi = diag(2) * c(.99, .99),
    mSigs = c(.1,.1),
    dsig = .15,
    vd = matrix(c(.4,.4), nrow = 2),
    alpha = 0
  )
  binder <- array(runif(500)-.5, c(5, 500, 1))
  d_sim$mx <- abind::abind(d_sim$mx, binder, along = 3)
  d_sim$mb <- cbind(d_sim$mb, 0)
  d_sim
}

d_sim <- simular()
ind_train <- 1:432
x_train <- d_sim$mx[,ind_train,]
y_train <- d_sim$vy[,ind_train]
x_test <- d_sim$mx[,-ind_train,]
y_test <- d_sim$vy[,-ind_train]

res <- ltm_mcmc(x_train, y_train, burnin = 2000, iter = 8000, K = 3)

# previsoes na base de treino
pred_post_train <- ltm_pred(res, x_train)
pred_mean_train <- apply(pred_post_train, 2:3, mean)
pred_sd_train <- apply(pred_post_train, 2:3, sd)

# previsoes na base de teste
pred_post_test <- ltm_pred(res, x_test)
pred_mean_test <- apply(pred_post_test, 2:3, mean)
pred_sd_test <- apply(pred_post_test, 2:3, sd)

prepare <- function(arr, lab = "train") {
  arr %>%
    as.data.frame() %>%
    as_tibble() %>%
    rowid_to_column() %>%
    gather(tempo, y, -rowid) %>%
    mutate(nm = lab, tempo = as.numeric(tempo))
}

# comparacoes
bind_rows(
  prepare(y_train, "data_train"),
  prepare(y_test, "data_test"),
  prepare(pred_mean_train, "pred_train"),
  prepare(pred_mean_test, "pred_test")
) %>%
  mutate(
    teste = str_detect(nm, "test"),
    pred = !str_detect(nm, "data"),
    tempo = if_else(teste & pred, tempo + 432, tempo)
  ) %>%
  filter(teste) %>%
  ggplot(aes(x = tempo, y = y, linetype = pred, colour = teste)) +
  geom_line() +
  theme_minimal(14) +
  facet_wrap(~factor(rowid))

(rmse_train <- sqrt(mean((pred_mean_train - y_train)^2)))
(rmse_test <- sqrt(sum((pred_mean_test - y_test)^2)))





