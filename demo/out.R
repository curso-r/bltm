# vamos partir do objeto "result", que e uma lista com
# 4 elementos, sendo cada elemento uma chain do modelo.

# Esse objeto contem todas as chains empilhadas
a <- do.call(rbind, result)

# tabelas ----------------------------------------------------------------------
tabelas <- a[,!str_detect(colnames(a), "beta\\[")] %>%
  as.data.frame() %>%
  tibble::as_tibble() %>%
  tibble::rowid_to_column() %>%
  tidyr::gather(key, val, -rowid) %>%
  dplyr::group_by(key) %>%
  dplyr::summarise(
    mediana = median(val),
    sd = sd(val),
    q025 = quantile(val, 0.025),
    q975 = quantile(val, 0.975)
  ) %>%
  dplyr::mutate(
    p = as.numeric(stringr::str_extract(key, "[0-9]+(?=\\])")),
    key = stringr::str_extract(key, "[a-z_]+"),
    is_alpha = stringr::str_detect(key, "alpha|sig$")
  ) %>%
  dplyr::arrange(p) %>%
  dplyr::select(key, p, dplyr::everything()) %>%
  dplyr::group_split(is_alpha, keep = FALSE)

alpha <- tabelas[[2]] %>%
  filter(key == "alpha") %>%
  group_by(key) %>%
  summarise(
    media_mediana = mean(mediana),
    media_sd = mean(sd),
    media_zero = mean((0 > q025)&(0 < q975))
  )

sigma <- tabelas[[2]] %>%
  filter(key != "alpha")

writexl::write_xlsx(
  list(parm = tabelas[[1]], alpha = alpha, sigma = sigma),
  "demo/tabela_resultsF.xlsx"
)

# grafico betas ----------------------------------------------------------------

D <- c(as.Date("1986/1/1"),
       as.Date("1994/4/1"),
       as.Date("2002/8/1"),
       as.Date("2010/12/1"))

p <- purrr::map(1:3, ~{
  plot_betas(a, .x, real_values = NULL)+
    scale_x_continuous(breaks=0:3*100, labels = D)+
    geom_line(linetype=1)+
    theme_bw()
})

dir.create("demo/graficosf")
purrr::iwalk(p, ~{
  ggsave(sprintf("demo/graficosf/grafico_%02d.png", .y), .x)
})

# vetor betas com training -----------------------------------------------------
betas_training <- a[,str_detect(colnames(a), "beta\\[")] %>%
  as.data.frame() %>%
  tibble::as_tibble() %>%
  tibble::rowid_to_column("iter") %>%
  tidyr::gather(key, val, -iter) %>%
  dplyr::mutate(
    time = as.numeric(stringr::str_extract(key, "(?<=\\[)[0-9]+")),
    p = as.numeric(stringr::str_extract(key, "(?<=\\[[0-9]{1,10},)[0-9]+"))
  ) %>%
  dplyr::select(iter, p, time, val) %>%
  tidyr::spread(p, val, sep = "") %>%
  # aqui eu tiro as medias a posteriori
  dplyr::group_by(time) %>%
  dplyr::summarise_at(dplyr::vars(dplyr::starts_with("p")), mean) %>%
  dplyr::ungroup()

# rhat -------------------------------------------------------------------------

# note que aqui ele usa como base a lista de arrays, nao o array empilhado
a_array <- result %>%
  abind::abind(along = 3) %>%
  aperm(c(3, 1, 2))

nm_all <- dimnames(a_array)[[3]]

# tabelas dos rhat parametros
parms <- c("d", "sig_eta", "mu", "phi")
purrr::map(parms, ~{
  nm <- nm_all[str_detect(nm_all, paste0(.x, "\\["))]
  rh <- a_array[,,nm] %>%
    `class<-`("mcmcarray") %>%
    mcmcr::rhat("term") %>%
    setNames(nm)
  rh
}) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  mutate(parm = parms) %>%
  select(parm, everything()) %>%
  set_names(c("Parm", paste0("P=", 1:(ncol(.)-1))))

# grafico histograma dos rhat dos betas
nm <- nm_all[str_detect(nm_all, paste0("beta", "\\["))]
rh <- a_array[,,nm] %>%
  `class<-`("mcmcarray") %>%
  mcmcr::rhat("term") %>%
  setNames(nm)
tibble(beta = rh) %>%
  ggplot(aes(x = beta)) +
  geom_histogram(fill = "lightblue", colour = "black", bins = 30) +
  geom_vline(xintercept = 1.1, colour = "red", linetype = 2) +
  labs(x = expression(hat(R)), y = "Count") +
  theme_bw(14)

# grafico histograma dos alphas
nm <- nm_all[str_detect(nm_all, paste0("alpha", "\\["))]
rh <- a_array[,,nm] %>%
  `class<-`("mcmcarray") %>%
  mcmcr::rhat("term") %>%
  setNames(nm)
tibble(alpha = rh) %>%
  ggplot(aes(x = alpha)) +
  geom_histogram(fill = "lightblue", colour = "black", bins = 30) +
  geom_vline(xintercept = 1.1, colour = "red", linetype = 2) +
  labs(x = expression(hat(R)), y = "Count") +
  theme_bw(14)

# predicao ---------------------------------------------------------------------

# previsoes na base de treino (demora)
pred_post_train <- ltm_pred(a, x_train)
# previsoes na base de teste  (demora)
pred_post_test <- ltm_pred(a, x_test)

# estatisticas
pred_mean_train <- apply(pred_post_train, 2:3, mean)
pred_median_train <- apply(pred_post_train, 2:3, median)
pred_sd_train <- apply(pred_post_train, 2:3, sd)

pred_mean_test <- apply(pred_post_test, 2:3, mean)
pred_median_test <- apply(pred_post_test, 2:3, median)
pred_sd_test <- apply(pred_post_test, 2:3, sd)

# RMSE
(rmse_train_mean <- sqrt(mean((pred_mean_train - y_train)^2)))
(rmse_test_mean <- sqrt(sum((pred_mean_test - y_test)^2)))
(rmse_train_median <- sqrt(mean((pred_median_train - y_train)^2)))
(rmse_test_median <- sqrt(sum((pred_median_test - y_test)^2)))

# R2
(r2_train_mean <- 1 - sum((pred_mean_train - y_train)^2) / sum(y_train^2))
(r2_test_mean <- 1 - sum((pred_mean_test - y_test)^2) / sum(y_test^2))
(r2_train_median <- 1 - sum((pred_median_train - y_train)^2) / sum(y_train^2))
(r2_test_median <- 1 - sum((pred_median_test - y_test)^2) / sum(y_test^2))

