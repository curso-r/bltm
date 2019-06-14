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

summary_table <- function(i, variaveis, nk, result) {
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
  purrr::map(seq_len(nk), summary_table, variaveis = variaveis, nk = nk, result = result)
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

grafico_betas <- function(a, betas = "all", scales = "free_y", real_values = NULL) {
  pegar_stat <- function(prob = .5, a, betas = "all") {
    nm <- stringr::str_subset(dimnames(a)[[2]], "beta\\[[^,]+,[0-9]")
    if (betas[1] == "all") {
      betas <- as.numeric(unique(stringr::str_extract(nm, "[0-9]+(?=\\])")))
    }
    a[,nm] %>%
      apply(2, quantile, prob) %>%
      as.numeric() %>%
      purrr::set_names(nm) %>%
      tibble::enframe() %>%
      dplyr::mutate(p = as.numeric(stringr::str_extract(name, "[0-9]+(?=\\])"))) %>%
      dplyr::filter(p %in% betas) %>%
      dplyr::mutate(id = as.numeric(stringr::str_extract(name, "(?<=\\[)[0-9]+")))
  }
  d <- a[,stringr::str_detect(colnames(a), "d\\[")] %>%
    apply(2, mean) %>%
    tibble::enframe() %>%
    dplyr::mutate(p = as.numeric(stringr::str_extract(name, "[0-9]+(?=\\])"))) %>%
    dplyr::filter(p %in% betas) %>%
    dplyr::select(p, d_mean = value)
  d_betas <- purrr::map_dfr(c(.1, .5, .9), pegar_stat, a, betas = betas, .id = "prob") %>%
    tidyr::spread(prob, value, sep = "_") %>%
    dplyr::inner_join(d, "p")


  p <- d_betas %>%
    ggplot(aes(x = id, y = prob_2)) +
    geom_line(linetype = 2) +
    geom_ribbon(aes(ymin = prob_1, ymax = prob_3), alpha = .2) +
    geom_hline(aes(yintercept = c(-1,1) * d_mean), linetype = 2) +
    facet_wrap(~factor(p), scales = scales)

  if (!is.null(real_values)) {
    betas_reais <- real_values$mb %>%
      as.data.frame() %>%
      tibble::as_tibble() %>%
      tibble::rowid_to_column("id") %>%
      tidyr::gather(p, prob_2, -id) %>%
      dplyr::mutate(p = readr::parse_number(p), d_mean = .4) %>%
      dplyr::filter(p %in% betas)
    p <- p +
      geom_line(data = betas_reais) +
      geom_hline(aes(yintercept = c(-1,1) * d_mean), data = betas_reais)
  }
  p
}
