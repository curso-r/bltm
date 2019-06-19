
#' Simulate LTM model
#'
#' Simulate LTM model
#'
#' @param ns number of times
#' @param nk number of covariates
#' @param ni number of series
#' @param vmu vector mu
#' @param mPhi phi matrix
#' @param mSigs sigma eta vector
#' @param dsig general sigma
#' @param vd threshold parameter
#' @param alpha intercept
#'
#' @return list containing the generated y, x, beta and thresholded beta
#'
#' @export
ltm_sim <- function(ns, nk, ni, vmu, mPhi, mSigs, dsig, vd, alpha) {
  # variavel preditora
  mx <- array(stats::runif(ns*nk*ni)-.5, dim = c(ni, ns, nk))
  # betas variando no tempo
  mb <- matrix(0, nrow = ns, ncol = nk)
  # valor do primeiro beta
  mb[1,] <- t(vmu + (chol(solve(diag(nk) - mPhi^2)) * mSigs) %*% stats::rnorm(nk))
  # equação do beta
  for (i in seq_len(ns-1)) {
    mb[i+1,] <- t(vmu + mPhi %*% (t(mb[i,,drop=F]) - vmu) + mSigs * stats::rnorm(nk))
  }
  # beta zerado
  mb_zerado <- mb * t(apply(abs(mb), 1, function(x) x > vd))
  # resposta

  xbs <- plyr::aaply(mx, 1, function(x) x * mb_zerado, .drop = FALSE)

  vy <- alpha + plyr::aaply(xbs, 1:2, sum) + dsig * stats::rnorm(ns*ni)

  # resultados
  list(vy = vy, mx = mx, mb = mb, mb_zerado = mb_zerado)
}


ltm_pred <- function(post, newdata) {
  n_iter <- nrow(post)
  nm <- colnames(post)
  ni <- sum(grepl("alpha", nm))
  nk <- sum(grepl("phi", nm))
  purrr::map(seq_len(n_iter), ltm_pred_iter, post = post,
             nm = nm, ni = ni, nk = nk, newdata = newdata) %>%
    abind::abind(along = 3) %>%
    aperm(c(3, 1, 2))
}

ltm_pred_iter <- function(iter, post, nm, ni, nk, newdata) {

  # recuperando parametros
  post_i <- post[iter, ]
  mb <- matrix(post_i[grepl("beta", nm)], ncol = nk)
  mPhi <- post_i[grepl("phi", nm)] * diag(nk)
  mSigs <- post_i[grepl("sig_eta", nm)]
  vd <- post_i[grepl("d\\[", nm)]
  dsig <- post_i[grepl("sig\\[", nm)]
  alpha <- post_i[grepl("alpha\\[", nm)]
  vmu <- post_i[grepl("mu\\[", nm)]

  ns_future <- dim(newdata)[2]
  # betas variando no tempo
  mb_new <- matrix(0, nrow = ns_future, ncol = nk)
  # valor do primeiro beta
  mb_new[1,] <- t(vmu + (chol(solve(diag(nk) - mPhi^2)) * mSigs) %*% stats::rnorm(nk))
  # equação do beta
  for (i in seq_len(ns_future-1)) {
    mb_new[i+1,] <- t(vmu + mPhi %*% (t(mb_new[i,,drop=F]) - vmu) + mSigs * stats::rnorm(nk))
  }
  # beta zerado
  mb_zerado <- mb_new * t(apply(abs(mb_new), 1, function(x) x > vd))
  # resposta
  xbs <- plyr::aaply(newdata, 1, function(x) x * mb_zerado, .drop = FALSE)
  vy <- alpha + plyr::aaply(xbs, 1:2, sum) + dsig * stats::rnorm(ns_future*ni)
  vy
}






