ltm_sim <- function(ns, nk, ni, vmu, mPhi, mSigs, dsig, vd) {
  # variavel preditora
  mx <- array(runif(ns*nk*ni)-.5, dim = c(ni, ns, nk))
  # betas variando no tempo
  mb <- matrix(0, nrow = ns, ncol = nk)
  # valor do primeiro beta
  mb[1,] <- t(vmu + (chol(solve(diag(nk) - mPhi^2)) * mSigs) %*% rnorm(nk))
  # equação do beta
  for (i in seq_len(ns-1)) {
    mb[i+1,] <- t(vmu + mPhi %*% (t(mb[i,,drop=F]) - vmu) + mSigs * rnorm(nk))
  }
  # beta zerado
  mb_zerado <- mb * t(apply(abs(mb), 1, function(x) x > vd))
  # resposta

  xbs <- plyr::aaply(mx, 1, function(x) x * mb_zerado, .drop = FALSE)

  vy <- plyr::aaply(xbs, 1:2, sum) + dsig * rnorm(ns*ni)

  # resultados
  list(vy = vy, mx = mx, mb = mb, mb_zerado = mb_zerado)
}
