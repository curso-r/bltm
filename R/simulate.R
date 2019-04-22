ltm_sim <- function(ns, nk, vmu, mPhi, mSigs, dsig, vd) {
  # variavel preditora
  mx <- matrix(runif(ns*nk)-.5, ncol = nk)
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
  vy <- apply(mx * mb, 1, sum) + dsig * rnorm(ns)

  # resultados
  list(vy = vy, mx = mx, mb = mb, mb_zerado = mb_zerado)
}
