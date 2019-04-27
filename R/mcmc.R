ltm_mcmc <- function(x, y, burnin = 2000, iter = 8000, K = 3) {

  # variables -----
  if (length(dim(x)) == 2) x <- array(x, c(1, dim(x)[1], dim(x)[2]))
  if (is.null(dim(y))) y <- matrix(y, nrow = 1)

  ni <- dim(x)[1] # number of series
  ns <- dim(x)[2] # number of sample
  nk <- dim(x)[3] # number of variables

  # variaveis p/ salvar
  betas_f <- list()
  mu_f <- list()
  phi_f <- list()
  d_f <- list()
  sig_f <- list()
  sig_eta_f <- list()

  # initial values
  dsig <- 0.1
  mSigs <- rep(0.01, nk)
  vmu <- matrix(0, nrow = nk)
  mPhi <- 0.9 * diag(nk)
  vd <- matrix(0, nk)
  betas <- matrix(.1, ncol = nk, nrow = ns)
  mu <- matrix(0, nrow = nk)

  # priors -----

  # sig2 ~ IG(n0/2, S0/2)
  n0 <- 6; S0 <- 0.06
  # sig_eta ~ IG(v0/2, V0/2)
  v0 <- 6; V0 <- 0.06
  # mu ~ N(m0, s0^2)
  m0 <- 0; s0 <- 1
  # (phi+1)/2 ~ Beta(a0, b0)
  a0 <- 20; b0 <- 1.5
  # d < |mu| + K * v

  # mcmc loop ------

  for (j in seq(-burnin, iter)) {

    if (j %% 10 == 0) cat("iter", j, ": ")

    # betas
    for (t in 1:ns) {
      betas[t,] <- sample_beta_t(
        t, betas[t,], vd, dsig, matrix(x[,t,], nrow = ni), y[,t],
        betas[t-(t!=1),], betas[t+(t!=ns),],
        mu, ns, mPhi, mSigs
      )
    }

    # sigma_eta, phi, mu
    for (i in 1:nk) {
      mSigs[i] <- sample_sig_eta(v0, V0, ns, betas[,i], mu[i,], mPhi[i,i], vd[i,], K)
      mPhi[i,i] <- sample_phi(betas[,i], mu[i,], ns, mSigs[i], vd[i,], a0, b0, mPhi[i,i], K)
      mu[i,] <- sample_mu(mPhi[i,i], mSigs[i], betas[,i], vd[i,], s0, ns, m0, K)
    }

    # sig, threshold
    dsig <- sample_sig(n0, S0, ns, betas[-(ns+1),], y, x)
    vd <- sample_d(mu, K, mSigs, mPhi, vd, x, y, betas, dsig)

    # store values
    if (j > 0) {
      betas_f[[j]] <- betas
      mu_f[[j]] <- mu
      phi_f[[j]] <- mPhi
      d_f[[j]] <- vd
      sig_f[[j]] <- dsig
      sig_eta_f[[j]] <- mSigs
    }

    if (j %% 10 == 0) {
      m <- "mu=%.2f, phi=%.2f, d=%.2f, sig_eta=%.2f, sig=%.2f"
      message(sprintf(m, mu[1,1], mPhi[1,1], vd[1,1], mSigs[1], dsig))
    }
  }

  list(beta = betas_f, mu = mu_f, phi = phi_f,
       d = d_f, sig = sig_f, sig_eta = sig_eta_f)
}



