#' simulacao LTM especial
#'
#' 6 betas, sendo
#' (1,2) gerados normalmente
#' (3,4) zerados na metade
#' (5,6) totalmente zerados
#'
set.seed(10)

# precisa carregar o pacote ltm, ou colocar esse script na pasta
devtools::load_all()

# parte 1 ---------
parte1 <- ltm_sim(400, 2, 1, c(.5,.5), diag(2) * .99, c(.1,.1), .15, c(.4,.4))
# parte 2 ---------
## primeiras 200 obs
parte2 <- ltm_sim(200, 2, 1, c(.6,.6), diag(2) * .99, c(.1,.1), .15, c(.4,.4))
## valor de x com beta zerado, t > 200
x_200 <- array(runif(400)-.5, c(1, 200, 2))
## beta zerado, t > 200
zero_200 <- matrix(0, 200, 2)
## montando parte 2
parte2$mx <- abind::abind(parte2$mx, x_200, along = 2)
parte2$mb <- rbind(parte2$mb, zero_200)
parte2$mb_zerado <- rbind(parte2$mb_zerado, zero_200)
# parte 3 --------
## valor de x com beta zerado
x_400 <- array(runif(2*400)-.5, c(1, 400, 2))
## todos os betas sao zero
zero_400 <- matrix(0, 400, 2)
## montando parte 3
parte3 <- list(
  mx = x_400,
  mb = zero_400,
  mb_zerado = zero_400
)
# juntando os inputs --------
d_sim <- list(
  mx = abind::abind(parte1$mx, parte2$mx, parte3$mx, along = 3),
  mb = cbind(parte1$mb, parte2$mb, parte3$mb),
  mb_zerado = cbind(parte1$mb_zerado, parte2$mb_zerado, parte3$mb_zerado)
)
# gerando o y -----------
xbs <- plyr::aaply(d_sim$mx, 1, function(x) x * d_sim$mb_zerado, .drop = FALSE)
d_sim$vy <- plyr::aaply(xbs, 1:2, sum) + .15 * rnorm(400)

readr::write_rds(d_sim, "d_sim_new.rds")
save(d_sim, file = "d_sim_new.RData")



