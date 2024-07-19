################
# Simulator
################

rmvnorm_prec <-
function( mu, # estimated fixed and random effects
          prec, # estimated joint precision
          n.sims) {

  require(Matrix)
  # Simulate values
  z0 = matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
  # Q = t(P) * L * t(L) * P
  L = Cholesky(prec, super=TRUE)
  # Calcualte t(P) * solve(t(L)) * z0 in two steps
  z = solve(L, z0, system = "Lt") # z = Lt^-1 * z
  z = solve(L, z, system = "Pt") # z = Pt    * z
  return(mu + as.matrix(z))
}

## Simulate a 2D AR1 Poisson process
sim_spde_dat <- function(n, sparse=TRUE, map=NULL){
  require(RTMB)
  set.seed(n)
  t0 <- Sys.time()
  n_x = n_y = n
  D_xx = abs(outer(1:n_x, 1:n_x, FUN="-"))
  A_xx = Matrix( ifelse(D_xx==1, 1, 0) )
  Q_xx = -0.4 * A_xx
  diag(Q_xx) = 1 + 0.4^2
  Q_zz = kronecker( Q_xx, Q_xx )
  z = rmvnorm_prec( n=1, mu=rep(0,nrow(Q_zz)), prec=Q_zz )
  lambda = exp( 2 + as.vector(z) )
  ## Simulate nuissance parameter z from oscillatory (day-night) process
  Data = data.frame( expand.grid(x=1:n_x, y=1:n_y), z=as.vector(z), lambda=lambda )
  Data$n = rpois( nrow(Data), lambda=Data$lambda )
  ## make mesh
  mesh = fmesher::fm_mesh_2d( Data[,c('x','y')] )
  spde = fmesher::fm_fem( mesh, refine=FALSE )
  ## #############
  ## RTMB objects
  ## ##############
  data <- list( n = Data$n, meshidxloc = mesh$idx$loc )
  parameters <- list(beta0=2, log_tau=-1.75, log_kappa=.3, x=rnorm(mesh$n))
 # map <- NULL#list(beta0=factor(1), log_tau=factor(NA), log_kappa=factor(NA))
  ## Objective function
  Q_spde <- function(spde, kappa) {
    kappa_pow2 <- kappa * kappa
    kappa_pow4 <- kappa_pow2 * kappa_pow2
    kappa_pow4 * spde$M0 + 2 * kappa_pow2 * spde$M1 + spde$M2    ## M0=G0, M1=G1, M2=G2
  }
  if(sparse){
    data$spde <- list( "M0" = spde$c0, "M1" = spde$g1, "M2" = spde$g2)
    f <- function(parms) {
      require(RTMB)
      getAll(parms, data)
      tau <- exp(log_tau)
      kappa <- exp(log_kappa)
      Q <- Q_spde(spde, kappa)
      nll <- -dgmrf(x, 0, Q, log=TRUE)        ## Negative log likelihood
      eta <- beta0 + x[meshidxloc] / tau
      nll <- nll - sum(dpois( n, exp(eta), log=TRUE ))
      return(nll)
    }
  } else {
    data$spde <- list( "M0" = as.matrix(spde$c0), "M1" = as.matrix(spde$g1), "M2" = as.matrix(spde$g2))
    f <- function(parms) {
      require(RTMB)
      require(Matrix)
      getAll(parms, data)
      tau <- exp(log_tau)
      kappa <- exp(log_kappa)
      Q <- Q_spde(spde, kappa)
      ## hack to force Q to be dense but retain the sparseMatrix
      ## type expected by dgmrf
      denseQ <- Matrix::sparseMatrix(i=row(Q), j=col(Q), x=as.numeric(Matrix::as.matrix(Q)))
      ## log.det.Q <- determinant(Q)$modulus
      ## Mu <- rep(0,length(x))
      ## nll <- as.numeric((x-Mu) %*% Q %*% (x - Mu)/2 - log.det.Q/2+log(2*pi)*length(x)/2)
      nll <- -dgmrf(x, 0, denseQ, log=TRUE)        ## Negative log likelihood
      eta <- beta0 + x[meshidxloc] / tau
      nll <- nll - sum(dpois( n, exp(eta), log=TRUE ))
      return(nll)
    }
  }
  obj <- RTMB::MakeADFun(f, parameters, random="x", silent=TRUE, map=map)
  return(obj)
}
