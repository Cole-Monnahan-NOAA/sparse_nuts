##TMB::runExample("simple",framework="TMBad",CPPFLAGS="-DTMBAD_INDEX_TYPE=uint64_t")
TMB::runExample("ar1xar1", clean = FALSE, framework="TMBad",CPPFLAGS="-DTMBAD_INDEX_TYPE=uint64_t")
## Approximating gaussian: N(xhat, Q^-1)
Q <- sdreport(obj, getJointPrecision=TRUE)$jointPrecision
xhat <- obj$env$last.par.best

obj$gr()
## rebuil without random effects
obj2 <- TMB::MakeADFun(data = obj$env$data, parameters = obj$env$parList(),
                       map = obj$env$map, random = NULL, silent = TRUE,
                       DLL = obj$env$DLL)

## rotate it with Q. x.cur is the rotated xhat value
rsparse <- adnuts:::.rotate_posterior('sparse', obj2$fn, obj2$gr, Q=Q, Qinv=solve(Q), y.cur=xhat)
grads_sparse <- rsparse$gr2(rsparse$x.cur)
rdense <- adnuts:::.rotate_posterior('dense', obj2$fn, obj2$gr, Q=Q, Qinv=solve(Q), y.cur=xhat)
grads_dense <- rdense$gr2(rdense$x.cur)

# math matches up
max(abs(grads_sparse-grads_dense))

## now try Kasper's fancy way
detach(package:TMB)
library(RTMB)
library(Matrix)
library(tmbstan)
chd <- Matrix::Cholesky(Q, super=TRUE, perm=TRUE)
L <- as(chd, "sparseMatrix")
perm <- chd@perm + 1L
iperm <- Matrix::invPerm(perm)
F <- RTMB::GetTape(obj)
## Q[perm,perm] = L %*% t(L)
## Q = L[iperm, ] %*% t(L[iperm, ])
obj3 <- RTMB::MakeADFun(function(u) {
  x <- solve(t(L), u)[iperm] + xhat
  REPORT(x)
  F(x)
}, numeric(length(xhat)))

## I thought this would match above? Same point and rotated gradient function?
max(abs(as.numeric(obj3$gr(xhat))-grads_sparse))

# RTMB:::vectorize(obj) ## Minor optimization - try with / without
# system.time(qw<-tmbstan(obj,chains=2, cores=2, seed=1))

