
## OLD DO NOT USE. Does not include the latest efficient perm rotation

source("code/startup.R")
source("code/load_tmb_objects.R")

objs <- list(obj.causal, obj.dlm, obj.gp_pois_regr, obj.petrel, obj.pollock, obj.salamanders, obj.sam,
     obj.sdmTMB, obj.simple, obj.swallows, obj.wham, obj.wildf)[2]
btmb <- lapply(objs, benchmark_metrics) |> bind_rows()

source("code/load_rtmb_objects.R")

bind_rows(benchmark_metrics(obj.causal, model='causal'),
          benchmark_metrics(obj.diamonds, model='diamonds')
          )


if ("RTMB" %in% .packages()) {
  detach("package:RTMB", unload=TRUE)
}
rm(list=ls())
## A TMB model that is big in dimesion and highly sparse
TMB::runExample("ar1xar1",framework="TMBad",CPPFLAGS="-DTMBAD_INDEX_TYPE=uint64_t")
#TMB::runExample("randomregression",framework="TMBad",CPPFLAGS="-DTMBAD_INDEX_TYPE=uint64_t")
##TMB::runExample("simple",framework="TMBad",CPPFLAGS="-DTMBAD_INDEX_TYPE=uint64_t", clean=TRUE)


opt <- with(obj, nlminb(par,fn,gr))
Q <- sdreport(obj, getJointPrecision = TRUE)$jointPrecision
# image(Cholesky(Q, LDL=FALSE, perm=FALSE))
# image(Cholesky(Q, LDL=FALSE, perm=TRUE))
# Rebuild without random effects
obj2 <-  TMB::MakeADFun(data=obj$env$data, parameters=obj$env$parList(),
                        map=obj$env$map,
                        random=NULL, silent=TRUE,
                        DLL=obj$env$DLL)

if ("TMB" %in% .packages()) {
  detach("package:TMB", unload=TRUE)
}
library(RTMB)
library(Matrix)
chd <- Matrix::Cholesky(Q, super=TRUE, perm=TRUE)
L <- as(chd, "sparseMatrix")
perm <- chd@perm + 1L
iperm <- Matrix::invPerm(perm)
F <- RTMB::GetTape(obj)
Lt <- Matrix::t(L)
## Q[perm,perm] = L %*% t(L)
## Q = L[iperm, ] %*% t(L[iperm, ])
obj3 <- RTMB::MakeADFun(function(y) {
  x <- Matrix::solve(Lt, y)[iperm] #+ xhat
  REPORT(x)
  F(x)
}, numeric(nrow(Q)), silent = TRUE)

# an arbitrary vector in original parameter space
x0 <- as.numeric(obj$env$last.par.best)+.121
# convert to transformed space
y0 <- as.vector(t(L) %*% x0[perm])
# convert back to original space as check (should be zero)
max(abs(Matrix::solve(t(L), y0)[iperm]-x0))
# these are the log-posterior functions Stan expects using RTMB's
# tape approach
f.tape <- function(x) -obj3$fn(x)
g.tape <- function(x) -as.numeric(obj3$gr(x))

# Try to do it without retaping, using easier math but less
# efficient code first
P <- as.matrix(0*Q)
for(i in 1:length(iperm))P[i,iperm[i]] <- 1
P <- as(P, 'sparseMatrix')
Lt <- t(L)
Pt <- Matrix::t(P)
f.inefficient <- function(y) -obj2$fn(as.numeric(P%*%solve(Lt)%*%y))
g.inefficient <- function(y)
  -as.numeric(obj2$gr(as.numeric(P%*% solve(Lt)%*%y)) %*% (P%*%solve(Lt)))
# checks out
f.tape(y0)-f.inefficient(y0)
max(abs(g.tape(y0)-g.inefficient(y0)))

# now for more efficient "solve" statements, in theory we don't
# have to invert L
#
# step 1: re-correlate y
y <- y0
p1 <- solve(Lt,y, system='L')
# how to get around transposing L here?
# p1 <- solve(L,y, system='Lt')
max(abs(solve(Lt)%*%y-p1))

# step 2: gradient of inverse-permuted vector
p2 <- obj2$gr(as.numeric(solve(Pt,p1, system='P')))
max(abs(obj2$gr(as.numeric(P%*%p1)) - p2))
# is it faster to do p1[iperm] instead?
max(abs(obj2$gr(as.numeric(p1[iperm])) - p2))

#step 3: multiply gradient vector by Jacobian
p3 <- solve(L, as.numeric(solve(P,t(p2))), system='Pt')
max(abs(as.numeric(p2 %*% P %*%solve(Lt)-p3)))

# put it all together
f.efficient <- function(y) {
  -obj2$fn(as.numeric(solve(Pt,solve(Lt,y, system='L'), system='P')))
}
g.efficient <- function(y){
  -solve(L, as.numeric(solve(P,t(obj2$gr(as.numeric(solve(Pt,solve(Lt,y, system='L'), system='P')))))), system='Pt')
}
f.efficient(y0)-f.tape(y0)
max(abs(g.efficient(y0)-g.tape(y0)))

# really naieve "dense" way of using the dense covariance and not
# sparse precision
Ldense <- as.matrix(t(chol(solve(Q))))               # lower triangular Cholesky decomp.
ydense <- as.numeric(solve(Ldense) %*% x0)
f.dense <- function(y) obj2$fn(Ldense %*% y)
g.dense <- function(y) {as.vector( obj2$gr(as.numeric(Ldense %*% y)) %*% Ldense )}

# Jim's sparse using "J" to reorder but ***not permuted for efficiency***
J <- Matrix::sparseMatrix( i=1:nrow(Q), j=nrow(Q):1 )
Lsparse <- Matrix::Cholesky(J%*%Q%*%J, super=TRUE, perm=FALSE)
Lsparseinv <- solve(as(Lsparse, "Matrix"))
ysparse <- as.numeric(J%*%chol(J%*%Q%*%J) %*% J%*%x0)
Linv_times_x <- function(L,x){
  as.numeric(J%*% Matrix::solve(L, Matrix::solve(L, J%*%x, system="Lt"), system="Pt"))
}
x_times_Linv <- function(L,x){
  as.numeric(J%*%Matrix::solve(L, Matrix::solve(L, Matrix::t(x%*%J), system="L"), system="Pt"))
}
f.sparse <- function(y){obj2$fn(Linv_times_x(Lsparse, y))}
g.sparse <- function(y) {
  as.vector(  x_times_Linv(Lsparse, obj2$gr(Linv_times_x(Lsparse, y))) )
}

library(microbenchmark)

n <- length(y0)
bench <- microbenchmark(original=obj2$gr(rnorm(n)),
                        dense=g.dense(rnorm(n)),
                        Jim_sparse=g.sparse(rnorm(n)),
                        RTMBtape=g.tape(rnorm(n)),
                        #inefficient=g.inefficient(rnorm(n)),
                        efficient=g.efficient(rnorm(n)),
                        unit='ms', times=50)
summary(bench) |> dplyr::arrange(mean) |> print(digits=1)
