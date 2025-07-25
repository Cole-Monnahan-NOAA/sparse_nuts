source(here('code/load_tmb_objects.R'))

# test initial values
obj <- obj.pollock

opt <- with(obj, nlminb(par,fn,gr))
Q <- sdreport(obj, getJointPrecision = TRUE)$jointPrecision
M <- as.matrix(solve(Q))
draws <- rmvnorm(n=1000, mu=obj$env$last.par.best, Sigma=M)

obj2 <- MakeADFun(data=obj$env$data, parameters=obj$env$parList(), map=obj$env$map, random=NULL, DLL=obj$env$DLL)

test <- apply(draws, 1, \(x) obj2$fn(x))

mcmc <- sample_sparse_tmb(obj, iter=10, chains=20, init='random', metric='unit')
obj2$fn(mcmc)
mcmc
# # chol decomp by hand from GMRF book
# Q <- matrix(rnorm(9), nrow=3)
# Q <- Q %*% t(Q)
# L0 <- t(chol(Q))
#
# n <- nrow(Q)
# v <- rep(NA, n)
# L <- 0*Q
# for(j in 1:n){
#   v[j:n] =Q[j:n,j]
#   for(k in 1:(j-1)){
#     if(k<1) break
#     v[j:n] <- v[j:n]- L[j:n,k]*L[j,k]
#   }
#   L[j:n,j]=v[j:n]/sqrt(v[j])
#   cat(j,k,"\n")
#   print(L)
# }

cov.all <- matrix(0, nrow=10, ncol=10)
Huu <- numDeriv::hessian(f=obj$env$f, x=obj$env$last.par.best)[1:8, 1:8]
ff <- function(theta){
  obj$fn(theta)
  obj$env$last.par
}
J <- numDeriv::jacobian(ff, x=opt$par)
V <- sdreport(obj)$cov.fixed
cov.all[1:8, 1:8] <- solve(Huu)
cov.all <- cov.all + J %*% V %*% t(J)
cov.all - solve(Q)
solve(cov.all)

library(Matrix)
r <- obj$env$random
nonr <- setdiff(seq_along(par), r)
H_Bhat <- optimHess(opt$par, obj$fn, obj$gr)
H_AA <- obj$env$spHess(obj$env$last.par.best, random = TRUE)
H_AB <- obj$env$f(par, order = 1, type = "ADGrad", keepx=nonr, keepy=r) ## TMBad only !!!
H_BA <- t(H_AB)
H_BB <- H_BA %*% solve(t(H_AA)) %*% H_AB + H_Bhat
Q <- Matrix(0, nrow = 10, ncol=10, sparse=TRUE)
Q[r,r] <- H_AA
Q[r,nonr] <- H_AB
Q[nonr,r] <- H_BA
Q[nonr, nonr] <- H_BB
max(abs(Q-sdrep$jointPrecision))

H_BB2 <- H_BA %*% solve(H_AA,H_AB) + hessian.fixed
H_BB-H_BB2

 <- rbind2(M1, M2)
M <- forceSymmetric(M, uplo = "L")
p <- invPerm(c(r, (1:length(par))[-r]))
# the joint precision Q
Q <- M[p, p]

library(Matrix)
Matrix::image(as(J, 'sparseMatrix'))
Matrix::image(as(V, 'sparseMatrix'))
Matrix::image(as(t(J), 'sparseMatrix'))
image(cov.all)



library(posteriordb)
library(tidyr); library(dplyr); library(ggplot2)
#my_pdb <- pdb_local(path='C:/Users/cole.monnahan/posteriordb-r/')
my_pdb <- pdb_github()
pnms <- posterior_names(my_pdb)
po <- posterior(pnms[i], my_pdb)

stan_code('mcycle_gp-accel_gp')
stan_data('hmm_example')

out <- lapply(pnms, function(nm){
  ref <- tryCatch(reference_posterior_draws(nm), error=function(e) 'none')
  if(is.character(ref)) return(NULL)
  df <- lapply(ref[1:length(ref)], \(x) as.data.frame(x)) |> bind_rows()
  cor <- cor(df)
  data.frame(name=nm, npars=nrow(cor), cor=cor[upper.tri(cor, diag=FALSE)])
})
out <- do.call(rbind, out)
ggplot(out, aes(name, y=abs(cor))) + geom_violin() + coord_flip()

outsum <- out |> group_by(name) |> summarize(npars=npars[1], maxcor=max(abs(cor)))
outsum |> arrange(desc(maxcor))
outsum |> arrange(desc(npars))
ggplot(outsum, aes(npars, maxcor)) + geom_point()





fits <- readRDS('results/swallows_fits.RDS')
mean(fits[[1]]$mle$Q==0)
tmp <- fits[[1]]$mle$Q
tmp[tmp!=0] <- 1
Matrix::image(tmp)
tmp <- as.matrix(tmp)

mean(tmp==0)
test <- sample_sparse_tmb(obj.swallows, iter=5000, warmup=1000,
                          metric='dense', chains=6, cores=6,
                          control=list(adapt_delta=.99))
pars <- c('sigmayearphi', 'yeareffphi_raw[1]','yeareffphi_raw[2]','yeareffphi_raw[3]','yeareffphi_raw[4]')[1:3]
pairs_admb(fits[[3]],  pars=pars)





## validate the model by running it on known distributions and
## checking for the correct properties
library(RTMB)
set.seed(23)
cov <- rWishart(n=1, df=2, Sigma=diag(2))[,,1]
cov2cor(cov)
pars <- list(x=c(0,0))
dat <- list(Sigma=cov)
f <- function(pars){
  getAll(pars,dat)
  -dmvnorm(x=x, mu=c(1,1), Sigma = Sigma, log=TRUE)
}
f <- function(pars){
  getAll(pars)
  -sum(dnorm(x=x, mean=c(1,1), sd = c(1,1), log=TRUE))
}
obj <- MakeADFun(f, pars, silent=TRUE)
obj$fn()
obj$gr()
opt <- with(obj, nlminb(par,gr))
sdrep <- sdreport(obj)
sdrep$cov.fixed

fits.validate <- fit_models(obj, iter=2000000, replicates=1, thin=10,
                            warmup=1000, cores=3, chains=3, globals=list(dat=dat),
                            metrics=c('unit', 'diag', 'dense')[1])
posts <- lapply(fits.validate, function(x) cbind(metric=x$metric, as.data.frame(x))) |>
  bind_rows() %>% setNames(c('metric', 'x1', 'x2')) %>%
  group_by(metric) %>%
  mutate(iter=1:n(), x1mean=cummean(x1)-1, x2mean=cummean(x2)-1) |> ungroup()

ggplot(filter(posts, iter>1000), aes(iter, abs(x1mean), color=metric)) + geom_line() + scale_y_log10() +
  scale_x_log10()

ggplot(filter(posts, iter>1000), aes(iter, x1mean, color=metric)) + geom_line() +
  scale_x_log10()


group_by(posts, metric) %>% filter(iter==max(iter))




library(ggplot2)
library(dplyr)
theme_set(theme_bw())
library(TMB)

TMB::runExample("simple",framework="TMBad",CPPFLAGS="-DTMBAD_INDEX_TYPE=uint64_t")
#TMB::runExample("ar1xar1", clean = FALSE, framework="TMBad",CPPFLAGS="-DTMBAD_INDEX_TYPE=uint64_t")
## Approximating gaussian: N(xhat, Q^-1)
Q <- sdreport(obj, getJointPrecision=TRUE)$jointPrecision
xhat <- obj$env$last.par.best

obj$gr()
## rebuil without random effects
obj2 <- TMB::MakeADFun(data = obj$env$data, parameters = obj$env$parList(),
                       map = obj$env$map, random = NULL, silent = TRUE,
                       DLL = obj$env$DLL)

## rotate it with Q. x.cur is the rotated xhat value
rsparse <- adnuts:::.rotate_posterior('sparse-naive', obj2$fn, obj2$gr, Q=Q, Qinv=solve(Q), y.cur=xhat)
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
  x <- solve(t(L), u)[iperm] #+ xhat
  REPORT(x)
  F(x)
}, numeric(nrow(Q)))
obj3$gr(obj$env$last.par.best)
rsparse_Jnoperm$gr2(rsparse_Jnoperm$x.cur)
## I thought this would match above? Same point and rotated gradient function?
max(abs(as.numeric(obj3$gr(xhat))-grads_sparse))

# RTMB:::vectorize(obj) ## Minor optimization - try with / without
# system.time(qw<-tmbstan(obj,chains=2, cores=2, seed=1))






library(TMB)
library(Matrix)

runExample('simple')
Q <- sdreport(obj, getJointPrecision = TRUE)$jointPrecision
image(Q)
ff <- function(Q, perm, super){
  L <-Matrix::Cholesky(Q, perm=perm, super=super)
  L2 <- as(L, 'sparseMatrix')
  pct.sparsity <- round(100*mean(L2[lower.tri(L2)] == 0),2)
  txt <- paste('perm=', perm, '; super=', super, ';\nsparsity=',pct.sparsity,'%')
  Matrix::image(L, main=txt)
}

ff(Q, perm=TRUE, super=TRUE)
ff(Q, perm=FALSE, super=TRUE)
ff(Q, perm=TRUE, super=FALSE)
ff(Q, perm=FALSE, super=FALSE)

mle <- obj$env$last.par.best
obj <- MakeADFun(data=obj$env$data, parameters = obj$env$parList(), DLL='simple')
rsparse <- adnuts:::.rotate_posterior('sparse-naive', obj$fn, obj$gr, Q=Q, Qinv=solve(Q), y.cur=mle)
rdense <- adnuts:::.rotate_posterior('dense', obj$fn, obj$gr, Q=Q, Qinv=solve(Q), y.cur=mle)

ind <- 112:118
ind <- 1:118
Qtest <- Q[ind, ind]
xtest <- as.numeric(mle[ind])# rnorm(nrow(Qtest))
M <- solve(Qtest)
Lm <- t(chol(M))
Lminv <- solve(Lm)
ytest <- as.numeric(Lminv %*% xtest)
J <- Matrix::sparseMatrix( i=1:nrow(Qtest), j=nrow(Qtest):1 )


plot_sparse_matrix <- function(x, main='', size=1, pch=16, plot=FALSE) {
  x <- t(as.matrix(as(x, 'sparseMatrix')))
  #x[x!=0] <- 1
  colnames(x)  <- nrow(x):1
  rownames(x) <- 1:nrow(x)
  g <- reshape2::melt(x) |>
    filter(value!=0) |>
    #ggplot(aes(x=Var1,y=Var2, fill=value,color=value))+
    #geom_tile()+xlab('x')+ylab('col')+
    ggplot(aes(x=Var1,y=Var2)) +
    geom_point(size=size, pch=pch) +
    labs(x=NULL, y=NULL, title=main) +
    #theme(legend.position = "none") +
    guides(x = "none", y = "none", color="none") +
    scale_fill_gradient2() +
    scale_color_gradient2()
  g
  if(plot) print(g)
  return(invisible(g))
}

size <- .8
pch <- 16

Lm <- chol(solve(Qtest)) |> t()

# Try without permutations but doing simpler sparse calcs
chd0 <- Matrix::Cholesky(J%*%Qtest%*%J, super=FALSE, perm=FALSE)
L0 <- as(chd0, "sparseMatrix")
L0inv <- solve(L0)
Ltilde0 <- t(J%*%L0%*%J)
Ltilde0inv <- solve(Ltilde0)
max(abs((L0%*%L0inv - diag(length(ind)))))
# can recover Q from LL'
max(abs(J%*%Qtest%*%J - (L0 %*% t(L0))))
# from x to y matches
ytest0 <- as.numeric(Ltilde0%*%xtest)
max(abs(ytest - ytest0))
max(abs(rsparse$x.cur-ytest0))
# from y to x matches
max(abs(xtest - as.numeric(Ltilde0inv %*% ytest0)))
# LP matches
rsparse$fn2(ytest)
rdense$fn2(ytest)
obj$fn(as.numeric(Ltilde0inv %*% ytest))
# gr matches
head(rsparse$gr2(ytest))
head(rdense$gr2(ytest))
gr0 <- function(y) as.numeric(obj$gr(as.numeric(Ltilde0inv %*% y)) %*% Ltilde0inv)
head(gr0(ytest))
g01 <- plot_sparse_matrix(Q, main='No perm: Q', size=size, pch=pch)
g02 <- plot_sparse_matrix(J%*%Q%*%J, main='No perm: JQJ', size=size, pch=pch)
g03 <- plot_sparse_matrix(L0, main='No perm: L=chol(JQJ)', size=size, pch=pch)
g04 <- plot_sparse_matrix(Ltilde0, main='No perm: Ltilde=(JLJ)\'', size=size, pch=pch)
g05 <- plot_sparse_matrix(Ltilde0inv, main='No perm: Ltildeinv', size=size, pch=pch)
cowplot::plot_grid(g01,g02,g03, g04,g05, ncol=1)


# Now with permutations
Qtilde <- J%*%Qtest%*%J
chd1 <- Matrix::Cholesky(Qtilde, super=TRUE, perm=TRUE)
# this is the optimal permutation for Qtilde
perm <- chd1@perm+1
iperm <- invPerm(perm)
# how to convert from perm to perm for Q?
P <- as.matrix(0*Qtest)
for(i in 1:length(perm))P[i,perm[i]] <- 1
P <- as(P, 'sparseMatrix')
A <- solve(J) %*% P %*% J
Ainv <- solve(A)
image(A)
image(P%*%J%*%Qtest %*%J%*%t(P))
image(J%*%A%*%Qtest%*%t(A)%*%J)

# this is the perm of Q
perm2 <- apply(as.matrix(A),1, \(x) which(x==1))
iperm2 <- invPerm(perm2)
Qperm <- Qtest[perm2, perm2]
xperm <- as.numeric(xtest[perm2])
as.numeric((Qperm %*% xperm)[iperm2])
as.numeric(Qtest %*% xtest)
image(Qperm)

# pre-order Qtest so perm is not needed after applying J
chd1perm <- Matrix::Cholesky(J%*%Qperm%*%J, super=TRUE, perm=FALSE)
#chd1perm <- Matrix::Cholesky(J%*%Qperm%*%J, super=FALSE, perm=FALSE)
image(chd1)
image(chd1perm)
chd1 <- chd1perm
L1 <- as(chd1, "sparseMatrix")
L1inv <- solve(L1)
Ltilde1 <- t(J%*%L1%*%J)
Ltilde1inv <- solve(Ltilde1)
image(Ltilde1)
image(Ltilde1inv)
max(abs((L1%*%L1inv - diag(length(ind)))))
max(abs((Ltilde1%*%Ltilde1inv - diag(length(ind)))))
# can recover Q from LL'
max(abs(J%*%Qperm%*%J - (L1 %*% t(L1))))

image(Ainv %*% Ltilde1inv)

# from x to y matches
yperm1 <- as.numeric(Ltilde1%*%xperm)
ytest1 <- yperm1[perm2]
max(abs(ytest - yperm1))
max(abs(rsparse$x.cur-ytest1))
# from y to x matches
max(abs(xperm-xperm1))
xtest1 <- as.numeric(solve(A) %*%as.numeric(Ltilde1inv %*% yperm1))
max(abs(xtest1-xtest))
# LP matches
rsparse$fn2(ytest)
rdense$fn2(ytest)
obj$fn(as.numeric(Ainv%*%Ltilde1inv %*% yperm1))
# gr matches
head(rsparse$gr2(ytest))
gr1 <- function(y)
  as.numeric(obj$gr(as.numeric(Ainv %*%Ltilde1inv %*% y)) %*% (Ainv %*%Ltilde1inv))
head(gr1(yperm1))
head(rdense$gr2(ytest))

g11 <- plot_sparse_matrix(Qperm, main='Perm: Qperm=AQA\'', size=size, pch=pch)
g12 <- plot_sparse_matrix(J%*%Qperm%*%J, main='Perm: JQpermJ', size=size, pch=pch)
g13 <- plot_sparse_matrix(L1, main='Perm: Lperm=chol(JQpermJ)', size=size, pch=pch)
g14 <- plot_sparse_matrix(Ltilde1, main='Perm: Ltilde=(JLpermJ)\'', size=size, pch=pch)
g15 <- plot_sparse_matrix(Ainv %*% Ltilde1inv, main='Perm: Ainv * Ltildeinv', size=size, pch=pch)
cowplot::plot_grid(g01,g03, g04,g05,g11,g13, g14,g15, ncol=2, byrow = FALSE)


# Now with permutations and effiency
chd2 <- Matrix::Cholesky(J%*%Qperm%*%J, super=TRUE, perm=FALSE)
# x to y
yperm2 <- as.numeric(J%*%chol(J%*%Qperm%*%J) %*% J%*%xperm)
yperm1-yperm2
# try to gmatch with Matrix fancy solve methods
xperm2 <- as.numeric(J%*% Matrix::solve(chd2, Matrix::solve(chd2, J%*%yperm2, system="Lt"), system="Pt"))
xperm2-xperm

Matrix::solve(chd2, J%*%yperm2, system="Lt")
L2 <- as(chd2, 'sparseMatrix')
(L2)%*%J%*%yperm2

Ltilde2 <- Matrix::t(J%*%chd2%*%J)
Ltilde1inv <- solve(Ltilde1)
yperm1 <- as.numeric(Ltilde1%*%xperm)
ytest1 <- yperm1[perm2]
# from y to x matches
xperm1 <- as.numeric(Ltilde1inv%*%yperm1)
as.numeric(Matrix::solve(chd1, Matrix::solve(chd, x, system="Lt"), system="Pt"))
max(abs(xperm-xperm1))
xtest1 <- as.numeric(Ainv %*%as.numeric(Ltilde1inv %*% yperm1))
max(abs(xtest1-xtest))
# LP matches
rsparse$fn2(ytest)
rdense$fn2(ytest)
obj$fn(as.numeric(Ainv%*%Ltilde1inv %*% yperm1))
# gr matches
head(rsparse$gr2(ytest))
gr1 <- function(y)
  as.numeric(obj$gr(as.numeric(Ainv %*%Ltilde1inv %*% y)) %*% (Ainv %*%Ltilde1inv))
head(gr1(yperm1))
head(rdense$gr2(ytest))

g11 <- plot_sparse_matrix(Qperm, main='Perm: Qperm=AQA\'', size=size, pch=pch)
g12 <- plot_sparse_matrix(J%*%Qperm%*%J, main='Perm: JQpermJ', size=size, pch=pch)
g13 <- plot_sparse_matrix(L1, main='Perm: Lperm=chol(JQpermJ)', size=size, pch=pch)
g14 <- plot_sparse_matrix(Ltilde1, main='Perm: Ltilde=(JLpermJ)\'', size=size, pch=pch)
g15 <- plot_sparse_matrix(Ainv %*% Ltilde1inv, main='Perm: Ainv * Ltildeinv', size=size, pch=pch)
cowplot::plot_grid(g01,g03, g04,g05,g11,g13, g14,g15, ncol=2, byrow = FALSE)




## check that Q and M give same initial trajectories for a simple
## model

## Demo with simple
message("starting chains to compare M vs Q and how long they match")
TMB::runExample('simple')
## TMB::runExample('ar1xar1')
thin <- 1
iter <- 500*thin
warmup <- 200
fitM <- sample_sparse_tmb(obj, iter=iter, warmup=warmup,
                          chains=1, cores=1, metric='dense', seed=124)
fitQ <- sample_sparse_tmb(obj, iter=iter, warmup=warmup,
                          chains=1, cores=1, metric='sparse', seed=124)
pM <- extract_samples(fitM, inc_warmup=TRUE, inc_lp=TRUE)
pQ <- extract_samples(fitQ, inc_warmup=TRUE, inc_lp=TRUE)
diff <- ((pM-pQ))
diff.long <- mutate(diff, iter=1:n()) %>%
  pivot_longer(-iter) %>% mutate(model='simple')
## filter(diff.long, abs(value)>.001)


pk <- readRDS('models/pollock/fit.RDS')
compile('models/pollock/goa_pk_tmb.cpp')
dyn.load('models/pollock/goa_pk_tmb.dll')
obj2 <- MakeADFun(data=pk$input$dat, parameters=pk$input$pars,
                  random=pk$input$random, DLL='goa_pk_tmb',
                  map=pk$input$map, silent=TRUE)
opt <- with(obj2, nlminb(par,fn,gr))
sdr <- sdreport(obj2, getJointPrecision=TRUE)
Q <- sdr$jointPrecision
Qinv <- as.matrix(solve(Q))
thin <- 1
iter <- 500*thin
warmup <- 200
td <- 3
fitM2 <- sample_sparse_tmb(obj2, iter=iter, warmup=warmup,
                           chains=1, cores=1, metric='dense',
                           seed=124, max_treedepth=td,
                           skip_optimization=TRUE,
                           Q=Q, Qinv=Qinv)
fitQ2 <- sample_sparse_tmb(obj2, iter=iter, warmup=warmup,
                           chains=1, cores=1, metric='sparse',
                           seed=124, max_treedepth=td,
                           skip_optimization=TRUE,
                           Q=Q, Qinv=Qinv)
pM2 <- extract_samples(fitM2, inc_warmup=TRUE, inc_lp=TRUE)
pQ2 <- extract_samples(fitQ2, inc_warmup=TRUE, inc_lp=TRUE)
diff2 <- ((pM2-pQ2))
diff.long2 <- mutate(diff2, iter=1:n()) %>%
  pivot_longer(-iter) %>% mutate(model='pollock')
## filter(diff.long2, abs(value)>.001)
xx <- rbind(diff.long, diff.long2)## %>% filter( (iter<150 & model=='simple') | (iter<20 & model=='pollock'))
g <- ggplot(xx, aes(iter, value, group=name)) + geom_line()+
  facet_wrap('model', nrow=2, scales='free') + geom_vline(xintercept=warmup, col=2) +
  ##  coord_cartesian(ylim=c(-1,1)) +
  ylab("Difference of pars in sparse vs dense chain")
ggsave('plots/checks_consistency_all_pars.png', g, width=8, height=4, units='in')
g <- ggplot(filter(xx, name=='lp__'), aes(iter, value, group=name)) + geom_line()+
  facet_wrap('model', nrow=2, scales='free_y') +
  geom_vline(xintercept=warmup, col=2) +
  coord_cartesian(ylim=c(-1,1)/20) +
  ylab("Difference of NLL of sparse vs dense chain")
ggsave('plots/checks_consistency_nll.png', g, width=8, height=4, units='in')


## run super long chains with thinning to ensure the same posteriors
TMB::runExample('simple')
thin <- 10
iter <- 10000*thin
warmup <- iter/10
fits <- fit_models(obj, iter=iter, warmup=warmup, chains=7,
                   cores=7, thin=thin)


get_diffs <- function(fits){
  xx <- lapply(fits, function(x) cbind(metric=x$metric, extract_samples(x, inc_lp=TRUE))) %>%
    bind_rows() %>%  pivot_longer(-metric) %>%
    group_by(metric, name)
  yy <- xx %>%
    summarize(mean=mean(value),
              median=median(value),
              upr=as.numeric(quantile(value, probs=.8)),
              lwr=as.numeric(quantile(value, probs=.2))) %>%
    ungroup() %>%
    pivot_longer(c(-metric, -name), names_to='quantity') %>%
    group_by(name, quantity) %>%
    mutate(absdiff=value-value[metric=='unit'],
           reldiff=(value-value[metric=='unit'])/value[metric=='unit']) %>%
    ungroup
  ## zz <- group_by(zz, name) %>% mutate(maxval=max(abs(value[quantity=='mean']))) %>% filter(maxval>.5)%>%ungroup

  pars <- filter(yy, quantity=='mean' & metric=='unit') %>%
    arrange(abs(value)) %>% pull(name)
  yy <- mutate(yy, name=factor(name, levels=pars))
  yy <- arrange(yy, name, metric, quantity)
  return(yy)
}

yy <- get_diffs(fits)
filter(yy, quantity=='mean') %>% tail(n=20)
g <- ggplot(yy, aes(factor(name), absdiff, group=metric, color=metric)) +
  geom_line()+ geom_point() + facet_wrap('quantity', nrow=4) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x='parameter (sorted by |mean|)', y='Absolute diff from unit metric')
ggsave('plots/checks_absdiffs_pars.png', g, width=9, height=6, unit='in')
g <- ggplot(yy, aes(factor(name), abs(reldiff), group=metric, color=metric)) +
  geom_line()+ geom_point() + facet_wrap('quantity', nrow=4) +
  coord_cartesian(ylim=c(0,1)/2)+
  theme(axis.text.x = element_text(angle = 90))+
  labs(x='parameter (sorted by |mean|)', y='Absolute Relative diff from unit metric')
ggsave('plots/checks_reldiffs_pars.png', g, width=9, height=6, unit='in')


filter(xx, name=='lp__') %>%
  ggplot(aes(metric, y=value)) + geom_violin()
