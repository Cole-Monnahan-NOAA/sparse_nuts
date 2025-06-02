#
# # remotes::install_github("fishfollower/SAM/stockassessment")
#
# setwd(tempdir())
# filestoget <- c("cn.dat", "cw.dat", "dw.dat", "lf.dat", "lw.dat",
#                 "mo.dat", "nm.dat", "pf.dat", "pm.dat", "sw.dat",
#                 "survey.dat")
# url <- "https://raw.githubusercontent.com/fishfollower/SAM/master/stockassessment/tests/nsher/"
# d <- lapply(filestoget, function(f)download.file(paste(url,f,sep=""), f))
#
# library(stockassessment)
# cn <- read.ices("cn.dat")
# cw <- read.ices("cw.dat")
# dw <- read.ices("dw.dat")
# lf <- read.ices("lf.dat")
# lw <- read.ices("lw.dat")
# mo <- read.ices("mo.dat")
# nm <- read.ices("nm.dat")
# pf <- read.ices("pf.dat")
# pm <- read.ices("pm.dat")
# sw <- read.ices("sw.dat")
# surveys <- read.ices("survey.dat")
#
# dat <- setup.sam.data(surveys=surveys,
#                       residual.fleet=cn,
#                       prop.mature=mo,
#                       stock.mean.weight=sw,
#                       catch.mean.weight=cw,
#                       dis.mean.weight=dw,
#                       land.mean.weight=lw,
#                       prop.f=pf,
#                       prop.m=pm,
#                       natural.mortality=nm,
#                       land.frac=lf)
#
# conf <- defcon(dat)
#
# conf$fbarRange <- c(2,6)
# conf$corFlag <- 1
# conf$keyLogFpar <- matrix(c(
#   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
#   -1,    0,    1,    2,    3,    4,    5,    6,   -1,
#   -1,    7,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
#   8,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1), nrow=4, byrow=TRUE)
#
# par <- defpar(dat,conf)
#
# par$logFpar <- rep(0,9)
#
# obj <- sam.fit(dat,conf,par, run=FALSE)$obj
# saveRDS(obj, 'obj.sam.RDS')


library(stockassessment)
library(TMB)
obj <- readRDS('obj.sam.RDS')
obj$retape()
obj$env$beSilent()

opt <- with(obj, nlminb(par,fn,gr))#TMBhelper::fit_tmb(obj, getsd=FALSE, newtonsteps = 0)
sdrep <- sdreport(obj, getJointPrecision=TRUE)
Q <- sdrep$jointPrecision
M <- as.matrix(solve(Q))

library(adnuts)
fit <- sample_sparse_tmb(obj, Q=Q, Qinv=M, skip_optimization = TRUE,
                         iter=2300, warmup=300, seed=1, chains=5)

pairs_admb(fit, pars=1:5, order='slow')
pairs_admb(fit, pars=1:5, order='fast')
pairs_admb(fit, pars=1:5, order='mismatch')
plot_uncertainties(fit)
plot_Q(fit)

obj$par |> length()
obj$env$par |> length() - obj$par |> length()
opt <- with(obj, nlminb(par,fn,gr))
Q <- fit$mle$Q
M <- solve(Q)
max(abs(cov2cor(M)[lower.tri(M)]))
hist(abs(cov2cor(M)[lower.tri(M)]))

post <- as.data.frame(fit)
cors <- cor(post)
max(abs(cors[lower.tri(cors)]))


post.tmb <- as.data.frame(fit)
minsd <- apply(post.tmb, 2, sd) |> min()
maxsd <- apply(post.tmb, 2, sd) |> max()
maxsd/minsd
minsd <- min(sqrt(diag(M)))
maxsd <- max(sqrt(diag(M)))
maxsd/minsd

setwd(here::here())
