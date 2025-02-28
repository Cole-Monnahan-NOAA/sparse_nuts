
## https://james-thorson-noaa.github.io/dsem/articles/vignette.html#multi-causal-ecosystem-synthesis
# Requires loading RTMB
library(dsem)
library(RTMB)

setwd(here('models/causal'))

# Load data
data(bering_sea)

#
Z = ts( bering_sea )
family = rep( "fixed", ncol(bering_sea) )

# Specify model
sem = "
    log_seaice -> log_CP, 0, seaice_to_CP
    log_CP -> log_Cfall, 0, CP_to_Cfall
    log_CP -> log_Esummer, 0, CP_to_E
    log_PercentEuph -> log_RperS, 0, Seuph_to_RperS
    log_PercentCop -> log_RperS, 0, Scop_to_RperS
    log_Esummer -> log_PercentEuph, 0, Esummer_to_Suph
    log_Cfall -> log_PercentCop, 0, Cfall_to_Scop
    SSB -> log_RperS, 0, SSB_to_RperS

    log_seaice -> log_seaice, 1, AR1, 0.001
    log_CP -> log_CP, 1,  AR2, 0.001
    log_Cfall -> log_Cfall, 1, AR4, 0.001
    log_Esummer -> log_Esummer, 1, AR5, 0.001
    SSB -> SSB, 1, AR6, 0.001
    log_RperS ->  log_RperS, 1, AR7, 0.001
    log_PercentEuph -> log_PercentEuph, 1, AR8, 0.001
    log_PercentCop -> log_PercentCop, 1, AR9, 0.001
  "

# Using fitRTMB
log_prior = function(p){
  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")
  sum(dnorm( p$beta_z[9:16], mean=0, sd=0.25, log=TRUE)) +
    sum(dnorm( p$beta_z[20], mean=0, sd=0.5, log=TRUE))

}
fitRTMB = dsemRTMB( sem = sem,
                    tsdata = Z,
                    family = family,
                    log_prior = log_prior,
                    control = dsem_control( use_REML = FALSE) )
inputs <- list(sem=sem, tsdata=Z,
               log_prior=log_prior, family=family)
saveRDS(inputs, 'inputs.RDS')

fit <- sample_sparse_tmb(obj, iter=2000, warmup=1000,
                          init='random', seed=1,
                          control=list(adapt_delta=.8))
pairs_admb(fit, order='mismatch', pars=1:5)
pairs_admb(fit, order='slow', pars=1:5)


obj$par |> length()
obj$env$par |> length()
opt <- with(obj, nlminb(par,fn,gr))
Q <- sdreport(obj, getJointPrecision=TRUE)$jointPrecision
M <- solve(Q)
max(abs(cov2cor(M)[lower.tri(M)]))

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
