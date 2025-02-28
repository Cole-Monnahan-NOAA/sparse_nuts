library(GOApollock)
## the 2024 assessment with a small tweak to the prior for selex,
## transf_rho, and map off DM_pars for srv6 since it wandered to
## Inf for some reason?

setwd(here::here('models/pollock'))
input <- prepare_pk_input(path=getwd(), datfile='pk24_12.txt',
                          version='23d: 2024 final',
                          modfile = 'pollock',
                          complike = 'D-M')
# age 1 and 2 Shelikof indices turned off so don't estimate catchabilities
input$map$log_q4 <- input$map$log_q5 <- factor(NA)
input$map$log_DM_pars <- factor(c(1:4,NA))
input$par$log_DM_pars[5] <- 1

# bump up so log_DM pars doesn't hit bounds. Needs to be moved into dat file later
input$dat$multN_srv1 <- input$dat$multN_srv1*2
input$dat$multN_srv3 <- input$dat$multN_srv3*2
input$dat$multN_srv6 <- input$dat$multN_srv6*2

# fit the model with defaults, since not doing any Newton steps
# during optimization Q and M are not positive definite. So
# cheating here by providing obj$par with the MLE. This breaks
# runtime for opt for this model.
opt <- fit_pk(input, filename = 'pollockfit.RDS', do.fit=TRUE, save.sdrep = FALSE, getsd = FALSE)
obj <- opt$obj
obj$par <- opt$opt$par
saveRDS(obj, 'obj.pollock.RDS')

## run longer chains
library(adnuts)
library(StanEstimators)
fit <- sample_sparse_tmb(obj, iter=4000, warmup=200, chains=4,
                         cores=4)
## compare correlations and marginal sds
post <- as.data.frame(fit)
cors <- cor(post)
max(abs(cors[lower.tri(cors)]))

obj$par |> length()
obj$env$par |> length() -obj$par |> length()
opt <- with(obj, nlminb(par,fn,gr))
Q <- sdreport(obj, getJointPrecision=TRUE)$jointPrecision
M <- solve(Q)
max(abs(cov2cor(M)[lower.tri(M)]))

minsd <- apply(post, 2, sd) |> min()
maxsd <- apply(post, 2, sd) |> max()
maxsd/minsd
minsd <- min(sqrt(diag(M)))
maxsd <- max(sqrt(diag(M)))
maxsd/minsd



pairs_admb(fit, pars=1:4, order='slow')
pairs_admb(fit, pars=1:4, order='mismatch')
plot_uncertainties(fit)
plot_sampler_params(fit)



setwd(here::here())
