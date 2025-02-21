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

# fit the model with defaults
fit <- fit_pk(input, filename = 'pollockfit.RDS', save.sdrep = FALSE, getsd = FALSE)

library(adnuts)
mcmc <- sample_sparse_tmb(fit$obj, iter=1000, warmup=200, chains=4, seed=1, cores=4,
                          init='random')
pairs_admb(mcmc, order='mismatch', pars=1:5)
#launch_shinyadmb(mcmc)

setwd(here::here())
