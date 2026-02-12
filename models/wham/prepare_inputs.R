

library(wham)
library(TMB)
library(ggplot2)
theme_set(theme_bw())
setwd(here::here('models/wham'))

# # example from the WHAM examples on NAA, pkcing the best AIC model (m13)
# # #by default do not perform bias-correction
# if(!exists("basic_info")) basic_info <- NULL
# library(ggplot2)
# library(tidyr)
# library(dplyr)
# #
# #
# wham.dir <- find.package("wham")
# path_to_examples <- system.file("extdata", package="wham")
# asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))
# env.dat <- read.csv(file.path(path_to_examples,"GSI.csv"), header=T)
# df.mods <- data.frame(NAA_cor = c('---','iid','ar1_y','iid','ar1_a','ar1_y','2dar1','iid','ar1_y','iid','ar1_a','ar1_y','2dar1'),
#                       NAA_sigma = c('---',rep("rec",2),rep("rec+1",4),rep("rec",2),rep("rec+1",4)),
#                       R_how = paste0(c(rep("none",7),rep("limiting-lag-1-linear",6))), stringsAsFactors=FALSE)
# n.mods <- 1
# mods <- vector("list",n.mods)
# m <- 13
# NAA_list <- list(cor=df.mods[m,"NAA_cor"], sigma=df.mods[m,"NAA_sigma"],
#                  decouple_recruitment = FALSE)
# if(NAA_list$sigma == '---') NAA_list = NULL
# ecov <- list(
#   label = "GSI",
#   mean = as.matrix(env.dat$GSI),
#   logsigma = 'est_1', #estimate obs sigma, 1 value shared across years
#   year = env.dat$year,
#   use_obs = matrix(1, ncol=1, nrow=dim(env.dat)[1]),# use all obs (=1)
#   process_model = 'ar1', #"rw" or "ar1"
#   recruitment_how = matrix(df.mods$R_how[m]))# n_Ecov x n_stocks
# input <- suppressWarnings(prepare_wham_input(asap3, recruit_model = 3, #Bev Holt recruitment
#                                              model_name = "Ex 6: Numbers-at-age",
#                                              selectivity=list(model=rep("age-specific",3), re=c("none","none","none"),
#                                                               initial_pars=list(c(0.1,0.5,0.5,1,1,1),
#                                                                                 c(0.5,0.5,0.5,1,0.5,0.5),
#                                                                                 c(0.5,0.5,1,1,1,1)),
#                                                               fix_pars=list(4:6,4,3:6)),
#                                              NAA_re = NAA_list,
#                                              ecov=ecov,
#                                              basic_info = basic_info,
#                                              age_comp = "logistic-normal-miss0")) #logistic
# obj0 <- fit_wham(input, do.retro=FALSE, do.osa=FALSE, do.fit = FALSE)
# saveRDS(obj0, file='obj0.wham.RDS')
#
# library(wham)
# library(TMB)
# obj0 <- readRDS('obj0.wham.RDS')
# obj0$retape()
# obj0$env$beSilent()
# opt0 <- TMBhelper::fit_tmb(obj0, getsd=FALSE)
# sdrep <- sdreport(obj0, getJointPrecision=TRUE)
# Q <- sdrep$jointPrecision
# M <- as.matrix(solve(Q))
# #pilot chains to see if anything totally bad
# library(adnuts)
# fit0 <- sample_sparse_tmb(obj0, Q=Q, Qinv=M, skip_optimization = TRUE,
#                          iter=700, warmup=300, seed=1, chains=6,
#                          control=list(max_treedepth=5))
# pairs_admb(fit0, order='slow', pars=1:6)
# #
# input <- obj0$input
# input$map$logit_selpars <- factor(NA*input$map$logit_selpars)
# input$map$Ecov_process_pars <- factor(c(NA,1,2))
# input$par$Ecov_process_pars[1] <- obj0$env$parList(opt0$par)$Ecov_process_pars[1]
# input$map$mean_rec_pars <- factor(c(1,NA))
# input$par$mean_rec_pars[2] <- obj0$env$parList(opt0$par)$mean_rec_pars[2]
# input$map$Ecov_obs_logsigma <- factor(NA*input$map$Ecov_obs_logsigma)
# input$par$Ecov_obs_logsigma <- obj0$env$parList(opt0$par)$Ecov_obs_logsigma
# input$par$logit_selpars <- obj0$env$parList(opt0$par)$logit_selpars
# obj <- MakeADFun(data=input$data, parameters = input$par, map=input$map,
#                  random=input$random, silent=TRUE, DLL=obj0$env$DLL)
# obj$gr()
# saveRDS(obj, file='obj.wham.RDS')

obj <- readRDS('obj.wham.RDS')
obj$retape()

opt <- TMBhelper::fit_tmb(obj, getsd=FALSE, newtonsteps = 0)
sdrep <- sdreport(obj, getJointPrecision=TRUE)
Q <- sdrep$jointPrecision
M <- as.matrix(solve(Q))
library(SparseNUTS)
fit <- sample_snuts(obj, Q=Q, Qinv=M, skip_optimization = TRUE,
                         iter=2300, warmup=300, seed=1, chains=5)
save.image()
pairs(fit, pars=1:5, order='slow')
pairs(fit, pars=1:5, order='fast')
pairs(fit, pars=1:5, order='mismatch')
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
