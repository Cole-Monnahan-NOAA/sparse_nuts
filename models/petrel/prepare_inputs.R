## modified from a file sent by Theo Michelot

# devtools::install_github("TheoMichelot/hmmTMB")
library(hmmTMB)
library(ggplot2)
theme_set(theme_bw())
setwd(here::here('models/petrel'))

##################
## Prepare data ##
##################
# # Data from: https://doi.org/10.5441/001/1.q4gn4q56
# data <- read.csv(paste0("At-sea distribution Antarctic Petrel, Antarctica 2012",
#                         " (data from Descamps et al. 2016)-gps.csv")) |>
#     subset(select = c("individual.local.identifier",
#                       "location.long",
#                       "location.lat",
#                       "timestamp")) |>
#     setNames(c("ID", "x", "y", "time")) |>
#     moveHMM::prepData() |>
#     moveHMM::splitAtGaps(maxGap = 60, shortestTrack = 60*24)
# ind_zeros <- which(data$step == 0)
# data$step[ind_zeros] <- runif(length(ind_zeros),
#                               min = 0,
#                               max = min(data$step[-ind_zeros], na.rm = TRUE))
# write.csv(x = data, file = "movement_data.csv", row.names = FALSE)

# data <- read.csv("movement_data.csv")
# saveRDS(data, file='dat.RDS')

data <- readRDS('dat.RDS')
##################
## Define model ##
##################
# Hidden state process
hid <- MarkovChain$new(data = data,
                       n_states = 2,
                       formula = ~s(ID, bs = "re"),
                       initial_state = 2)

# Observation model
dists <- list(step = "gamma2")
par0 <- list(step = list(mean = c(2, 10), sd = c(2, 10)))
obs <- Observation$new(data = data,
                       dists = dists,
                       par = par0,
                       n_states = 2)

# HMM
hmm <- HMM$new(obs = obs, hid = hid)

################
## Set priors ##
################
priors <- list(coeff_fe_obs = matrix(c(log(2), 0.5,
                                       log(10), 0.5,
                                       log(2), 0.5,
                                       log(10), 0.5),
                                     ncol = 2, byrow = TRUE),
               coeff_fe_hid = matrix(c(-2, 0.5,
                                       -2, 0.5),
                                     ncol = 2, byrow = TRUE),
               log_lambda_hid = matrix(c(log(2), 0.5,
                                         log(250), 2),
                                       ncol = 2, byrow = TRUE))
hmm$set_priors(new_priors = priors)

library(adnuts)
hmm$setup()
obj <- hmm$tmb_obj()
saveRDS(obj, file='obj.petrel.RDS')

fit <- sample_sparse_tmb(obj, iter=2000, warmup=500, seed=1)
pairs_admb(fit, pars=1:5, order='slow')
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

#
# ###############
# ## Fit model ##
# ###############
# # ~30 min for 1 chain and 2000 iterations
# hmm$fit_stan(chains = 4, cores=4)
#
# # Histograms of posterior samples
# hmm$iters() |>
#     as.data.frame.table() |>
#     setNames(c("iter", "parameter", "value")) |>
#     ggplot(aes(value)) +
#     geom_histogram(bins = 50) +
#     facet_wrap("parameter", scales = "free", nrow = 2)














