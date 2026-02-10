
## library(SparseNUTS)
library(StanEstimators)
#library(rstan)
library(shinystan)
##library(TMB)
library(TMBhelper)
library(dplyr)
library(ggplot2)
library(tidyr)
library(dsem)
library(cowplot)
library(mvtnorm)
library(microbenchmark)
#library(bench)
library(future)
library(fmesher)
library(SparseNUTS)
library(here)
library(Matrix)
theme_set(theme_bw())

#' Convert metric to algorithm
algf <- function(metric){
  factor(ifelse(metric=='unit' | metric=='Stan default', 'NUTS','SNUTS'))
}

#' Make a factor of model names sorted alphabetically (for plotting)
modelf <- function(x){
  labs <- sort(unique(as.character(x)))
  factor(x, levels=labs)
}

metricf <- function(x){
  lvl <- c('unit', 'diag', 'dense', 'sparse', 'stan', 'tmbstan')
  labs <- lvl
  labs[5] <- 'Stan default'
  factor(x, levels=lvl, labels=labs)
}

#' Wrapper function to run warmup tests.
fit_warmups <- function(obj, cpus, reps, ...){
  fit_models(obj=obj, num_samples=250,
             num_warmup=1000, chains=1,
             adapt_metric = TRUE,
             outpath='results/warmups',
             metrics='auto', plot=FALSE,
             cpus=cpus, replicates=reps,
             ...)
}



#' Fit a model with different metrics. Returns a list of fits)
fit_models<- function(obj,  num_samples=NULL, num_warmup=NULL, chains=4,
                      cores=chains,
                      replicates=1:3, cpus=1, thin=1,
                      adapt_metric=NULL,
                      init='random',
                      metrics=c('unit', 'auto'),
                      outpath='results',
                      do.tmbstan=FALSE,
                      control=NULL, model=NULL, plot=TRUE, refresh=500,
                      ...){

  library(snowfall)
  isRTMB <- if(obj$env$DLL=='RTMB') TRUE else FALSE
  if(is.null(model)) model <- obj$env$DLL #otherwise all RTMB models will have same name
  fp <- gsub('.dll','', unlist(getLoadedDLLs()[[obj$env$DLL]])$path)
  sfInit(parallel=cpus>1, cpus=cpus)
  sfExportAll()
  sfLibrary(SparseNUTS)
  if(isRTMB){
    sfLibrary(RTMB)
  } else {
    sfLibrary(TMB)
  }

  fits <- sfLapply(replicates, function(i){
    dyn.load(TMB::dynlib(fp))
    obj$retape()
    fits <- list()
    k <- 1
    for(metric in metrics){
      message("--- Starting model ", model, " replicate ", i, " for metric ", metric, ' ---')
      fit <- sample_snuts(obj, num_samples=num_samples,
                          num_warmup=num_warmup,
                          chains=chains, cores=cores,
                          seed=i, metric=metric, thin=thin,
                          control=control, init=init,
                          refresh=refresh, model_name = model,
                          adapt_stan_metric = adapt_metric,
                          ...)
      fit$par_type <- ifelse(seq_along(fit$par_names) %in% obj$env$random, 'random', 'fixed')
      fit$replicate <- i; fit$model <- model
      fits[[k]] <- fit
      k <- k+1
      # png(paste0('plots/', model, '_', metric, '_pairs_slow.png'), width=7, height=5, units='in', res = 300)
      # pairs_admb(fit, order='slow', pars=1:5)
      # dev.off()
      # png(paste0('plots/', model, '_', metric, '_pairs_mismatch.png'), width=7, height=5, units='in', res = 300)
      # pairs_admb(fit, order='mismatch', pars=1:5)
      # dev.off()
    }
    if(do.tmbstan){
      message("Starting tmbstan run")
      set.seed(i)
      inits <- sample_inits(fit, chains=chains)
      tmp <- tmbstan::tmbstan(obj, cores=cores, chains=chains, init=inits,
                              seed=i, refresh=500)
      fit.tmbstan <- list(samples=rstan::extract(tmp, permuted=FALSE, inc_warmup=TRUE),
                          sampler_params=rstan::get_sampler_params(tmp, inc_warmup=TRUE),
                          # confirmed that this appears to drop warmup samples automatically
                          monitor=posterior::summarise_draws(tmp),
                          model=model, metric='tmbstan', warmup=tmp@stan_args[[1]]$warmup,
                          iter=tmp@stan_args[[1]]$iter-tmp@stan_args[[1]]$warmup,
                          thin=tmp@stan_args[[1]]$thin, replicate=i,
                          time.warmup=as.numeric(rstan::get_elapsed_time(tmp)[,1]),
                          time.sampling=as.numeric(rstan::get_elapsed_time(tmp)[,2]),
                          time.total=sum(rstan::get_elapsed_time(tmp)),
                          time.Q=0, time.Qinv=0, time.opt=0)
      fits[[k]] <- fit.tmbstan
      k <- k +1
    }
    return(invisible(fits))
  })

  fits <- do.call(rbind, fits)
  if(!is.null(outpath))  saveRDS(fits, file=paste0(outpath,'/',model, '_fits.RDS'))
  if(plot) plot_output(fits)
  return(invisible(fits))
}

plot_output <- function(fits, do.pairs=FALSE){
  model <- fits[[1]]$model
  pdf(paste0('plots/case_studies/',model,'.pdf'), width=6, height=9)
  plot_Q(fits[[1]])
  g <- plot_stats(fits)
  print(g)
  #ggsave(paste0('plots/',model, '_stats.png'), g, width=5, height=6)
  g4 <- plot_maxcors(fits)
  #ggsave(paste0('plots/',model, '_maxcors.png'), g, width=5, height=3.5, units='in')
  g5 <- plot_vars(fits)
  xx <- plot_mon(fits)
  #ggsave(paste0('plots/',model, '_monitor.png'), g, width=5, height=3.5, units='in')
  yy <- plot_timings(fits)
  #ggsave(paste0('plots/',model, '_timing.png'), g, width=5, height=4.5, units='in')
  g <- cowplot::plot_grid(xx[[1]], xx[[2]], xx[[3]], g4,g5, yy[[1]], ncol=2, byrow = FALSE)
  #ggsave(paste0('plots/', model, '_results.png'), g, width=7, height=7.5, units='in')
  print(g)
  if(do.pairs){
    for(ii in 1:length(fits)){
      fittmp <- fits[[ii]]
      if(fittmp$replicate==1){
        for(jj in c('slow', 'mismatch')){
          npar <- min(5,length(fittmp$par_names))
          if(jj=='slow') npar=npar+1
          pairs_admb(fittmp, order=jj, pars=1:npar)
          mtext(text=paste0(fittmp$model, " | ", fittmp$metric, ' | ', jj), 3,
                line=-2, outer=TRUE, cex=1.2)
        }
      }
    }
  }
  dev.off()
}

#' Fit a model with all 4 variants of NUTS. Returns a list of fits)
fit_models_old <- function(obj,  model, iter, warmup=NULL, chains=4, cores=4,
                      replicates=1, adapt_metric=TRUE, init='random',
                      metrics=c('unit', 'diag', 'dense', 'sparse'),
                      control=NULL, ...){
  fits <- list()
  k <- 1
  for(i in 1:replicates){
    for(metric in metrics){
      control2 <- control; wm <- warmup
      if(is.null(warmup)) wm <- ifelse(metric=='unit', floor(iter/2), 150)
      if(!adapt_metric & metric!='unit'){
        control2$metric <- 'unit_e'
        control2$adapt_init_buffer <- 25
        control2$adapt_term_buffer <- 75
        control2$adapt_window <- 50
      }
      if(is.null(wm)) browser()
      fit <- sample_snuts(obj, iter=iter, warmup=wm,
                               chains=chains, cores=cores,
                               seed=i, metric=metric,
                               control=control2, init=init,
                               ...)
      fit$replicate <- i; fit$model <- model
      fits[[k]] <- fit
      k <- k+1
    }
  }
  return(invisible(fits))
}


get_maxcors <- function(fits){
  x <- lapply(fits, function(fit){
    if(fit$metric %in% c('stan', 'tmbstan')) return(NULL)
    post <- as.data.frame(fit)
    post.cor <- cor(post)
    diag(post.cor) <- 0 # zero out so can take max along rows
    max.cors <- sapply(1:ncol(post), function(i) post.cor[i,which.max(abs(post.cor[i,]))])
    ess <- fit$monitor$ess_bulk[1:ncol(post)]
    eff <- ess/sum(fit$time.total + fit$time.Q + fit$time.Qinv)
    data.frame(model=fit$model, replicate=fit$replicate,
               metric=fit$metric, par=fit$par_names,
               ess=ess, max.cor=max.cors, eff=eff)
  })
  x <- do.call(rbind,x) %>% mutate(metric=metricf(metric))
  x
}


get_cors <- function(fits){
  x <- lapply(fits, function(fit){
    if(fit$metric %in% c('stan', 'tmbstan')) return(NULL)
    post <- as.data.frame(fit)
    post.cor <- cor(post)
    post.cor <- as.numeric(post.cor[lower.tri(post.cor)])
    mle.cor <- fit$mle$cor
    mle.cor <- as.numeric(mle.cor[lower.tri(mle.cor)])
    data.frame(model=fit$model, replicate=fit$replicate,
               post.cor=post.cor, mle.cor=mle.cor,
               metric=fit$metric)
  })
  x <- do.call(rbind,x) %>% mutate(metric=metricf(metric))
  x
}


get_vars <- function(fits){
  x <- lapply(fits, function(fit){
    if(fit$metric %in% c('stan', 'tmbstan')) return(NULL)
    post <- as.data.frame(fit)
    post.sd <- apply(post, 2, sd)
    mle.sd <- fit$mle$se
    ess <- fit$monitor$ess_bulk[1:ncol(post)]
    eff <- ess/sum(fit$time.total + fit$time.Q + fit$time.Qinv)
    data.frame(model=fit$model, replicate=fit$replicate,
               metric=fit$metric, par=fit$par_names,
               ess=ess, marginal.sd=mle.sd, post.sd=post.sd, eff=eff)
  })
  x <- do.call(rbind,x) %>% mutate(metric=metricf(metric))
  x
}

#' Helper function to extract key outputs from a fitted model
get_stats <- function(fits){
  x <- lapply(fits, function(fit) {
    ess <- min(fit$monitor$ess_bulk)
    rhat <- max(fit$monitor$rhat)
    time.total <- sum(fit$time.total + fit$time.Q + fit$time.Qinv)
    eff <- ess/time.total
    sp <- extract_sampler_params(fit, inc_warmup = FALSE)
    data.frame(model=fit$model, metric=fit$metric, replicate=fit$replicate, ess=ess,
               rhat=rhat, eff=eff, #time.total=time.total,
               time.Qall=fit$time.Q +fit$time.Qinv + fit$time.opt,
               time.total=time.total,
               pct.divs=100*mean(sp$divergent__),
               n.leapfrog=mean(sp$n_leapfrog__),
               time.warmup=sum(fit$time.warmup),
               time.sampling=sum(fit$time.sampling))
  })
  do.call(rbind,x)
}

get_timings <- function(fits){
  x <- lapply(fits, function(fit) {
    with(fit, data.frame(model=model, metric=metricf(metric),
               replicate=replicate, chain=1:dim(fit$samples)[2],
               Q=time.Q,
               ##gr2=time.gr2,
               Qinv=time.Qinv,
               opt=time.opt,
              ## gr=time.gr,
               warmup=time.warmup,
               sampling=time.sampling))
 })
 x <- do.call(rbind,x)
 x
}

plot_timings <- function(fits){
  x <- get_timings(fits) %>%
    pivot_longer(-c(model, metric, replicate, chain),
                 names_to='operation', values_to='time') %>%
    group_by(model, metric, replicate,chain) %>%
    mutate(pct.time=100*time/sum(time),
           operation=factor(operation,
                            levels=c('opt', 'Q', 'Qinv', 'warmup', 'sampling'),
                            labels=c('optimization', 'calculate Q', 'invert Q', 'NUTS warmup', 'NUTS sampling')))%>%
    filter(time>0)
  g1 <- ggplot(x, aes(metric, pct.time, color=operation)) + geom_jitter(width=.05)+
    ylim(0,100) + theme(legend.position='top') +
    labs(x=NULL, y='% time', color=NULL)+
    guides(color=guide_legend(ncol=3))
  g2 <- ggplot(x, aes(metric, time, color=operation)) +
    geom_jitter( width=.05) + scale_y_log10() +
    theme(legend.position='none') + labs(y='Time (s)')
  # g <- cowplot::plot_grid(g1,g2, rel_heights=c(1.25,1), ncol=1)
  # return(invisible(g))
  return(list(g1,g2))
}

plot_stats <- function(fits){
  stats <- get_stats(fits) %>% select(-time.Qall, -time.sampling) %>%
    pivot_longer(c(-replicate, -model, -metric)) %>%
    mutate(name=factor(name,
        levels=c('rhat', 'time.warmup', 'time.total', 'pct.divs', 'ess', 'eff'),
        labels=c('rhat', 'NUTS warmup (s)',
                 'Total time (s)', '% divergent', 'ESS',
                 'Efficiency (ESS/time)')))
  g <- ggplot(stats, aes(x=metric, y=value)) +
    geom_jitter(width=.03, pch=1) +
    facet_wrap('name', scales='free', ncol=2) + ylim(0,NA) +
    labs(x=NULL, y=NULL)
  return(invisible(g))
}

plot_maxcors <- function(fits, logy=TRUE, abscor=FALSE, summarize=FALSE){
  mc <- get_maxcors(fits)
  if(summarize)
    mc <- group_by(model, metric, par) %>%
      summarize(max.cor=mean(max.cor), ess=mean(ess))
  if(abscor){
    mc$max.cor <- abs(mc$max.cor)
    xlim <- xlim(0,1)
    xlab <- 'Absolute value of maximum correlation'
  } else {
    xlim <- xlim(-1,1)
    xlab <- 'Maximum correlation'
  }
  g <- ggplot(mc, aes(max.cor, eff, color=metric)) + geom_point(alpha=.5) +
    theme(legend.position='top')+
    #stat_smooth(method='loess', formula=y~x) +
    xlim +labs(x=xlab, y="Efficiency (ESS/time)", color=NULL)
  if(logy) g <- g + scale_y_log10()
  if(!logy) g <- g+ coord_cartesian(ylim=c(0,NA))
  return(invisible(g))
}

plot_vars <- function(fits, logy=TRUE, summarize=FALSE){
  vars <- get_vars(fits)
  if(summarize)
    vars <- vars %>% group_by(model, metric, par) %>%
      summarize(marginal.sd=mean(marginal.sd), ess=mean(ess))
  g <- ggplot(vars, aes(marginal.sd, eff, color=metric)) +
    geom_point(alpha=.5) + scale_x_log10()+
   # stat_smooth(method='loess', formula=y~x) +
    labs(x='Marginal MLE SD', y="Efficiency (ESS/time)", color=NULL) +
    theme(legend.position='top')
  if(logy) g <- g + scale_y_log10()
  if(!logy) g <- g+ coord_cartesian(ylim=c(0,NA))
  return(invisible(g))
}


get_mon <- function(fits){
  lapply(fits, function(fit){
    divs <- sum(extract_sampler_params(fit)$divergent___)
  data.frame(model=fit$model, metric=fit$metric,
             replicate=fit$replicate, par=fit$monitor$variable,
             partype=factor(c(fit$par_type, 'lp__'),
                            levels=c('random', 'fixed', 'lp__')),
             ESS=fit$monitor$ess_bulk,
             rhat=fit$monitor$rhat, divergences=divs)
  }) %>% bind_rows %>% mutate(metric=metricf(metric)) %>%
    arrange(partype)
}

plot_mon <- function(fits){
  mon <- get_mon(fits) %>%
    pivot_longer(cols=c('ESS', 'rhat', 'divergences'))
 g1 <- ggplot(filter(mon, name=='ESS'),
              aes(metric, y=value, color=partype,
                  group=interaction(par,replicate))) +
    geom_line(alpha=.25) + geom_jitter(width=.05)+
   ylim(0,NA) + theme(legend.position = 'top')+
   labs(x=NULL, y='Effective sample size', color=NULL)

 g2 <- ggplot(filter(mon, name=='rhat'),
              aes(metric, y=value, color=partype,
                  group=interaction(par,replicate))) +
   geom_line(alpha=.25) + geom_jitter(width=.05)+
   labs(x=NULL, y='rhat', color=NULL) + ylim(1,NA) +
   theme(legend.position = 'top')
 g3 <- ggplot(filter(mon, name=='divergences'),
              aes(metric, y=value)) +
  geom_point()+
   labs(x=NULL, y='No. of divergences') + ylim(0,NA)

 # facet_wrap('name', scales='free_y', ncol=1) +
   #labs(color=NULL, y=NULL,x=NULL)
  return(invisible(list(g1,g2,g3)))
}


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
sim_spde_dat <- function(n, sparse=TRUE, map=NULL, seed=NULL){
  require(RTMB)
  if(is.null(seed)) set.seed(n)
  t0 <- Sys.time()
  n_x = n_y = n
  # make banded sparse to represent SAR (simultaneous autoregression) model
  A_xx = Matrix::bandSparse(n=n, k=c(1,-1)) * 1
  Q_xx = -0.4 * A_xx
  diag(Q_xx) = 1 + 0.4^2
  Q_zz = kronecker( Q_xx, Q_xx )
  z = RTMB:::rgmrf0(n=1, Q=Q_zz )
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
      # broad prior to stabilize when estimated
      nll <- nll - dnorm(log_tau, mean=-1.75, sd=2, log=TRUE)
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
  # allows tmbstan to work
  environment(f) <- new.env()
  environment(f)$data <- data
  obj <- RTMB::MakeADFun(f, parameters, random="x", silent=TRUE, map=map)
  return(obj)
}




get_wasserstein <- function(reps, obj, post, model, savedraws=FALSE,
                            init=c('mode', 'post', 'unif-2', 'unif-1', '0')){
  # # use recursion to loop over multiple values (this fails!)
  # if(length(rep)>1)
  #   bind_rows(lapply(rep, \(i) get_wasserstein(rep=i, obj, post, model, init)) )
  isRTMB <- obj$env$DLL=='RTMB'
  # the only tricky part here is getting the inits right
  # especially when there's a map used. in that case I think only
  # default makes sense
  init <- match.arg(init)
  draws.post <- post |> select(-c(model, metric, replicate))
  out <- list()
  for(rep in reps){
    message("Starting w1d for model=", model, "; init=",init,"; rep=", rep, " at time=", Sys.time())
    set.seed(rep)
    if(init=='post')  inits <- as.numeric(draws.post[sample(1:nrow(draws.post), size=1),])
    if(init=='mode')  inits <- obj$env$last.par
    if(init=='unif-2') inits <- runif(ncol(draws.post), min=-2, max=2)
    if(init=='unif-1') inits <- runif(ncol(draws.post), min=-1, max=1)
    if(init=='0') inits <- 0*runif(ncol(draws.post), min=-1, max=1)
    # do joint first so can use parList() to build objm
    if(isRTMB){
      objj <- RTMB::MakeADFun(obj$env$data, parameters=obj$env$parList(),
                              map=obj$env$map, silent=TRUE)
      objm <- RTMB::MakeADFun(obj$env$data, parameters=objj$env$parList(inits),
                              map=obj$env$map, random=obj$env$random,
                              silent=TRUE)
    } else {
      objj <- TMB::MakeADFun(data=obj$env$data, parameters=obj$env$parList(),
                             map=obj$env$map, silent=TRUE, DLL=obj$env$DLL)
      objm <- TMB::MakeADFun(data=obj$env$data, parameters=objj$env$parList(inits),
                             map=obj$env$map, random=obj$env$random, DLL=obj$env$DLL,
                             silent=TRUE)
    }
    stopifnot(length(objj$env$last.par)==ncol(draws.post))

    # Get Q draws
    time0 <- Sys.time()
    opt <- with(objm, nlminb(par,fn,gr))
    opt <- with(objm, nlminb(opt$par, fn, gr))
    time.opt <- Sys.time()-time0
    mle <- objm$env$last.par.best
    sdrep <- sdreport(objm, getJointPrecision=TRUE)
    time.sdrep <- Sys.time()-time0- time.opt
    Q <- sdrep$jointPrecision
    if(is.null(Q)){
      M <- sdrep$cov.fixed
    } else {
      M <- solve(Q) |> as.matrix()
    }
    draws.Q <- mvtnorm::rmvnorm(n=1000, mean=mle, sigma=M) |>
      as.data.frame()
    time.Qdraws <- Sys.time()-time0-time.sdrep
    time.total <- Sys.time()-time0

    # run pathfinder with default of 1000 draws
    fn <- function(x) -objj$fn(x)
    grad_fun <- function(x) -objj$gr(x)
    time0 <- Sys.time()
    pf <- stan_pathfinder(fn=fn, grad_fun=grad_fun,
                          par_inits=inits, quiet=TRUE, refresh=0)
    time.pf <- Sys.time()-time0
    draws.pf <- pf@draws |> as.data.frame() |> select(-(1:2)) |>
      select(-c(.chain, .iteration, .draw, path__))
    if(savedraws & rep==1){
      saveRDS(draws.pf, file=paste0('results/pathfinder/', model, '_draws.RDS'))
    }

    ## calculate wasserstein distance for both
    library(transport)
    a.pf = wpp(draws.pf, mass = rep(1 / nrow(draws.pf), nrow(draws.pf)))
    a.Q = wpp(draws.Q, mass = rep(1 / nrow(draws.Q), nrow(draws.Q)))
    b = wpp(draws.post, mass = rep(1 / nrow(draws.post), nrow(draws.post)))
    w.pf <- wasserstein(a=a.pf, b=b, p = 1)
    w.Q <- wasserstein(a=a.Q, b=b, p = 1)

    df.Q <- data.frame(model=model, type='Q', w1d=w.Q, rep=rep,
                       time.total=as.numeric(time.total, 'secs'),
                       time.opt=as.numeric(time.opt, 'secs'),
                       time.sdrep=as.numeric(time.sdrep, 'secs'),
                       time.Qdraws=as.numeric(time.Qdraws, 'secs'))
    df.pf <- data.frame(model=model, type='Pathfinder', w1d=w.pf,
                        rep=rep,
                        time.total=as.numeric(time.pf, 'secs'))
    out <- bind_rows(out, df.Q, df.pf)
  }
  saveRDS(out, paste0('results/pathfinder/', model, '_pathfinder.RDS'))
  message("Done with wassersteine distance for model=", model)
  return(out)
}



# quick plot function for the gradient benchmarks
plot.bench <- function(bench) {
  bench_spde_gr <- bench |> group_by(metric,type, nrepars) |>
    summarize(time=median(time), size=median(size), .groups='drop') |>
    mutate(reltime=time/time[metric=='unit'],
           relsize=size/size[metric=='unit']) |>
    mutate(type=factor(type,
                       levels=c('simple', 'original'),
                       labels = c('Preconditioning', 'Preconditioning + gradient'))) |>
    filter(metric %in% c('dense', 'diag','sparse')) |>
    mutate(metric=metricf(metric)) |>
    pivot_longer(cols=c(reltime,relsize)) |>
    # some weird data massaging to get time and size into same long factor
    filter(! (name=='relsize' & type =='Simple')) |>
    mutate(type2=ifelse(name=='reltime', as.character(type), 'Memory size')) |>
    mutate(type2=factor(type2,
                        levels=c('Preconditioning',
                                 'Preconditioning + gradient',
                                 'Memory size')))
  g <- ggplot(bench_spde_gr, aes(x=nrepars, y=value,
                                 color=metric)) +
    geom_line(linewidth=1, alpha=.8) + geom_point() +
    labs(y='Value relative to original model (q)',
         x='Number of parameters',
         color=NULL) +
    facet_wrap('type2', nrow=1, scales='free_y') + scale_y_log10() +
    scale_x_log10() +   geom_hline(yintercept=1)+
    #guides(col=guide_legend(nrow=2)) +
    theme(legend.position='inside',
          legend.position.inside = c(.075,.77),
          legend.background = element_rect(fill = "transparent"))
  g
}
