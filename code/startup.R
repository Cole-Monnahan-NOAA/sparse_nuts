
## library(adnuts)
library(StanEstimators)
library(rstan)
library(shinystan)
library(TMB)
library(dplyr)
library(ggplot2)
library(tidyr)
library(dsem)
library(cowplot)
library(mvtnorm)
library(microbenchmark)
library(future)
theme_set(theme_bw())


#' Fit a model with all 4 variants of NUTS. Returns a list of fits)
fit_models<- function(obj,  iter, warmup=NULL, chains=3,
                      cores=chains,
                      replicates=1:3, cpus=3,
                      adapt_metric=FALSE, init='random',
                      metrics=c('unit', 'diag', 'dense', 'sparse'),
                      control=NULL, model=NULL, ...){

  library(snowfall)
  if(is.null(model)) model <- obj$env$DLL
  fp <- gsub('.dll','', unlist(getLoadedDLLs()[[obj$env$DLL]])$path)
  sfInit(parallel=cpus>1, cpus=cpus)
  sfExportAll()
  sfLibrary(adnuts)
  sfLibrary(StanEstimators)
  sfLibrary(TMB)

  fits <- sfLapply(replicates, function(i){
    dyn.load(TMB::dynlib(fp))
    obj$retape()
    fits <- list()
    k <- 1
    for(metric in metrics){
      control2 <- control; wm <- warmup
      if(is.null(warmup)) wm <- ifelse(metric=='unit', floor(iter/2), 150)
      if(!adapt_metric & metric!='unit'){
        control2$metric <- 'unit_e'
        control2$adapt_init_buffer <- 25
        control2$adapt_term_buffer <- 75
        control2$adapt_window <- 50
      }
      fit <- sample_sparse_tmb(obj, iter=iter, warmup=wm,
                               chains=chains, cores=cores,
                               seed=i, metric=metric,
                               control=control2, init=init,
                               ...)
      fit$par_type <- ifelse(seq_along(fit$par_names) %in% obj$env$random, 'random', 'fixed')
      fit$replicate <- i; fit$model <- model
      fits[[k]] <- fit
      k <- k+1
    }
    return(invisible(fits))
  })

  fits <- do.call(rbind, fits)
  saveRDS(fits, file=paste0('results/',model, '_fits.RDS'))
  plot_output(fits)
  return(invisible(fits))
}

plot_output <- function(fits){
  model <- fits[[1]]$model
  g <- plot_stats(fits)
  ggsave(paste0('plots/',model, '_stats.png'), g, width=5, height=6)
  g <- plot_maxcors(fits)
  ggsave(paste0('plots/',model, '_maxcors.png'), g, width=5, height=3.5, units='in')
  g <- plot_mon(fits)
  ggsave(paste0('plots/',model, '_monitor.png'), g, width=5, height=3.5, units='in')
  g <- plot_timings(fits)
  ggsave(paste0('plots/',model, '_timing.png'), g, width=5, height=4.5, units='in')
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
      fit <- sample_sparse_tmb(obj, iter=iter, warmup=wm,
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
    post <- as.data.frame(fit)
    post.cor <- cor(post)
    diag(post.cor) <- 0 # zero out so can take max along rows
    max.cors <- sapply(1:ncol(post), function(i) post.cor[i,which.max(abs(post.cor[i,]))])
    ess <- fit$monitor$n_eff[1:ncol(post)]
    data.frame(model=fit$model, replicate=fit$replicate,
               metric=fit$metric, par=fit$par_names,
               ess=ess, max.cor=max.cors)
  })
  x <- do.call(rbind,x) %>% mutate(metric=metricf(metric))
  x
}

#' Helper function to extract key outputs from a fitted model
get_stats <- function(fits){
  x <- lapply(fits, function(fit) {
    ess <- min(fit$monitor$n_eff)
    rhat <- max(fit$monitor$rhat)
    time.total <- sum(fit$time.total + fit$time.Q + fit$time.Qinv)
    eff <- ess/time.total
    sp <- extract_sampler_params(fit)
    data.frame(model=fit$model, metric=metricf(fit$metric), replicate=fit$replicate, ess=ess,
               rhat=rhat, eff=eff, #time.total=time.total,
               gr=fit$time.gr,
               pct.divs=100*mean(sp$divergent__),
               time.warmup=sum(fit$time.warmup),
               time.sampling=sum(fit$time.sampling),
               gr2=fit$time.gr2)
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
  g <- cowplot::plot_grid(g1,g2, rel_heights=c(1.25,1), ncol=1)
  return(invisible(g))
}

metricf <- function(x) factor(x, levels=c('unit', 'diag', 'dense', 'sparse'))
plot_stats <- function(fits){
  stats <- get_stats(fits) %>% select(-rhat, -gr) %>%
    pivot_longer(c(-replicate, -model, -metric)) %>%
    mutate(name=factor(name,
        levels=c('gr2', 'time.warmup', 'time.sampling', 'pct.divs', 'ess', 'eff'),
        labels=c('Gradient eval (s)', 'NUTS warmup (s)',
                 'NUTS sampling (s)', '% divergent', 'ESS',
                 'Efficiency (ESS/time)')))
  g <- ggplot(stats, aes(x=metric, y=value)) +
    geom_jitter(width=.03, pch=1) +
    facet_wrap('name', scales='free', ncol=2) + ylim(0,NA)
  return(invisible(g))
}

plot_maxcors <- function(fits, logy=FALSE, abscor=FALSE, summarize=FALSE){
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
  g <- ggplot(mc, aes(max.cor, ess, color=metric)) + geom_point(alpha=.5) +
    stat_smooth(method='loess', formula=y~x) + xlim +labs(x=xlab, y="Effective sample size")
  if(logy) g <- g + scale_y_log10()
  if(!logy) g <- g+ coord_cartesian(ylim=c(0,NA))
  return(invisible(g))
}

get_mon <- function(fits){
  lapply(fits, function(fit){
  data.frame(model=fit$model, metric=fit$metric,
             replicate=fit$replicate, par=fit$monitor$variable,
             partype=factor(c(fit$par_type, 'lp__'),
                            levels=c('random', 'fixed', 'lp__')),
             ESS=fit$monitor$n_eff,
             Rhat=fit$monitor$Rhat)
  }) %>% bind_rows %>% mutate(metric=metricf(metric)) %>%
    arrange(partype)
}

plot_mon <- function(fits){
  mon <- get_mon(fits) %>% pivot_longer(cols=c('ESS', 'Rhat'))
 g <- ggplot(mon, aes(metric, y=value, color=partype,
                  group=interaction(par,replicate))) +
    geom_line(alpha=.25) + geom_jitter(width=.05)+
   facet_wrap('name', scales='free_y', ncol=1) +
   labs(color=NULL, y=NULL,x=NULL)
  return(invisible(g))
}
