

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
