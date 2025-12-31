

setwd(here::here())
source("code/startup.R")
skip.models <- c('ratios', 'ratio', 'RTMB', 'simple', 'wildf2', 'wildf3', 'wildf4', 'cors', 'cor','VAR','glmmTMB')
skip.models <- c(skip.models, paste0(skip.models, '_ELA'))

width1 <- 180/25.4 # 180mm full width in inches
width2 <- 80/25.4 # half width

message("Processing simulated results...")
stats_sim <- bind_rows(readRDS('results/spde_stats.RDS'),
                       #readRDS('results/VAR_stats.RDS'),
                       readRDS("results/glmmTMB_stats.RDS")) |>
 # group_by(model, metric, nrepars) |> #summarize(n=n()) |> print(n=Inf)
  filter(! (model=='VAR' & nrepars <205)) |>
  ungroup()

perf_sim <- stats_sim |>
  mutate(alg=algf(metric),
         time=time.warmup+time.sampling) |>
  group_by(model, alg, nrepars)  |>
  # mean across replicates for each algorithm
  mutate(eff_avg=mean(eff),
         ess_avg=mean(ess),
         time_avg=mean(time)) |> ungroup() |>
  pivot_longer(cols=c('eff', 'ess', 'time')) |>
  # relative to unit metric
  group_by(model, nrepars) |>
  mutate(relval=value/value[alg=='NUTS']) |> ungroup() |>
  mutate(name2=factor(name, levels=c('time', 'ess', 'eff'),
                      labels=c('Total time (s)', 'minimum ESS', 'Efficiency (ESS/s)')))
# print some stats for the MS
filter(perf_sim, replicate==1 & name=='time' & alg=='SNUTS') |>
  select(nrepars, model, relval, name2, metric)
filter(perf_sim, replicate==1 & name=='ess' & alg=='SNUTS') |>
  select(nrepars, model, relval, name2)
filter(perf_sim, replicate==1 & name=='eff' & alg=='SNUTS') |>
  select(nrepars, model, relval, name2)



g <- ggplot(perf_sim, aes(x=nrepars, y=value, color=alg, group=interaction(alg,replicate))) +
  facet_grid(name2~model, scales='free') + geom_point() +
  geom_line()+
  scale_x_log10() + scale_y_log10() +
  labs(x='Number of random effects',
       y=NULL, color=NULL) +
  guides(col=guide_legend(nrow=1))+
   theme(legend.position = 'inside', legend.position.inside = c(.25,.4))
ggsave('plots/perf_sim.png', g, width=5, height=4.2)


group_by(perf_sim, model) |>
  filter(nrepars==max(nrepars) & replicate==1 & alg=='SNUTS') |>
  select(model, nrepars, alg, name2, relval)




message("Starting case studies...")
message("Reading in large case study fit files...")
xx <- list.files('results/case_studies',  pattern='_fits.RDS', full.names = TRUE)
xx <- xx[!xx %in% paste0('results/case_studies/',skip.models,'_fits.RDS')]
results <- xx |> lapply(readRDS)
mods <- sapply(results, \(x) x[[1]]$model)
ela <- grepl('_ELA', mods)
results_ela <- results[ela]
results <- results[!ela]
mods <- sapply(results, \(x) x[[1]]$model)
names(results) <- mods

sp1 <- extract_sampler_params(results$diamonds[[1]])
sp2 <- extract_sampler_params(results$diamonds[[4]])
sp1$n_leapfrog__ |> mean()
sp2$n_leapfrog__ |> mean()


sp1 <- extract_sampler_params(results$wham[[1]])
sp2 <- extract_sampler_params(results$wham[[4]])
sp1$n_leapfrog__ |> mean()
sp2$n_leapfrog__ |> mean()

message("Making MLE vs posterior correlation and variance plots (slow)...")
cors <- lapply(results, get_cors) |> bind_rows() |>filter(replicate==1 & metric=='Stan default')
vars <- lapply(results, get_vars) |> bind_rows() |>filter(replicate==1 & metric=='Stan default')
g1 <- ggplot(vars, aes(post.sd, marginal.sd, color=model)) +
  geom_abline(intercept=0, slope=1)+geom_point(size=1) +
  scale_y_log10()+scale_x_log10() +
  theme(legend.position='none')  +
  labs(x='Posterior SD', y='Approximate SD (Q)')
#cors2 <- group_by(cors, model) |> slice_sample(n=1)
g2 <- ggplot(cors, aes(post.cor, mle.cor, color=model)) +
  geom_abline(intercept=0, slope=1)+geom_point(size=1) +
  labs(x='Posterior correlation', y='Approximate correlation (Q)', color=NULL) +
    guides(color = guide_legend(ncol = 2)) +
  theme(legend.key.spacing.y = unit(0, "cm"))+
  xlim(-1,1) + ylim(-1,1)
g <- cowplot::plot_grid(g1,g2, rel_widths = c(1,1.9),
                        labels=c('(a)', '(b)'),
                        label_size = 12)
ggsave('plots/post_vs_Q.png', g, width=8, height=2.75, units='in', dpi=300)

message("making ESS change plots")
# get ESS changes with NUTS to SNUTS
ess_results <- lapply(results,
        function(x) {
          mon <- get_mon(x)
          stats <- get_stats(x) |> mutate(metric=metricf(metric))
          out <- merge(x=mon, y=stats, by=c('model', 'replicate', 'metric'), all.x=TRUE)
          select(mutate(out, eff=ESS/time.total), model, metric=metric,
               replicate, par, partype, ESS, time.total, eff)})|>
                    bind_rows() |>
 mutate(alg=algf(metric))
g <- ggplot(ess_results,
            aes(alg, y=eff, color=partype,
                group=interaction(par,replicate))) +
  geom_line(alpha=.25) + geom_jitter(width=.08, size=1)+
  theme(legend.position = 'top')+
  facet_wrap('model', scales='free_y') +scale_y_log10()+
  labs(x=NULL, y='Efficiency (ESS/t)', color=NULL)
ggsave("plots/eff_nuts_snuts_all.png", g, width=8, height=7)



message("Making case study plots...")
stats <- lapply(results, get_stats) |> bind_rows() |>
  filter(metric!='tmbstan') |>
  group_by(model) |>
  mutate(alg=algf(metric),
         rel_eff=eff/mean(eff[metric=='unit']),
         rel_ess=ess/mean(ess[metric=='unit']),
         rel_leapfrog=n.leapfrog/mean(n.leapfrog[metric=='unit']),
         rel_time=time.total/mean(time.total[metric=='unit']),
         #gr.ratio=gr2/gr,
         mean_rel_eff=mean(rel_eff)) |>
  ungroup() |>
  arrange(desc(mean_rel_eff)) |>
  mutate(model=factor(model, levels=unique(model)),
         metricf=metricf(metric))

perf_long <- stats |> pivot_longer(cols=c('rel_ess', 'rel_time', 'rel_leapfrog', 'rel_eff')) |>
  mutate(pre_metric=metricf(metric),
         name2=factor(name, levels=c( 'rel_leapfrog','rel_time', 'rel_ess', 'rel_eff'),
                      labels=c('Trajectory length','Time', 'min ESS', 'Efficiency (ESS/s)')))
g <- ggplot(perf_long, aes(x=model, y=value, color=alg)) +
  geom_hline(yintercept=1) +
  geom_jitter(width=.1, height=0, alpha=.5) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=.5)) +
  #theme(legend.position='top')+
  theme(legend.position='inside', legend.position.inside = c(.9,.35))+
  scale_y_log10() + #coord_flip()+
  facet_wrap('name2', ncol=2, scales='free_y') +
  labs(x=NULL, y='Relative value\n(SNUTS / NUTS)',
       color=NULL)
ggsave('plots/perf_mods.pdf', g, width=7, height=5)


filter(perf_long, model=='wham') |> select(mean_rel_eff)
filter(perf_long, model=='wham' & replicate==1) |> select(mean_rel_eff)
filter(perf_long, model=='diamonds' & replicate==2) |> select(mean_rel_eff)

stats_long <- stats |> mutate(gr.ratio=log10(gr.ratio)) |>
  select(model, metric, pct.divs, rhat, gr.ratio) |>
    pivot_longer(c('pct.divs', 'gr.ratio', 'rhat')) |>
  mutate(alg=algf(metric),
    name2=factor(name, levels=c('gr.ratio', 'pct.divs', 'rhat'),
    labels=c('log10 ratio gradient', '% Divergences', 'max Rhat'))) |>
  select(model, alg, metric, value, name2)
g <- ggplot(stats_long, aes(model, y=value, color=alg)) +
  geom_jitter(width=.1) +
  facet_wrap('name2', ncol=1, scales='free_y') +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=.5))
ggsave('plots/case_studies_stats.png', g, width=4, height=6)

message("Making table of case study timings...")
grad_timings <- stats |>
  #each metric has it's own gr and gr2, so don't need unit here
  filter(metric!='unit') |>
  group_by(model) |>
  summarize(metric=metric[1],
            # gr=1000*mean(gr),
            # gr2=1000*mean(gr2),
            # gr.ratio=mean(gr.ratio),
            pct.Qall=100*mean(time.Qall/time.total)) |>
  arrange(model)
bench <- read.csv('results/bench_casestudies.csv')

# remove unneccesary models
mods_timings <- merge(bench, grad_timings, by='model', all=TRUE) |>
  filter(!model=='wildf_adapted') |>
  select(model, metric=metric.x, npar, pct.sparsity, rel_gr, pct.Qall)
stopifnot(sum(is.na(mods_timings))==0)
write.csv(mods_timings, file='results/timings_casestudies_table.csv', row.names=FALSE)


message("Making plot of wildf issues...")
fitswildf2 <- readRDS('results/wildf2_fits.RDS')
fitswildf3 <- readRDS('results/wildf3_fits.RDS')
post2 <- lapply(fitswildf2, \(x) cbind(adapt_delta=.9, metric=x$metric, as.data.frame(x))) |>
  bind_rows()
post3 <- lapply(fitswildf3, \(x) cbind(adapt_delta=.99, metric=x$metric, as.data.frame(x))) |>
  bind_rows()

# pull out the mode and cor for two worst params
ind <- which(fitswildf3[[1]]$mle$parnames %in% c('plantSlopeSD', 'slope'))
mle <- fitswildf3[[1]]$mle$est[ind]
cors <- fitswildf3[[1]]$mle$cor[ind,ind]
se <- fitswildf3[[1]]$mle$se[ind]
cov <- cors*(se %*% t(se))

# Rotate the posterior to see sampling space better
L <- t(chol(fitswildf3[[1]]$mle$Qinv))
Linv <- solve(L)
post3=post2
post3_rotated <- cbind(post3[,1:2], t(Linv %*% t(as.matrix(post3[, -(1:2)]))))
mle_rotated <- Linv %*% as.numeric(fitswildf3[[1]]$mle$est)
mle_rotated <- mle_rotated[ind]

png('plots/wildf_pairs.png', width=6, height=3, units='in', res=300)
xlim1 <- range(c(post3$plantSlopeSD))
ylim1 <- range(c(post3$slope))
xlim2 <- range(c(post3_rotated$plantSlopeSD))
ylim2 <- range(c(post3_rotated$slope))
par(mfrow=c(1,2), mar=c(2.5,2.5,1.5,.5), mgp=c(1.25,.25,0), tck=-.01,
    cex.axis=.8, col.axis=gray(.4), cex.lab=.9)
with(post3, plot(plantSlopeSD, slope, col=ifelse(metric=='unit', 1,2), cex=.5,
                 xlab='log plantSlopeSD', xlim=xlim1, ylim=ylim1))
mtext('Model space', line=.25)
lines(ellipse::ellipse(cors[1,2], centre=mle, scale=se), col='blue')
points(mle[1], mle[2], pch=17, cex=1, col='blue')
with(post3_rotated, plot(plantSlopeSD, slope, col=ifelse(metric=='unit', 1,2), cex=.5,
                         xlab='log plantSlopeSD',
                         xlim=xlim2, ylim=ylim2))
mtext('Transformed via Q', line=.25)
lines(ellipse::ellipse(0, centre=mle_rotated, scale=c(1,1)), col='blue')
points(mle_rotated[1], mle_rotated[2], pch=17, cex=1, col='blue')
box(col=gray(.4))
legend('bottomright', legend=c('unit', 'sparse'), pch=1, cex=.8,  col=1:2, bty='n')
dev.off()




message("Making ELA plots...")
# embedded laplace approximation (ELA) results compared to 'auto'
# for full integration
stats_ela <- lapply(results_ela, get_stats) |> bind_rows() |>
  mutate(model=gsub('_ELA', '', model)) |>
  # stupid typo in model name for ELA, quick fix here for now
  mutate(model=ifelse(model=='sdmbTMB', 'sdmTMB', model)) |>
  filter(!model %in% c('glmmTMB', 'VAR', 'spde', 'wildf2', 'wildf3')) |>
  group_by(model) |>
  mutate(rel_eff=eff/mean(eff[metric=='unit']),
         rel_ess=ess/mean(ess[metric=='unit']),
         rel_time=time.total/mean(time.total[metric=='unit']),
        # gr.ratio=gr2/gr,
         mean_rel_eff=mean(rel_eff)) |>
  ungroup() |>
  arrange(desc(mean_rel_eff)) |>
  mutate(model=factor(model, levels=unique(model)))
# #### old, overly complicated way to do it
# # very carefully massage output to plot for two scenariors: NUTS and SNUTS
# stats_tmp <- bind_rows(mutate(stats_ela, ELA=TRUE), mutate(stats, ELA=FALSE)) |>
#   filter(metric!='tmbstan') |>
#   select(model, metric, eff, ELA) |>
#   arrange(model, metric, ELA) |> group_by(model) |>
#   mutate(eff_rel=eff/mean(eff[metric=='unit' & ELA==FALSE]),
#     eff_rel_unit=eff/mean(eff[metric=='unit' & ELA==TRUE]),
#          eff_rel_full=eff/mean(eff[metric!='unit' & ELA==FALSE])) |>
#   filter(ELA==TRUE) |> ungroup()
# # this one is efficiency of ELA:full with both using 'unit'
# # metric so comparable to previous studies
# x1 <- filter(stats_tmp, metric=='unit') |> mutate(y=eff_rel, type='NUTS')
# # now look at ELA:full using SNUTS 'auto' for both
# x2 <- filter(stats_tmp, metric!='unit') |> mutate(y=eff_rel_full, type='SNUTS')
# xx <- bind_rows(x1,x2) |> ungroup() |> mutate(modelf=modelf(model), metricf=metricf(metric))
# g <- ggplot(xx, aes(modelf, y=y, color=metricf, shape=type)) +
#   geom_hline(yintercept=1) +
#   geom_jitter(width=.1, height=0, alpha=.7) +
#   theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=.5)) +
#   #scale_y_log10() + #coord_flip()+
#   labs(x=NULL, y='Relative Efficiency\n(ELA / full NUTS)',
#        color=NULL, shape=NULL)+
#   theme(legend.position='right')
# g
#

stats_tmp <- bind_rows(mutate(stats_ela, ELA=TRUE), mutate(stats, ELA=FALSE)) |>
  filter(metric!='tmbstan') |>
  mutate(alg=algf(metric), model=modelf(model)) |>
  select(model, metric, eff, ELA, alg) |>
  arrange(model, metric, ELA) |>
  group_by(model) |>
  mutate(eff_rel=eff/mean(eff[alg=='SNUTS' & ELA==FALSE])) |>
  ungroup() |> filter(ELA==TRUE, is.finite(eff_rel)) |>
  mutate(alg2=factor(alg, levels=c('NUTS','SNUTS'), labels=c('ELA-NUTS', 'ELA-SNUTS')))

g <- ggplot(stats_tmp, aes(x=model, shape=alg2, color=alg2, y=eff_rel)) +
 # facet_wrap('alg', ncol=1) +
  geom_hline(yintercept=1) +
  geom_jitter(width=.1, height=0, alpha=.7) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=.5)) +
  scale_y_log10() + #coord_flip()+
  labs(x=NULL, y='Efficiency relative to full SNUTS',
       color=NULL, shape=NULL) +
  #  theme(legend.position='bottomleft') +
  guides(col=guide_legend(nrow=2))+
  theme(legend.position='inside', legend.position.inside = c(.21,.1),
        legend.title = element_text(size = 0), # Adjust title font size
        legend.text = element_text(size = 8),   # Adjust label font size
        axis.title.y=element_text(size=10),
    #    axis.text=element_blank(), axis.ticks=element_blank(),
        legend.key.size = unit(0.2, "cm") # Adjust the size of the key box/symbol area
  )
g
ggsave('plots/perf_ELA.pdf', g, width=width2, height=3, units='in')

if(TRUE){
  message("skipping reference posterior lists since fits did not change..")
  message("!!! make sure to rerun pathfinder if model fits change !!!")
} else {

  message("Creating reference posterior list for pathfinder...")
  ref_posts_list <- lapply(results, \(mod) {
    #print(mod[[1]]$model)
    df=lapply(mod, \(x)
              cbind(model=x$model, metric=x$metric, replicate=x$replicate,
                    as.data.frame(x))) |> bind_rows()
    if('sparse' %in% df$metric){
      df = df |> filter(metric=='sparse')
    } else if('dense' %in% df$metric){
      df = df |> filter(metric=='dense')
    } else if('diag' %in% df$metric){
      df = df |> filter(metric=='diag')
    } else {stop("no Q metric found in results file for mod", x$model)}
    # thin down to 10k draws
    if(nrow(df)>10000) df=df[seq(from=1, to=nrow(df), len=10000),]
    df
  }
  )
  # convert to named list for easier access later
  names(ref_posts_list) <- lapply(ref_posts_list, \(x) x$model[1]) |> unlist()
  saveRDS(ref_posts_list, 'results/ref_posts_list.RDS')
}

message("Making Pathfinder plots...")
warmups <- lapply(list.files('results/warmups/', full.names = TRUE), readRDS)
sp <- lapply(warmups, function(y){
  lapply(y, function(x){
    cat(x$model, x$metric, x$replicate, "\n")
    cbind(model=x$model, metric=x$metric, replicate=x$replicate, plot_sampler_params(x, plot=FALSE)$data )
  }
  )
}
) |> bind_rows()
eps <- filter(sp, variable=='log_stepsize' & iteration < 1000)
g <- ggplot(eps, aes(iteration, y=value,
                     group=interaction(chain,replicate, metric))) +
  geom_line(alpha=.1) +
  facet_wrap('model', ncol=3,scale='free_y') + labs(y='log step-size')
ggsave('plots/warmup.png',g, width=8, height=8)

message("Done processing results!")
