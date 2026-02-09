## this file processes the output and makes plots and tables and figures (for the MS)


setwd(here::here())
source("code/startup.R")
skip.models <- c('ratios', 'ratio', 'RTMB', 'simple', 'wildf2', 'wildf3', 'wildf4', 'cors', 'cor','VAR','glmmTMB')
skip.models <- c(skip.models, paste0(skip.models, '_ELA'))

# for a Wiley journal the figure widths are:
width1 <- 180/25.4 # 180mm full width in inches
width2 <- 80/25.4 # half width

message("Making simple timing figure 1")
fits_ratio <- readRDS('results/ratio_fits.RDS')
stats.sd <- readRDS('results/ratio_stats.RDS') |>
  group_by(model, sd, metric) |>
  summarize(time.warmup=mean(time.warmup), time.sampling=mean(time.sampling), .groups='drop') |>
  mutate(time=time.warmup+time.sampling, metric=metricf(metric)) |>
  mutate(metric=case_when(metric=='diag'~'Descale', metric=='dense'~'Descale + decorrelate', metric=='Stan default'~'Stan default'))
stats.cor <- readRDS('results/cor_stats.RDS') |>
  group_by(model, cor, metric) |>
  summarize(time.warmup=mean(time.warmup), time.sampling=mean(time.sampling), .groups='drop') |>
  mutate(time=time.warmup+time.sampling, metric=metricf(metric)) |>
  mutate(metric=case_when(metric=='diag'~'Descale', metric=='dense'~'Descale + decorrelate', metric=='Stan default'~'Stan default'))
timing.sd <- select(stats.sd, sd, metric, time.warmup, time.total=time) %>%
  pivot_longer(-c(sd,metric), names_to='sampling', values_to='time') %>%
  mutate(sampling=gsub('time.','',sampling))
timing.cor <- select(stats.cor, cor, metric, time.warmup, time.total=time) %>%
  pivot_longer(-c(cor,metric), names_to='sampling', values_to='time') %>%
  mutate(sampling=gsub('time.','',sampling))
g1 <- ggplot(timing.sd, aes(sd, y=time, color=metric, lty=sampling)) +
  geom_line() +
  scale_x_log10() + scale_y_log10() +
  labs(x=expression(paste('Ratio of marginal SDs (', tau, ')')),
       y='Time (s)', color=NULL, lty=NULL) +
  theme(legend.position='inside', legend.position.inside = c(.35,.68)) +
  theme(legend.text = element_text(size = 9),   # Adjust label font size
              legend.key.size = unit(.4, "cm"))
cors <- c(0,.5,.9,.99,.999, .9999)
unique(stats.cor)
cors.lab <- -log(2/(cors+1)-1)
g2 <- #ggplot(timing.cor, aes(log(cor)-log(1-cor), y=time, color=metric, lty=sampling)) +
  ggplot(timing.cor, aes(x=-log(2/(cor+1)-1), y=time, color=metric, lty=sampling)) +
  geom_line() +
  scale_y_log10() +
  scale_x_continuous(breaks=cors.lab, labels=cors) +
  labs(x=expression(paste('Correlation (', rho, ')')),
       y='Time (s)', color=NULL)+
  #theme(legend.position='inside', legend.position.inside = c(.2,.65))
  theme(legend.position='none')
g <- cowplot::plot_grid(g1,g2, labels=c('(a)', '(b)'), label_size=10, ncol=2)
ggsave("plots/fig1_example_timings.pdf", g, width=width1, height=3, units='in')



message("Making Q vs correlation vs L figure 2 from sdmTMB model")
tmp <- readRDS('results/case_studies/sdmTMB_fits.RDS')[[1]]
# really hack this together to get a correlation and sparse
# precision matrix to show up on the same figure and be
# interpretable
library(Matrix)
Q <- tmp$mle$Q
Q <- as(Q, 'sparseMatrix')
dimnames(Q) <- NULL
M <- solve(as.matrix(Q))
Q2 <- as.matrix(Q)
Q2[Q2==0] <- NA
Q2[Q2!=0] <- 1.1
Lsparse <- as(Cholesky(Q, LDL=FALSE, super=TRUE, perm=FALSE), 'sparseMatrix')
Lsparse <- as.matrix(Lsparse)
Lsparse[Lsparse==0] <- NA
Lsparse[Lsparse!=0] <- 1.1
Lperm <- as(Cholesky(Q, LDL=FALSE, super=TRUE, perm=TRUE), 'sparseMatrix')
Lperm <- as.matrix(Lperm)
Lperm[Lperm==0] <- NA
Lperm[Lperm!=0] <- 1.1
test <- bind_rows(cbind(x='Q',reshape2::melt(Q2)),
                  cbind(x='Corr', reshape2::melt(cov2cor(M))),
                  cbind(x='Lsparse', reshape2::melt(Lsparse)),
                  cbind(x='Lperm', reshape2::melt(Lperm))) |>
  filter(!is.na(value)) |>
  # reorder to look like matrices
  mutate(Var1=max(Var1)-Var1)
test$x <- factor(test$x, levels=c('Q', 'Corr', 'Lsparse', 'Lperm'),
                 labels=c('Precision~(Q)',
                          'Correlation', expression(Chol(Q)),
                          expression(Chol(PQP^T))))
                 #labels=c('Q','Correlation', 'Chol(Q)', "Chol(PQP')"))
g <- ggplot(test, aes(Var2, Var1, color=value)) + geom_point(pch='.') +
  facet_wrap('x', labeller=label_parsed) + labs(x=NULL, y=NULL) +
  theme(axis.text=element_blank(), axis.ticks=element_blank()) +
  scale_color_gradient2(
    high = "blue",      # Color for strong negative correlations (-1)
    mid = 'white',     # Color for zero correlation (0)
    low = "red",      # Color for strong positive correlations (1)
    midpoint = 0,      # Ensures the midpoint color is at 0
    limit = c(-1, 1),  # Sets the full range of correlation values
    name = NULL)+
  theme(axis.text=element_blank(), axis.ticks=element_blank()) +
  theme(legend.position = 'inside', legend.position.inside = c(.92,.8))+
  theme(legend.key.size = unit(0.35, "cm")) # Adjust the size of the key box/symbol area
ggsave('plots/fig2_sparsity_patterns.pdf', g, width=width1, height=5.5, units='in', dpi=400)

message("Making figure 3 benchmark metrics")
bench <- readRDS('results/bench_spde_gr.RDS')
g <- plot.bench(bench)
ggsave('plots/fig3_spde_gradient_bechmark.pdf', g, width=width1, height=2.5, units='in')
#pct.sparsity <- round(100*mean(as.matrix(Q)[lower.tri(Q)] == 0),2)
#xx <- as.matrix(cov2cor(Qinv))
#maxcor <- max(abs(xx[lower.tri(xx)]))
#hist(abs(xx[lower.tri(xx)]))

message("making figure 4 Pathfinder")
w.all <- readRDS('results/w1d_all.RDS')
w.all$type <- gsub(x=w.all$type, replacement = 'TMB', pattern = 'Q')
w.means <- group_by(w.all, model, type) |>
  summarize(w1d_rel=mean(w1d_rel), time_rel=mean(time_rel), .groups='drop')
# order by worst to best for Q
lvls <- filter(w.means, type=='Pathfinder') |> arrange(w1d_rel) |> pull(model)
w.all <- mutate(w.all, model=factor(model, levels=lvls))
w.means <- mutate(w.means, model=factor(model, levels=lvls))
w.all.long <- pivot_longer(w.all, cols=c('w1d_rel', 'time_rel')) |>
  mutate(name2=factor(name, levels=c('time_rel', 'w1d_rel'),
                      labels=c('Time', '1D Wasserstein distance')))
w.means <- group_by(w.all.long, model, type, name2) |>
  summarize(value=mean(value), .groups='drop')
g <- ggplot(w.all.long, aes(x=model, y=value, color=type)) +
  geom_hline(yintercept=1, color=gray(.5))+
  geom_jitter(width=.1, alpha=.5, pch=1, size=.8)+ facet_wrap('name2', ncol=2) +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=.5)) +
  theme(legend.position = 'inside', legend.position.inside = c(.22, .35))+
  labs(y='Value relative to TMB', color=NULL, x=NULL) +
  geom_point(data=w.means, pch=17, size=2)
ggsave('plots/fig4_pathfinder_all.pdf', g, width=width2, height=5, dpi=300)



message("Processing simulated results to make figure 5...")
stats_sim <- bind_rows(readRDS('results/spde_stats.RDS'),
                       #readRDS('results/VAR_stats.RDS'),
                       readRDS("results/glmmTMB_stats.RDS")) |>
 # group_by(model, metric, nrepars) |> #summarize(n=n()) |> print(n=Inf)
  filter(! (model=='VAR' & nrepars <205)) |>
  ungroup()

perf_sim <- stats_sim |>
  mutate(alg=algf(metric),
         time=time.warmup+time.sampling) |>
  #group_by(model, alg, nrepars) |>
  # mutate(eff_avg=mean(eff),
  #        L_avg=mean(n.leapfrog),
  #        ess_avg=mean(ess),
  #        time_avg=mean(time)) |> ungroup() |>
   pivot_longer(cols=c('eff', 'n.leapfrog', 'ess', 'time')[-2]) |>
  # relative to unit metric
  group_by(model, nrepars) |>
  mutate(relval=value/value[alg=='NUTS']) |> ungroup() |>
  mutate(name2=factor(name, levels=c('n.leapfrog', 'time', 'ess', 'eff'),
                      labels=c('Trajectory length', 'Total time (s)', 'Minimum ESS', 'Efficiency (ESS/s)')))
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
  labs(x='Number of parameters',
       y=NULL, color=NULL) +
  guides(col=guide_legend(nrow=1))+
   theme(legend.position = 'inside', legend.position.inside = c(.25,.4))
ggsave('plots/fig5_perf_sim.pdf', g, width=width1, height=5)


group_by(perf_sim, model) |>
  filter(nrepars==max(nrepars) & replicate==1 & alg=='SNUTS') |>
  select(model, nrepars, alg, name2, relval)




## need to filter out the tmbstan runs
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


message("Making case study plots...")
stats <- lapply(results, get_stats) |> bind_rows() |>
#  filter(metric!='tmbstan') |>
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
  mutate(metric=metricf(metric),
         name2=factor(name, levels=c( 'rel_leapfrog','rel_time', 'rel_ess', 'rel_eff'),
                      labels=c('Mean Trajectory length','Total Time', 'Minimum ESS', 'Efficiency (ESS/s)')))
g <- ggplot(filter(perf_long, metric!='tmbstan'),
            aes(x=model, y=value, color=alg)) +
  geom_hline(yintercept=1) +
  geom_jitter(width=.1, height=0, alpha=.5) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=.5)) +
  #theme(legend.position='top')+
  theme(legend.position='inside', legend.position.inside = c(.9,.35))+
  scale_y_log10() + #coord_flip()+
  facet_wrap('name2', ncol=2, scales='free_y') +
  labs(x=NULL, y='Relative value\n(SNUTS / NUTS)',
       color=NULL)
ggsave('plots/fig6_perf_mods.pdf', g, width=width1, height=5)

# tmbstan comparison to unit SNUTS, shows overhead of R and StanEstimators
g <- ggplot(filter(perf_long, metric %in%c('unit', 'tmbstan')),
            aes(x=model, y=value, color=metric)) +
  geom_hline(yintercept=1) +
  geom_jitter(width=.1, height=0, alpha=.5) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=.5)) +
  #theme(legend.position='top')+
  theme(legend.position='inside', legend.position.inside = c(.9,.35))+
  scale_y_log10() + #coord_flip()+
  facet_wrap('name2', ncol=2, scales='free_y') +
  labs(x=NULL, y='Relative value\n(SNUTS / NUTS)',
       color=NULL)

g

# stats to report for the MS
stats.means <- group_by(stats, model, alg) |>
  summarize(mean_leapfrog=mean(n.leapfrog),
            mean_time=mean(time.total),
            mean_ess=mean(ess),
            mean_eff=mean(eff))
stats.rel <- filter(stats, alg=='SNUTS') |> group_by(model) |>
  summarize(mean_rel_leapfrog=mean(rel_leapfrog),
            mean_rel_time=mean(rel_time),
            mean_rel_ess=mean(rel_ess),
            mean_rel_eff=mean(rel_eff))
filter(stats.rel, model=='wham')
filter(stats.rel, model=='diamonds')
filter(stats.means, model=='diamonds')


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

# get some stats to report in the MS
stats_ela |> filter(model %in% c('causal', 'sam')) |>
  group_by(model, metric) |>
  summarize_all(mean)


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
ggsave('plots/fig7_perf_ELA.pdf', g, width=width2, height=3, units='in')

# get stats to put in MS

# how much faster is SNUTS vs NUTS when doing ELA?
tmp <- stats_tmp |> group_by(model) |>
  summarize(xx=mean(eff_rel[alg=='SNUTS']/mean(eff_rel[alg=='NUTS']))) |>
  arrange(xx)
tmp
mean(tmp$xx)
# which models outperform full SNUTS?
filter(stats_tmp) |> group_by(model, metric,alg2) |>
  summarize(mean_eff_rel=mean(eff_rel)) |>
  arrange(desc(mean_eff_rel))


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
#stopifnot(sum(is.na(mods_timings))==0)
write.csv(mods_timings, file='results/timings_casestudies_table.csv', row.names=FALSE)


message("Making MLE vs posterior correlation and variance plots (slow)...")
cors <- lapply(results, get_cors) |> bind_rows() |>filter(replicate==1 & metric!='unit')
vars <- lapply(results, get_vars) |> bind_rows() |>filter(replicate==1 & metric!='unit')
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
ggsave('plots/post_vs_Q.png', g, width=width1, height=2.75, units='in', dpi=300)



message("making ESS change plots")
# get ESS changes with NUTS to SNUTS
ess_results <- lapply(results,
                      function(x) {
                                                mon <- get_mon(x)
                        stats <- get_stats(x) |> mutate(metric=metricf(metric))
                        out <- merge(x=mon, y=stats, by=c('model', 'replicate', 'metric'), all.x=TRUE)
                        select(mutate(out, eff=ESS/time.total), model, metric=metric,
                               replicate, par, partype, ESS, time.total, eff)})|>
  bind_rows() |> filter(metric!='tmbstan') |>
  mutate(alg=algf(metric))
g <- ggplot(ess_results,
            aes(alg, y=eff, color=partype,
                group=interaction(par,replicate))) +
  geom_line(alpha=.25) + geom_jitter(width=.08, size=1)+
  theme(legend.position = 'top')+
  facet_wrap('model', scales='free_y') +scale_y_log10()+
  labs(x=NULL, y='Efficiency (ESS/t)', color=NULL)
ggsave("plots/eff_nuts_snuts_all.png", g, width=width1, height=7)
g <- ggplot(ess_results,
            aes(alg, y=ESS, color=partype,
                group=interaction(par,replicate))) +
  geom_line(alpha=.25) + geom_jitter(width=.08, size=1)+
  theme(legend.position = 'top')+
  facet_wrap('model', scales='free_y') +scale_y_log10()+
  labs(x=NULL, y='ESS', color=NULL)
ggsave("plots/ess_nuts_snuts_all.png", g, width=width1, height=7)




message("Making plot of wildf issues...")
# this is the extra long run to get more stable estimates
fitswildf2 <- readRDS('results/case_studies/wildf_fits2.RDS')
xx <- plot_uncertainties(fitswildf2[[3]],plot=FALSE)
g1 <- ggplot(xx, aes(sd.post, sd.mle)) + geom_point(pch=1) +
  labs(x='Posterior marginal SD', y='Asymptotic marginal SD (Q)') +
  geom_abline(slope=1, intercept=0)


ess_wildf2 <- lapply(list(fitswildf2),
                      function(x) {
                        mon <- get_mon(x)
                        stats <- get_stats(x) |> mutate(metric=metricf(metric))
                        out <- merge(x=mon, y=stats, by=c('model', 'replicate', 'metric'), all.x=TRUE)
                        select(mutate(out, eff=ESS/time.total), model, metric=metric,
                               replicate, par, partype, ESS, time.total, eff)})|>
  bind_rows() |>
  mutate(alg=algf(metric))
ess_wildf2$model |> unique()
ess_wildf2$alg2 <- as.character(ess_wildf2$alg)
ess_wildf2$alg2[ess_wildf2$model=='wildf2_adapt'] <- 'SNUTS\n+ adaptation'
g2 <- ggplot(ess_wildf2,
            aes(alg2, y=ESS, color=partype,
                group=interaction(par,replicate))) +
  geom_line(alpha=.25) + geom_jitter(width=.01, size=1)+
  theme(legend.position = 'inside', legend.position.inside = c(.4,.93))+
  labs(x=NULL, y='ESS', color=NULL) +
  theme(legend.text = element_text(size = 8),   # Adjust label font size
        legend.key.size = unit(.1, "cm")) +
  guides(col=guide_legend(nrow=1)) + scale_y_log10(limits=c(NA,65000))


ind <- which(fitswildf2[[1]]$par_names %in% c('plantSlopeSD', 'slope'))
mle <- fitswildf2[[1]]$mle$est[ind]
cors <- fitswildf2[[1]]$mle$cor[ind,ind]
se <- fitswildf2[[1]]$mle$se[ind]
cov <- cors*(se %*% t(se))
tmp <- ellipse::ellipse(cors[1,2], centre=mle, scale=se) |> as.data.frame()
n <- nrow(tmp)
ellipsedf <- data.frame(x=tmp$x[-n], y=tmp$y[-n], xend=tmp$x[-1], yend=tmp$y[-1])
wildf2post <- as.data.frame(fitswildf2)
#pfwildf <- readRDS('results/pathfinder/wildf_draws.RDS')
g3 <- ggplot(wildf2post, aes(plantSlopeSD,slope)) + geom_point(alpha=.3, size=.15) +
  labs(x='log plantSlopeSD', y='log slope') +
   geom_point(data=data.frame(x=mle[1], y=mle[2]), mapping=aes(x,y), color='red') +
  geom_segment(data=ellipsedf, mapping=aes(x=x,y=y, xend=xend,yend=yend),
               color='red', lwd=1)
g <- cowplot::plot_grid(g1,g2,g3, labels = c('(a)', '(b)', '(c)'), label_size = 12,
                        ncol = 3)
ggsave("plots/wildf_mismatch.png", g, width=width1*1.3, height=3)
# ## This is the old way of doing it. I don't think it makes sense to plot qprime
# ## since that's just a linear combination of other parameters and not directly
# ## interetable as the origina pars
# pairs(fits.wildf3[[2]], order='slow')
# post3 <- lapply(fitswildf3, \(x) cbind(adapt_delta=.999, metric=x$metric, as.data.frame(x))) |>
#   bind_rows()
# post3 <- filter(post3, metric=='sparse') |> slice_sample(n=500)
# # pull out the mode and cor for two worst params
# ind <- which(fitswildf3[[1]]$par_names %in% c('plantSlopeSD', 'slope'))
# mle <- fitswildf3[[1]]$mle$est[ind]
# cors <- fitswildf3[[1]]$mle$cor[ind,ind]
# se <- fitswildf3[[1]]$mle$se[ind]
# cov <- cors*(se %*% t(se))
# # Rotate the posterior via Sigma to see sampling space better
# L <- t(chol(solve(as.matrix(fitswildf3[[1]]$mle$Q))))
# Linv <- solve(L)
# post3_rotated <- cbind(post3[,1:2], t(Linv %*% t(as.matrix(post3[, -(1:2)]))))
# mle_rotated <- Ltinv %*% as.numeric(fitswildf3[[1]]$mle$est)
# mle_rotated <- mle_rotated[ind]
# # Rotate the posterior via Q to see sampling space better
# Lt <- chol(as.matrix(fitswildf3[[1]]$mle$Q)) # upper triangle
# post3_rotated <- cbind(post3[,1:2], t(Lt %*% t(as.matrix(post3[, -(1:2)]))))
# mle_rotated <- Lt %*% as.numeric(fitswildf3[[1]]$mle$est)
# mle_rotated <- mle_rotated[ind]
# png('plots/wildf_pairs.png', width=6, height=3, units='in', res=300)
# xlim1 <- range(c(post3$plantSlopeSD))
# ylim1 <- range(c(post3$slope))
# xlim2 <- range(c(post3_rotated$plantSlopeSD))
# ylim2 <- range(c(post3_rotated$slope))
# par(mfrow=c(1,2), mar=c(2.5,2.5,1.5,.5), mgp=c(1.25,.25,0), tck=-.01,
#     cex.axis=.8, col.axis=gray(.4), cex.lab=.9)
# with(post3, plot(plantSlopeSD, slope, col=ifelse(metric=='unit', 1,2), cex=.5,
#                  xlab='log plantSlopeSD', xlim=xlim1, ylim=ylim1))
# mtext('Original', line=.25)
# lines(ellipse::ellipse(cors[1,2], centre=mle, scale=se), col='blue')
# points(mle[1], mle[2], pch=17, cex=1, col='blue')
# with(post3_rotated, plot(plantSlopeSD, slope, col=ifelse(metric=='unit', 1,2), cex=.5,
#                          xlab='log plantSlopeSD',
#                          xlim=xlim2, ylim=ylim2))
# mtext('Decorrelated via Q', line=.25)
# lines(ellipse::ellipse(0, centre=mle_rotated, scale=c(1,1)), col='blue')
# points(mle_rotated[1], mle_rotated[2], pch=17, cex=1, col='blue')
# box(col=gray(.4))
# legend('bottomright', legend=c('Stan default', 'sparse'), pch=1, cex=.8,  col=1:2, bty='n')
# dev.off()
#

message("Making NUTS trajectory plots")
source("code/make_NUTS_trajectories.R")


if(TRUE){
  warning("skipping reference posterior lists since fits did not change..")
  warning("!!! make sure to rerun pathfinder if model fits change !!!")
} else {

  message("Creating reference posterior list for pathfinder...")
  ref_posts_list <- lapply(results, \(mod) {
    df=lapply(mod, \(x){
      if(x$metric=='tmbstan') return(NULL)
              message(x$model,' ', x$metric)
              cbind(model=x$model, metric=x$metric, replicate=x$replicate,
                    as.data.frame(x))}) |> bind_rows()
    if('sparse' %in% df$metric){
      df = df |> filter(metric=='sparse')
    } else if('dense' %in% df$metric){
      df = df |> filter(metric=='dense')
    } else if('diag' %in% df$metric){
      df = df |> filter(metric=='diag')
    } else {stop("no Q metric found in results file for mod ", x$model)}
    # thin down to 10k draws
    if(nrow(df)>10000) df=df[seq(from=1, to=nrow(df), len=10000),]
    df
  }
  )
  # convert to named list for easier access later
  names(ref_posts_list) <- lapply(ref_posts_list, \(x) x$model[1]) |> unlist()
  saveRDS(ref_posts_list, 'results/pathfinder/ref_posts_list.RDS')
}

message("Making warmup plots...")
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
