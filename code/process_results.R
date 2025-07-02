

setwd(here::here())
source("code/startup.R")

skip.models <- c('ratios', 'ratio', 'RTMB', 'simple', 'wildf2', 'wildf3', 'wildf4', 'cors', 'cor')
skip.models <- c(skip.models, paste0(skip.models, '_ELA'))

stats_sim <- bind_rows(readRDS('results/spde_stats.RDS'),
                       readRDS('results/VAR_stats.RDS'),
                       readRDS("results/glmmTMB_stats.RDS")) |>
 # group_by(model, metric, nrepars) |> #summarize(n=n()) |> print(n=Inf)
  filter(! (model=='VAR' & nrepars <205)) |>
  ungroup()

perf_sim <- stats_sim |>
  group_by(model, metric, nrepars)  |>
  # mean across replicates for each metric
  mutate(eff_avg=mean(eff),
         ess_avg=mean(ess),
         time.total=time.warmup+time.sampling,
         time_avg=mean(time.total)) |> ungroup() |>
  pivot_longer(cols=c('eff_avg', 'ess_avg', 'time_avg')) |>
  # relative to unit metric
  group_by(model, nrepars) |>
  mutate(relval=value/value[as.character(metric)=='unit']) |> ungroup() |>
  mutate(metricf=metricf(metric),
         name2=factor(name, levels=c('time_avg', 'ess_avg', 'eff_avg'),
                      labels=c('Total time', 'min(ESS)', 'Efficiency (ESS/s)')))

# #, rel_eff=eff/eff_avg[as.character(metric)=='unit'])
# g <- ggplot(filter(perf_sim, metric!='unit'),
#             aes(x=nrepars, y=relval, color=metric)) +
#   facet_grid(name2~model, scales='free') + geom_point() + geom_line()+
#   scale_x_log10() + scale_y_log10() +
#   labs(x='Number of random effects',
#        y="Mean value relative to metric 'unit'",
#        color=NULL) +
#   geom_hline(yintercept = 1)
# g
# ggsave('plots/perf_sim_rel.png', g, width=7, height=5)

g <- ggplot(perf_sim, aes(x=nrepars, y=value, color=metricf)) +
  facet_grid(name2~model, scales='free') + geom_point() + geom_line()+
  scale_x_log10() + scale_y_log10() +
  labs(x='Number of random effects',
       y=NULL, color=NULL)
ggsave('plots/perf_sim.png', g, width=7, height=4)


stats_sim |> str()



message("Reading in large case study fit files...")
xx <- list.files('results',  pattern='_fits.RDS', full.names = TRUE)
xx <- xx[!xx %in% paste0('results/',skip.models,'_fits.RDS')]
results <- xx |> lapply(readRDS)
mods <- sapply(results, \(x) x[[1]]$model)

ela <- grepl('_ELA', mods)
results_ela <- results[ela]
results <- results[!ela]

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

# get ESS changes with NUTS to SNUTS
ess_results <- lapply(results,
        function(x) {
          mon <- get_mon(x)
          stats <- get_stats(x) |> mutate(metric=metricf(metric))
          out <- merge(x=mon, y=stats, by=c('model', 'replicate', 'metric'), all.x=TRUE)
          select(mutate(out, eff=ESS/time.total), model, metric=metric,
               replicate, par, partype, ESS, time.total, eff)})|>
                    bind_rows() |>
 mutate(SNUTS=ifelse(as.character(metric)=='Stan default', 'NUTS', 'SNUTS'))
g <- ggplot(ess_results,
             aes(SNUTS, y=eff, color=partype,
                 group=interaction(par,replicate))) +
  geom_line(alpha=.25) + geom_jitter(width=.05)+
  theme(legend.position = 'top')+
  facet_wrap('model') + scale_y_log10()+
  labs(x=NULL, y='Efficiency (ESS/t)', color=NULL)
ggsave("plots/essf_nuts_snuts_all.png", g, width=6, height=6)
g <- ggplot(ess_results,
            aes(SNUTS, y=ESS, color=partype,
                group=interaction(par,replicate))) +
  geom_line(alpha=.25) + geom_jitter(width=.05)+
  theme(legend.position = 'top')+
  facet_wrap('model') + scale_y_log10()+
  labs(x=NULL, y='Effective sample size (ESS)', color=NULL)
ggsave("plots/ess_nuts_snuts_all.png", g, width=6, height=6)

message("Creating reference posterior list...")
ref_posts_list <- lapply(results, \(mod) {
  print(mod[[1]]$model)
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

stats <- lapply(results, get_stats) |> bind_rows() |>
  group_by(model) |>
  mutate(rel_eff=eff/mean(eff[metric=='unit']),
         rel_ess=ess/mean(ess[metric=='unit']),
         rel_time=time.total/mean(time.total[metric=='unit']),
         gr.ratio=gr2/gr,
         mean_rel_eff=mean(rel_eff)) |>
  ungroup() |>
  arrange(desc(mean_rel_eff)) |>
  mutate(model=factor(model, levels=unique(model)),
         metricf=metricf(metric))
# g1 <- ggplot(stats, aes(x=model, y=rel_time, color=metric)) +
#   geom_jitter(width=.1, height=0) +
#   theme(axis.text.x = element_text(angle = 90)) +
#   scale_y_log10() + geom_hline(yintercept=1) +
#   labs(x=NULL, y='Wall time relative to unit metric',
#        color=NULL)
# g2 <- ggplot(stats, aes(x=model, y=rel_ess, color=metric)) +
#   geom_jitter(width=.1, height=0) +
#   theme(axis.text.x = element_text(angle = 90)) +
#   scale_y_log10() + geom_hline(yintercept=1) +
#   labs(x=NULL, y='ESS relative to unit metric',
#        color=NULL)
# g3 <- ggplot(stats, aes(x=model, y=rel_eff, color=metric)) +
#   geom_jitter(width=.1, height=0) +
#   theme(axis.text.x = element_text(angle = 90)) +
#   scale_y_log10() + geom_hline(yintercept=1) +
#   labs(x=NULL, y='Efficiency (minESS/time)\nrelative to unit metric',
#        color=NULL)
#ggsave('plots/perf_rel.png', g3, width=6, height=4)


perf_long <- stats |> pivot_longer(cols=c('rel_ess', 'rel_time', 'rel_eff')) |>
  mutate(pre_metric=metricf(metric),
         name2=factor(name, levels=c('rel_time', 'rel_ess', 'rel_eff'),
                      labels=c('Time', 'min ESS', 'Efficiency (ESS/s)')))
g <- ggplot(perf_long, aes(x=model, y=value, color=pre_metric)) +
  geom_hline(yintercept=1) +
  geom_jitter(width=.1, height=0, alpha=.5) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=.5)) +
  theme(legend.position='top')+
  scale_y_log10() + #coord_flip()+
  facet_wrap('name2', ncol=3, scales='free_y') +
  labs(x=NULL, y='Value relative to Stan defaults',
       color=NULL)
ggsave('plots/perf_mods.png', g, width=7, height=3)

# g <- ggplot(filter(perf_long, name=='rel_eff'),
#             aes(x=model, y=value, color=pre_metric)) +
#   geom_hline(yintercept=1) +
#   geom_jitter(width=.1, height=0, alpha=.5) +
#   theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=.5),
#         legend.position = 'inside', legend.position.inside = c(.9,.7)) +
#   scale_y_log10() + #coord_flip()+
#   facet_wrap('name2', ncol=1, scales='free_y') +
#   labs(x=NULL, y='Value relative to no pre-metric',
#        color=NULL)
# g
# ggsave('plots/case_studies_eff_rel.png', g, width=5, height=3.5)


filter(perf_long, model=='wildf_adapted' & replicate==1) |> select(mean_rel_eff)

stats_long <- stats |> mutate(gr.ratio=log10(gr.ratio)) |>
  select(model, metric, pct.divs, rhat, gr.ratio) |>
    pivot_longer(c('pct.divs', 'gr.ratio', 'rhat')) |>
  mutate(name2=factor(name, levels=c('gr.ratio', 'pct.divs', 'rhat'),
    labels=c('log10 ratio gradient', '% Divergences', 'max Rhat'))) |>
  select(model, metric, value, name2)
g <- ggplot(stats_long, aes(model, y=value, color=metric)) +
  geom_jitter(width=.1) +
  facet_wrap('name2', ncol=1, scales='free_y') +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=.5))
ggsave('plots/case_studies_stats.png', g, width=4, height=6)

grad_timings <- stats |>
  #each metric has it's own gr and gr2, so don't need unit here
  filter(metric!='unit') |>
  group_by(model) |>
  summarize(metric=metric[1],
            gr=1000*mean(gr),
            gr2=1000*mean(gr2),
            gr.ratio=mean(gr.ratio),
            pct.Qall=100*mean(time.Qall/time.total)) |>
arrange(gr.ratio)
grad_timings
write.csv(grad_timings, file='results/grad_timings.csv', row.names=FALSE)


fitswildf2 <- readRDS('results/wildf2_fits.RDS')
fitswildf3 <- readRDS('results/wildf3_fits.RDS')

# wildf2 <- lapply(fitswildf2, \(x) cbind(adapt_delta=.9, metric=x$metric, as.data.frame(x))) |>
#   bind_rows()
# parsall <- g3[[1]]$data |> arrange(value) |> pull(par)
# p <- parsall[1:4]
# pairs_admb(fitswildf2[[1]], pars=p)
# pairs_admb(fitswildf2[[2]], pars=p)
# pairs_admb(fitswildf2[[2]], pars=c('plantSlopeSD', 'slope'))
# pairs_admb(fitswildf3[[1]], pars=p)
# pairs_admb(fitswildf3[[2]], pars=p)
post2 <- lapply(fitswildf2, \(x) cbind(adapt_delta=.9, metric=x$metric, as.data.frame(x))) |>
  bind_rows()
post3 <- lapply(fitswildf3, \(x) cbind(adapt_delta=.99, metric=x$metric, as.data.frame(x))) |>
  bind_rows()
# postwildf <- bind_rows(post2, post3) |> select(all_of(c('metric','adapt_delta',p)))
# ggplot(postwildf, aes(plantSlopeSD, slope, color=metric)) + geom_point(alpha=.5) +
#   facet_wrap('adapt_delta')

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
mtext('Rotated via Q', line=.25)
lines(ellipse::ellipse(0, centre=mle_rotated, scale=c(1,1)), col='blue')
points(mle_rotated[1], mle_rotated[2], pch=17, cex=1, col='blue')
box(col=gray(.4))
legend('bottomright', legend=c('unit', 'sparse'), pch=1, cex=.8,  col=1:2, bty='n')
dev.off()


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
         gr.ratio=gr2/gr,
         mean_rel_eff=mean(rel_eff)) |>
  ungroup() |>
  arrange(desc(mean_rel_eff)) |>
  mutate(model=factor(model, levels=unique(model)))


# very carefully massage output to plot for two scenariors: NUTS and SNUTS
stats_tmp <- bind_rows(mutate(stats_ela, ELA=TRUE), mutate(stats, ELA=FALSE)) |>
  select(model, metric, eff, ELA) |>
  arrange(model, metric, ELA) |> group_by(model) |>
  mutate(eff_rel=eff/mean(eff[metric=='unit' & ELA==FALSE]),
    eff_rel_unit=eff/mean(eff[metric=='unit' & ELA==TRUE]),
         eff_rel_full=eff/mean(eff[metric!='unit' & ELA==FALSE])) |>
  filter(ELA==TRUE) |> ungroup()
# this one is efficiency of ELA:full with both using 'unit'
# metric so comparable to previous studies
x1 <- filter(stats_tmp, metric=='unit') |> mutate(y=eff_rel, type='NUTS')
# now look at ELA:full using SNUTS 'auto' for both
x2 <- filter(stats_tmp, metric!='unit') |> mutate(y=eff_rel_full, type='SNUTS')
xx <- bind_rows(x1,x2) |> ungroup() |> mutate(modelf=modelf(model), metricf=metricf(metric))
g <- ggplot(xx, aes(modelf, y=y, color=metricf, shape=type)) +
  geom_hline(yintercept=1) +
  geom_jitter(width=.1, height=0, alpha=.7) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=.5)) +
  scale_y_log10() + #coord_flip()+
  labs(x=NULL, y='Relative Efficiency (ELA/full NUTS)',
       color=NULL, shape=NULL)+
  theme(legend.position='top')
ggsave('plots/perf_ELA.png', g, width=4.5, height=3)



perf_long |> filter(model=='wham') |> group_by(metric) |> summarize(ess=mean(ess), time=mean(time.total)/(60*60)/3, eff=mean(eff))
