

setwd(here::here())
source("code/startup.R")


stats_sim <- bind_rows(readRDS('results/spde_stats.RDS'),
                       readRDS('results/VAR_stats.RDS'),
                       readRDS("results/glmmTMB_stats.RDS")) |>
  group_by(model, metric, nrepars) |> #summarize(n=n()) |> print(n=Inf)
  filter(! (model=='VAR' & nrepars <205)) |> mutate(metric2=metricf2(metric)) |>
  ungroup()

perf_sim <- stats_sim |>
  mutate(eff_avg=mean(eff),
         ess_avg=mean(ess),
         time.total=time.warmup+time.sampling,
         time_avg=mean(time.total)) |> ungroup() |>
  pivot_longer(cols=c('eff_avg', 'ess_avg', 'time_avg')) |>
  group_by(model, nrepars) |>
  mutate(relval=value/value[as.character(metric)=='unit']) |> ungroup() |>
  mutate(name2=factor(name, levels=c('time_avg', 'ess_avg', 'eff_avg'),
                      labels=c('Total time', 'min(ESS)', 'Efficiency (ESS/s)')))

#, rel_eff=eff/eff_avg[as.character(metric)=='unit'])
g <- ggplot(filter(perf_sim, metric!='unit'),
            aes(x=nrepars, y=relval, color=metric)) +
  facet_grid(name2~model, scales='free') + geom_point() + geom_line()+
  scale_x_log10() + scale_y_log10() +
  labs(x='Number of random effects',
       y="Mean value relative to metric 'unit'",
       color=NULL) +
  geom_hline(yintercept = 1)
g
ggsave('plots/perf_sim_rel.png', g, width=7, height=5)

g <- ggplot(perf_sim, aes(x=nrepars, y=value, color=metric)) +
  facet_grid(name2~model, scales='free') + geom_point() + geom_line()+
  scale_x_log10() + scale_y_log10() +
  labs(x='Number of random effects',
       y=NULL)
ggsave('plots/perf_sim.png', g, width=7, height=4)


stats_sim |> str()


message("Reading in large case study fit files...")
results <- list.files('results',  pattern='_fits.RDS', full.names = TRUE) |> lapply(readRDS)
mods <- sapply(results, \(x) x[[1]]$model)
ela <- grepl('_ELA', mods)
results_ela <- results[ela]
results <- results[!ela]

message("Creating reference posterior list...")
ref_posts_list <- lapply(results, \(mod) {
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
  filter(!model %in% c('glmmTMB', 'VAR', 'spde', 'wildf2', 'wildf3')) |>
  group_by(model) |>
  mutate(rel_eff=eff/mean(eff[metric=='unit']),
         rel_ess=ess/mean(ess[metric=='unit']),
         rel_time=time.total/mean(time.total[metric=='unit']),
         gr.ratio=gr2/gr,
         mean_rel_eff=mean(rel_eff)) |>
  ungroup() |>
  arrange(desc(mean_rel_eff)) |>
  mutate(model=factor(model, levels=unique(model)),
         metric2=metricf2(metric))
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
  mutate(name2=factor(name, levels=c('rel_time', 'rel_ess', 'rel_eff'),
                      labels=c('Time', 'min ESS', 'Efficiency (ESS/s)')))
g <- ggplot(perf_long, aes(x=model, y=value, color=metric)) +
  geom_hline(yintercept=1) +
  geom_jitter(width=.1, height=0, alpha=.5) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=.5)) +
  scale_y_log10() + #coord_flip()+
  facet_wrap('name2', ncol=1, scales='free_y') +
  labs(x=NULL, y='Value relative to unit metric',
       color=NULL)
ggsave('plots/case_studies_perf_rel.png', g, width=4, height=6)

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
  filter(!model %in% c('glmmTMB', 'VAR', 'spde', 'wildf2', 'wildf3')) |>
  group_by(model) |>
  mutate(rel_eff=eff/mean(eff[metric=='unit']),
         rel_ess=ess/mean(ess[metric=='unit']),
         rel_time=time.total/mean(time.total[metric=='unit']),
         gr.ratio=gr2/gr,
         mean_rel_eff=mean(rel_eff)) |>
  ungroup() |>
  arrange(desc(mean_rel_eff)) |>
  mutate(model=factor(model, levels=unique(model)),
         metric2=metricf2(metric))


stats_tmp <- bind_rows(mutate(stats_ela, ELA=TRUE), mutate(stats, ELA=FALSE)) |>
  select(model, metric, eff, ELA) |>
  arrange(model, metric, ELA) |> group_by(model) |>
  mutate(eff_rel_unit=eff/mean(eff[metric=='unit' & ELA==TRUE]),
         eff_rel_full=eff/mean(eff[metric!='unit' & ELA==FALSE]))
g <- ggplot(filter(stats_tmp, metric!='unit' & ELA==TRUE), aes(model, y=eff_rel_full, color=metric)) +
  geom_hline(yintercept=1) +
  geom_jitter(width=.1, height=0, alpha=.5) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=.5)) +
  scale_y_log10() + #coord_flip()+
    labs(x=NULL, y='Efficiency relative to full SNUTS',
       color=NULL)
ggsave('plots/case_studies_perf_ELA.png', g, width=4, height=3)
