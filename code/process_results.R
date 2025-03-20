

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
  filter(!model %in% c('glmmTMB', 'VAR', 'spde')) |>
  group_by(model) |>
  mutate(rel_eff=eff/mean(eff[metric=='unit']),
         rel_ess=ess/mean(ess[metric=='unit']),
         rel_time=time.total/mean(time.total[metric=='unit']),
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

stats_long <- stats |> pivot_longer(c('pct.divs', 'gr', 'gr2', 'rhat')) |>
  select(model, metric, value, name)
ggplot(stats_long, aes(model, y=value, color=metric)) + geom_point() +
  facet_wrap('name', ncol=1, scales='free_y') +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=.5)) +
  scale_y_log10()


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
