# Script to run Pathfinder through StanEstimators and compare to
# approximate samples derived by samples from Q around the mode.

# reference posteriors come from previously fitted models (auto
# metric, see code/process_results.R)
posts <- readRDS('results/ref_posts_list.RDS')

# which replicates to run, where each is a single optimization
# path for Q and 4 for multiPathfinder in serial

# still need to figure out how to do the timing more fairly!!

reps <- 1:20
source("code/load_rtmb_objects.R")
# check sensitivity to initial values for a simple model like
# schools
w.inits <-
  lapply(c('post', 'unif-1', 'unif-2', 'mode', '0'),
       \(init) cbind(init=init, get_wasserstein(reps, obj.schools, posts[['schools']], 'schools', init=init))) |> bind_rows() |>
  group_by(init) |>
  mutate(w1d_rel=w1d/mean(w1d[type=='Q'])) |>
  mutate(time_rel=time.total/mean(time.total[type=='Q']))
saveRDS(w.inits, file='results/w1d_inits.RDS')

w.inits <- readRDS('results/w1d_inits.RDS')
w.means <- group_by(w.inits, init, type) |>
  summarize(w1d_rel=mean(w1d_rel), time_rel=mean(time_rel))
g1 <- ggplot(w.inits, aes(x=init, y=w1d, color=type)) +
  geom_jitter(width=.1, alpha=.5) +
  coord_flip() + labs(x=NULL, y='1D Wasserstein distance')+
  scale_y_log10()
g2 <- ggplot(w.inits, aes(x=init, y=w1d_rel, color=type)) +
  geom_jitter(width=.1, alpha=.5, pch=1, size=1) +
  coord_flip()+
  labs(x=NULL, y='1D Wasserstein distance relative to Q')+
  scale_y_log10() + geom_point(data=w.means, pch=17, size=2)
g3 <- ggplot(w.inits, aes(x=init, y=time_rel, color=type)) +
  geom_jitter(width=.1, alpha=.5, pch=1, size=1) +
  coord_flip() + scale_y_log10() +
  labs(x='Model', y='Time relative to Q') +
  theme(legend.position='none') +geom_point(data=w.means, pch=17, size=2)
g <- plot_grid(g3,  g2, nrow=1, labels = c('(a)', '(b)', '(c)')[-3],
               label_size = 10,
               rel_widths=c(1, 1, 1.35)[-2])
ggsave('plots/pathfinder_inits.png', width=7, height=2, dpi=200)



w.schools <- get_wasserstein(reps, obj.schools, posts[['schools']], 'schools', init='post')
w.radon <- get_wasserstein(reps, obj.radon, posts[['radon']], 'radon', init='post')
w.diamonds <- get_wasserstein(reps, obj.diamonds, posts[['diamonds']], 'diamonds', init='post')
w.causal <- get_wasserstein(reps, obj.causal, posts[['causal']], 'causal', init='default')
w.kilpisjarvi <- get_wasserstein(reps, obj.kilpisjarvi, posts[['kilpisjarvi']], 'kilpisjarvi', init='post')


source("code/load_tmb_objects.R")
w.dlm <- get_wasserstein(reps, obj.dlm, posts[['dlm']], 'dlm', init='post')
w.wildf <- get_wasserstein(reps, obj.wildf, posts[['wildf']], 'wildf', init='post')
w.swallows <- get_wasserstein(reps, obj.swallows, posts[['swallows']], 'swallows', init='post')
w.simple <- get_wasserstein(reps, obj.simple, posts[['simple']], 'simple', init='post')
# this one fails for rep=5 and stalls out the R session so start it from the mode
w.sdmTMB <- get_wasserstein(reps, obj.sdmTMB, posts[['sdmTMB']], 'sdmTMB', init='default')
w.pollock <- get_wasserstein(reps, obj.pollock, posts[['pollock']], 'pollock', init='post')



# quick plots
w.all <- rbind(w.schools, w.radon, w.diamonds, w.causal, w.kilpisjarvi,
               w.dlm, w.wildf, w.swallows, w.simple, w.sdmTMB, w.pollock) |>
  group_by(model) |>
  mutate(w1d_rel=w1d/mean(w1d[type=='Q'])) |>
  mutate(time_rel=time.total/mean(time.total[type=='Q']))

saveRDS(w.all, file='results/w1d_all.RDS')

w.all <- readRDS('results/w1d_all.RDS')
w.means <- group_by(w.all, model, type) |>
  summarize(w1d_rel=mean(w1d_rel), time_rel=mean(time_rel))
g1 <- ggplot(w.all, aes(x=model, y=w1d, color=type)) +
  geom_jitter(width=.1, alpha=.5) +
  coord_flip() + labs(x=NULL, y='1D Wasserstein distance')+
  scale_y_log10()
g2 <- ggplot(w.all, aes(x=model, y=w1d_rel, color=type)) +
  geom_jitter(width=.1, alpha=.5, pch=1, size=1) +
  coord_flip()+
  labs(x=NULL, y='1D Wasserstein distance relative to Q')+
  scale_y_log10() + geom_point(data=w.means, pch=17, size=2)
g3 <- ggplot(w.all, aes(x=model, y=time_rel, color=type)) +
  geom_jitter(width=.1, alpha=.5, pch=1, size=1) +
  coord_flip() + scale_y_log10() +
  labs(x='Model', y='Time relative to Q') +
  theme(legend.position='none') +geom_point(data=w.means, pch=17, size=2)
g <- plot_grid(g3,  g2, nrow=1, labels = c('(a)', '(b)', '(c)')[-3],
               label_size = 10,
               rel_widths=c(1, 1, 1.35)[-2])
ggsave('plots/pathfinder_all.png', width=7, height=2, dpi=200)



