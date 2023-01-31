################################################################################
#
#                    Extract data from MHP database 
#
#                        Eli Strauss Dec 2022
################################################################################


options(stringsAsFactors = FALSE)
library(dplyr)
library(here)
library(rethinking)
library(ggplot2)

plot.dir <- '~/Dropbox/Documents/Research/Full_projects/2023 Inheritancy_mobility/plots/'
load('0_hyena_data.RData')

ars <- ars %>%
  group_by(clan, year) %>%
  mutate(stan_rank = 1 - (rank-1)/(max(rank)-1),
  repro_events = sum(ars),
  clan_size = length(id))


ars$clan_num <- as.numeric(as.factor(ars$clan))

### Number of reproductive events experienced for guiding markov chain simulation
svg(filename = paste0(plot.dir, '/repro_events_per_year.svg'), width = 6, height = 4)
ars %>%
  group_by(id) %>%
  mutate(repro_scaled = repro_events/clan_size) %>%
  summarize(total_repro_scaled = sum(repro_scaled),
            repro_longevity = length(repro_events))%>%
  ggplot(aes(x = repro_longevity, y = total_repro_scaled))+
  geom_point() + 
  geom_smooth(color = 'black', method = 'loess')+
  theme_classic()+
  ylab('Recruitment events experienced\n(scaled by clan size)')+
  xlab('Number of years observed reproducing')+
  geom_segment(x = 10, xend = 10, y = -2, yend = 4, lty = 2)+
  geom_segment(x = -2, xend = 10, y = 4, yend = 4, lty = 2)
dev.off()
  
ars_dat <- list(ars = ars$ars, 
                stan_rank = ars$stan_rank, 
                clan = ars$clan_num, 
                repro_events = ars$repro_events)


ars.binom <- ulam(
  alist(
    ars ~ dbinom(repro_events, p),
    logit(p) <- u + a_clan[clan] + B_rank * stan_rank,
    
    a_clan[clan] ~ dnorm(0, sigma),
    u ~ dnorm(0, 1.5),
    sigma ~ dexp(1),
    B_rank ~ dnorm(0, 1.5)
  ), data = ars_dat, chains = 3, iter = 6000, cores = 2, log_lik = T)

precis(ars.binom, depth = 2)

model_fit <- precis(ars.binom)
u_fit = model_fit['u','mean']
B_rank_fit = model_fit['B_rank','mean']
sr <-  seq(from = 1, to = 0, length.out = 30)

repro_function <- function(stan_rank, 
                      B_rank = B_rank_fit, 
                      u = u_fit,
                      rank.effect = 1){
  x = u + stan_rank*B_rank * rank.effect
  p <- 1/(1 + exp(-x))
  return(p)
}

save(u_fit, B_rank_fit, repro_function, file = '2_repro_function.RData')

samps <- extract.samples(ars.binom,n = 300)

rf_post <- function(u, B_rank){
  x = u + seq(from = 0, to = 1, length.out = 30)*B_rank
  p <- 1/(1 + exp(-x))
  return(p)
}

p.preds.list <- list()
for(i in 1:length(samps[[2]])){
  p.preds.list[[length(p.preds.list)+1]] <- 
    data.frame(probs = rf_post(samps$u[i], samps$B_rank[i]),
               sample = i,
               category = 'post',
               rank = 1:30)
}
p.preds <- do.call(rbind, p.preds.list)

p.preds <- rbind(p.preds,
                 data.frame(probs = repro_function(seq(from = 0, to = 1, length.out = 30)),
                            sample = -1,
                            category = 'mean',
                            rank = 1:30))

svg(paste0(plot.dir, 'prob_reproduce.svg'), width = 6, height = 4)
ggplot(data = filter(p.preds, category == 'post'), aes(x = rank, y = probs, group = sample))+
  geom_line(size = 0.1, alpha = 0.2)+
  geom_line(data = filter(p.preds, category == 'mean'), size = 2)+
  theme_classic() +
  xlab('Rank')+
  ylab('Probablity of being selected to reproduce')+
  scale_x_continuous(breaks = c(1,30), labels = c('lowest', 'highest'))
dev.off()





