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

load('0_hyena_data.RData')


ars <- ars %>%
  group_by(clan, year) %>%
  mutate(stan_rank = 1 - (rank-1)/(max(rank)-1),
  repro_events = sum(ars),
  clan_size = length(id))

ars$clan_num <- as.numeric(as.factor(ars$clan))

### Number of reproductive events experienced by individual
ars %>%
  group_by(id) %>%
  mutate(repro_scaled = repro_events/clan_size) %>%
  summarize(total_repro_scaled = sum(repro_scaled))%>%
  ggplot(aes(x = total_repro_scaled))+
  geom_density() + 
  theme_classic()
  
ars_dat <- list(ars = ars$ars, 
                stan_rank = ars$stan_rank, 
                clan = ars$clan_num, 
                repro_events = ars$repro_events,
                clan_size = scale(ars$clan_size))



ars.binom <- ulam(
  alist(
    ars ~ dbinom(repro_events, p),
    logit(p) <- a_clan[clan] + B_rank * stan_rank + B_clan_size * clan_size,
    
    a_clan[clan] ~ dnorm(u, sigma),
    u ~ dnorm(0, 1.5),
    sigma ~ dexp(1),
    B_rank ~ dnorm(0, 1.5),
    B_clan_size ~ dnorm(0, 1.5)
  ), data = ars_dat, chains = 3, iter = 3000, cores = 3, log_lik = T)

precis(ars.binom, depth = 2)

model_fit <- precis(ars.binom)
u_fit = model_fit['u','mean']
B_rank_fit = model_fit['B_rank','mean']
B_clan_size_fit = model_fit['B_clan_size','mean']
sr <-  seq(from = 1, to = 0, length.out = 25)

repro_function <- function(stan_rank, clan_size = 0, 
                      B_clan_size = B_clan_size_fit,
                      B_rank = B_rank_fit, 
                      u = u_fit,
                      rank.effect = 1){
  x = u + clan_size * B_clan_size + stan_rank*B_rank * rank.effect
  p <- 1/(1 + exp(-x))
  return(p)
}

save(u_fit, B_rank_fit, B_clan_size_fit, repro_function, file = '2_repro_function.RData')









ps <- link(ars.binom)
ps_sim <- link(ars.binom, data = list(stan_rank = seq(from = 1, to = 0, length.out = 25),
                                      clan = rep(4, 25),
                                      clan_size = rep(0, 25)),
               replace = list(clan = rep(0,25)))

p_repro <- data.frame(stan_rank = ars$stan_rank,
                      p = apply(ps, mean, MARGIN = 2))
p_sim <- data.frame(stan_rank = seq(from = 1, to = 0, length.out = 25),
                    p = apply(ps_sim, mean, MARGIN = 2))

ggplot(p_repro, aes(x = stan_rank, y = p))+
  geom_point()+
  geom_smooth(color = 'black') +
  #geom_line(data = p_sim, color = 'dodgerblue')+
  geom_line(data = data.frame(stan_rank = sr,
                              p = inv_logit(u + 0 * B_clan_size + sr*B_rank)),
            size = 1, color = 'dodgerblue')+
  theme_classic()


# ### Reproductive events by clan size
# ggplot(data  = ars, aes(x = clan_size, y = repro_events))+
#   geom_point() +
#   geom_smooth (method = 'lm') +
#   theme_classic() +
#   geom_vline(aes(xintercept = 25))


# ars.mod <- ulam(
#   alist(
#     ars ~ dpois(l),
#     log(l) <- a_clan[clan] + B_rank * stan_rank,
#     
#     a_clan[clan] ~ dnorm(u, sigma),
#     u ~ dnorm(0, 1.5),
#     sigma ~ dexp(1),
#     B_rank ~ dnorm(0, 1.5)
#   ), data = ars_dat, chains = 4, iter = 1000)
# 

# summary(ars.mod)
# ars$predicted <- apply(X = sim(ars.mod), mean, MARGIN = 2)
# 
# for(i in 1:30){
#   offspring <- rpois(length(sr), lambda = exp(-2.36 + sr * 1.2))
# cat(sum(offspring), ', ')
# }
# 
# ggplot(ars, aes(x = stan_rank, y = ars))+ 
#   geom_smooth(color=  'black') +
#   geom_smooth(aes(y = predicted), color = 'dodgerblue', lty = 2, se = F)+
#   theme_classic()+
#   coord_cartesian(ylim = c(0, 0.8))
# 
# rpois(length(sr), lambda = exp(-2.36 + sr * 1.2)/2)
# mean(exp(-2.36 + sr * 1.2))
