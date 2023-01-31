################################################################################
#
#  Markov models of state transitions through hierarchy in hyenas
#
# 
#
################################################################################

library(markovchain)
library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)
library(patchwork)
library(DynaRankR)
library(ggridges)
library(expm)
source('1_define_functions.R')
set.seed(1989)

plot.dir <- '~/Dropbox/Documents/Research/Full_projects/2023 Inheritancy_mobility/plots/'

load('0_hyena_data.RData')

ranks <- rbind(get_dynamics(ranks[ranks$clan == 'talek',], type = 'rank'),
                        get_dynamics(ranks[ranks$clan == 'serena.n',], type = 'rank'),
                        get_dynamics(ranks[ranks$clan == 'serena.s',], type = 'rank'),
                        get_dynamics(ranks[ranks$clan == 'happy.zebra',], type = 'rank'))


#### Extract rank changes from ranks data

ranks <- ranks %>%
  group_by(id) %>%
  mutate(delta.active.decile = c(NA, (-delta.active[-1])/(max.rank[-length(id)]-1)),
         delta.passive.decile = c(NA, (rank[-1]-1)/(max.rank[-1]-1) -
                                    (rank[-length(id)]-1 - delta.active[-1])/(max.rank[-length(id)]-1)),
         delta.decile = c(NA, (rank[-1]-1)/(max.rank[-1]-1) - 
                            (rank[-length(id)]-1)/(max.rank[-length(id)]-1)))


ranks.full.life <- ranks %>% 
  group_by(id) %>%
  mutate(first.year = first(period),
         last.year = last(period)) %>% 
  ungroup() %>%
  group_by(clan) %>%
  mutate(first.year.clan = first(period), 
         last.year.clan = last(period)) %>% 
  ungroup () %>%
  filter(first.year > first.year.clan,
         last.year < last.year.clan)
length(unique(ranks.full.life$id))

dynamics <- ranks.full.life %>% 
  group_by(id) %>%
  summarize(total.passive = abs(sum(delta.passive.decile, na.rm = T)),
            total.active = abs(sum(delta.active.decile, na.rm = T))) %>%
  ungroup() %>%
  select(total.passive, total.active) %>%
  gather(key = 'dynamics', value = 'cumulative.change')

dynamics$dynamics <- factor(dynamics$dynamics, levels = c('total.active', 'total.passive'), labels = c('Active\ndynamics', 'Passive\ndynamics'))

dynamics.test <- t.test(abs(filter(dynamics, dynamics == 'Passive\ndynamics')$cumulative.change),
                        abs(filter(dynamics, dynamics == 'Active\ndynamics')$cumulative.change),
                        paired = T)


## Mean change due to active dynamics
mean(filter(dynamics, dynamics == 'Active\ndynamics')$cumulative.change) * 100
## Mean change due to passive dynamics
mean(filter(dynamics, dynamics == 'Passive\ndynamics')$cumulative.change) * 100


#pdf(file = 'dynamics.pdf', width = 4, height = 2.5)
dynamics.plot <- ggplot(dynamics)+
  geom_density_ridges(aes(x = cumulative.change, y = dynamics, fill = dynamics, height = ..density..), alpha = 0.8,
                      stat = 'density', trim = T) +
  theme_classic()+
  theme(axis.title.y = element_blank(), legend.position = 'none', 
        axis.text.y = element_text(hjust = 0.5, angle = 90))+
  xlab('Magnitude of net displacement from starting rank')+
  scale_fill_manual(values = c('black', 'gold2'))+
  scale_y_discrete(expand = c(0.1,0))+
  labs(tag = 'A')
#dev.off()



### Define states
states <- letters[1:10]

rank.state.list <- list()
rank.state.passive.list <- list()
rank.state.active.list <- list()
cumulative.changes.list <- list()
for(tid in unique(ranks$id)){
  ## Skip individuals who are only around for 1 year
  if(nrow(ranks[ranks$id == tid,]) == 1)
    next
  start.rank <- ranks[ranks$id == tid,]$rank.decile[1]
  changes.passive <- ranks[ranks$id == tid,]$delta.passive.decile[-1]
  changes.active <- ranks[ranks$id == tid,]$delta.active.decile[-1]
  
  rank.traj.passive <- start.rank + cumsum(changes.passive)
  rank.traj.active <- start.rank + cumsum(changes.active)
  rank.traj <- round(start.rank + cumsum(changes.passive) + cumsum(changes.active), 6)
  
  rank.state.list[[tid]] <-  cut(rank.traj, 
                                 breaks = c(-0.001, 0.1, 0.2, 0.3, 0.4,
                                            0.5, 0.6, 0.7, 0.8, 0.9, 1.001), 
                                 labels = states)
  rank.state.passive.list[[tid]] <- cut(rank.traj.passive, 
                                        breaks = c(-0.001, 0.1, 0.2, 0.3, 0.4,
                                                   0.5, 0.6, 0.7, 0.8, 0.9, 1.001), 
                                        labels = states)
  rank.state.active.list[[tid]] <- cut(rank.traj.active, 
                                       breaks = c(-0.001, 0.1, 0.2, 0.3, 0.4,
                                                  0.5, 0.6, 0.7, 0.8, 0.9, 1.001), 
                                       labels = states)
  
  cumulative.changes.list[[tid]] <- data.frame(id = tid, lifetime = length(rank.traj), 
                                               starting.rank = start.rank,
                                               total.change = cumsum(changes.passive) + cumsum(changes.active),
                                               passive.change = cumsum(changes.passive),
                                               active.change = cumsum(changes.active))
}

cumulative.changes <- do.call(rbind, cumulative.changes.list)
  
#### Fit markov models

## Total changes
transitions.total <- markovchainFit(data = rank.state.list, method = 'map', confint = T)
transitions.total.estimate <- transitions.total$estimate@transitionMatrix
limit.hyena.total <- data.frame(probs = find_limit(transitions.total.estimate),
                                states)

transitions.passive <- markovchainFit(data = rank.state.passive.list, method = 'map', confint = T)
transitions.passive.estimate <- transitions.passive$estimate@transitionMatrix
limit.hyena.passive <- data.frame(probs = find_limit(transitions.passive.estimate),
                                states)

transitions.active <- markovchainFit(data = rank.state.active.list, method = 'map', confint = T)
transitions.active.estimate <- transitions.active$estimate@transitionMatrix
limit.hyena.active <- data.frame(probs = find_limit(transitions.active.estimate),
                                  states)


### Set up markov chain

total.mc <- new('markovchain', states = states,
                  transitionMatrix = transitions.total.estimate)

passive.mc <- new('markovchain', states = states,
                           transitionMatrix = transitions.passive.estimate)
active.mc <- new('markovchain', states = states,
                            transitionMatrix = transitions.active.estimate)

active <- active.mc@transitionMatrix %^% 10000


lifespan <- 10
reps = 10000
for(condition in c('total.', 'passive.', 'active.')){
  predicted.lifetime.ranks <- list()
  mrf <- get(paste0(condition, 'mc'))
  for(start.rank in states){
    predicted.lifetime.ranks[[start.rank]] <- data.frame(
      start = rep(start.rank, lifespan+1),
      time = seq(1:(lifespan+1))
    )
    replicates = replicate(reps, markovchainSequence(lifespan, mrf, t0 = start.rank, include.t0 = T))
    replicates = matrix(match(replicates, states), nrow = lifespan+1, ncol = reps)
    predicted.lifetime.ranks[[start.rank]]$rank = apply(replicates, median, MARGIN = 1)
    predicted.lifetime.ranks[[start.rank]]$rank.upper = apply(replicates, quantile, MARGIN = 1, 0.25, type = 1)
    predicted.lifetime.ranks[[start.rank]]$rank.lower = apply(replicates, quantile, MARGIN = 1, 0.75, type = 1)
    predicted.lifetime.ranks[[start.rank]]$rank.mean = apply(replicates, mean, MARGIN = 1)
    predicted.lifetime.ranks[[start.rank]]$rank.mean.sd = apply(replicates, sd, MARGIN = 1)
  }
  assign(paste0(condition, 'predicted.ranks'), do.call(rbind, predicted.lifetime.ranks))
}

#pdf('plots/total_change_hyena.pdf', width = 4, height = 4)
hyena.total <- ggplot(data = total.predicted.ranks, aes(x = time, y = rank.mean, col = start, fill = start))+
  geom_ribbon(aes(ymin = rank.mean + rank.mean.sd, ymax = rank.mean - rank.mean.sd), alpha = 0.3, color = NA)+
  geom_line(size = 1)+
  theme_classic()+
  scale_color_manual(values = colorRampPalette(colors = c('black', 'gold'))(10))+
  scale_fill_manual(values = colorRampPalette(colors = c('black', 'gold'))(10))+
  theme(legend.position = 'none', axis.text.x = element_blank(), 
        axis.ticks = element_blank(), axis.text.y = element_text(angle = 90))+
  scale_y_continuous(breaks = c(10,1), labels = c('low', 'high'))+
  xlab('Time')+
  ylab('Rank')+
  coord_cartesian(ylim = c(10,1))+
  labs(tag = 'B')+
  ggtitle('Total dynamics')

hyena.total.limit <- ggplot(data = limit.hyena.total, 
                             aes(y = states[10:1], x = probs, fill = states[10:1])) +
  geom_bar(stat= 'identity') +
  scale_fill_manual(values = colorRampPalette(c('gold', 'black'))(10))+
  theme_void()+
  theme(legend.position = 'none',
        axis.line.x = element_line(), axis.text.x = element_text(),
        axis.title.x = element_text())+
  scale_x_continuous(breaks = c(0, 0.5),
                     limits = c(0, 0.5), name = 'Steady state')
#total
#dev.off()

#pdf('plots/passive_change_hyena.pdf', width = 4, height = 4)
hyena.passive <- ggplot(data = passive.predicted.ranks, aes(x = time, y = rank.mean, col = start, fill = start))+
  geom_ribbon(aes(ymin = rank.mean + rank.mean.sd, ymax = rank.mean - rank.mean.sd), alpha = 0.3, color = NA)+
  geom_line(size = 1)+
  theme_classic()+
  scale_color_manual(values = colorRampPalette(colors = c('black', 'gold'))(10))+
  scale_fill_manual(values = colorRampPalette(colors = c('black', 'gold'))(10))+
  theme(legend.position = 'none', axis.text.x = element_blank(), 
        axis.ticks = element_blank(), axis.text.y = element_text(angle = 90))+
  scale_y_continuous(breaks = c(10,1), labels = c('low', 'high'))+
  xlab('Time')+
  ylab('Rank')+
  coord_cartesian(ylim = c(10,1)) + 
  labs(tag = 'C')+
  ggtitle('Passive dynamics')


hyena.passive.limit <- ggplot(data = limit.hyena.passive, 
                            aes(y = states[10:1], x = probs, fill = states[10:1])) +
  geom_bar(stat= 'identity') +
  scale_fill_manual(values = colorRampPalette(c('gold', 'black'))(10))+
  theme_void()+
  theme(legend.position = 'none',
        axis.line.x = element_line(), axis.text.x = element_text(),
        axis.title.x = element_text())+
  scale_x_continuous(breaks = c(0, 0.5),
                     limits = c(0, 0.51), name = 'Steady state')
#passive
#dev.off()

#pdf('plots/active_change_hyena.pdf', width = 4, height = 4)
hyena.active <- ggplot(data = active.predicted.ranks, aes(x = time, y = rank.mean, col = start, fill = start))+
  geom_ribbon(aes(ymin = rank.mean + rank.mean.sd, ymax = rank.mean - rank.mean.sd), alpha = 0.3, color = NA)+
  geom_line(size = 1)+
  theme_classic()+
  scale_color_manual(values = colorRampPalette(colors = c('black', 'gold'))(10))+
  scale_fill_manual(values = colorRampPalette(colors = c('black', 'gold'))(10))+
  theme(legend.position = 'none', axis.text.x = element_blank(), 
        axis.ticks = element_blank(), axis.text.y = element_text(angle = 90))+
  scale_y_continuous(breaks = c(10,1), labels = c('low', 'high'))+
  xlab('Time')+
  ylab('Rank')+
  coord_cartesian(ylim = c(10,1))+
  labs(tag = 'D')+
  ggtitle('Active dynamics')

hyena.active.limit <- ggplot(data = limit.hyena.active, 
                              aes(y = states[10:1], x = probs, fill = states[10:1])) +
  geom_bar(stat= 'identity') +
  scale_fill_manual(values = colorRampPalette(c('gold', 'black'))(10))+
  theme_void()+
  theme(legend.position = 'none',
        axis.line.x = element_line(), axis.text.x = element_text(),
        axis.title.x = element_text())+
  scale_x_continuous(breaks = c(0, 0.5),
                     limits = c(0, 0.5), name = 'Steady state')
#active
#dev.off()

svg(paste0(plot.dir, '/hyena_dynamics.svg'), height = 4.8, width = 9)
dynamics.plot + hyena.total + hyena.total.limit + 
  hyena.passive + hyena.passive.limit + 
  hyena.active + hyena.active.limit + 
  plot_layout(design = 
                'AAAAABBBCC
              DDDEEFFFGG')
dev.off()

# Median hierarchy position decile (10 = highest 10% s, 1 = lowest 10%) expected rank in the long run
seq(from = 10, to = 1)[median(rep(1:length(states), times = round(limit.hyena.total$probs * 1000000)))]
seq(from = 10, to = 1)[median(rep(1:length(states), times = round(limit.hyena.passive$probs * 1000000)))]
seq(from = 10, to = 1)[median(rep(1:length(states), times = round(limit.hyena.active$probs * 1000000)))]





