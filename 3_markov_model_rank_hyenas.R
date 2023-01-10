################################################################################
#
#  Markov models of state transitions through hierarchy in hyenas
#
# 
#
################################################################################

library(markovchain)
library(dplyr)
library(purrr)
library(ggplot2)
library(patchwork)
library(DynaRankR)

load('0_hyena_data.RData')

ranks <- rbind(get_dynamics(ranks[ranks$clan == 'talek',], type = 'rank'),
                        get_dynamics(ranks[ranks$clan == 'serena.n',], type = 'rank'),
                        get_dynamics(ranks[ranks$clan == 'serena.s',], type = 'rank'),
                        get_dynamics(ranks[ranks$clan == 'happy.zebra',], type = 'rank'))


#### Extract rank changes from ranks data

## delta.active.decile = (old.rank + overtakes)/old.group.size - old.rank/old.group.size = overtakes/old.group.size
## delta.passive.decile = (new.rank/new.group.size) - (old.rank+overtakes)/old.group size

ranks <- ranks %>%
  group_by(id) %>%
  mutate(delta.active.decile = c(NA, (-delta.active[-1])/(max.rank[-length(id)]-1)),
         delta.passive.decile = c(NA, (rank[-1]-1)/(max.rank[-1]-1) -
                                    (rank[-length(id)]-1 - delta.active[-1])/(max.rank[-length(id)]-1)),
         delta.decile = c(NA, (rank[-1]-1)/(max.rank[-1]-1) - 
                            (rank[-length(id)]-1)/(max.rank[-length(id)]-1)))


### Define states
states <- letters[1:11]

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
  
  rank.state.list[[tid]] <- states[ceiling(rank.traj * 10) + 1]
  rank.state.passive.list[[tid]] <- states[ceiling(rank.traj.passive * 10)+1]
  rank.state.active.list[[tid]] <- states[ceiling(rank.traj.active * 10)+1]
  
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

transitions.passive <- markovchainFit(data = rank.state.passive.list, method = 'map', confint = T)
transitions.passive.estimate <- transitions.passive$estimate@transitionMatrix

transitions.active <- markovchainFit(data = rank.state.active.list, method = 'map', confint = T)
transitions.active.estimate <- transitions.active$estimate@transitionMatrix


### Set up markov chain monte carlo

total.mcmc <- new('markovchain', states = states,
                  transitionMatrix = transitions.total.estimate)

passive.mcmc <- new('markovchain', states = states,
                           transitionMatrix = transitions.passive.estimate)
active.mcmc <- new('markovchain', states = states,
                            transitionMatrix = transitions.active.estimate)


lifespan <- 25
reps = 1000
for(condition in c('total.', 'passive.', 'active.')){
  predicted.lifetime.ranks <- list()
  mrf <- get(paste0(condition, 'mcmc'))
  for(start.rank in states[c(1, 3, 5, 7, 9, 11)]){
    predicted.lifetime.ranks[[start.rank]] <- data.frame(
      start = rep(start.rank, lifespan+1),
      time = seq(1:(lifespan+1))
    )
    replicates = replicate(reps, markovchainSequence(lifespan, mrf, t0 = start.rank, include.t0 = T))
    replicates = matrix(match(replicates, states), nrow = lifespan+1, ncol = reps)
    predicted.lifetime.ranks[[start.rank]]$rank = apply(replicates, median, MARGIN = 1)
    predicted.lifetime.ranks[[start.rank]]$rank.upper = apply(replicates, quantile, MARGIN = 1, 0.25, type = 1)
    predicted.lifetime.ranks[[start.rank]]$rank.lower = apply(replicates, quantile, MARGIN = 1, 0.75, type = 1)
  }
  assign(paste0(condition, 'predicted.ranks'), do.call(rbind, predicted.lifetime.ranks))
}

pdf('plots/total_change_hyena.pdf', width = 4, height = 4)
total <- ggplot(data = total.predicted.ranks, aes(x = time, y = rank, col = start, fill = start))+
  geom_ribbon(aes(ymin = rank.lower + 0.25, ymax = rank.upper - 0.25), alpha = 0.3, color = NA)+
  geom_line(size = 1)+
  ylim(11.25,0.75)+
  theme_classic()+
  scale_color_manual(values = viridis::magma(7))+
  scale_fill_manual(values = viridis::magma(7))+
  theme(legend.position = 'none', axis.text = element_blank(), 
        axis.ticks = element_blank())+
  xlab('Time')+
  ylab('Rank')+
  ggtitle('Total dynamics')
total
dev.off()

pdf('plots/passive_change_hyena.pdf', width = 4, height = 4)
passive <- ggplot(data = passive.predicted.ranks, aes(x = time, y = rank, col = start, fill = start))+
  geom_ribbon(aes(ymin = rank.lower + 0.25, ymax = rank.upper - 0.25), alpha = 0.3, color = NA)+
  geom_line(size = 1)+
  ylim(11.25,0.75)+
  theme_classic()+
  scale_color_manual(values = viridis::magma(7))+
  scale_fill_manual(values = viridis::magma(7))+
  theme(legend.position = 'none', axis.text = element_blank(), 
        axis.ticks = element_blank())+
  xlab('Time')+
  ylab('Rank')+
  ggtitle('Passive dynamics')
passive
dev.off()

pdf('plots/active_change_hyena.pdf', width = 4, height = 4)
active <- ggplot(data = active.predicted.ranks, aes(x = time, y = rank, col = start, fill = start))+
  geom_ribbon(aes(ymin = rank.lower + 0.25, ymax = rank.upper - 0.25), alpha = 0.3, color = NA)+
  geom_line(size = 1)+
  ylim(11.25,0.75)+
  theme_classic()+
  scale_color_manual(values = viridis::magma(7))+
  scale_fill_manual(values = viridis::magma(7))+
  theme(legend.position = 'none', axis.text = element_blank(), 
        axis.ticks = element_blank())+
  xlab('Time')+
  ylab('Rank')+
  ggtitle('Active dynamics')
active
dev.off()

total + passive + active
