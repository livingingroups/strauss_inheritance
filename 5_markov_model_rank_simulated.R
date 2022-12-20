################################################################################

# Calculate rank trajectories of simulated societies#


#Eli Strauss, Sep 2022#

################################################################################


library(dplyr)
library(here)
library(markovchain)
library(ggplot2)
library(DynaRankR)

rm(list = ls())
source('1_define_functions.R')
load('rank_data_simulated.Rdata')
set.seed(1989)

################################################################################
### Extract dynamics and standardize to decimal ranks

sim.ranks.list <- list()
for(i in unique(output$inheritance)){
  for(re in unique(output$rank.effect)){
    for(r in unique(output$replicate)){
      ranks <- data.frame(filter(output, replicate == r, inheritance == i, rank.effect == re))
      ranks <- ranks %>% 
        rename(period = generation) %>%
        select(id, rank, period, rank.effect, inheritance, replicate) %>% 
        group_by(period) %>%
        mutate(max.rank = max(rank), 
               rank.decile = (rank-1)/(max.rank-1))
      ranks <- get_dynamics(ranks, type = 'rank')
      sim.ranks.list[[length(sim.ranks.list) +1]] <- ranks
    }
  }
}
sim.ranks <- do.call(rbind, sim.ranks.list)

sim.ranks <- sim.ranks %>%
  group_by(inheritance, rank.effect, replicate, id) %>%
  mutate(delta.active.decile = c(NA, (-delta.active[-1])/(max.rank[-length(id)]-1)),
         delta.passive.decile = c(NA, (rank[-1]-1)/(max.rank[-1]-1) -
                                    (rank[-length(id)]-1 - delta.active[-1])/(max.rank[-length(id)]-1)),
         delta.decile = c(NA, (rank[-1]-1)/(max.rank[-1]-1) - 
                            (rank[-length(id)]-1)/(max.rank[-length(id)]-1)))


dev.off()
par(mfrow = c(3,3))
inh <- 'mri_youngest'
re <- 0
for(i in sample(1:max(sim.ranks$replicate), 9)){
  plot_ranks(as.data.frame(filter(sim.ranks, inheritance == inh, rank.effect == re, replicate == i)))
}


#states <- c('alpha', letters[2:11])
states <- letters[1:10]

rank.state.list <- list()
cumulative.changes.list <- list()
for(i in unique(sim.ranks$inheritance)){
  for(re in unique(output$rank.effect)){
    for(r in unique(output$replicate)){
      sim.ranks.treatment <- filter(sim.ranks, inheritance == i, 
                                    rank.effect == re,
                                    replicate == r)
      for(tid in unique(sim.ranks.treatment$id)){
        ## Skip individuals who are only around for 1 year
        if(nrow(sim.ranks.treatment[sim.ranks.treatment$id == tid,]) == 1)
          next
        
        start.rank <- sim.ranks.treatment[sim.ranks.treatment$id == tid,]$rank.decile[1]
        changes <- sim.ranks.treatment[sim.ranks.treatment$id == tid,]$delta.passive.decile[-1]
        rank.traj <- round(start.rank + cumsum(changes), 6)
        
       # rank.state.list[[paste0('re', re)]][[i]][[tid]] <- states[ceiling(rank.traj * 10)]
        rank.state.list[[paste0('re', re)]][[i]][[tid]] <- as.character(cut(rank.traj, 
                                                                            breaks = c(-0.001, 0.1, 0.2, 0.3, 0.4, 
                                                                                       0.5, 0.6, 0.7, 0.8, 0.9, 1.001), 
                         labels = letters[1:10]))
        
        cumulative.changes.list[[paste0('re', re)]][[i]][[tid]] <- 
          data.frame(id = tid, lifetime = length(rank.traj), 
                     starting.rank = start.rank,
                     total.change = cumsum(changes))
      }
    }
  }
}
  
cumulative.changes <- do.call(rbind, cumulative.changes.list)




## Total changes
transitions.total.mri.youngest <- markovchainFit(data = rank.state.list[["re1"]][["mri_youngest"]], method = 'map', confint = T)
transitions.total.mri.youngest.estimate <- transitions.total.mri.youngest$estimate@transitionMatrix

transitions.total.mri.youngest.norank <- markovchainFit(data = rank.state.list[["re0"]][["mri_youngest"]], method = 'map', confint = T)
transitions.total.mri.youngest.norank.estimate <- transitions.total.mri.youngest.norank$estimate@transitionMatrix

transitions.total.mri.oldest <- markovchainFit(data = rank.state.list[["re1"]][["mri_oldest"]], method = 'map', confint = T)
transitions.total.mri.oldest.estimate <- transitions.total.mri.oldest$estimate@transitionMatrix

transitions.total.mocor.weak <- markovchainFit(data = rank.state.list[["re1"]][["parent_offspring_cor0"]], method = 'map', confint = T)
transitions.total.mocor.weak.estimate <- transitions.total.mocor.weak$estimate@transitionMatrix

transitions.total.mocor.weak.norank <- markovchainFit(data = rank.state.list[["re0"]][["parent_offspring_cor0"]], method = 'map', confint = T)
transitions.total.mocor.weak.norank.estimate <- transitions.total.mocor.weak.norank$estimate@transitionMatrix

transitions.total.mocor.medium <- markovchainFit(data = rank.state.list[["re1"]][["parent_offspring_cor2"]], method = 'map', confint = T)
transitions.total.mocor.medium.estimate <- transitions.total.mocor.medium$estimate@transitionMatrix

transitions.total.mocor.strong <- markovchainFit(data = rank.state.list[["re1"]][["parent_offspring_cor4"]], method = 'map', confint = T)
transitions.total.mocor.strong.estimate <- transitions.total.mocor.strong$estimate@transitionMatrix


### Set up monte carlo

mcmc.mri.youngest <- new('markovchain', states = states,
                  transitionMatrix = transitions.total.mri.youngest.estimate)

mcmc.mri.youngest.norank <- new('markovchain', states = states,
                         transitionMatrix = transitions.total.mri.youngest.norank.estimate)

mcmc.mri.oldest <- new('markovchain', states = states,
                         transitionMatrix = transitions.total.mri.oldest.estimate)

mcmc.mocor.weak <- new('markovchain', states = states,
                       transitionMatrix = transitions.total.mocor.weak.estimate)

mcmc.mocor.weak.norank <- new('markovchain', states = states,
                       transitionMatrix = transitions.total.mocor.weak.norank.estimate)

mcmc.mocor.medium <- new('markovchain', states = states,
                           transitionMatrix = transitions.total.mocor.medium.estimate)

mcmc.mocor.strong <- new('markovchain', states = states,
                             transitionMatrix = transitions.total.mocor.strong.estimate)

## Run monte carlo simulation

lifespan <- 100
reps = 1000
for(condition in c('mri.youngest', 'mri.youngest.norank', 'mri.oldest', 'mocor.weak', 'mocor.weak.norank', 'mocor.medium', 'mocor.strong')){
  predicted.lifetime.ranks <- list()
  mrf <- get(paste0('mcmc.', condition))
  for(start.rank in states[c(1, 2, 4, 6, 8, 10)]){
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
  assign(paste0(condition, '.predicted.ranks'), do.call(rbind, predicted.lifetime.ranks))
}


pdf('plots/mri_youngest_sim.pdf', width = 4, height = 4)
mri.youngest <- ggplot(data = mri.youngest.predicted.ranks, aes(x = time, y = rank, col = start, fill = start))+
  geom_ribbon(aes(ymin = rank.lower + 0.25, ymax = rank.upper - 0.25), alpha = 0.3, color = NA)+
  geom_line(size = 1)+
  ylim(10.25,0.75)+
  theme_classic()+
  scale_color_manual(values = viridis::magma(7))+
  scale_fill_manual(values = viridis::magma(7))+
  theme(legend.position = 'none', axis.text = element_blank(), 
        axis.ticks = element_blank())+
  xlab('Time')+
  ylab('Rank')+
  ggtitle('MRI youngest ascendancy')
mri.youngest
dev.off()

pdf('plots/mri_youngest_norank_sim.pdf', width = 4, height = 4)
mri.youngest.norank <- ggplot(data = mri.youngest.norank.predicted.ranks, aes(x = time, y = rank, col = start, fill = start))+
  geom_ribbon(aes(ymin = rank.lower + 0.25, ymax = rank.upper - 0.25), alpha = 0.3, color = NA)+
  geom_line(size = 1)+
  ylim(10.25,0.75)+
  theme_classic()+
  scale_color_manual(values = viridis::magma(7))+
  scale_fill_manual(values = viridis::magma(7))+
  theme(legend.position = 'none', axis.text = element_blank(), 
        axis.ticks = element_blank())+
  xlab('Time')+
  ylab('Rank')+
  ggtitle('MRI youngest ascendancy norank')
mri.youngest.norank
dev.off()

pdf('plots/mri_oldest_sim.pdf', width = 4, height = 4)
mri.oldest <- ggplot(data = mri.oldest.predicted.ranks, aes(x = time, y = rank, col = start, fill = start))+
  geom_ribbon(aes(ymin = rank.lower + 0.25, ymax = rank.upper - 0.25), alpha = 0.3, color = NA)+
  geom_line(size = 1)+
  ylim(10.25,0.75)+
  theme_classic()+
  scale_color_manual(values = viridis::magma(7))+
  scale_fill_manual(values = viridis::magma(7))+
  theme(legend.position = 'none', axis.text = element_blank(), 
        axis.ticks = element_blank())+
  xlab('Time')+
  ylab('Rank')+
  ggtitle('MRI oldest ascendancy')
mri.oldest
dev.off()

pdf('plots/mri_mocor_weak_sim.pdf', width = 4, height = 4)
mocor.weak <- ggplot(data = mocor.weak.predicted.ranks, aes(x = time, y = rank, col = start, fill = start))+
  geom_ribbon(aes(ymin = rank.lower + 0.25, ymax = rank.upper - 0.25), alpha = 0.3, color = NA)+
  geom_line(size = 1)+
  ylim(10.25,0.75)+
  theme_classic()+
  scale_color_manual(values = viridis::magma(7))+
  scale_fill_manual(values = viridis::magma(7))+
  theme(legend.position = 'none', axis.text = element_blank(), 
        axis.ticks = element_blank())+
  xlab('Time')+
  ylab('Rank')+
  ggtitle('Mother-offspring correlation (weak)')
mocor.weak
dev.off()

pdf('plots/mri_mocor_weak_norank_sim.pdf', width = 4, height = 4)
mocor.weak.norank <- ggplot(data = mocor.weak.norank.predicted.ranks, aes(x = time, y = rank, col = start, fill = start))+
  geom_ribbon(aes(ymin = rank.lower + 0.25, ymax = rank.upper - 0.25), alpha = 0.3, color = NA)+
  geom_line(size = 1)+
  ylim(10.25,0.75)+
  theme_classic()+
  scale_color_manual(values = viridis::magma(7))+
  scale_fill_manual(values = viridis::magma(7))+
  theme(legend.position = 'none', axis.text = element_blank(), 
        axis.ticks = element_blank())+
  xlab('Time')+
  ylab('Rank')+
  ggtitle('Mother-offspring correlation (weak), norank')
mocor.weak.norank
dev.off()

pdf('plots/mri_mocor_medium_sim.pdf', width = 4, height = 4)
mocor.medium <- ggplot(data = mocor.medium.predicted.ranks, aes(x = time, y = rank, col = start, fill = start))+
  geom_ribbon(aes(ymin = rank.lower + 0.25, ymax = rank.upper - 0.25), alpha = 0.3, color = NA)+
  geom_line(size = 1)+
  ylim(10.25,0.75)+
  theme_classic()+
  scale_color_manual(values = viridis::magma(7))+
  scale_fill_manual(values = viridis::magma(7))+
  theme(legend.position = 'none', axis.text = element_blank(), 
        axis.ticks = element_blank())+
  xlab('Time')+
  ylab('Rank')+
  ggtitle('Mother-offspring correlation (medium)')
mocor.medium
dev.off()

pdf('plots/mri_mocor_strong_sim.pdf', width = 4, height = 4)
mocor.strong <- ggplot(data = mocor.strong.predicted.ranks, aes(x = time, y = rank, col = start, fill = start))+
  geom_ribbon(aes(ymin = rank.lower + 0.25, ymax = rank.upper - 0.25), alpha = 0.3, color = NA)+
  geom_line(size = 1)+
  ylim(10.25,0.75)+
  theme_classic()+
  scale_color_manual(values = viridis::magma(7))+
  scale_fill_manual(values = viridis::magma(7))+
  theme(legend.position = 'none', axis.text = element_blank(), 
        axis.ticks = element_blank())+
  xlab('Time')+
  ylab('Rank')+
  ggtitle('Mother-offspring correlation (strong)')
mocor.strong
dev.off()







