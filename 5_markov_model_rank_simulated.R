################################################################################

# Calculate rank trajectories of simulated societies#


#Eli Strauss, Sep 2022#

################################################################################


library(dplyr)
library(here)
library(markovchain)
library(ggplot2)
# library(DynaRankR)
library(patchwork )
library(gridExtra)

rm(list = ls())
source('1_define_functions.R')
load('rank_data_simulated.Rdata')
set.seed(1989)

################################################################################
### Extract dynamics and standardize to decimal ranks

sim.ranks <- output %>% 
  rename(period = generation) %>%
  select(id, rank, period, rank.effect, inheritance, replicate)

states <- letters[1:11]

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
        rank.seq <- sim.ranks.treatment[sim.ranks.treatment$id == tid,]$rank
        state.seq <- as.character(cut((rank.seq - 1)/(group.size - 1), 
                                      breaks = c(-0.001, 0.001, 0.1, 0.2, 0.3, 0.4, 
                                                 0.5, 0.6, 0.7, 0.8, 0.9, 1.001), 
                                      labels = letters[1:11]))
        rank.state.list[[paste0('re', re)]][[i]][[length(rank.state.list[[paste0('re', re)]][[i]]) + 1]] <- state.seq

      }
    }
  }
}


## Total changes
transitions.total.mri.youngest <- markovchainFit(data = rank.state.list[["re1"]][["mri_youngest"]], method = 'map', confint = F)
transitions.total.mri.youngest.estimate <- transitions.total.mri.youngest$estimate@transitionMatrix

transitions.total.mri.youngest.norank <- markovchainFit(data = rank.state.list[["re0"]][["mri_youngest"]], method = 'map', confint = T)
transitions.total.mri.youngest.norank.estimate <- transitions.total.mri.youngest.norank$estimate@transitionMatrix

transitions.total.mri.oldest <- markovchainFit(data = rank.state.list[["re1"]][["mri_oldest"]], method = 'map', confint = T)
transitions.total.mri.oldest.estimate <- transitions.total.mri.oldest$estimate@transitionMatrix

transitions.total.mocor.weak <- markovchainFit(data = rank.state.list[["re1"]][["parent_offspring_cor0"]], method = 'map', confint = T)
transitions.total.mocor.weak.estimate <- transitions.total.mocor.weak$estimate@transitionMatrix

transitions.total.mocor.weak.norank <- markovchainFit(data = rank.state.list[["re0"]][["parent_offspring_cor0"]], method = 'map', confint = T)
transitions.total.mocor.weak.norank.estimate <- transitions.total.mocor.weak.norank$estimate@transitionMatrix

transitions.total.mocor.medium <- markovchainFit(data = rank.state.list[["re1"]][["parent_offspring_cor4"]], method = 'map', confint = T)
transitions.total.mocor.medium.estimate <- transitions.total.mocor.medium$estimate@transitionMatrix

transitions.total.mocor.strong <- markovchainFit(data = rank.state.list[["re1"]][["parent_offspring_cor2"]], method = 'map', confint = T)
transitions.total.mocor.strong.estimate <- transitions.total.mocor.strong$estimate@transitionMatrix


### Set up markov chain

mc.mri.youngest <- new('markovchain', states = states,
                  transitionMatrix = transitions.total.mri.youngest.estimate)

mc.mri.youngest.norank <- new('markovchain', states = states,
                         transitionMatrix = transitions.total.mri.youngest.norank.estimate)

mc.mri.oldest <- new('markovchain', states = states,
                         transitionMatrix = transitions.total.mri.oldest.estimate)

mc.mocor.weak <- new('markovchain', states = states,
                       transitionMatrix = transitions.total.mocor.weak.estimate)

mc.mocor.weak.norank <- new('markovchain',
                       transitionMatrix = transitions.total.mocor.weak.norank.estimate)

mc.mocor.medium <- new('markovchain', states = states,
                           transitionMatrix = transitions.total.mocor.medium.estimate)

mc.mocor.strong <- new('markovchain', states = states,
                             transitionMatrix = transitions.total.mocor.strong.estimate)

## Run markov chain simulation

repro_events <- 120 ## Based on the expected number of reproductive events experience by a moderately long-lived female
reps = 10000
for(condition in c('mri.youngest', 'mri.youngest.norank', 'mri.oldest', 'mocor.weak', 'mocor.weak.norank', 'mocor.medium', 'mocor.strong')){
  predicted.lifetime.ranks <- list()
  mrf <- get(paste0('mc.', condition))
  for(start.rank in states[c(1,3,5,7,9,11)]){
    predicted.lifetime.ranks[[start.rank]] <- data.frame(
      start = rep(start.rank, repro_events+1),
      time = seq(1:(repro_events+1))
    )
    replicates = replicate(reps, markovchainSequence(repro_events, mrf, t0 = start.rank, include.t0 = T))
    replicates = matrix(match(replicates, states), nrow = repro_events+1, ncol = reps)
    predicted.lifetime.ranks[[start.rank]]$rank = apply(replicates, median, MARGIN = 1)
    predicted.lifetime.ranks[[start.rank]]$rank.upper = apply(replicates, quantile, MARGIN = 1, 0.25, type = 1)
    predicted.lifetime.ranks[[start.rank]]$rank.lower = apply(replicates, quantile, MARGIN = 1, 0.75, type = 1)
  }
  assign(paste0(condition, '.predicted.ranks'), do.call(rbind, predicted.lifetime.ranks))
}

### Four by four comparison of combinations of inheritance and rank effects
mri.youngest <- ggplot(data = mri.youngest.predicted.ranks, aes(x = time, y = rank, col = start, fill = start))+
  geom_ribbon(aes(ymin = rank.lower + 0.25, ymax = rank.upper - 0.25), alpha = 0.3, color = NA)+
  geom_line(size = 1)+
  ylim(11.25,0.75)+
  theme_classic()+
  scale_color_manual(values = viridis::magma(7))+
  scale_fill_manual(values = viridis::magma(7))+
  theme(legend.position = 'none', axis.text = element_blank(), 
        axis.ticks = element_blank(), axis.title.x = element_blank())+
  ylab('Rank')
mri.youngest 
hyenalike <- mri.youngest + theme(legend.position = 'none', axis.text = element_blank(), 
                                    axis.ticks = element_blank(), axis.title.x = element_text())+
  ylab('Rank') + 
  xlab('Time')

mri.youngest.norank <- ggplot(data = mri.youngest.norank.predicted.ranks, aes(x = time, y = rank, col = start, fill = start))+
  geom_ribbon(aes(ymin = rank.lower + 0.25, ymax = rank.upper - 0.25), alpha = 0.3, color = NA)+
  geom_line(size = 1)+
  ylim(11.25,0.75)+
  theme_classic()+
  scale_color_manual(values = viridis::magma(7))+
  scale_fill_manual(values = viridis::magma(7))+
  theme(legend.position = 'none', axis.text = element_blank(), 
        axis.ticks = element_blank(), axis.title = element_blank())
mri.youngest.norank

mocor.weak.norank <- ggplot(data = mocor.weak.norank.predicted.ranks, aes(x = time, y = rank, col = start, fill = start))+
  geom_ribbon(aes(ymin = rank.lower + 0.25, ymax = rank.upper - 0.25), alpha = 0.3, color = NA)+
  geom_line(size = 1)+
  ylim(11.25,0.75)+
  theme_classic()+
  scale_color_manual(values = viridis::magma(7))+
  scale_fill_manual(values = viridis::magma(7))+
  theme(legend.position = 'none', axis.text = element_blank(), 
        axis.ticks = element_blank(), axis.title.y = element_blank())+
  xlab('Time')
mocor.weak.norank


mocor.weak <- ggplot(data = mocor.weak.predicted.ranks, aes(x = time, y = rank, col = start, fill = start))+
  geom_ribbon(aes(ymin = rank.lower + 0.25, ymax = rank.upper - 0.25), alpha = 0.3, color = NA)+
  geom_line(size = 1)+
  ylim(11.25,0.75)+
  theme_classic()+
  scale_color_manual(values = viridis::magma(7))+
  scale_fill_manual(values = viridis::magma(7))+
  theme(legend.position = 'none', axis.text = element_blank(), 
        axis.ticks = element_blank(), axis.text.x = element_blank())+
  ylab('Rank')+
  xlab('Time')
mocor.weak


pdf('plot_grid.pdf', width = 5, height = 5)
(mri.youngest + mri.youngest.norank) /
  ( mocor.weak + mocor.weak.norank) + plot_annotation(tag_levels = 'A') 
dev.off()


mri.oldest <- ggplot(data = mri.oldest.predicted.ranks, aes(x = time, y = rank, col = start, fill = start))+
  geom_ribbon(aes(ymin = rank.lower + 0.25, ymax = rank.upper - 0.25), alpha = 0.3, color = NA)+
  geom_line(size = 1)+
  ylim(11.25,0.75)+
  theme_classic()+
  scale_color_manual(values = viridis::magma(7))+
  scale_fill_manual(values = viridis::magma(7))+
  theme(legend.position = 'none', axis.text = element_blank(), 
        axis.ticks = element_blank())+
  xlab('Time')+
  ylab('Rank')
mri.oldest


### Not included in paper, but intermediate as expected ### 
pdf('plots/mri_mocor_medium_sim.pdf', width = 4, height = 4)
mocor.medium <- ggplot(data = mocor.medium.predicted.ranks, aes(x = time, y = rank, col = start, fill = start))+
  geom_ribbon(aes(ymin = rank.lower + 0.25, ymax = rank.upper - 0.25), alpha = 0.3, color = NA)+
  geom_line(size = 1)+
  ylim(11.25,0.75)+
  theme_classic()+
  scale_color_manual(values = viridis::magma(7))+
  scale_fill_manual(values = viridis::magma(7))+
  theme(legend.position = 'none', axis.text = element_blank(), 
        axis.ticks = element_blank())+
  xlab('Time')+
  ylab('Rank')
mocor.medium
dev.off()


mocor.strong <- ggplot(data = mocor.strong.predicted.ranks, aes(x = time, y = rank, col = start, fill = start))+
  geom_ribbon(aes(ymin = rank.lower + 0.25, ymax = rank.upper - 0.25), alpha = 0.3, color = NA)+
  geom_line(size = 1)+
  ylim(11.25,0.75)+
  theme_classic()+
  scale_color_manual(values = viridis::magma(7))+
  scale_fill_manual(values = viridis::magma(7))+
  theme(legend.position = 'none', axis.text = element_blank(), 
        axis.ticks = element_blank())+
  xlab('Time')+
  ylab('Rank')

mocor.strong

pdf('different_inheritance_styles.pdf', height = 2.5, width = 7.5)
hyenalike + mri.oldest + mocor.strong + plot_annotation(tag_levels = 'A')
dev.off()



