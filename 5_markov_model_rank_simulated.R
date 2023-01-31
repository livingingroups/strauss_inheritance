################################################################################

# Calculate rank trajectories of simulated societies#


#Eli Strauss, Sep 2022#

################################################################################


library(dplyr)
library(here)
library(markovchain)
library(ggplot2)
library(patchwork )
library(gridExtra)
library(expm)

rm(list = ls())
plot.dir <- '~/Dropbox/Documents/Research/Full_projects/2023 Inheritancy_mobility/plots/'
source('1_define_functions.R')
load('4_rank_data_simulated.Rdata')
load('2_repro_function.RData')
set.seed(1989)

states <- letters[1:10]
output.seq <- output %>% 
  group_by(inheritance, rank.effect, replicate,id) %>%
  mutate(state = as.character(cut((rank - 1)/(group.size - 1), 
                                  # breaks = c(-0.001, 0.001, 0.1, 0.2, 0.3, 0.4, 
                                  #            0.5, 0.6, 0.7, 0.8, 0.9, 1.001), 
                                  breaks = c(-0.001, 0.1, 0.2, 0.3, 0.4,
                                             0.5, 0.6, 0.7, 0.8, 0.9, 1.001),
                                  labels = states))) %>%
  summarize(sequence = c(state, NA))

## Total changes
transitions.total.mri.youngest <- markovchainFit(data = filter(output.seq, 
                                                               rank.effect == 1,
                                                               inheritance == 'mri_youngest')$sequence,
                                                 method = 'map', confint = F)
transitions.total.mri.youngest.estimate <- transitions.total.mri.youngest$estimate@transitionMatrix
limit.mri.youngest <- data.frame(probs = find_limit(transitions.total.mri.youngest.estimate),
                                        state = states)

transitions.total.mri.youngest.norank <- markovchainFit(data = filter(output.seq, 
                                                                      rank.effect == 0,
                                                                      inheritance == 'mri_youngest')$sequence,
                                                        method = 'map', confint = F)
transitions.total.mri.youngest.norank.estimate <- transitions.total.mri.youngest.norank$estimate@transitionMatrix
limit.mri.youngest.norank <- data.frame(probs = find_limit(transitions.total.mri.youngest.norank.estimate),
                                        state = states)

transitions.total.mri.oldest <- markovchainFit(data = filter(output.seq, 
                                                             rank.effect == 1,
                                                             inheritance == 'mri_oldest')$sequence,
                                               method = 'map', confint = F)
transitions.total.mri.oldest.estimate <- transitions.total.mri.oldest$estimate@transitionMatrix
limit.oldest <- data.frame(probs = find_limit(transitions.total.mri.oldest.estimate),
                           states)

transitions.total.none <- markovchainFit(data = filter(output.seq, 
                                                       rank.effect == 1,
                                                       inheritance == 'none')$sequence,
                                         method = 'map', confint = F, possibleStates = states)
transitions.total.none.estimate <- transitions.total.none$estimate@transitionMatrix
limit.none <- data.frame(probs = find_limit(transitions.total.none.estimate),
                         states)

transitions.total.none.norank <- markovchainFit(data = filter(output.seq, 
                                                              rank.effect == 0,
                                                              inheritance == 'none')$sequence,
                                                method = 'map', confint = F)
transitions.total.none.norank.estimate <- transitions.total.none.norank$estimate@transitionMatrix
limit.none.norank <- data.frame(probs = find_limit(transitions.total.none.norank.estimate),
                                states)

transitions.total.mocor <- markovchainFit(data = filter(output.seq, 
                                                        rank.effect == 1,
                                                        inheritance == 'parent_offspring_cor')$sequence,
                                          method = 'map', confint = F)
transitions.total.mocor.estimate <- transitions.total.mocor$estimate@transitionMatrix
limit.mocor <- data.frame(probs = find_limit(transitions.total.mocor.estimate),
                          states)

### Set up markov chain

mc.mri.youngest <- new('markovchain', states = states,
                  transitionMatrix = transitions.total.mri.youngest.estimate)

mc.mri.youngest.norank <- new('markovchain', states = states,
                         transitionMatrix = transitions.total.mri.youngest.norank.estimate)

mc.mri.oldest <- new('markovchain', states = states,
                         transitionMatrix = transitions.total.mri.oldest.estimate)

mc.none <- new('markovchain', states = states,
                       transitionMatrix = transitions.total.none.estimate)

mc.none.norank <- new('markovchain',
                       transitionMatrix = transitions.total.none.norank.estimate)

mc.mocor <- new('markovchain', states = states,
                             transitionMatrix = transitions.total.mocor.estimate)



## Run markov chain simulation

repro_events <- 120 ## Based on the expected number of reproductive events experience by a moderately long-lived female
reps = 10000
conditions <- c('mri.youngest', 'mri.youngest.norank', 'mri.oldest', 'none', 'none.norank', 'mocor')
predicted.trajectories <- array(NA, dim = c(length(states), repro_events+1, reps, length(conditions)), 
                                dimnames = list(states, 1:(repro_events+1), 1:10000, conditions))
for(condition in conditions){
  predicted.lifetime.ranks <- list()
  mrf <- get(paste0('mc.', condition))
  for(start.rank in states){
    predicted.lifetime.ranks[[start.rank]] <- data.frame(
      start = rep(start.rank, repro_events+1),
      time = seq(1:(repro_events+1))
    )
    replicates = replicate(reps, markovchainSequence(repro_events, mrf, t0 = start.rank, include.t0 = T))
    replicates = matrix(match(replicates, states), nrow = repro_events+1, ncol = reps)
    
    predicted.trajectories[start.rank, , ,condition] <- replicates
    
    predicted.lifetime.ranks[[start.rank]]$rank = apply(replicates, median, MARGIN = 1)
    predicted.lifetime.ranks[[start.rank]]$rank.mean = apply(replicates, mean, MARGIN = 1)
    predicted.lifetime.ranks[[start.rank]]$rank.mean.sd = apply(replicates, sd, MARGIN = 1)
    predicted.lifetime.ranks[[start.rank]]$rank.upper = apply(replicates, quantile, MARGIN = 1, 0.25, type = 1)
    predicted.lifetime.ranks[[start.rank]]$rank.lower = apply(replicates, quantile, MARGIN = 1, 0.75, type = 1)
  }
  assign(paste0(condition, '.predicted.ranks'), do.call(rbind, predicted.lifetime.ranks))
}

save(mc.mri.youngest, mc.mri.oldest, mc.none, mc.mocor,
     mc.mri.youngest.norank, mc.none.norank,
     predicted.trajectories, mri.youngest.predicted.ranks, 
     mri.oldest.predicted.ranks, none.predicted.ranks, mocor.predicted.ranks,
     states,
     file = '5_sim_markov_chains.Rdata')


### Four by four comparison of combinations of inheritance and rank effects
mri.youngest <- ggplot(data = filter(mri.youngest.predicted.ranks, start %in% states), 
                       aes(x = time, y = rank.mean, col = start, fill = start))+
  geom_ribbon(aes(ymin = rank.mean - rank.mean.sd, ymax = rank.mean + rank.mean.sd), alpha = 0.3, color = NA)+
  geom_line(size = 1)+
  theme_classic()+
  scale_color_manual(values = colorRampPalette(c('black', 'dodgerblue'))(10))+
  scale_fill_manual(values = colorRampPalette(c('black', 'dodgerblue'))(10))+
  ylab('Rank')+
  xlab('Time')+
  coord_cartesian(ylim = c(10,1))+
  labs(tag = 'A')+
  theme(legend.position = 'none', axis.text.x = element_blank(), 
        axis.ticks = element_blank(), axis.text.y = element_text(angle = 90))+
  scale_y_continuous(breaks = c(10,1), labels = c('low', 'high'))

mri.youngest.limit <- ggplot(data = limit.mri.youngest, 
       aes(y = states[10:1], x = probs, fill = states[10:1])) +
  geom_bar(stat= 'identity') +
  scale_fill_manual(values = colorRampPalette(c('dodgerblue2', 'black'))(10))+
  theme_void()+
  theme(legend.position = 'none',
        axis.line.x = element_line(), axis.text.x = element_text(),
        axis.title.x = element_text())+
  scale_x_continuous(breaks = c(0, 0.5),
                  limits = c(0, 0.52), name = 'Steady state')

#mri.youngest + mri.youngest.limit + plot_layout(design = 'AAB')
#mri.youngest.raw.data + mri.youngest.limit + mri.youngest

mri.youngest.norank <- ggplot(data = filter(mri.youngest.norank.predicted.ranks, start %in% states), 
                              aes(x = time, y = rank.mean, col = start, fill = start))+
  geom_ribbon(aes(ymin = rank.mean - rank.mean.sd, ymax = rank.mean + rank.mean.sd), alpha = 0.2, color = NA)+
  geom_line(size = 1)+
  theme_classic()+
  scale_color_manual(values = colorRampPalette(c('black', 'dodgerblue'))(10))+
  scale_fill_manual(values = colorRampPalette(c('black', 'dodgerblue'))(10))+
  ylab('Rank')+
  xlab('Time')+
  coord_cartesian(ylim = c(10,1))+
  labs(tag = 'B')+
  theme(legend.position = 'none', axis.text.x = element_blank(), 
        axis.ticks = element_blank(), axis.text.y = element_text(angle = 90))+
  scale_y_continuous(breaks = c(10,1), labels = c('low', 'high'))
#mri.youngest.norank

mri.youngest.norank.limit <- ggplot(data = limit.mri.youngest.norank, 
                             aes(y = states[10:1], x = probs, fill = states[10:1])) +
  geom_bar(stat= 'identity') +
  scale_fill_manual(values = colorRampPalette(c('dodgerblue2', 'black'))(10))+
  theme_void()+
  theme(legend.position = 'none',
        axis.line.x = element_line(), axis.text.x = element_text(),
        axis.title.x = element_text())+
  scale_x_continuous(breaks = c(0, 0.5),
                     limits = c(0, 0.5), name = 'Steady state')

# mri.youngest.norank + mri.youngest.norank.limit + plot_layout(design = 'AAB')
# mri.youngest.norank.raw.data + mri.youngest.norank.limit  + mri.youngest.norank

none.norank <- ggplot(data = filter(none.norank.predicted.ranks,start %in% states), 
                            aes(x = time, y = rank.mean, col = start, fill = start))+
  geom_ribbon(aes(ymin = rank.mean - rank.mean.sd, ymax = rank.mean + rank.mean.sd), alpha = 0.2, color = NA)+
  geom_line(size = 1)+
  theme_classic()+
  scale_color_manual(values = colorRampPalette(c('black', 'dodgerblue'))(10))+
  scale_fill_manual(values = colorRampPalette(c('black', 'dodgerblue'))(10))+
  theme(legend.position = 'none', axis.text.x = element_blank(), 
        axis.ticks = element_blank(), axis.text.y = element_text(angle = 90))+
  scale_y_continuous(breaks = c(10,1), labels = c('low', 'high'))+
  ylab('Rank')+
  xlab('Time')+
  coord_cartesian(ylim = c(10,1))+
  labs(tag= 'D')
#none.norank

none.norank.limit <- ggplot(data = limit.none.norank, 
                                    aes(y = states[10:1], x = probs, fill = states[10:1])) +
  geom_bar(stat= 'identity') +
  scale_fill_manual(values = colorRampPalette(c('dodgerblue2', 'black'))(10))+
  theme_void()+
  theme(legend.position = 'none',
        axis.line.x = element_line(), axis.text.x = element_text(),
        axis.title.x = element_text())+
  scale_x_continuous(breaks = c(0, 0.5),
                     limits = c(0, 0.5), name = 'Steady state')


#none.norank.rawdata + none.norank.limit + none.norank 

none <- ggplot(data = filter(none.predicted.ranks, start %in% states),
                     aes(x = time, y = rank.mean, col = start, fill = start))+
  geom_ribbon(aes(ymin = rank.mean - rank.mean.sd, ymax = rank.mean + rank.mean.sd), alpha = 0.2, color = NA)+
  geom_line(size = 1)+
  theme_classic()+
  scale_color_manual(values = colorRampPalette(c('black', 'dodgerblue'))(10))+
  scale_fill_manual(values = colorRampPalette(c('black', 'dodgerblue'))(10))+
  theme(legend.position = 'none', axis.text.x = element_blank(), 
        axis.ticks = element_blank(), axis.text.y = element_text(angle = 90))+
  scale_y_continuous(breaks = c(10,1), labels = c('low', 'high'))+
  ylab('Rank')+
  xlab('Time')+
  coord_cartesian(ylim = c(10,1)) + 
  labs(tag = 'C')
#none

none.limit <- ggplot(data = limit.none, 
       aes(y = states[10:1], x = probs, fill = states[10:1])) +
  geom_bar(stat= 'identity') +
  scale_fill_manual(values = colorRampPalette(c('dodgerblue2', 'black'))(10))+
  theme_void()+
  theme(legend.position = 'none',
        axis.line.x = element_line(), axis.text.x = element_text(),
        axis.title.x = element_text())+
  scale_x_continuous(breaks = c(0, 0.5),
                     limits = c(0, 0.52), name = 'Steady state')

#none.raw.data + none.limit + none


pdf(paste(plot.dir, 'simulations.pdf'), width = 8, height = 4)
mri.youngest + mri.youngest.limit +
  mri.youngest.norank + mri.youngest.norank.limit +
  none + none.limit +
  none.norank + none.norank.limit + 
  plot_layout(design = '
  AAABBCCCDD
  EEEFFGGGHH
  ')
dev.off()

# Median hierarchy position decile (10 = highest 10% s, 1 = lowest 10%) expected rank in the long run
seq(from = 10, to = 1)[median(rep(1:length(states), times = round(limit.mri.youngest$probs * 1000000)))]
seq(from = 10, to = 1)[median(rep(1:length(states), times = round(limit.mri.youngest.norank$probs * 1000000)))]
seq(from = 10, to = 1)[median(rep(1:length(states), times = round(limit.none$probs * 1000000)))]
seq(from = 10, to = 1)[median(rep(1:length(states), times = round(limit.none.norank$probs * 1000000)))]



