

library(dplyr)
library(here)
library(ggplot2)
library(patchwork )
library(gridExtra)
library(ineq)
library(tidyr)
library(markovchain)

rm(list = ls())
source('1_define_functions.R')
load('4_rank_data_simulated.Rdata')
load('2_repro_function.RData')
load('5_sim_markov_chains.Rdata')
set.seed(1989)

plot.dir <- '~/../Dropbox/Documents/Research/Full_projects/2023 Inheritancy_mobility/plots/'


#### Queen effects

queen.states <- letters[1:11]
queen.output.seq <- output %>% 
  group_by(inheritance, rank.effect, replicate,id) %>%
  mutate(state = as.character(cut((rank - 1)/(group.size - 1), 
                                  breaks = c(-0.001, 0.001, 0.1, 0.2, 0.3, 0.4,
                                             0.5, 0.6, 0.7, 0.8, 0.9, 1.001),
                                  labels = queen.states))) %>%
  summarize(sequence = c(state, NA))



queen.transitions.total.mri.youngest <- markovchainFit(data = filter(queen.output.seq, 
                                                               rank.effect == 1,
                                                               inheritance == 'mri_youngest')$sequence,
                                                 method = 'map', confint = F)
queen.transitions.total.mri.youngest.estimate <- queen.transitions.total.mri.youngest$estimate@transitionMatrix

queen.transitions.total.mri.youngest.norank <- markovchainFit(data = filter(queen.output.seq, 
                                                                      rank.effect == 0,
                                                                      inheritance == 'mri_youngest')$sequence,
                                                        method = 'map', confint = F)
queen.transitions.total.mri.youngest.norank.estimate <- queen.transitions.total.mri.youngest.norank$estimate@transitionMatrix


queen.transitions.total.mri.oldest <- markovchainFit(data = filter(queen.output.seq, 
                                                             rank.effect == 1,
                                                             inheritance == 'mri_oldest')$sequence,
                                               method = 'map', confint = F)
queen.transitions.total.mri.oldest.estimate <- queen.transitions.total.mri.oldest$estimate@transitionMatrix

queen.transitions.total.none <- markovchainFit(data = filter(queen.output.seq, 
                                                       rank.effect == 1,
                                                       inheritance == 'none')$sequence,
                                         method = 'map', confint = F)
queen.transitions.total.none.estimate <- queen.transitions.total.none$estimate@transitionMatrix


queen.transitions.total.none.norank <- markovchainFit(data = filter(queen.output.seq, 
                                                              rank.effect == 0,
                                                              inheritance == 'none')$sequence,
                                                method = 'map', confint = F)
queen.transitions.total.none.norank.estimate <- queen.transitions.total.none.norank$estimate@transitionMatrix

queen.transitions.total.mocor <- markovchainFit(data = filter(queen.output.seq, 
                                                        rank.effect == 1,
                                                        inheritance == 'parent_offspring_cor')$sequence,
                                          method = 'map', confint = F)
queen.transitions.total.mocor.estimate <- queen.transitions.total.mocor$estimate@transitionMatrix

### Set up markov chain

queen.mc.mri.youngest <- new('markovchain', states = queen.states,
                       transitionMatrix = queen.transitions.total.mri.youngest.estimate)

queen.mc.mri.youngest.norank <- new('markovchain', states = queen.states,
                              transitionMatrix = queen.transitions.total.mri.youngest.norank.estimate)

queen.mc.mri.oldest <- new('markovchain', states = queen.states,
                     transitionMatrix = queen.transitions.total.mri.oldest.estimate)

queen.mc.none <- new('markovchain', states = queen.states,
               transitionMatrix = queen.transitions.total.none.estimate)

queen.mc.none.norank <- new('markovchain', states = queen.states,
                      transitionMatrix = queen.transitions.total.none.norank.estimate)

queen.mc.mocor <- new('markovchain', states = queen.states,
                transitionMatrix = queen.transitions.total.mocor.estimate)


### Rerun markov models with 11th state for queen
repro_events <- 120 ## Based on the expected number of reproductive events experience by a moderately long-lived female
reps = 10000
conditions <- c('mri.youngest', 'mri.youngest.norank', 'mri.oldest', 'none', 'none.norank', 'mocor')
queen.predicted.trajectories <- array(NA, dim = c(length(queen.states), repro_events+1, reps, length(conditions)), 
                                      dimnames = list(queen.states, 1:(repro_events+1), 1:10000, conditions))
for(condition in conditions){
  queen.predicted.lifetime.ranks <- list()
  mrf <- get(paste0('queen.mc.', condition))
  for(start.rank in queen.states){
    queen.predicted.lifetime.ranks[[start.rank]] <- data.frame(
      start = rep(start.rank, repro_events+1),
      time = seq(1:(repro_events+1))
    )
    replicates = replicate(reps, markovchainSequence(repro_events, mrf, t0 = start.rank, include.t0 = T))
    replicates = matrix(match(replicates, queen.states), nrow = repro_events+1, ncol = reps)
    
    queen.predicted.trajectories[start.rank, , ,condition] <- replicates
    
    queen.predicted.lifetime.ranks[[start.rank]]$rank = apply(replicates, median, MARGIN = 1)
    queen.predicted.lifetime.ranks[[start.rank]]$rank.mean = apply(replicates, mean, MARGIN = 1)
    queen.predicted.lifetime.ranks[[start.rank]]$rank.mean.sd = apply(replicates, sd, MARGIN = 1)
    queen.predicted.lifetime.ranks[[start.rank]]$rank.upper = apply(replicates, quantile, MARGIN = 1, 0.25, type = 1)
    queen.predicted.lifetime.ranks[[start.rank]]$rank.lower = apply(replicates, quantile, MARGIN = 1, 0.75, type = 1)
  }
  assign(paste0(condition, 'queen.predicted.ranks'), do.call(rbind, queen.predicted.lifetime.ranks))
}




#### New types of inheritance  

none.queen <- ggplot(data = filter(nonequeen.predicted.ranks, start %in% queen.states),
               aes(x = time, y = rank.mean, col = start, fill = start))+
  geom_ribbon(aes(ymin = rank.mean - rank.mean.sd, ymax = rank.mean + rank.mean.sd), alpha = 0.2, color = NA)+
  geom_line(size = 1)+
  theme_classic()+
  scale_color_manual(values = colorRampPalette(c('black', 'dodgerblue'))(11))+
  scale_fill_manual(values = colorRampPalette(c('black', 'dodgerblue'))(11))+
  theme(legend.position = 'none', axis.text.x = element_blank(), 
        axis.ticks = element_blank(), axis.text.y = element_text(angle = 90))+
  scale_y_continuous(breaks = c(11,1), labels = c('low', 'high'))+
  ylab('Rank')+
  xlab('Time')+
  coord_cartesian(ylim = c(11,1)) 

none.norank.queen <- ggplot(data = filter(none.norankqueen.predicted.ranks, start %in% queen.states),
                     aes(x = time, y = rank.mean, col = start, fill = start))+
  geom_ribbon(aes(ymin = rank.mean - rank.mean.sd, ymax = rank.mean + rank.mean.sd), alpha = 0.2, color = NA)+
  geom_line(size = 1)+
  theme_classic()+
  scale_color_manual(values = colorRampPalette(c('black', 'dodgerblue'))(11))+
  scale_fill_manual(values = colorRampPalette(c('black', 'dodgerblue'))(11))+
  theme(legend.position = 'none', axis.text.x = element_blank(), 
        axis.ticks = element_blank(), axis.text.y = element_text(angle = 90))+
  scale_y_continuous(breaks = c(11,1), labels = c('low', 'high'))+
  ylab('Rank')+
  xlab('Time')+
  coord_cartesian(ylim = c(11,1)) 

mri.youngest.queen <- ggplot(data = filter(mri.youngestqueen.predicted.ranks, start %in% queen.states), 
                       aes(x = time, y = rank.mean, col = start, fill = start))+
  geom_ribbon(aes(ymin = rank.mean - rank.mean.sd, ymax = rank.mean + rank.mean.sd), alpha = 0.3, color = NA)+
  geom_line(size = 1)+
  theme_classic()+
  scale_color_manual(values = colorRampPalette(c('black', 'dodgerblue'))(11))+
  scale_fill_manual(values = colorRampPalette(c('black', 'dodgerblue'))(11))+
  theme(legend.position = 'none', axis.text.x = element_blank(), 
        axis.ticks = element_blank(), axis.text.y = element_text(angle = 90))+
  scale_y_continuous(breaks = c(11,1), labels = c('low', 'high'))+
  ylab('Rank')+
  xlab('Time')+
  coord_cartesian(ylim = c(11,1))

mri.youngest.norank.queen <- ggplot(data = filter(mri.youngest.norankqueen.predicted.ranks, start %in% states), 
                              aes(x = time, y = rank.mean, col = start, fill = start))+
  geom_ribbon(aes(ymin = rank.mean - rank.mean.sd, ymax = rank.mean + rank.mean.sd), alpha = 0.2, color = NA)+
  geom_line(size = 1)+
  theme_classic()+
  scale_color_manual(values = colorRampPalette(c('black', 'dodgerblue'))(10))+
  scale_fill_manual(values = colorRampPalette(c('black', 'dodgerblue'))(10))+
  theme(legend.position = 'none', axis.text.x = element_blank(), 
        axis.ticks = element_blank(), axis.text.y = element_text(angle = 90))+
  scale_y_continuous(breaks = c(11,1), labels = c('low', 'high'))+
  ylab('Rank')+
  xlab('Time')+
  coord_cartesian(ylim = c(11,1))

mri.oldest.queen <- ggplot(data = filter(mri.oldestqueen.predicted.ranks, start %in% queen.states),
                     aes(x = time, y = rank.mean, col = start, fill = start))+
  geom_ribbon(aes(ymin = rank.mean + rank.mean.sd, ymax = rank.mean - rank.mean.sd), alpha = 0.3, color = NA)+
  geom_line(size = 1)+
  theme_classic()+
  scale_color_manual(values = colorRampPalette(c('black', 'dodgerblue'))(11))+
  scale_fill_manual(values = colorRampPalette(c('black', 'dodgerblue'))(11))+
  theme(legend.position = 'none', axis.text.x = element_blank(), 
        axis.ticks = element_blank(), axis.text.y = element_text(angle = 90))+
  scale_y_continuous(breaks = c(11,1), labels = c('low', 'high'))+
  xlab('Time')+
  ylab('Rank')+
  coord_cartesian(ylim = c(11, 1))
mri.oldest.queen


mocor.queen <- ggplot(data = filter(mocorqueen.predicted.ranks, start %in% queen.states),
                aes(x = time, y = rank.mean, col = start, fill = start))+
  geom_ribbon(aes(ymin = rank.mean - rank.mean.sd, ymax = rank.mean + rank.mean.sd), alpha = 0.3, color = NA)+
  geom_line(size = 1)+
  theme_classic()+
  scale_color_manual(values = colorRampPalette(c('black', 'dodgerblue'))(11))+
  scale_fill_manual(values = colorRampPalette(c('black', 'dodgerblue'))(11))+
  theme(legend.position = 'none', axis.text.x = element_blank(), 
        axis.ticks = element_blank(), axis.text.y = element_text(angle = 90))+
  scale_y_continuous(breaks = c(10,1), labels = c('low', 'high'))+
  ylab('Rank')+
  xlab('Time')+
  coord_cartesian(ylim = c(11,1))

mocor.queen

mri.youngest + 
  geom_line(data = filter(mri.youngestqueen.predicted.ranks, start == queen.states[1]), col = 'coral', lty = 2, size = 1)

none + 
  geom_line(data = filter(nonequeen.predicted.ranks, start == queen.states[1]), col = 'coral', lty = 2, size = 1)


pdf(paste0(plot.dir, 'simulations.pdf'), width = 8, height = 4)
mri.youngest + geom_line(data = filter(mri.youngestqueen.predicted.ranks, start == queen.states[1]), col = 'coral', lty = 2, size = 1) + 
mri.youngest.limit +
  mri.youngest.norank + geom_line(data = filter(mri.youngest.norankqueen.predicted.ranks, start == queen.states[1]), col = 'coral', lty = 2, size = 1) +  mri.youngest.norank.limit +
  none + geom_line(data = filter(nonequeen.predicted.ranks, start == queen.states[1]), col = 'coral', lty = 2, size = 1) +
  none.limit +
  none.norank + geom_line(data = filter(none.norankqueen.predicted.ranks, start == queen.states[1]), col = 'coral', lty = 2, size = 1) +
  none.norank.limit + 
  plot_layout(design = '
  AAABBCCCDD
  EEEFFGGGHH
  ')
dev.off()



#### Predicted reproduction over lifetime

repro_events <- 120 ## Based on the expected number of reproductive events experience by a moderately long-lived female
reps = 10000
repro.conditions <- c('none', 'mri.youngest', 'mri.oldest', 'mocor')
predicted.repro <- array(NA, dim = c(length(queen.states), reps, length(repro.conditions)), dimnames = list(queen.states, 1:10000, repro.conditions))

for(condition in repro.conditions){
  for(start.rank in queen.states){
    
    replicates <- queen.predicted.trajectories[start.rank, , , condition]
    predicted.repro[start.rank, , condition] <- apply(X = replicates, MARGIN = 2, FUN = function(x)(sum(repro_function(stan_rank = 1 - (x-1)/(length(queen.states)-1)))))
  }
}


### Predicted reproduction
reproductive.inequality.list <- list()

for(condition in repro.conditions){
  reproductive.inequality.list[[length(reproductive.inequality.list) + 1]] <- 
    data.frame(gini = apply(X = predicted.repro[, ,condition], MARGIN = 2, FUN = Gini),
               condition)
}
reproductive.inequality <- do.call(rbind, reproductive.inequality.list)


## Directly compare means instead of doing null hypothesis test (White et al. 2014)
mean(filter(reproductive.inequality, condition == 'none')$gini)
mean(filter(reproductive.inequality, condition == 'mocor')$gini)
mean(filter(reproductive.inequality, condition == 'mri.youngest')$gini)
mean(filter(reproductive.inequality, condition == 'mri.oldest')$gini)
Gini(repro_function(stan_rank = seq(0, 1, length.out = 11)))



repro.plot <- reproductive.inequality %>%
  group_by(condition) %>%
  summarize(gini.u = mean(gini), gini.sd = sd(gini),
            gini.upper = quantile(gini, 0.975), gini.lower = quantile(gini, 0.025))

repro.plot$condition <- factor(repro.plot$condition,
                               levels = c('none', 'mri.youngest', 'mri.oldest', 'mocor'),
                               labels = c('No\ninheritance', 'MRI:\nYoungest ascendency', 
                                          'MRI:\nPrimogeniture', 'Parent-offspring\n correlation')) 

#pdf('reproductive_inequality.pdf', width = 7, height = 4)
reproduction <- ggplot(repro.plot, aes(x = condition, y = gini.u))+
  geom_point(size = 3) + 
  geom_errorbar(aes(x = condition, ymin = gini.u - gini.sd, ymax = gini.u + gini.sd), width = 0.1, size = 0.5)+
  theme_classic()+
  ylab("Inequality in expected lifetime reproduction")+
  theme(axis.text = element_text(color = 'black'),
        axis.title.x = element_blank())+
  geom_hline(yintercept = Gini(x = repro_function(seq(from = 0, to = 1, length.out = 11))),
             lty =2)+
  ylim(c(0.12, 0.26))+
  scale_x_discrete(position = 'top')
#dev.off()


pdf(paste0(plot.dir, 'reproductive_inequality.pdf'), width = 9, height = 5)       
reproduction + 
  inset_element(none.queen, 0.08, 0.6, 0.30, 0.9, align_to = 'full') +
  inset_element(mri.youngest.queen, 0.30, 0.6, 0.52, 0.9, align_to = 'full')+ 
  inset_element(mri.oldest.queen, 0.52, 0.6, 0.74, 0.9, align_to = 'full')+
  inset_element(mocor.queen, 0.74, 0.6, 0.96, 0.9, align_to = 'full')
dev.off()
