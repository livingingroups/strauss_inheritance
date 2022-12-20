################################################################################

#Simulate impact of inheritance on social structure#


#Eli Strauss, Sep 2022#

################################################################################

library(dplyr)
library(here)
library(viridis)
library(lhs)
library(stringr)
library(plot.matrix)
library(purrr)
library(markovchain)
library(ggplot2)
library(DynaRankR)
library(ineq)
library(here)
library(patchwork)

rm(list = ls())
source(here('1_define_functions.R'))
load('2_repro_function.RData')
set.seed(1989)


################################################################################
# Set invariant global parameters
group.size <- 30
turnover = 1 ## 10% of group size
# reproduction <- function(rank, rank.effect = 8){
#   rank <- (rank-1)/(max(rank)-1)
#   1/(10*(1+exp((rank-mean(rank))*rank.effect)))
# }




# mortality <- function(rank, rank.effect = 8){
#   rank <- (rank-1)/(max(rank)-1)
#   1/(10*(1+exp((-rank+mean(rank))*rank.effect)))
# }
# reproduction <- function(rank, rank.effect = 8){
#   rank <- (rank-1)/(max(rank)-1)
#   shape <- (0.1 * exp(rank.effect/2 * -rank))
#   shape/mean(shape) * 0.05
# }
# mortality <- function(rank, rank.effect = 0){
#   rank <- (rank-1)/(max(rank)-1)
#   rank <- rank[length(rank):1]
#   1/(10*(1+exp((rank-mean(rank))*rank.effect)))
# }

## Flat mortality
mortality <- function(stan_rank, rank.effect){
  return(rep(1/length(stan_rank), length(stan_rank)))
}


letters.long <- c(letters,as.vector(t(outer(X = letters, Y = letters, paste0))))

sr <- seq(0, 1, length.out = 30)
plot(sr,repro_function(stan_rank = sr), type = 'l')
points(sr, repro_function(stan_rank = sr, rank.effect = 0))
points(sr, repro_function(stan_rank = sr, rank.effect = 2))
lines(sr, mortality(sr, rank.effect = 0))


################################################################################
### Set up group object
group <- data.frame(id = paste0(letters.long[1:group.size],0), 
                    rank = 1:group.size,
                    matriline = paste0(letters.long[1:group.size],0),
                    generation = 0,
                    rank.inher.precis = NA,
                    rank.effect = NA,
                    replicate = NA,
                    mother = NA,
                    inheritance = NA)


################################################################################
### Set up data objects
mor.list <- list()
rank.diff.list <- list()
output.list <- list()

################################################################################
### Run simulations
rank.inher.precis.list <- c(0, 2, 4)
rank.effect.list <- c(0, 1, 2)
nreplicates <- 20
generations = 15

for(replicate in 1:nreplicates){
  group$replicate <- replicate
  rank.inher.precis <- NA
  for(rank.effect in rank.effect.list){
    group$rank.effect <- rank.effect
    output.list[[length(output.list)+1]] <- run_generations(group, generations, 0, 
                                                            f_repro = repro_function,
                                                            f_mort = mortality,
                                                            f_inherit = inherit_mri_youngest,
                                                            inheritance = 'mri_youngest',
                                                            rank.inher.precis = rank.inher.precis,
                                                            rank.effect = rank.effect)
    
    output.list[[length(output.list)+1]] <- run_generations(group, generations, 0, 
                                                            f_repro = repro_function,
                                                            f_mort = mortality,
                                                            f_inherit = inherit_mri_oldest,
                                                            rank.inher.precis = rank.inher.precis,
                                                            inheritance = 'mri_oldest',
                                                            rank.effect = rank.effect)
    for(rank.inher.precis in rank.inher.precis.list){
      group$rank.inher.precis <- rank.inher.precis
      output.list[[length(output.list)+1]] <- run_generations(group, generations, 0, 
                                                              f_repro = repro_function,
                                                              f_mort = mortality,
                                                              f_inherit = inherit_parent_offspring_cor,
                                                              rank.inher.precis = rank.inher.precis,
                                                              inheritance = paste0('parent_offspring_cor', 
                                                                                   rank.inher.precis),
                                                              rank.effect = rank.effect)
    }
  }
}

mor <- do.call(rbind, mor.list)
rank.diff <- do.call(rbind, rank.diff.list)
output <- do.call(rbind, output.list)

par(mfrow = c(2,3))
for(i in unique(mor$inheritance)){
  plot(mor[mor$inheritance == i,]$mother.rank, mor[mor$inheritance == i,]$offspring.rank, 
       col = 'dodgerblue')
  lines(x = 1:50, y = 1:50)
}



save(output, file = 'rank_data_simulated.Rdata')

# 
# 
# ################################################################################
# ### Visualize output in grid
png(filename= here('plots/mri_grid.png'), width = 5, height = 5, units = 'in',
    res = 300)
par(mfrow = c(3,3),
    mar = c(0.2,0.2,0,0),
    yaxt = 'n',
    xaxt = 'n',
    oma = c(1,3.5,3,1))
for(re in rank.effect.list[3:1]){
  for(ris in rank.inher.precis.list[3:1]){
    output.example <- filter(output, replicate == 1, rank.inher.precis == ris, rank.effect == re)
    output.matrix <- matrix(data = NA, nrow = max(output.example$rank), ncol = max(output.example$generation+1))
    output.matrix[cbind(output.example$rank, output.example$generation+1)] <- factor(output.example$matriline, 
                                                                                        levels = unique(output.example$matriline))
                                                                                      
    plot(output.matrix, col = magma(30), border = NA, key = NULL, xlab = '', ylab = '',
         main = '')
    box(col = 'black', bty = 'L')
    if(which(ris == rank.inher.precis.list) == 3)
      mtext('Rank', side = 2, line = 0, cex = 0.7)
    if(which(re == rank.effect.list) == 1)
      mtext('Time', side = 1, line = 0, cex = 0.7)
    if(which(ris == rank.inher.precis.list) == 3 & which(re == rank.effect.list) == 2){
      arrows(-1, 45, -1, -20, xpd = NA, length = 0.05, lwd = 2)
      mtext('Decreasing effect of rank on fitness', side = 2, cex = 1.3, line = 2)
    }

    if(which(ris == rank.inher.precis.list) == 2 & which(re == rank.effect.list) == 3){
      arrows(-10, 32, 30, 32, xpd = NA, length = 0.05, lwd = 2)
      mtext('Decreasing fidelity of rank inheritance', side = 3, cex = 1.3, line = 1)
    }
  }
}
dev.off()


### What on average are the ranks of added or removed individuals? 

reproranks <- output %>%
  filter(!id %in% output[1:30,'id'], rank.effect == 0) %>%
  group_by(replicate, inheritance, rank.inher.precis) %>%
  filter(!duplicated(id)) %>%
  ungroup()

mortranks <- output %>%
  filter(!id %in% output[1:30,'id'], rank.effect == 0) %>%
  group_by(replicate, inheritance, rank.inher.precis) %>%
  filter(!duplicated(id, fromLast = T)) %>%
  ungroup()

repro <- ggplot(data = reproranks, aes(x = rank, color = inheritance)) +
  geom_histogram(bins = 60)+
  theme_classic()+
  facet_wrap(~ inheritance)+
  ggtitle('Reproduction')

mort <- ggplot(data = mortranks, aes(x = rank, color = inheritance)) +
  geom_histogram(bins = 60)+
  theme_classic()+
  facet_wrap(~ inheritance)+
  ggtitle("Mortality")

repro + mort

