################################################################################

#Simulate impact of inheritance on social structure#


#Eli Strauss, Sep 2022#

################################################################################

library(dplyr)
library(here)
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


################################################################################
# Set invariant global parameters
group.size <- 30
nreplicates <- 100
generations <- 4 * group.size
rank.effect.list <- c(0, 1)

letters.long <- c(letters,as.vector(t(outer(X = letters, Y = letters, paste0))))


################################################################################
### Set up group object
group <- data.frame(id = paste0(letters.long[1:group.size],0), 
                    rank = 1:group.size,
                    matriline = paste0(letters.long[1:group.size],0),
                    generation = 0,
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

for(replicate in 1:nreplicates){
  group$replicate <- replicate
  for(re in rank.effect.list){
    group$rank.effect <- re
    output.list[[length(output.list)+1]] <- run_generations(group, generations, 0, 
                                                            f_repro = repro_function,
                                                            f_mort = mortality,
                                                            f_inherit = inherit_mri_youngest,
                                                            inheritance = 'mri_youngest',
                                                            rank.effect = re)
    
    output.list[[length(output.list)+1]] <- run_generations(group, generations, 0, 
                                                            f_repro = repro_function,
                                                            f_mort = mortality,
                                                            f_inherit = inherit_mri_oldest,
                                                            inheritance = 'mri_oldest',
                                                            rank.effect = re)
    
    output.list[[length(output.list)+1]] <- run_generations(group, generations, 0, 
                                                            f_repro = repro_function,
                                                            f_mort = mortality,
                                                            f_inherit = inherit_parent_offspring_cor,
                                                            inheritance = 'parent_offspring_cor',
                                                            rank.effect = re)
    
    output.list[[length(output.list)+1]] <- run_generations(group, generations, 0, 
                                                            f_repro = repro_function,
                                                            f_mort = mortality,
                                                            f_inherit = inherit_random,
                                                            inheritance = 'none',
                                                            rank.effect = re)
  }
}

mor <- do.call(rbind, mor.list)
rank.diff <- do.call(rbind, rank.diff.list)
output <- do.call(rbind, output.list)

par(mfrow = c(2,3))
for(i in unique(mor$inheritance)){
  plot(mor[mor$inheritance == i,]$mother.rank, mor[mor$inheritance == i,]$offspring.rank, 
       col = alpha('dodgerblue', 0.01))
  lines(x = 1:50, y = 1:50)
}



save(output, group.size, mor, file = '4_rank_data_simulated.Rdata')


### What on average are the ranks of added or removed individuals? 
## These are visual checks to confirm that the simulation is working as expected and to inspect outpu

reproranks <- output %>%
  filter(!id %in% output[1:30,'id'], rank.effect == 0) %>%
  group_by(replicate, inheritance) %>%
  filter(!duplicated(id)) %>%
  ungroup()

mortranks <- output %>%
  filter(!id %in% output[1:30,'id'], rank.effect == 0) %>%
  group_by(replicate, inheritance) %>%
  filter(!duplicated(id, fromLast = T)) %>%
  ungroup()

repro <- ggplot(data = reproranks, aes(x = rank, fill = inheritance)) +
  geom_histogram(bins = 31)+
  theme_classic()+
  facet_wrap(~ inheritance)+
  ggtitle('Reproduction')

mort <- ggplot(data = mortranks, aes(x = rank, fill = inheritance)) +
  geom_histogram(bins = 30)+
  theme_classic()+
  facet_wrap(~ inheritance)+
  ggtitle("Mortality")

repro + mort

#### Check average ranks of starting ids
start.ids <- output[1:group.size,]$id
nothing <- filter(output, rank.effect == 0, inheritance == 'none')
none.norank.rawdata <- nothing %>% 
  mutate(rank = as.character(cut((rank - 1)/(group.size - 1),
                                 # breaks = c(-0.001, 0.001, 0.1, 0.2, 0.3, 0.4,
                                 #            0.5, 0.6, 0.7, 0.8, 0.9, 1.001),
                                 breaks = c(-0.001, 0.1, 0.2, 0.3, 0.4,
                                            0.5, 0.6, 0.7, 0.8, 0.9, 1.001),
                                 labels = letters[1:10]))) %>%
  mutate(rank = match(rank, letters[1:10])) %>%
  group_by(replicate, id) %>%
  mutate(first_rank = first(rank),
         age = generation - first(generation)) %>%
  ungroup() %>%
  group_by(first_rank, age) %>% 
  summarize(rank = mean(rank)) %>%
  ggplot(aes(x = age, y = rank, group = first_rank, color = first_rank)) + 
  geom_line() + 
  theme_classic() + 
  ylim(c(10,1))+
  theme(legend.position = 'none')

random <- filter(output, rank.effect == 1, inheritance == 'none')
none.raw.data <- random %>% 
  mutate(rank = as.character(cut((rank - 1)/(group.size - 1),
                                 # breaks = c(-0.001, 0.001, 0.1, 0.2, 0.3, 0.4,
                                 #            0.5, 0.6, 0.7, 0.8, 0.9, 1.001),
                                 breaks = c(-0.001, 0.1, 0.2, 0.3, 0.4,
                                            0.5, 0.6, 0.7, 0.8, 0.9, 1.001),
                                 labels = letters[1:10]))) %>%
  mutate(rank = match(rank, letters[1:10])) %>%
  group_by(replicate, id) %>%
  mutate(first_rank = first(rank),
         age = generation - first(generation)) %>%
  ungroup() %>%
  group_by(first_rank, age) %>% 
  summarize(rank.sd = sd(rank), rank = mean(rank)) %>%
  ggplot(aes(x = age, y = rank, letters[1:10], group = first_rank, color = first_rank)) +
  geom_line() + 
  ylim(c(10, 1))+
  theme_classic()+
  theme(legend.position = 'none')

mri <- filter(output, inheritance == "mri_youngest", rank.effect == 1)
mri.youngest.raw.data <- mri %>% 
  mutate(rank = as.character(cut((rank - 1)/(group.size - 1),
                                 # breaks = c(-0.001, 0.001, 0.1, 0.2, 0.3, 0.4,
                                 #            0.5, 0.6, 0.7, 0.8, 0.9, 1.001),
                                 breaks = c(-0.001, 0.1, 0.2, 0.3, 0.4,
                                            0.5, 0.6, 0.7, 0.8, 0.9, 1.001),
                                 labels = letters[1:10]))) %>%
  mutate(rank = match(rank, letters[1:10])) %>%
  group_by(replicate, id) %>%
  mutate(first_rank = first(rank),
         age = generation - first(generation)) %>%
  ungroup() %>%
  group_by(first_rank, age) %>% 
  summarize(rank = mean(rank)) %>%
  ggplot(aes(x = age, y = rank, group = first_rank, color = first_rank)) + 
  geom_line() + 
  ylim(c(10, 1))+
  theme_classic()+
  theme(legend.position = 'none')

mri.norank <- filter(output, inheritance == "mri_youngest", rank.effect == 0)
mri.youngest.norank.raw.data <- mri.norank %>% 
  mutate(rank = as.character(cut((rank - 1)/(group.size - 1),
                                 # breaks = c(-0.001, 0.001, 0.1, 0.2, 0.3, 0.4,
                                 #            0.5, 0.6, 0.7, 0.8, 0.9, 1.001),
                                 breaks = c(-0.001, 0.1, 0.2, 0.3, 0.4,
                                            0.5, 0.6, 0.7, 0.8, 0.9, 1.001),
                                 labels = letters[1:10]))) %>%
  mutate(rank = match(rank, letters[1:10])) %>%
  group_by(replicate, id) %>%
  mutate(first_rank = first(rank),
         age = generation - first(generation)) %>% 
  ungroup() %>%
  group_by(first_rank, age) %>%
  summarize(rank = mean(rank)) %>%
  ungroup() %>% 
  ggplot(aes(x = age, y = rank, group = first_rank, color = first_rank)) +
  geom_line() + 
  ylim(c(10, 1))+
  theme_classic()+
  theme(legend.position = 'none')

mri.norank %>% 
  filter(id %in% mri.norank$id[1:group.size]) %>% 
  group_by(id, generation) %>%
  summarize(rank.u = mean(rank)) %>% 
  ungroup %>% 
  ggplot(aes(x = generation, y = rank.u, group = id, color = id)) +
  geom_line() + 
  ylim(c(30, 1))+
  theme_classic()+
  theme(legend.position = 'none')
  
mri.oldest <- filter(output, inheritance == "mri_oldest", rank.effect == 1)
mri.oldest%>% 
  group_by(replicate, id) %>%
  mutate(first_rank = first(rank),
         age = generation - first(generation)) %>%
  ungroup() %>%
  group_by(first_rank, age) %>% 
  summarize(rank = mean(rank)) %>%
  ggplot(aes(x = age, y = rank, group = first_rank, color = first_rank)) + 
  geom_line() + 
  theme_classic()

  