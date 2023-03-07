################################################################################
#                       Get ranks using DynaRank                               #
#                                                                              #
#                                                                              #
#                           By Eli Strauss                                     #
#                                                                              #
#                            March 2023                                        #
################################################################################


rm(list = ls())
options(stringsAsFactors = FALSE)
setwd('C:\\Users\\strau\\Documents\\code\\inheritance\\0_create_ranks')

library(DynaRankR)
library(dplyr)

load('1.hyena_intx_data.RData')

set.seed(1989)


################################################################################
#
#    Talek
#
talek_ranks <- informed_matreorder(contestants = contestants.talek,
                                   convention = 'mri',
                                   n = 50,
                                   shuffles = 30,
                                   require.corroboration = TRUE,
                                   initial.ranks = initial.ranks.talek,
                                   interactions = interactions.talek)

plot_ranks(talek_ranks)



################################################################################
#
#    North
#

north_ranks <- informed_matreorder(contestants = contestants.serena.n,
                        convention = 'mri',
                        n = 100,
                        shuffles = 30,
                        require.corroboration = TRUE,
                        initial.ranks = initial.ranks.serena.n,
                        interactions = interactions.serena.n)


plot_ranks(north_ranks)

################################################################################
#
#    South
#

south_ranks <- informed_matreorder(contestants = contestants.serena.s,
                        convention = 'mri',
                        n = 100,
                        shuffles = 30,
                        require.corroboration = TRUE,
                        initial.ranks = initial.ranks.serena.s,
                        interactions = interactions.serena.s)

plot_ranks(south_ranks)

################################################################################
#
#    Happy Zebra
#

hz_ranks <- informed_matreorder(contestants = contestants.happy.zebra,
                     convention = 'mri',
                     n = 100,
                     shuffles = 30,
                     require.corroboration = TRUE,
                     initial.ranks = initial.ranks.happy.zebra,
                     interactions = interactions.happy.zebra)

plot_ranks(hz_ranks)




## Combine data
talek_ranks$clan <- 'talek'
north_ranks$clan <- 'serena.n'
south_ranks$clan <- 'serena.s'
hz_ranks$clan <- 'happy.zebra'

#filter to 2017
ranks <- rbind(filter(talek_ranks, period <= 2017),
               filter(north_ranks, period <= 2017),
               filter(south_ranks, period <= 2017),
               filter(hz_ranks,period <= 2017))

tblFemaleRanks <- ranks %>%
  rename(year = period) %>% 
  select('clan', 'year', 'id', 'rank')

tblFemaleRanks$year  <- as.numeric(tblFemaleRanks$year)

save(tblFemaleRanks, file = 'female_ranks.RData')


