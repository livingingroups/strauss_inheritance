################################################################################
#
#                    Extract data from MHP database 
#
#                        Eli Strauss Dec 2022
################################################################################

library(hyenadata) ### Version 1.2.92
library(dplyr)
library(here)

load('0_create_ranks/female_ranks.RData')
data(tblLifeHistory.wide)
data(tblHyenas)
tblFemaleRanks <- left_join(tblFemaleRanks, tblHyenas[,c('id', 'sex')])

### Assemble annual reproductive success data
reproduction <- tblLifeHistory.wide %>%
  filter(dob_event_data %in% c('talek', 'serena.s', 
                               'serena.n', 'happy.zebra'),
         !is.na(dob),
         is.na(disappeared) | ((disappeared - dob)/365.25) >= 2) %>% ## Offspring must survive to 2
  mutate(year = as.numeric(format(dob, '%Y'))) %>%
  select(id, dob, year, dob_event_data) %>%
  rename(clan = dob_event_data) %>% 
  left_join(tblHyenas[c('id', 'mom', 'sex')]) %>% 
  filter(!is.na(mom)) %>% 
  group_by(mom, year) %>%
  summarize(ars = length(id))

ars <- left_join(tblFemaleRanks, reproduction, by = c('id' = 'mom', 'year'))
ars <- ars %>% filter(sex == 'f') %>% select(clan, year, id, rank, ars)
ars[is.na(ars$ars),]$ars <- 0


### Assemble rank data
tblLifeHistory.wide <- left_join(tblLifeHistory.wide, tblHyenas[,c('id', 'sex')])

ranks <- tblFemaleRanks[c('clan', 'year', 'id', 'rank')]

ranks <- ranks %>%
  rename(period = year) %>%
  group_by(clan, period) %>%
  mutate(max.rank = max(rank), 
         rank.decile = (rank-1)/(max.rank-1))


save(ars, ranks, file = '0_hyena_data.RData')
