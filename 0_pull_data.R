################################################################################
#
#                    Extract data from MHP database 
#
#                        Eli Strauss Dec 2022
################################################################################

library(hyenadata)
library(dplyr)
library(here)

data(tblFemaleRanks)
data(tblLifeHistory.wide)
data(tblHyenas)


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
ars[is.na(ars$ars),]$ars <- 0


### Assemble rank data
tblLifeHistory.wide <- left_join(tblLifeHistory.wide, tblHyenas[,c('id', 'sex')])
lifespans <- filter(tblLifeHistory.wide, sex == 'f', !is.na(disappeared))
lifespans$lifespan <- as.numeric(lifespans$disappeared - lifespans$dob)/365.25
mean(filter(lifespans, lifespan >= 3)$lifespan)


ranks <- tblFemaleRanks[c('clan', 'year', 'id', 'rank')]
ranks <- ranks %>%
  rename(period = year) %>%
  group_by(clan, period) %>%
  mutate(max.rank = max(rank), 
         rank.decile = (rank-1)/(max.rank-1))


save(ars, ranks, file = '0_hyena_data.RData')
