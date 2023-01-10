################################################################################

#Define functions for analysis of effects of inheritance on social structure#


#Eli Strauss, Sep 2022#

################################################################################


################################################################################
trunc_norm <- function(mean, n, sd, max, min){
  x <- rnorm(mean = mean, n = 1000*n, sd = sd)
  x <- x[x <= max & x >= min]
  return(sample(x, n))
}
################################################################################

################################################################################
## Flat mortality
mortality <- function(stan_rank, rank.effect){
  return(rep(1/length(stan_rank), length(stan_rank)))
}
################################################################################

################################################################################
standardize_rank <- function(rank){
  return(1 -(rank-1)/(max(rank)-1))
}
################################################################################

################################################################################
# pick_reproducer <- function(g, f_repro){
#   return(seq(1, nrow(g), 1)[as.logical(rbinom(nrow(g), 1, f_repro(nrow(g),
#                                                                         seq(1, 1-1, length.out = nrow(g)))))])
# }


# pick_reproducer <- function(g, f_repro, rank.effect){
#   return(seq(1, nrow(g), 1)[as.logical(rbinom(nrow(g), 1, 
#                                               f_repro(rank = 1:nrow(g), 
#                                                       rank.effect = rank.effect)))])
# }

pick_reproducer <- function(g, f_repro, rank.effect){
  
  sample(1:nrow(g), size = 1, 
         prob = f_repro(stan_rank = standardize_rank(g$rank),
                               rank.effect = rank.effect))
}

################################################################################

################################################################################
# pick_mortality <- function(g, f_mort, rank.effect){
#   return(seq(1, nrow(g), 1)[as.logical(rbinom(nrow(g), 1, 
#                                               f_mort(rank = 1:nrow(g), 
#                                                       rank.effect = rank.effect)))])
# }
pick_mortality <- function(g, f_mort, rank.efect){
  sample(1:nrow(g), size = 1, 
         prob = f_mort(stan_rank = standardize_rank(g$rank),
                        rank.effect = rank.effect))
}
################################################################################

################################################################################
inherit_mri_youngest <- function(g, new.id, mother, rank.inher.precis){
  new.id$rank <- g$rank[match(mother$id, g$id)]+0.1
  return(new.id)
}
################################################################################

################################################################################
get_descendants <- function(g, id){
  kin <- g[g$mother %in% id,]$id
  if(length(kin)){
    return(unique(c(kin, get_descendants(g, kin))))
  }else{return(NULL)}
}
################################################################################

################################################################################
inherit_mri_oldest <- function(g, new.id, mother, rank.inher.precis){
  
  sibs <- g[!is.na(g$mother) & g$mother == new.id$mother,]
  if(nrow(sibs)){
    sibs.birth.gen <- as.numeric(do.call(rbind,(strsplit(sibs$id, split = '\\.')))[,2])
    oldest.sib <- sibs[which.min(sibs.birth.gen),]
    sibs.matriline <- rbind(oldest.sib, g[g$id %in% get_descendants(g, oldest.sib$id),])
    new.id$rank <- max(sibs.matriline$rank)+0.1
  }else{
    new.id$rank  <- g$rank[g$id == new.id$mother] +0.1
  }
  return(new.id)
}
################################################################################

################################################################################
inherit_parent_offspring_cor <- function(g, new.id, mother, rank.inher.precis){
  
  min.change <-  -mother$rank 
  max.change <- max(g$rank) - mother$rank + 1
  
  if(rank.inher.precis == 0){
    place.above <- sample(c(g$rank, max(g$rank) + 1), 1)
    new.id$rank <- place.above - 0.1
  }else{
    new.id$rank <- mother$rank + trunc_norm(mean = 0, n = 1, 
                                            sd = max(g$rank)/2^rank.inher.precis,
                                            max = max.change,
                                            min = min.change)
  }
  
  return(new.id)
}
################################################################################

################################################################################
add_new_id <- function(g, mother, reproduce, die, f_inherit, inheritance, rank.inher.precis){
  
  reproducer <- g$id[reproduce]
  n.id <- paste0(str_extract(reproducer,
                             '[:alpha:]+'), 
                 as.numeric(str_extract(reproducer,
                                        '[:digit:]+'))+1, 
                 '.',
                 unique(g$generation+1))
  
  new.id <- data.frame(id = n.id,
                       rank = NA,
                       matriline = g$matriline[reproduce],
                       generation = unique(g$generation),
                       rank.inher.precis = rank.inher.precis,
                       rank.effect = rank.effect,
                       replicate = replicate,
                       mother = mother$id,
                       inheritance = inheritance)
  
  new.id <- f_inherit(g, new.id, mother, rank.inher.precis)
  
  g <- arrange(rbind(g[-die,], new.id), rank)
  g$rank <- 1:nrow(g)
  
  return(g)
}



################################################################################

################################################################################
generation_step <- function(g, f_inherit, f_repro, f_mort, inheritance, rank.inher.precis, rank.effect){
  
  initial.group <- g
  reproduce <- pick_reproducer(g, f_repro, rank.effect)
  mother <- g[reproduce,]
  die <- pick_mortality(g, f_mort, rank.effect)

  g <- add_new_id(g, mother, reproduce, die, f_inherit, inheritance, rank.inher.precis)
  
  g$generation <- g$generation + 1
  
  #### Record data into objects in global environment
  .GlobalEnv$mor.list[[length(.GlobalEnv$mor.list)+1]] <- data.frame(mother.rank = g[match(g[!g$id %in% initial.group$id,]$mother, g$id),]$rank,
                                                                     offspring.rank = g[!g$id %in% initial.group$id,]$rank,
                                                                     rank.inher.precis = rank.inher.precis,
                                                                     rank.effect = rank.effect,
                                                                     replicate = replicate,
                                                                     mother = g[match(g[!g$id %in% initial.group$id,]$mother, g$id),]$id,
                                                                     inheritance)
  
  
  mean.ranks <- by(g,
                   INDICES = g$matriline,
                   FUN = function(x)(mean(x$rank)))
  mean.ranks <- data.frame(matriline = names(mean.ranks),
                           mean.rank = mean.ranks[], row.names = NULL)
  
  mean.ranks <- as.data.frame(arrange(mean.ranks, matriline))
  
  .GlobalEnv$rank.diff.list[[length(.GlobalEnv$rank.diff.list)+1]] <- 
    data.frame(matrilines = paste(mean.ranks[1:(nrow(mean.ranks)-1),'matriline'], 
                                  mean.ranks[2:nrow(mean.ranks),'matriline'],
                                  sep = '-'), 
               rank.diff = mean.ranks[2:nrow(mean.ranks),'mean.rank'] - 
                 mean.ranks[1:(nrow(mean.ranks)-1),'mean.rank'],
               generation = unique(g$generation),
               rank.inher.precis = rank.inher.precis,
               rank.effect = rank.effect,
               replicate = replicate)
  
  return(g)
}

################################################################################
run_generations <- function(g, num.gens, generation, f_inherit, f_repro, f_mort, 
                            inheritance, rank.inher.precis, rank.effect){
  g$inheritance <- inheritance
  if(num.gens == generation){
    return(g)
    cat('Finished inheritance treatment: ', inheritance, '\n')
  }else{
    return(rbind(g, run_generations(generation_step(g, f_inherit, f_repro, f_mort, 
                                                    inheritance, rank.inher.precis, rank.effect), 
                                    num.gens, generation+1, f_inherit = f_inherit, 
                                    f_repro = f_repro, f_mort = f_mort, inheritance, 
                                    rank.inher.precis, rank.effect)))
  }
}

################################################################################
