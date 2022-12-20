################################################################################
mutate_color <- function(color){
  color <- t(col2rgb(color))
  nudge <- randomLHS(100, 3)
  nudge <- qnorm(nudge, mean = 0, sd = nudge.distance)
  color <- sweep(nudge, MARGIN = 2, color, "+")
  color <- color[apply(color, MARGIN = 1, function(x)(all(x >=0 & x < 255))),]
  return(rgb(color[1,1], color[1,2], color[1,3], maxColorValue = 255))
}

################################################################################
generation_step <- function(g){
  initial.group <- g
  reproduce <- seq(1, nrow(g), 1)[as.logical(rbinom(nrow(g), 1, reproduction.sigmoid(nrow(g),
                                                                                     seq(1, 1-1, length.out = nrow(g)))))]
  mothers <- g[reproduce,]
  die <- seq(1, nrow(g), 1)[as.logical(rbinom(nrow(g), 1, mortality.sigmoid(nrow(g))))]
  
  if(sum(reproduce)){
    new.ids <- data.frame(id = paste0(str_extract(g$id[reproduce],
                                                 '[:alpha:]'), 
                                      as.numeric(str_extract(g$id[reproduce],
                                                             '[:digit:]+'))+1, 
                                      '.',
                                      unique(g$generation+1)),
                          rank = g$rank[reproduce]+0.1,
                          matriline = g$matriline[reproduce],
                          color = sapply(g$color[reproduce], mutate_color, USE.NAMES = FALSE),
                          generation = unique(g$generation),
                          rank.inher.sd = rank.inher.sd,
                          rank.effect = rank.effect,
                          replicate = replicate,
                          mother = mothers$id)
    
    if(rank.inher.sd > 0){
      new.ranks <- rnorm(mean = new.ids$rank, n = length(new.ids$rank), sd = rank.inher.sd)
      while(any(new.ranks < 0.5 | new.ranks > (nrow(g) + nrow(new.ids) - length(die)))){
        new.ranks <- rnorm(mean = new.ids$rank, n = length(new.ids$rank), sd = rank.inher.sd)
      }
      new.ids$rank <- new.ranks
    }
    
    if(sum(die)){
      g <- arrange(rbind(g[-die,], new.ids), rank)
    }else{
      g <- arrange(rbind(g, new.ids), rank)
    }
  }else{
    if(sum(die)){
      g <- arrange(g[-die,], rank)
    }else{
      g <- g
    }
  }
  
  g$rank <- 1:nrow(g)
  g$generation <- g$generation + 1
  
  #### Record data into objects in global environment
  
  if(nrow(mothers))
    .GlobalEnv$mor.list[[length(.GlobalEnv$mor.list)+1]] <- data.frame(mother.rank = mothers$rank,
                                                                       offspring.rank = g[g$id %in% new.ids$id,]$rank,
                                                                       rank.inher.sd = rank.inher.sd,
                                                                       rank.effect = rank.effect,
                                                                       replicate = replicate,
                                                                       mother = mothers$id)
  
  
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
               rank.inher.sd = rank.inher.sd,
               rank.effect = rank.effect,
               replicate = replicate)
  
  return(g)
}

################################################################################
run_generations <- function(g, num.gens, generation){
  if(num.gens == generation){
    return(g)
  }else{
    return(rbind(g, run_generations(generation_step(g), num.gens, generation+1)))
  }
}

################################################################################
calc_genetic_contribution <- function(id){
  return(0.5 ^ as.numeric(str_extract(id, '[:digit:]+')))
}

calc_reproductive_value <- function(mother){
  
}


