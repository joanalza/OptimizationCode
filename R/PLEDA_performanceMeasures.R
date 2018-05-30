library(gtools)

#### PERFORMANCE MEASURE FUNCTIONS ####


#' @description Gets the average of the best-of-generation of an algorithm after several runs.
#' @param Data: data.frame with the best fitness of each run as column.
#' @note Drawback: the result is not normalized.
#'
get_best.of.gen <- function(data){
  best.of.runs <- apply(data, MARGIN = 1, function(row){
    return(min(row))
  })
  return(best.of.runs)
}

#' @description Provides convergence behaviour of algorithms, range(0,1).
#' @param Data: data that show the changes and best since last change.
#' @param Optimum: the global optimum at every change.
#' @note Take into account the generations and the changes.
#'
get_ARR <- function(data, optimum){
  pi <- nrow(data)
  fbest.1 <- data$bestPerChange[1]
  dividend <- sapply(data$bestPerChange, FUN = function(fij){
    return(fij - fbest.1)
  })
  return(sum(dividend)/(pi * (optimum - fbest.1)))
}

#' @description Based on the relative ratio between the best value gained by an algorithm and the global optimum value.
#' @param Data: data that show the changes and best since last change.
#' @param Optimum: global optimum of the search space.
#' @param Runs: total number of runs.
#' @param Num.changes: total number of changes.
#' @note Drawback: the result is not normalized.
#' @note Drawback: time of change.
#'
get_relative.ratio <- function(data, runs, num.changes, optimum){
  n.gen.per.change <- nrow(data)
  f.last <- tail(data,1)
  r.last <- abs(optimum / f.last$bestPerChange )
  sum.rij <- sapply(data$bestFit, FUN = function(fi.best){
    r.s <- abs(optimum / fi.best)
    return((1 - r.s)/n.gen.per.change)
  })
  r.ij <- r.last / (1 + sum(sum.rij))
  return(r.ij/(num.changes * runs))
}

#' @description Diversity of the population
#' @param Data: data that show the changes and best since last change.
#' @param Runs: total number of runs.
#' @param Num.changes: total number of changes.
#' @note Drawback: the result is not normalized.
#' @note Drawback: time of change.
#'
get_distance.based.diversity<- function(data, n){
  n.Best <- unique(data$nBestInd)
  norm.diversity <- sapply(data$diversity, FUN = function(row){
    return(row/((n - 1) * (n.Best *(n.Best - 1)/2)))
  })
  return(mean(norm.diversity))
}

## BEST-OF-GENERATION ####
F.best.of.gen <- function(files){
  df.list <- lapply(files, FUN = function(run){
    df <- read.csv(run)
    return(df$bestFit)
  })
  runs <- sub("*-.csv","",sub("([a-zA-Z0-9//_//.][-]?)*-","",files))
  df.bog <- as.data.frame(df.list)
  colnames(df.bog) <- runs
  f.bog <- get_best.of.gen(df.bog)
  F.bog <- mean(f.bog)
  return(F.bog)
}

## ABSOLUTE RECOVERY RATE ####
F.ARR <- function(files, optimums){
  df.list <- sapply(files, FUN = function(run){
    instance <- unlist(strsplit(run,"-"))[grep("tai", unlist(strsplit(run,"-")))]
    optima <- optimums[which(instance == optimums$Instance),"TotalFlowTime"] 
    df <- read.csv(run)
    n.changes <- max(unique(df$change))
    if (n.changes!=1){
      ARR <- sapply(2:n.changes, FUN = function(i){
        sub.df <- df[df$change==i,]
        return(get_ARR(sub.df, optima))
      })
    }else{
      sub.df <- df[df$change==n.changes,]
      return(get_ARR(sub.df, optima))
    }
    ARR.mean <- mean(ARR) # mean of the formula
    return(ARR.mean)
  })
}

## RELATIVE RATIO ####
F.relative.ratio <- function(files, optimums){
  n.runs <- length(files)
  df.list <- sapply(files, FUN = function(run){
    instance <- unlist(strsplit(run,"-"))[grep("tai", unlist(strsplit(run,"-")))]
    optima <- optimums[which(instance == optimums$Instance),"TotalFlowTime"] 
    df <- read.csv(run)
    n.changes <- max(unique(df$change))
    relative.ratio <- sapply(1:n.changes, FUN = function(i){
      sub.df <- df[df$change==i,]
      return(get_relative.ratio(sub.df, n.runs, n.changes,optima))
    })
    return(sum(relative.ratio))
  })
  return(sum(df.list))
}

## DISTANCE-BASED MEASURE ####
F.distance.based.measure <- function(files, n){
  df.list <- sapply(files, FUN = function(run){
    df <- read.csv(run)
    P.diversity <- get_distance.based.diversity(df,n)
    return(P.diversity)
  })
  return(mean(df.list))
}

## MAIN ####
directory <- "C:/Project/Code/OptimizationCode/results/cluster/tai50/tai50_5/"
# directory <- "C:/Project/Code/RKEDAC++/results/cluster/SEED_1/tai100_10_0/Cayley 90/"
setwd(directory)

# Get files
files.progress <- list.files(directory, pattern = "*.csv")
#files.best <- list.files(directory, pattern = "best*")
#files.improvement <- list.files(directory, pattern = "improvement*")

# Get parameters
rates <- unlist(unique(sapply(strsplit(files.progress,"-"), FUN = function(x){
  return(x[grep("ratio",x)])
})))

instances <- unlist(unique(sapply(strsplit(files.progress,"-"), FUN = function(x){
  return(x[grep("tai",x)])
})))

changes <- unlist(unique(sapply(strsplit(files.progress,"-"), FUN = function(x){
  return(x[grep("(c[0-9]+|noChange)",x)])
})))

cayley <- unlist(unique(sapply(strsplit(files.progress,"-"), FUN = function(x){
  return(x[grep("Cayley",x)])
})))

optimums <- read.csv("C:/Project/Code/RKEDAC++/optimums/optimum.csv")

# Get running files

rate.files <- lapply(rates, FUN = function(x){
  files <- list.files(directory, pattern = x)
  return(files[grep(".csv",files)])
})

instance.files <- lapply(instances, FUN = function(x){
  files <- list.files(directory, pattern = x)
  return(files[grep(".csv",files)])
})

changes.files <- lapply(changes, FUN = function(x){
  files <- list.files(directory, pattern = x)
  return(files[grep(".csv",files)])
})

cayley.files <- lapply(cayley, FUN = function(x){
  files <- list.files(directory, pattern = x)
  return(files[grep(".csv",files)])
})

cayley.changes.files <- unlist(lapply(changes, FUN = function(x){
  files <- list.files(directory, pattern = x)
  if (x!="noChange"){
    new.files <- lapply(cayley, FUN = function(Cay){
      files <- files[grep(Cay,files)]
      return(files[grep(".csv",files)])
    })
  }else{
    new.files <- list(files[grep(".csv",files)])
  }
  return(new.files)
}), recursive = F)

cayley.changes.files <- cayley.changes.files[lengths(cayley.changes.files) > 0L]

names.cayley.changes <- mixedsort(apply(expand.grid(changes, cayley),
                                        1, paste, collapse = "-", sep = ""))[-grep("noChange*",sort(apply(expand.grid(changes, cayley), 1, paste, collapse = "-", sep = "")))]
names.cayley.changes <- append(names.cayley.changes,"noChange")


global.files <- unlist(lapply(instances, FUN = function(x){
  files <- list.files(directory, pattern = x)
  new.files <- unlist(lapply(rates, FUN = function(rat){
      files <- files[grep(rat,files)]
      new.new.files <- unlist(lapply(changes, FUN = function(ch){
          files <- files[grep(ch,files)]
          new.new.new.files <- lapply(cayley, FUN = function(Cay){
              files <- files[grep(Cay,files)]
              return(files[grep(".csv",files)])
          })
          return(new.new.new.files)
          browser()
      }), recursive = F)
      return(new.new.files)
  }), recursive = F)
  return(new.files)
}), recursive = F)

global.files <- global.files[lengths(global.files) > 0L]

names.global <- sapply(global.files, FUN = function (x){
  result.name <- x[1]
  return(sub("*-PL.*","",sub("progress-*","",result.name)))
})

# Run performance measures

# best.of.gen <- sapply(rate.files, FUN = function(x){
# best.of.gen <- sapply(cayley.changes.files, FUN = function(x){
best.of.gen <- sapply(global.files, FUN = function(x){
  return(F.best.of.gen(x))
})

rank.bog <- rank(best.of.gen)

# mean.ARR <- sapply(rate.files, FUN = function(x){
# mean.ARR <- sapply(cayley.changes.files, FUN = function(x){
mean.ARR <- sapply(global.files, FUN = function(x){
  return(mean(F.ARR(x, optimums)))
})

rank.ARR <- rank(-mean.ARR) # DESCENDING ORDER

# relative.ratio <- sapply(rate.files, FUN = function(x){
# relative.ratio <- sapply(cayley.changes.files, FUN = function(x){
relative.ratio <- sapply(global.files, FUN = function(x){
  return(F.relative.ratio(x,optimums))
})

rank.rr <- rank(-relative.ratio)

df.ranks <- data.frame(rank.bog,rank.ARR, rank.rr)

average.rank <- apply(df.ranks,1,mean)

# Get results as csv
df.pm <- data.frame(best.of.gen, rank.bog, mean.ARR, rank.ARR, relative.ratio, rank.rr, average.rank)

# Order by the average of the ranks
df.pm <- df.pm[order(df.pm$average.rank),]
# row.names(df.pm) <- rates[as.numeric(row.names(df.pm))]
# row.names(df.pm) <- names.cayley.changes[as.numeric(row.names(df.pm))]
row.names(df.pm) <- names.global[as.numeric(row.names(df.pm))]
dir.name <- substr(getwd(),nchar(dirname(getwd()))+2,nchar(getwd()))

# Write he data frame on a .csv file
write.csv(df.pm, paste0("C:/Project/Code/OptimizationCode/performanceMeasures/cluster/PL-",dir.name,".csv") ,row.names = TRUE)
closeAllConnections()
