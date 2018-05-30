
# Put instance results in the directory
directory <- "C:/Project/Code/OptimizationCode/results/cluster/ratio/tai100_10_0/"
setwd("C:/Project/Code/OptimizationCode/plots/cluster/ratio/tai100_10_0/")

# Get instance names
instances <- list.files(directory, pattern = "*.csv")

# Get optimums
optimums <- read.csv("C:/Project/Code/RKEDAC++/optimums/optimum.csv")

apply(as.array(instances), MARGIN = 1, FUN = function(instance){
  df <- read.csv(paste0(directory,"/",instance), header = TRUE)
  # browser()
  name <- unlist(strsplit(instance, "-"))
  name <- name[grep("tai",name)]
  optima <- optimums[which(name == optimums$Instance),"TotalFlowTime"] 
  #file.rename(from = paste0(directory,"/",instance), to = paste0(directory,"/ploted/",instance))
  
  pdf(paste0(sub('\\.csv$', '',instance),".pdf"))
  layout(matrix(c(1,2), ncol = 1))
  
  # Progress
  plot(x =df$gen, y = df$bestPerChange, 
       cex.lab=0.8, ylab = "Fitness", xlab = "Generation", type = "l", 
       #xlim = c(0,50000), 
       ylim = c(min(df$bestFound), max(df$avgFit)*1.1)
       )
  title(main = paste0("DPL-EDA over ", name))
  lines(x = df$gen, y = df$avgFit, lty=2, col = "red")
  lines(x = df$gen, y = df$bestFound, lty=4, col = "skyblue3")
  abline(h= optima, col="green")
  legend("topright", legend=c("Current fitness", "Average fitness", "Global optimum", "Best found"),
         col=c("black", "red", "green","skyblue3"), lty=c(1,2,1,4), cex=0.5)
  
  # COOLING SCHEME
  plot(x = df$gen, y = df$rate, 
       col = "blue", cex.lab=0.8, ylab = "Boltzmann distribution", xlab = "Generation", type = "l",
       #xlim = c(0,50000),
       ylim = c(min(df$rate), max(df$rate))
       )
  title(main = "Variance parameter (temperature) during the execution")
  legend("topright", legend= "Temperature",
         col= "blue",lty = 1, cex=0.5)
  dev.off()
  return ()
})