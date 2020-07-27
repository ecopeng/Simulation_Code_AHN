rm(list = ls())

Transmission <- function(AHN, N_agent, ProbStay, TimeSteps, pi, rho){
  library(igraph)
  SIR_His <- data.frame(ID_Agent = 1:N_agent, SIR_Initial = 0, stringsAsFactors = F)
  SIR_His$SIR_Initial[sample(1:N_agent, 1)] <- 1
  SIR_His$SIR_Temp <- SIR_His$SIR_Initial
  ifelse(is.weighted(AHN), 
         AHN_mat <- as_adjacency_matrix(AHN, sparse = F, attr = "weight"), 
         AHN_mat <- as_adjacency_matrix(AHN, sparse = F, attr = NULL))
  Move_His <- data.frame(ID_Agent = 1:N_agent, Site_Initial = sample(1:nrow(AHN_mat), N_agent, replace = T), stringsAsFactors = F)
  Move_His$Site_Temp <- Move_His$Site_Initial
  MoveProb <- sweep(AHN_mat, 1, (1-ProbStay)/rowSums(AHN_mat), "*")
  diag(MoveProb) <- ProbStay
  MoveProb[is.na(MoveProb)] <- 0
  for(u in 1:TimeSteps){
    Probs <- MoveProb[Move_His$Site_Temp, ]
    Move_His[, u+3] <- apply(Probs, 1, function(x){sample(1:nrow(AHN_mat), 1, prob = x)})
    names(Move_His)[u+3] <- u
    Move_His$Site_Temp <- Move_His[, u+3]
    Temp <- outer(Move_His[, u+3], Move_His[, u+3], "==")*1
    diag(Temp) <- 0
    Contacts <- graph_from_adjacency_matrix(Temp, mode = "undirected")
    ID_I <- which(SIR_His$SIR_Temp == 1)
    if(length(ID_I) != 0){
      for(v in 1:length(ID_I)){
        Nei_all <- neighbors(Contacts, ID_I[v])
        if(length(Nei_all) != 0){
          for(w in 1:length(Nei_all)){
            if(SIR_His$SIR_Temp[Nei_all[w]] == 0){
              SIR_His$SIR_Temp[Nei_all[w]] <- sample(c(1, 0), 1, prob = c(pi, 1 - pi))
            }
          }
        }
        SIR_His$SIR_Temp[ID_I[v]] <- sample(c(2, 1), 1, prob = c(rho, 1 - rho))
      }
      SIR_His[, u+3] <- SIR_His$SIR_Temp
      names(SIR_His)[u+3] <- u
    }else{
      SIR_His[, u+3] <- SIR_His$SIR_Temp
      names(SIR_His)[u+3] <- u
    }
  }
  Rate <- SIR_His[, -c(1, 3)]
  return(data.frame(Time = 0:(ncol(Rate) - 1), 
                    Proportion = unname(apply(Rate, 2, function(x){sum(x == 1)}))/nrow(Rate),
                    Density = graph.density(AHN),
                    Clustering = transitivity(AHN, type = "average", isolates = "zero"), 
                    Modularity = modularity(fastgreedy.community(AHN)), 
                    Diameter = diameter(AHN, directed = NULL, unconnected = T), 
                    stringsAsFactors = F))
}


N <- 20
L <- c(5, 50)
mu <- 5
la <- c(.001, 30)
rep <- 100
library(igraph)
library(AnimalHabitatNetwork)
ahn <- list()
h <- 1
for(i in 1:rep){
  for(j in 1:length(L)){
    xcood <- runif(N, 0, L[j])
    ycood <- runif(N, 0, 25/L[j])
    for(k in 1:length(la)){
      ahn[[h]] <- list(ahn = ahn_gen(N, L[j], mu = mu, lamda = la[k], 
                                     Connected = T, Weighted = T, eta = 1, 
                                     X = xcood, Y = ycood), 
                       metadata = data.frame(rep = i, N = N, L = L[j], mu = mu, lambda = la[k]))
      h <- h + 1
    }
  }
}

save(ahn, file = "ahn.rdata")

N_agent <- 100
ProbStay <- .95
TimeSteps <- 500
pi <- .05
rho <- .01

n_simu <- length(ahn)
library(doParallel)
nc <- detectCores() - 1
cl <- makeCluster(nc, outfile = "")
registerDoParallel(cl)

foreach (u = 1:n_simu)%dopar%{
  res <- data.frame()
  res <- Transmission(ahn[[u]]$ahn, N_agent = N_agent, ProbStay = ProbStay, TimeSteps = TimeSteps, pi = pi, rho = rho)
  res <- cbind(res, ahn[[u]]$metadata)
  save(res, file = paste("SIR", "_rep_", u, ".rdata", sep = ""))
  cat(u, "is Done!\n")
}
stopCluster(cl)