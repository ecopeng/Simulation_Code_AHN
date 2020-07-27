rm(list = ls())

library(igraph)
# approximation functions AHN
ahn_gen_unwei <- function(N, L, mu, lamda, A = 25, Connected = TRUE){
  dm <- as.matrix(dist(data.frame(x = runif(N, 0, L), y = runif(N, 0, A/L))), method = 'euclidean', diag = FALSE)
  dm_0 <- 1/dm
  dm_0[is.infinite(dm_0)] <- 0
  mat_bin <- dm_0
  mat_bin[mat_bin > 0] <- 1
  mat_bin[lower.tri(mat_bin, diag = TRUE)] <- NA
  tr <- which(!is.na(mat_bin))
  prob <- 1/(1 + exp(-lamda*(dm[tr] - mu)))
  for(t in 1:length(tr)){
    if(sample(c(1, 0), size = 1, prob = c(prob[t], 1 - prob[t]))){
      mat_bin[tr][t] <- 0
    }
  }
  mat_bin[lower.tri(mat_bin)] <- t(mat_bin)[lower.tri(mat_bin)]
  ahn <- graph_from_adjacency_matrix(mat_bin, mode = 'undirected', diag = FALSE, weighted = NULL)
  if(!is.connected(ahn) && Connected){
    memb <- unname(components(ahn)$membership)
    while(max(memb) > 1){
      r_memb <- sample(memb, 1)
      temp <- dm_0[which(memb == r_memb), which(memb != r_memb), drop = FALSE]
      rn <- as.numeric(rownames(temp)[which(temp == max(temp), arr.ind = T)[1]])
      cn <- as.numeric(colnames(temp)[which(temp == max(temp), arr.ind = T)[2]])
      mat_bin[rn, cn] <- 1
      mat_bin[cn, rn] <- 1
      ahn <- graph_from_adjacency_matrix(mat_bin, mode = 'undirected', diag = FALSE, weighted = NULL)
      memb <- unname(components(ahn)$membership)
    }
  }
  return(ahn)
}

# mod
ahn_approx_mod <- function(N, L = seq(5, 30, length.out = 6), mu = c(.1, 2, 5, 7, 10), lamda = c(.001, .1, .15, .35, .4, .75, 1.25, 5, 30), ExpectedValue){
  dif <- 2 # [-1, 1]
  for(u in L){
    for(v in mu){
      for(w in lamda){
        temp <- ahn_gen_unwei(N, u, v, w)
        mod <- modularity(fastgreedy.community(temp))
        d <- abs(mod - ExpectedValue)
        if(dif >= d){
          dif <- d
          ahn_temp <- temp
        }
      }
    }
  }
  return(ahn_temp)
}

# clu
ahn_approx_clu <- function(N, L = seq(5, 30, length.out = 6), mu = c(.1, 2, 5, 7, 10), lamda = c(.001, .1, .15, .35, .4, .75, 1.25, 5, 30), ExpectedValue){
  dif <- 1 # [0, 1]
  for(u in L){
    for(v in mu){
      for(w in lamda){
        temp <- ahn_gen_unwei(N, u, v, w)
        clu <- transitivity(temp, type = "average", isolates = "zero")
        d <- abs(clu - ExpectedValue)
        if(dif > d){
          dif <- d
          ahn_temp <- temp
        }
      }
    }
  }
  return(ahn_temp)
}

# dia
ahn_approx_dia <- function(N, L = seq(5, 30, length.out = 6), mu = c(.1, 2, 5, 7, 10), lamda = c(.001, .1, .15, .35, .4, .75, 1.25, 5, 30), ExpectedValue){
  dif <- N # [0, N-1]
  for(u in L){
    for(v in mu){
      for(w in lamda){
        temp <- ahn_gen_unwei(N, u, v, w)
        dia <- diameter(temp, directed = NULL, unconnected = T)
        d <- abs(dia - ExpectedValue)
        if(dif > d){
          dif <- d
          ahn_temp <- temp
        }
      }
    }
  }
  return(ahn_temp)
}

# simulations
rep <- 15
N_space <- 270 # length(L)*length(mu)*length(lamda) from approximation functions
load("Friesen_et_al.rdata")

library(doParallel)
nc <- detectCores() - 2
cl <- makeCluster(nc, outfile = "")
registerDoParallel(cl)
foreach (u = 1:length(Friesen_et_al), .combine = rbind, .packages = c("doParallel", "igraph")) %dopar% {
  network <- Friesen_et_al[[u]]
  nv <- vcount(network)
  EMP_tra <- transitivity(network, type = "average", isolates = "zero")
  EMP_mod <- modularity(fastgreedy.community(network))
  EMP_dia <- diameter(network, directed = NULL, unconnected = T)
  AHN <- list()
  for(i in 1:rep){
    ahn_net_tra <- ahn_approx_clu(N = nv, ExpectedValue = EMP_tra)
    AHN_tra <- transitivity(ahn_net_tra, type = "average", isolates = "zero")
    ahn_net_mod <- ahn_approx_mod(N = nv, ExpectedValue = EMP_mod)
    AHN_mod <- modularity(fastgreedy.community(ahn_net_mod))
    ahn_net_dia <- ahn_approx_dia(N = nv, ExpectedValue = EMP_dia)
    AHN_dia <- diameter(ahn_net_dia, directed = NULL, unconnected = T)
    dMetrics <- data.frame(iNet_EMP = u, Network = "AHN", dClu = AHN_tra - EMP_tra, dMod = AHN_mod - EMP_mod, dDia = AHN_dia - EMP_dia)
    Networks <- list(AHN_net_tra = ahn_net_tra, AHN_net_mod = ahn_net_mod, AHN_net_dia = ahn_net_dia)
    AHN[[i]] <- list(dMetrics = dMetrics, Networks = Networks)
  }
  save(AHN, file = paste("Tests_", u, ".rdata", sep = ""))
}
stopCluster(cl)