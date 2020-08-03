rm(list = ls())
library(igraph)
library(AnimalHabitatNetwork)
# functions
ahn_approx_mod <- function(N, 
                           L = seq(5, 30, length.out = 6), 
                           mu = c(.1, 2, 5, 7, 10), 
                           lamda = c(.001, .1, .15, .35, .4, .75, 1.25, 5, 30), 
                           ExpectedValue){
  dif <- 2 # [-1, 1]
  for(u in L){
    for(v in mu){
      for(w in lamda){
        temp <- ahn_gen(N = N, L = u, mu = v, lamda = w, Connected = T, Weighted = F)
        mod <- modularity(fastgreedy.community(temp))
        d <- abs(mod - ExpectedValue)
        if(dif >= d){
          dif <- d
          pars <- data.frame(N = N, L = u, mu = v, lamda = w)
        }
      }
    }
  }
  return(pars)
}


ahn_approx_clu <- function(N, 
                           L = seq(5, 30, length.out = 6), 
                           mu = c(.1, 2, 5, 7, 10), 
                           lamda = c(.001, .1, .15, .35, .4, .75, 1.25, 5, 30), 
                           ExpectedValue){
  dif <- 1 # [0, 1]
  for(u in L){
    for(v in mu){
      for(w in lamda){
        temp <- ahn_gen(N = N, L = u, mu = v, lamda = w, Connected = T, Weighted = F)
        clu <- transitivity(temp, type = "average", isolates = "zero")
        d <- abs(clu - ExpectedValue)
        if(dif >= d){
          dif <- d
          pars <- data.frame(N = N, L = u, mu = v, lamda = w)
        }
      }
    }
  }
  return(pars)
}


ahn_approx_dia <- function(N, 
                           L = seq(5, 30, length.out = 6), 
                           mu = c(.1, 2, 5, 7, 10), 
                           lamda = c(.001, .1, .15, .35, .4, .75, 1.25, 5, 30), 
                           ExpectedValue){
  dif <- N # [0, N-1]
  for(u in L){
    for(v in mu){
      for(w in lamda){
        temp <- ahn_gen(N = N, L = u, mu = v, lamda = w, Connected = T, Weighted = F)
        dia <- diameter(temp, directed = NULL, unconnected = T)
        d <- abs(dia - ExpectedValue)
        if(dif >= d){
          dif <- d
          pars <- data.frame(N = N, L = u, mu = v, lamda = w)
        }
      }
    }
  }
  return(pars)
}

# simulation
rep <- 15
load("Friesen_et_al.rdata") # 60 networks in the dataset and 58 are used

library(doParallel)
nc <- detectCores() - 1
cl <- makeCluster(nc, outfile = "")
registerDoParallel(cl)
foreach (u = 1:length(Friesen_et_al))%dopar%{
  library(igraph)
  library(AnimalHabitatNetwork)
  network <- Friesen_et_al[[u]]
  nv <- vcount(network)
  EMP_tra <- transitivity(network, type = "average", isolates = "zero")
  EMP_mod <- modularity(fastgreedy.community(network))
  EMP_dia <- diameter(network, directed = NULL, unconnected = T)
  AHN <- list()
  for(i in 1:rep){
    # clustering coef
    par_tra <- ahn_approx_clu(N = nv, ExpectedValue = EMP_tra)
    ahn_net_tra <- ahn_gen(N = par_tra[1, 1], L = par_tra[1, 2], mu = par_tra[1, 3], lamda = par_tra[1, 4], Connected = T, Weighted = F)
    AHN_tra <- transitivity(ahn_net_tra, type = "average", isolates = "zero")
    # modularity
    par_mod <- ahn_approx_mod(N = nv, ExpectedValue = EMP_mod)
    ahn_net_mod <- ahn_gen(N = par_mod[1, 1], L = par_mod[1, 2], mu = par_mod[1, 3], lamda = par_mod[1, 4], Connected = T, Weighted = F)
    AHN_mod <- modularity(fastgreedy.community(ahn_net_mod))
    # diameter
    par_dia <- ahn_approx_dia(N = nv, ExpectedValue = EMP_dia)
    ahn_net_dia <- ahn_gen(N = par_dia[1, 1], L = par_dia[1, 2], mu = par_dia[1, 3], lamda = par_dia[1, 4], Connected = T, Weighted = F)
    AHN_dia <- diameter(ahn_net_dia, directed = NULL, unconnected = T)
    
    dMetrics <- data.frame(iNet_EMP = u, rep = i, dClu = AHN_tra - EMP_tra, dMod = AHN_mod - EMP_mod, dDia = AHN_dia - EMP_dia)
    Networks <- list(iNet_EMP = u, rep = i, AHN_net_tra = ahn_net_tra, AHN_net_mod = ahn_net_mod, AHN_net_dia = ahn_net_dia)
    AHN[[i]] <- list(dMetrics = dMetrics, Networks = Networks)
  }
  save(AHN, file = paste("Tests_Best_Fit_Pars_", u, ".rdata", sep = ""))
  cat(u, "is Done!\n")
}
stopCluster(cl)