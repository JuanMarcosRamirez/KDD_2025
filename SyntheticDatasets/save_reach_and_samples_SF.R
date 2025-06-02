library(igraph)
library(Matrix)
library(pracma)

set.seed(8192)
setwd("~/Documents/Research/NSUM/SyntheticDatasets")

N   <- 1e6
d   <- 30
rho <- 30

num_nets <- 10
num_hpop <- 10

for (i in 1:num_nets){
  print(i)
  graph <- sample_fitness_pl(N, 0.5*N*d, 2.2, 2.3)
  GM    <- as_adjacency_matrix(graph)
  rm(graph)
  saveRDS(colSums(GM), paste0("scale_free/synthetic_data/rho_",as.character(rho),"/G",as.character(i),".rds"))
  for (j in 1:num_hpop){
    indices <- randperm(N,round((rho/100)*N))
    saveRDS(colSums(GM[indices,]), paste0("scale_free/synthetic_data/rho_",as.character(rho),"/H",as.character(i),"_",as.character(j),".rds"))
    rm(indices)
  }
  rm(GM)
}
