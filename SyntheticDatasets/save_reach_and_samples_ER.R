library(igraph)
library(Matrix)
library(pracma)

set.seed(8192)
setwd("~/Documents/Research/NSUM/SyntheticDatasets")

N   <- 1e6
d   <- 15
md  <- d/(N-1)
rho <- 1

num_nets <- 10
num_hpop <- 10

for (i in 1:num_nets){
  print(i)
  graph <- sample_gnp(N, md, directed=TRUE)
  GM    <- as_adjacency_matrix(graph)
  rm(graph)
  saveRDS(colSums(GM), paste0("erdos_renyi/synthetic_data/mean_degree_",as.character(d),"/rho_",as.character(rho),"/G",as.character(i),".rds"))
  for (j in 1:num_hpop){
    indices <- randperm(N,round((rho/100)*N))
    saveRDS(colSums(GM[indices,]), paste0("erdos_renyi/synthetic_data/mean_degree_",as.character(d),"/rho_",as.character(rho),"/H",as.character(i),"_",as.character(j),".rds"))
    rm(indices)
  }
  rm(GM)
}
