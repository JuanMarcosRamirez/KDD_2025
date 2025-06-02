library(igraph)
library(Matrix)
library(pracma)

set.seed(8192)
setwd("~/Documents/Research/NSUM/SyntheticDatasets")

N   <- 1e6
d   <- c(50)
md  <- d/(N-1)
rho <- 5

S <- c(100, 200, 500, 800, 1000, 2000, 5000, 8000, 10000, 20000, 50000, 80000, 100000)

for (i in 1:length(d))
{
	file_MoR <- paste0("erdos_renyi/results/PMoR_vs_epsilon_deg_",as.character(d[i]),"_rho_",as.character(rho),".rds")
	Pe_MoR <- readRDS(file_MoR)
	
	file_RoS <- paste0("erdos_renyi/results/PRoS_vs_epsilon_deg_",as.character(d[i]),"_rho_",as.character(rho),".rds")
	Pe_RoS <- readRDS(file_RoS)
	
	file_RoS <- paste0("erdos_renyi/results/PRoSB_vs_epsilon_deg_",as.character(d[i]),"_rho_",as.character(rho),".rds")
	Pe_RoSB <- readRDS(file_RoS)
		
	yQuery <- 0.05
	sp <- 10
	ep <- 103  # 103, 131, 187
	y1 <- matrix(0,ep-sp+1,1)
	x1 <- matrix(0,ep-sp+1,1)
	
	y2 <- matrix(0,ep-sp+1,1)
	x2 <- matrix(0,ep-sp+1,1)
	
	y3 <- matrix(0,ep-sp+1,1)
	x3 <- matrix(0,ep-sp+1,1)
	for (i in sp:ep) {
    	print(i)
  		Y <- Pe_MoR[, i]
  		unique_Y <- Y[!duplicated(Y)]   # 'stable' option in R to preserve order
  		unique_X <- S[!duplicated(Y)]
    	
  		# Interpolation using 'pchip' (piecewise cubic Hermite interpolation)
  		y1[i-(sp-1)] <- interp1(unique_Y, unique_X, yQuery, method = "pchip")
  		x1[i-(sp-1)] <- i / 1000
  		
  		Y1 <- Pe_RoS[, i]
  		unique_Y1 <- Y1[!duplicated(Y1)]   # 'stable' option in R to preserve order
  		unique_X1 <- S[!duplicated(Y1)]
    	    	
  		# Interpolation using 'pchip' (piecewise cubic Hermite interpolation)
  		y2[i-(sp-1)] <- interp1(unique_Y1, unique_X1, yQuery, method = "pchip")
  		x2[i-(sp-1)] <- i / 1000
  		
  		Y2 <- Pe_RoSB[, i]
  		unique_Y2 <- Y2[!duplicated(Y1)]   # 'stable' option in R to preserve order
  		unique_X2 <- S[!duplicated(Y1)]
    	    	
  		# Interpolation using 'pchip' (piecewise cubic Hermite interpolation)
  		y3[i-(sp-1)] <- interp1(unique_Y2, unique_X2, yQuery, method = "pchip")
  		x3[i-(sp-1)] <- i / 1000
	}
#	print(x2)
#	print(y2)
	
	step <- 20
	x11 <- matrix(0,length(seq(1,ep-sp+1,step)),1)
	y11 <- matrix(0,length(seq(1,ep-sp+1,step)),1)
	x21 <- matrix(0,length(seq(1,ep-sp+1,step)),1)
	y21 <- matrix(0,length(seq(1,ep-sp+1,step)),1)
	x31 <- matrix(0,length(seq(1,ep-sp+1,step)),1)
	y31 <- matrix(0,length(seq(1,ep-sp+1,step)),1)
	
	jj <- 1
	for (j in seq(1,ep-sp+1,step))
	{
		print(j)
		x11[jj] <- x1[j]
		y11[jj] <- y1[j]
		x21[jj] <- x2[j]
		y21[jj] <- y2[j]
		x31[jj] <- x3[j]
		y31[jj] <- y3[j]
		jj <- jj + 1
	}
	
	eMoR <- data.frame(matrix(ncol = 0, nrow = nrow(x11)))
	eMoR$size <- y11
	eMoR$epsilon <- x11
	write.table(eMoR, paste0("erdos_renyi/results/size_vs_epsilonMoR_deg_",as.character(d),"_rho_",as.character(rho),".dat"),sep=" ",quote=FALSE,row.names = FALSE)
	
	eRoS <- data.frame(matrix(ncol = 0, nrow = nrow(x21)))
	eRoS$size <- y21
	eRoS$epsilon <- x21
	write.table(eRoS, paste0("erdos_renyi/results/size_vs_epsilonRoS_deg_",as.character(d),"_rho_",as.character(rho),".dat"),sep=" ",quote=FALSE,row.names = FALSE)
	
	eRoSB <- data.frame(matrix(ncol = 0, nrow = nrow(x31)))
	eRoSB$size <- y31
	eRoSB$epsilon <- x31
	write.table(eRoSB, paste0("erdos_renyi/results/size_vs_epsilonRoSB_deg_",as.character(d),"_rho_",as.character(rho),".dat"),sep=" ",quote=FALSE,row.names = FALSE)
}
