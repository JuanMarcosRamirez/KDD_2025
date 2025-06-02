library(igraph)
library(Matrix)
library(pracma)

set.seed(8192)
setwd("~/Documents/Research/NSUM/SyntheticDatasets")

N   <- 1e6
d   <- 50
md  <- d/(N-1)

S <- c(100, 200, 500, 800, 1000, 2000, 5000, 8000, 10000, 20000, 50000, 80000, 100000)
rhos <- c(1, 2, 5, 10, 15, 20)


Pe_MoR = matrix(0,length(S),length(rhos))
Pe_MoRB = matrix(0,length(S),length(rhos))
Pe_RoS = matrix(0,length(S),length(rhos))
Pe_RoSB = matrix(0,length(S),length(rhos))

for (i in 1:length(rhos))
{
#	file_MoR <- paste0("scale_free/results/PMoR_vs_epsilon_rho_",as.character(rhos[i]),".rds")
#	Pe_MoRt <- readRDS(file_MoR)
#	Pe_MoR[,i] <- Pe_MoRt[,50]
	
#	file_MoRB <- paste0("scale_free/results/PMoRB_vs_epsilon_rho_",as.character(rhos[i]),".rds")
#	Pe_MoRBt <- readRDS(file_MoRB)
#	Pe_MoRB[,i] <- Pe_MoRBt[,50]
	
#	file_RoS <- paste0("scale_free/results/PRoS_vs_epsilon_rho_",as.character(rhos[i]),".rds")
#	Pe_RoSt <- readRDS(file_RoS)
#	Pe_RoS[,i] <- Pe_RoSt[,50]
#	
	file_RoSB <- paste0("scale_free/results/PRoSB_vs_epsilon_rho_",as.character(rhos[i]),".rds")
	Pe_RoSBt <- readRDS(file_RoSB)
	Pe_RoSB[,i] <- Pe_RoSBt[,50]
}

yQuery <- 0.05
sp <- 2
ep <- 6  # 103, 131, 187
y1 <- matrix(0,ep-sp+1,1)
x1 <- matrix(0,ep-sp+1,1)
	
#y2 <- matrix(0,ep-sp+1,1)
#x2 <- matrix(0,ep-sp+1,1)
#	
#y3 <- matrix(0,ep-sp+1,1)
#x3 <- matrix(0,ep-sp+1,1)
for (i in sp:ep) {
    print(i)
  	Y <- Pe_RoSB[, i]
  	unique_Y <- Y[!duplicated(Y)]   # 'stable' option in R to preserve order
  	unique_X <- S[!duplicated(Y)]
    	
  	# Interpolation using 'pchip' (piecewise cubic Hermite interpolation)
  	y1[i-(sp-1)] <- interp1(unique_Y, unique_X, yQuery, method = "pchip")
  	x1[i-(sp-1)] <- rhos[i] / 100
  		
#  	Y1 <- Pe_RoS[, i]
#  	unique_Y1 <- Y1[!duplicated(Y1)]   # 'stable' option in R to preserve order
#  	unique_X1 <- S[!duplicated(Y1)]
#    	    	
#  	# Interpolation using 'pchip' (piecewise cubic Hermite interpolation)
#  	y2[i-(sp-1)] <- interp1(unique_Y1, unique_X1, yQuery, method = "pchip")
#  	x2[i-(sp-1)] <- rhos[i] / 100
#  		
#  	Y2 <- Pe_RoSB[, i]
#  	unique_Y2 <- Y2[!duplicated(Y1)]   # 'stable' option in R to preserve order
#  	unique_X2 <- S[!duplicated(Y1)]
#    	    	
#  	# Interpolation using 'pchip' (piecewise cubic Hermite interpolation)
#  	y3[i-(sp-1)] <- interp1(unique_Y2, unique_X2, yQuery, method = "pchip")
#  	x3[i-(sp-1)] <- rhos[i] / 100
}

#eMoR <- data.frame(matrix(ncol = 0, nrow = nrow(x1)))
#eMoR$size <- y1
#eMoR$rho <- x1
#write.table(eMoR, paste0("scale_free/results/size_vs_rhoMoR.dat"),sep=" ",quote=FALSE,row.names = FALSE)

#eMoRB <- data.frame(matrix(ncol = 0, nrow = nrow(x1)))
#eMoRB$size <- y1
#eMoRB$rho <- x1
#write.table(eMoRB, paste0("scale_free/results/size_vs_rhoMoRB.dat"),sep=" ",quote=FALSE,row.names = FALSE)
	
#eRoS <- data.frame(matrix(ncol = 0, nrow = nrow(x1)))
#eRoS$size <- y1
#eRoS$rho <- x1
#write.table(eRoS, paste0("scale_free/results/size_vs_rhoRoS.dat"),sep=" ",quote=FALSE,row.names = FALSE)
#	
eRoSB <- data.frame(matrix(ncol = 0, nrow = nrow(x1)))
eRoSB$size <- y1
eRoSB$rho <- x1
write.table(eRoSB, paste0("scale_free/results/size_vs_rhoRoSB.dat"),sep=" ",quote=FALSE,row.names = FALSE)
