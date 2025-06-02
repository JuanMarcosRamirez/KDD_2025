library(pracma)

setwd("~/Documents/Research/NSUM/SyntheticDatasets")

N   <- 1e6
d   <- 50
md  <- d/(N-1)
rho <- 15

num_nets <- 10
num_hpop <- 10

trials <- 200
#S <- c(100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000)
S <- c(100, 200, 500, 800, 1000, 2000, 5000, 8000, 10000, 20000, 50000, 80000, 100000)

e <- data.frame(matrix(ncol = 0, nrow = length(S)))
Pe <- data.frame(matrix(ncol = 0, nrow = length(S)))

eMoR <- matrix(0,trials*num_hpop*num_nets,length(S))
eRoS <- matrix(0,trials*num_hpop*num_nets,length(S))
Rs <- matrix(0,trials*num_hpop*num_nets,length(S))

for (l in 1:length(S))
{ 
	for (i in 1:num_nets)
	{
		file_G <- paste0("erdos_renyi/synthetic_data/mean_degree_",as.character(d),"/rho_",as.character(rho),"/G", as.character(i),".rds")
		G <- readRDS(file_G)
		for (j in 1:num_hpop)
		{	
			print(paste0("Size:",as.character(S[l]),". Network: ",as.character(i),". Unknown Population: ", as.character(j)))
			file_H <- paste0("erdos_renyi/synthetic_data/mean_degree_",as.character(d),"/rho_",as.character(rho),"/H", as.character(i),"_",as.character(j),".rds")
			H <- readRDS(file_H)
			H1 <- H[G>0]
			G1 <- G[G>0]
			N <- length(H1)
			for (k in 1:trials)
			{
				indices <- randperm(length(G1),S[l])
				Hi = H1[indices]
				Gi = G1[indices]
				MoR = mean(Hi/Gi)
				RoS = sum(Hi) / sum(Gi)
				eMoR[(i-1)*num_hpop*trials + (j-1)*trials + k,l] = max(c(1,MoR/(rho/100)))
				eRoS[(i-1)*num_hpop*trials + (j-1)*trials + k,l] = max(c(1,RoS/(rho/100)))
				Rs[(i-1)*num_hpop*trials + (j-1)*trials + k,l] = sum(Gi)
			}
		}
	}
}

e$eMoR = colMeans(eMoR)
e$eRoS = colMeans(eRoS)
e$sizes = S

boxplot_eMoR <- apply(eMoR, 2, function(col) quantile(col, probs = c(0, 0.25, 0.5, 0.75, 1)))
boxplot_eRoS <- apply(eRoS, 2, function(col) quantile(col, probs = c(0, 0.25, 0.5, 0.75, 1)))

print(boxplot_eMoR)
print(boxplot_eRoS)

file_errors = paste0("erdos_renyi/results/mean_eMoR_deg_",as.character(d),"_rho_",as.character(rho),".dat")
write.table(e,file_errors,sep=" ",quote=FALSE,row.names = FALSE)

maxeps <- 200
alpha  <- 0.50

PeMoR <- matrix(0,length(S),maxeps)
PeRoS <- matrix(0,length(S),maxeps)
PeMoR_bound <- matrix(0,length(S),maxeps)
PeRoS_bound <- matrix(0,length(S),maxeps)
Sample_Size <- matrix(0,maxeps,1)

ones <- matrix(1,length(S),1)
for (epsilon in 1:maxeps)
{
	beta = 1 + epsilon/1000
	PeMoR[,epsilon] = t(colMeans(eMoR > 1 + epsilon/1000))
	PeRoS[,epsilon] = t(colMeans(eRoS > 1 + epsilon/1000))
	PeMoR_bound[,epsilon] = pmin(ones, (exp(beta-1)/(beta^beta))^(S*(rho/100)) + (exp((1/beta) - 1)/(beta^(-1/beta)))^(S*(rho/100)))
	Sample_Size[epsilon] = (log(2) + alpha * log(N)) / ((rho/100) * (1 - (1/beta)*(log(beta) + 1)))
	
	for (l in 1:length(S))
	{ 
#		ht<-histc(Rs[,l], seq(min(Rs[,l]), max(Rs[,l])+1, 1))
#		Rv_temp <- ht$bin
#		Rv = Rv_temp[1:length(Rv_temp)-1]
#		Prob <- ht$cnt/(num_nets*num_hpop*trials)
		
		ht<-hist(Rs[,l], seq(min(Rs[,l]), max(Rs[,l])+1, 1), plot=FALSE)
		Rv_temp <- ht$breaks
		Rv = Rv_temp[1:length(Rv_temp)-1]
		Prob <- ht$counts/(num_nets*num_hpop*trials)
#		print(length(Rv))
#		print(length(Prob))
		
		BRoS <- sum(((exp(beta-1)/(beta^beta))^(Rv*(rho/100)) + (exp((1/beta)-1)/(beta^(-1/beta)))^(Rv*(rho/100)) )*Prob)
		PeRoS_bound[l,epsilon] <- min(c(1,BRoS))
	}	
}

eps_search = 50
Pe$MoR = PeMoR[,eps_search]
Pe$RoS = PeRoS[,eps_search]
Pe$MoRB = PeMoR_bound[,eps_search]
Pe$RoSB = PeRoS_bound[,eps_search]
Pe$Samp = Sample_Size[eps_search] * matrix(1,length(S),1)
Pe$S = S
file_Perrors = paste0("erdos_renyi/results/Pe_deg_",as.character(d),"_rho_",as.character(rho),"_epsilon_",as.character(eps_search),".dat")
write.table(Pe,file_Perrors,sep=" ",quote=FALSE,row.names = FALSE)

saveRDS(PeMoR, paste0("erdos_renyi/results/PMoR_vs_epsilon_deg_",as.character(d),"_rho_",as.character(rho),".rds"))
saveRDS(PeRoS, paste0("erdos_renyi/results/PRoS_vs_epsilon_deg_",as.character(d),"_rho_",as.character(rho),".rds"))
saveRDS(PeMoR_bound, paste0("erdos_renyi/results/PMoRB_vs_epsilon_deg_",as.character(d),"_rho_",as.character(rho),".rds"))
saveRDS(PeRoS_bound, paste0("erdos_renyi/results/PRoSB_vs_epsilon_deg_",as.character(d),"_rho_",as.character(rho),".rds"))
