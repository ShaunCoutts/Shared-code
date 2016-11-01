# R script for pub-R sheffield
#Does a simple phylogenetic and spatially lagged model

#need a few libraries 
library(GEOmap) # gets the spatial distance 
library(ape) # reads in the tree
# libraies for runing jags and looking at output
library(coda)
library(R2WinBUGS)
#install.packages(pkgs = 'R2WinBUGS', repos = 'http://cran.us.r-project.org')
library(R2jags)
#install.packages(pkgs = 'R2jags', repos = 'http://cran.us.r-project.org')

setwd('/home/shauncoutts/Desktop/git_repos/sheffield_R')

### NEED TO MAKE A FEW BITS AND PIECIES FIRST
# get the spatial and phylogenetic distance matricies 

# First get the tree
Tree <- read.nexus("PMDbl260314A1.nex")
plot(Tree)

# Also get the data
shef_dat <- read.csv('sheffield_R_data.csv', header = TRUE, stringsAsFactors = FALSE)
summary(shef_dat)

#caclulate distance to every other population, both geographic and phylogenetic. 
geo_dist <- matrix(nrow = dim(shef_dat)[1], ncol = dim(shef_dat)[1] - 1) #minus 1 becasue distance of each pop to itself is 0
phy_dist <- matrix(nrow = dim(shef_dat)[1], ncol = dim(shef_dat)[1] - 1) #minus 1 becasue distance of each pop to itself is 0_dist <- matrix(nrow = length(sp_data$vr_elast_PC1), ncol = length(sp_data$vr_elast_PC1) - 1) #minus 1 becasue distance of each pop to itself is 0
TLCA_all <- cophenetic(Tree) #Time to last comon ancestor for all species, need to pick out each pair of species in the actual data set 

#first pass to make the distance and neighbour value matricies
for(y in seq_along(shef_dat$SpeciesAuthor)){ 
  count = 1
  for(x in seq_along(shef_dat$SpeciesAuthor)[-y]){
    #geographic distance
    geo_dist[y, count] <- Ellipsoidal.Distance(olat = shef_dat$lat_dd[y], olon = shef_dat$lon_dd[y], 
      tlat = shef_dat$lat_dd[x], tlon = shef_dat$lon_dd[x])$dist 
    #gets phylogenetic distance 
    phy_dist[y, count] <- as.numeric(TLCA_all[shef_dat$SpeciesAuthor[y], shef_dat$SpeciesAuthor[x]])
    count = count + 1
  }
} 

# make a matrix, same dim as the distance matricies, put the Y value for each neigbour in each coloumn
Y_neigh <- matrix(nrow = dim(geo_dist)[1], ncol = dim(geo_dist)[2]) 
  for(y in seq_along(shef_dat$damping_ratio)){ 
  count = 1
  for(x in seq_along(shef_dat$damping_ratio)[-y]){
    Y_neigh[y, count] <- shef_dat$damping_ratio[x]
    count = count + 1
  }
}

### NOW DEFINE THE MODEL IN JAGS

sink("nofixed_2ac_decay_model.jags")
  cat('
    model{
      #PRIORS
      tau.global ~ dgamma(0.0001, 0.0001)#overall percison
      sd.global <- 1/sqrt(tau.global)#overall sd

      a ~ dnorm(0, 0.0001)
      for(i in 1:2){
	b[i] ~ dnorm(0, 0.0001)
      }

      #auto correlation parameters
      # exponent on the phlyogenetic prediction model which is a weighted average of other species, 
      # alpha controls the decay of the weight with distance  
      alpha ~ dunif(min_alpha, max_alpha) 
      
      # exponent on the geographic prediction model which is a weighted average of other population, 
      # weighted by distance to between study sites, rho controls the decay of the weight with distance  
      rho ~ dunif(min_rho, max_rho) 

      #LIKLIEHOOD 
      for(n in 1:num_data){
	#fit the linear model
	Y[n] ~ dnorm(mu[n], tau.global) 
	mu[n] <- a + b[1]*phy_pred[n] + b[2]*geo_pred[n] 

	#fit the phylo and geo models
	for(k in 1:(num_data - 1)){
	  phy_weight[n, k] <- exp(-alpha*phy_dist[n, k])
	  geo_weight[n, k] <- exp(-rho*geo_dist[n, k]) 
	}
	#make the sum of weights in denominator > 0
	phy_pred[n] <- sum(Y_neigh_phy[n, 1:phy_max_ind[n]]*phy_weight[n, 1:phy_max_ind[n]]) / max(0.00000000000001, sum(phy_weight[n, 1:phy_max_ind[n]]))
	geo_pred[n] <- sum(Y_neigh_geo[n, 1:geo_max_ind[n]]*geo_weight[n, 1:geo_max_ind[n]]) / max(0.00000000000001, sum(geo_weight[n, 1:geo_max_ind[n]]))

      }

      #MODEL PERFORMANCE	
      for(n in 1:num_data){
	#model performance predictors to make sure the model is not compleatly rubbish
	residual[n] <- Y[n] - mu[n]#takes the residual out of the eq below so can be tracked
	sq[n] <- pow(residual[n], 2)#residual squared 
	Y.new[n] ~ dnorm(mu[n], tau.global)#new data
	sq.new[n] <- pow(Y.new[n] - mu[n], 2)#residual squared for the new data
      }
      fit <- sum(sq[])#sum of squares for data
      fit.new <- sum(sq.new[])#sum of squares for simulated data
      test <- step(fit.new - fit)
      bpvalue <- mean(test)
    }',fill = TRUE)
sink()


### BRING TOGETHER ALL THE DATA, MODEL AND DISTANCE MATRICIES AND FIT THE MODEL
### USING THE R2JAGS() INTERFACE

jags.data <- list(Y = log(shef_dat$damping_ratio), geo_dist = geo_dist, phy_dist = phy_dist, 
  Y_neigh_geo = log(Y_neigh), Y_neigh_phy = log(Y_neigh), num_data = length(shef_dat$damping_ratio), 
  geo_max_ind = rep(dim(Y_neigh)[2], dim(Y_neigh)[1]), phy_max_ind = rep(dim(Y_neigh)[2], dim(Y_neigh)[1]), 
  min_alpha = 0, max_alpha = 1, min_rho = 0, max_rho = 1)

#intial values
inital <- function() list(a = runif(1, -20, 20), b = runif(2, -20, 20), tau.global = runif(1, 0.0001, 100), 
  alpha = runif(1, 0, 1), rho = runif(1, 0, 1)) 

# Parameters monitored
monitored<- c('a', "b", "alpha", "rho", "sd.global", "bpvalue")

#mcmc params
nIter <- 500
nThin <- 1
nBurnIn <- 1
nChains <- 3

#call to JAGS
set.seed(123)
damp_ratio_nofixed_allpops <- jags(jags.data, inital, monitored, "nofixed_2ac_decay_model.jags", 
  n.chains = nChains, n.thin = nThin, n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

# Some simple diagnostic plots
jags_obj_mc <- as.mcmc(damp_ratio_nofixed_allpops)
plot(jags_obj_mc, ask = TRUE)
autocorr.plot(jags_obj_mc, ask = TRUE, lag.max = 10)
