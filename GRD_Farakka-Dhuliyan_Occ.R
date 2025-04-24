library(jagsUI)
library(rjags)
library(mcmcOutput)

df <- read.csv("Data/Farakka-Dhuliyan_occdata_3J.csv")
rivwidth_df <- read.csv("Data/Farakka-Dhuliyan_rivwidth_2020-2024_3J.csv")
boat_df <- read.csv("Data/Farakka-Dhuliyan_occ_sitecovs.csv")
  
y <- df[, 2:ncol(df)]
y_analysis <- data.frame(matrix(nrow = nrow(y), ncol = ncol(y))) 
log_rivwidth <- log10(rivwidth_df[, 2:ncol(rivwidth_df)]) 

## converting observed abundance to detection/non-detection
for(i in 1:nrow(y)){
  for(j in 1:ncol(y)){
    y_analysis[i,j] <- ifelse(y[i,j] > 0, 1, 0)
  }
}
M <- nrow(y_analysis) ## no. of sites
J <- 3 ## no. of secondary sampling surveys
T <- ncol(y_analysis) / J ## no. of primary time periods

## creating T matrices of dimensions M * J 
y_array <- array(unlist(y_analysis), dim = c(M, J, T))
rivwidth_array <- array(unlist(log_rivwidth), dim = c(M, J, T))

## colnames(y) <- NULL
data <- list(y = y_array, 
             M = M,
             J = J,
             T = T,
             n_smallboats = boat_df$n_smallboats,
             n_mediumboats = boat_df$n_mediumboats,
             rivwidth = rivwidth_array)

## model
cat("

model {

# Occupancy 
for(i in 1:M) { 
  z[i, 1] ~ dbern(psi1[i]) ## first year occupancy state
  logit(psi1[i]) <- beta0 + 
                    beta1 * n_smallboats[i] + 
                    beta2 * n_mediumboats[i]  
  # probability of occupancy in each year after the first year as a function of 
  ## previous occupancy status and either gamma (probability of colonization) or epsilon (probability of extinction)
    for(t in 2:T) {
    z[i, t] ~ dbern(z[i, t-1] * (1 - eps) + (1 - z[i, t-1]) * gam) 
    }
}

# Detection 
for(i in 1:M) { 
  for(j in 1:J) {
    for(t in 1:T) {
      y[i, j, t] ~ dbern(p[i, j, t] * z[i, t]) 
      logit(p[i, j, t]) <- alpha0 + 
                         alpha1 * rivwidth[i, j, t]
    }
  }
}

# Priors
beta0 ~ dnorm(0, 0.01)
beta1 ~ dnorm(0, 0.01)
beta2 ~ dnorm(0, 0.01)
alpha0 ~ dnorm(0, 0.01)
alpha1 ~ dnorm(0, 0.01)
psi1 ~ dunif(0, 1)
eps ~ dunif(0, 1)
gam ~ dunif(0, 1)
p ~ dunif(0, 1)
}
", fill = TRUE, file = "grd_farakka-dhuliyan_dynocc.txt")

# Define parameters to save and MCMC settings
params <- c("beta0", "beta1","beta2", 
            "alpha0", "alpha1", 
            "psi1", "eps", "gam", "p")
nc <- 3 ; ni <- 11000 ; nb <- 1000 ; nt <- 1

output <- jags(data, inits = NULL, params, 
               "grd_farakka-dhuliyan_dynocc.txt", 
               n.thin = nt, n.chains = nc, n.burnin = nb, n.iter = ni)

mco <- mcmcOutput(output)
diagPlot(mco)
plot(mco)