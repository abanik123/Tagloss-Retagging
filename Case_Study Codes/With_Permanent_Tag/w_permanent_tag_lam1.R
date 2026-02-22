library(nimble)
library(MCMCvis)

####################################################################
## ------------------Custom Distribution------------------------- ##
####################################################################

dcatt <- nimbleFunction(
  run = function(x = double(0),
                 prob = double(1),
                 log = logical(0)) {
    returnType(double(0))
    
    if (x < 1 | x > length(prob)) {
      if (log) 
        return(-Inf)
      else
        return(0)
    }
    
    if (prob[x] == 0) {
      if (log)
        return(0)
      else
        return(prob[x])
    }
    if (prob[x] != 0) {
      if (log) {
        return(log(prob[x]))
      } else {
        return(prob[x])
      }
    }
  }
)

rcatt <- nimbleFunction(
  
  run = function(n = integer(), prob = double(1)) {
    
    returnType(integer())
    
    return(rcat(1, prob))
    
})


####################################################################
## ---------------------------CODE------------------------------- ##
####################################################################

get_first_non_zero_index <- function(row) {
  non_zero_indices <- which(row != 0)
  if (length(non_zero_indices) > 0) {
    return(non_zero_indices[1])
  } else {
    return(NA)
  }
}

cap <- read.csv("Capture_matrix.csv")
alive <- read.csv("Alive_matrix.csv") # Partially Latent
tag <- read.csv("Tag_matrix_final.csv") # Partially Latent

retag_L <- read.csv("Left_retag_matrix.csv") # When retagging happen for tag 1 (left tag)
retag_R <- read.csv("Right_retag_matrix.csv") # When retagging happen for tag 2 (right tag)

s_indi_L <- read.csv("Indi_left.csv")
s_indi_R <- read.csv("Indi_right.csv")

code_m1 <- nimbleCode({
  # Priors
  phi ~ dbeta(1,1) # Survival Probability
  p ~ dbeta(1,1) # Capture Probability
  
  lam_L ~ dbeta(1,1) # Tag 1 Retention Probability
  lam_R ~ dbeta(1,1) # Tag 2 Retention Probability
  
  # Latent state process probabilities
  px[1, 1] <- phi
  px[1, 2] <- 1 - phi
  
  px[2, 1] <- 0
  px[2, 2] <- 1
  
  # Observation process probabilities
  po[1, 1] <- 1 - p
  po[1, 2] <- p
  
  po[2, 1] <- 1
  po[2, 2] <- 0
  
  # Tag state process probabilities
  # probabilities of state z(t+1) given z(t)
  for (i in 1:N) {
    for (t in seq(f_1[i] + 1, nyears, 
                  length.out = max(0, nyears - f_1[i]))) {
      
      pt[1,1,i,t] <- lam_L*lam_R
      pt[1,2,i,t] <- lam_L*(1-lam_R)
      pt[1,3,i,t] <- (1-lam_L)*lam_R
      pt[1,4,i,t] <- (1-lam_L)*(1-lam_R)
      
      pt[2,1,i,t] <- lam_L*((s_R[i,t-1]*lam_R)+
                        ((1-s_R[i,t-1])*r_R[i,t-1]*lam_R))
      
      pt[2,2,i,t] <- lam_L*((s_R[i,t-1]*(1-lam_R))+
                        ((1-s_R[i,t-1])*((1-r_R[i,t-1])+r_R[i,t-1]*(1-lam_R))))
      
      pt[2,3,i,t] <- (1-lam_L)*((s_R[i,t-1]*lam_R)+
                        ((1-s_R[i,t-1])*r_R[i,t-1]*lam_R))
      
      pt[2,4,i,t] <- (1-lam_L)*((s_R[i,t-1]*(1-lam_R))+
                        ((1-s_R[i,t-1])*((1-r_R[i,t-1])+r_R[i,t-1]*(1-lam_R))))
      
      pt[3,1,i,t] <- ((s_L[i,t-1]*lam_L)+
                      ((1-s_L[i,t-1])*r_L[i,t-1]*lam_L))*lam_R
      
      pt[3,2,i,t] <- ((s_L[i,t-1]*lam_L)+
                      ((1-s_L[i,t-1])*r_L[i,t-1]*lam_L))*(1-lam_R)
      
      pt[3,3,i,t] <- ((s_L[i,t-1]*(1-lam_L))+
                      ((1-s_L[i,t-1])*((1-r_L[i,t-1])+r_L[i,t-1]*(1-lam_L))))*lam_R
      
      pt[3,4,i,t] <- ((s_L[i,t-1]*(1-lam_L))+
                      ((1-s_L[i,t-1])*((1-r_L[i,t-1])+r_L[i,t-1]*(1-lam_L))))*(1-lam_R)
      
      pt[4,1,i,t] <- ((s_L[i,t-1]*lam_L)+((1-s_L[i,t-1])*r_L[i,t-1]*lam_L))*
                      ((s_R[i,t-1]*lam_R)+((1-s_R[i,t-1])*r_R[i,t-1]*lam_R))
      
      pt[4,2,i,t] <- ((s_L[i,t-1]*lam_L)+((1-s_L[i,t-1])*r_L[i,t-1]*lam_L))*
                      ((s_R[i,t-1]*(1-lam_R))+((1-s_R[i,t-1])*
                                    ((1-r_R[i,t-1])+r_R[i,t-1]*(1-lam_R))))
      
      pt[4,3,i,t] <- ((s_L[i,t-1]*(1-lam_L))+((1-s_L[i,t-1])*
                              ((1-r_L[i,t-1])+r_L[i,t-1]*(1-lam_L))))*
                      ((s_R[i,t-1]*lam_R)+((1-s_R[i,t-1])*r_R[i,t-1]*lam_R))
      
      pt[4,4,i,t] <- ((s_L[i,t-1]*(1-lam_L))+((1-s_L[i,t-1])*
                              ((1-r_L[i,t-1])+r_L[i,t-1]*(1-lam_L))))*
                      ((s_R[i,t-1]*(1-lam_R))+((1-s_R[i,t-1])*
                              ((1-r_R[i,t-1])+r_R[i,t-1]*(1-lam_R))))
    }
  }
  
  #Likelihood 2 (from 2nd year of capture)
  for (i in 1:N){
    a[i, f_1[i]] ~ dcat(inprob_a[1:2])
    g[i, f_1[i]] ~ dcat(inprob_g[1:4])
    
    for (t in seq(f_1[i] + 1, nyears, 
                  length.out = max(0, nyears - f_1[i]))) { 
      
      a[i,t] ~ dcat(px[a[i, t - 1], 1:2])
      g[i,t] ~ dcatt(pt[g[i, t - 1], 1:4,i,t])
      h[i,t] ~ dcat(po[a[i, t], 1:2])
      
    }
  }
})

my.data <- list(h = cap + 1, a = alive, g = tag,
                s_L = s_indi_L, s_R = s_indi_R,
                r_L = retag_L, r_R = retag_R)

N = nrow(cap)
nyears = ncol(cap)

a_inits <- alive
a_inits[is.na(alive)] <- 1
a_inits[!is.na(alive)] <- NA

g_inits <- tag
g_inits[is.na(tag)] <- 1
#g_inits[is.na(tag)] <- sample(1:4, sum(is.na(tag)), replace = TRUE)
g_inits[!is.na(tag)] <- NA

f_1 <- apply(cap, 1, get_first_non_zero_index)

for (i in 1:N) {
  if (f_1[i] > 1) a_inits[i, 1:(f_1[i]-1)] <- NA
  if (f_1[i] > 1) g_inits[i, 1:(f_1[i]-1)] <- NA
  
}

my.constants <- list(N = N, nyears = nyears, 
                     f_1 = f_1, inprob_a=c(1, 0), 
                     inprob_g=c(0.25, 0.25, 0.25, 0.25))

initial.values <- list(phi = runif(1,0,1), 
                       p = runif(1,0,1), 
                       lam_L = runif(1,0,1),
                       lam_R = runif(1,0,1),
                       a = a_inits,
                       g = g_inits)

parameters.to.save <- c("phi", "p", "lam_L", "lam_R")

#parameters.to.save <- c("phi", "p", "lam_L", "lam_R", "g_1", "g_2", "h")

n.iter <- 50000
n.burnin <- 2000
n.chains <- 3

# m <- nimbleModel(code = code_m1,
#                  constants = my.constants,
#                  data = my.data,
#                  inits = initial.values)
# #m$calculate()

#m$g
#m$logProb_g

mcmc.multistate <- nimbleMCMC(code = code_m1, 
                              constants = my.constants,
                              data = my.data,              
                              inits = initial.values,
                              monitors = parameters.to.save,
                              niter = n.iter,
                              nburnin = n.burnin, 
                              nchains = n.chains,
                              summary = TRUE, WAIC = TRUE)

#Model Summary
samples<- mcmc.multistate$samples

pdf(file = "ms_wpt_mp_lam1.pdf")
MCMCplot(samples, HPD = T)
dev.off()

s <- MCMCsummary(samples, round = 5)
MCMCtrace(samples,pdf = T,open_pdf = F,filename = "ms_wpt_lam1", ind = TRUE,
          Rhat = FALSE, n.eff = FALSE)
write.csv(s, file = "ms_wpt_lam1_sum.csv")