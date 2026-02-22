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
## ----------------Custom Nimble Function------------------------ ##
####################################################################

ind_1 <- nimbleFunction(
  run = function(x = double(0)) {
    out_1 <- 0
    if (x < 4) {
      out_1 <- 1
    }
    returnType(double(0))
    return(out_1)
  }
)

ind_2 <- nimbleFunction(
  run = function(x = double(0), y=double(0)) {
    out_2 <- 2
    if (x*y == 1) {
      out_2 <- 1
    }
    returnType(double(0))
    return(out_2)
  }
)

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

retag_Lf <- read.csv("Left_retag_f_matrix.csv")
retag_Ln <- read.csv("Left_retag_n_matrix.csv")

retag_Rf <- read.csv("Right_retag_f_matrix.csv")
retag_Rn <- read.csv("Right_retag_n_matrix.csv")

s_indi_L <- read.csv("Indi_left.csv")
s_indi_R <- read.csv("Indi_right.csv")

code_m1 <- nimbleCode({
  # Priors
  phi ~ dbeta(1,1) # Survival Probability
  p ~ dbeta(1,1) # Capture Probability
  
  lam_Lf ~ dbeta(1,1) # Tag 1 Retention Probability First Time
  lam_Ln ~ dbeta(1,1) # Tag 1 Retention Probability Subseq Time
  
  lam_Rf ~ dbeta(1,1) # Tag 2 Retention Probability First Time
  lam_Rn ~ dbeta(1,1) # Tag 2 Retention Probability Subseq Time
  
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
      
      pt[1,1,i,t] <- ((s_L[i,t-1]*lam_Lf)+((1-s_L[i,t-1])*r_Lf[i,t-1]*lam_Ln))*
        ((s_R[i,t-1]*lam_Rf)+((1-s_R[i,t-1])*r_Rf[i,t-1]*lam_Rn))
      
      pt[1,2,i,t] <- ((s_L[i,t-1]*lam_Lf)+((1-s_L[i,t-1])*r_Lf[i,t-1]*lam_Ln))*
        (1-((s_R[i,t-1]*lam_Rf)+((1-s_R[i,t-1])*r_Rf[i,t-1]*lam_Rn)))
      
      pt[1,3,i,t] <- (1-((s_L[i,t-1]*lam_Lf)+((1-s_L[i,t-1])*r_Lf[i,t-1]*lam_Ln)))*
        ((s_R[i,t-1]*lam_Rf)+((1-s_R[i,t-1])*r_Rf[i,t-1]*lam_Rn))
      
      pt[1,4,i,t] <- (1-((s_L[i,t-1]*lam_Lf)+((1-s_L[i,t-1])*r_Lf[i,t-1]*lam_Ln)))*
        (1-((s_R[i,t-1]*lam_Rf)+((1-s_R[i,t-1])*r_Rf[i,t-1]*lam_Rn)))
      
      
      pt[2,1,i,t] <- ((s_L[i,t-1]*lam_Lf)+((1-s_L[i,t-1])*r_Lf[i,t-1]*lam_Ln))*
        ((s_R[i,t-1]*lam_Rf)+((1-s_R[i,t-1])*r_Rf[i,t-1]*r_Rn[i,t-1]*lam_Rn))
      
      pt[2,2,i,t] <- ((s_L[i,t-1]*lam_Lf)+((1-s_L[i,t-1])*r_Lf[i,t-1]*lam_Ln))*
        (1-((s_R[i,t-1]*lam_Rf)+((1-s_R[i,t-1])*r_Rf[i,t-1]*r_Rn[i,t-1]*lam_Rn)))
      
      pt[2,3,i,t] <- (1-((s_L[i,t-1]*lam_Lf)+((1-s_L[i,t-1])*r_Lf[i,t-1]*lam_Ln)))*
        ((s_R[i,t-1]*lam_Rf)+((1-s_R[i,t-1])*r_Rf[i,t-1]*r_Rn[i,t-1]*lam_Rn))
      
      pt[2,4,i,t] <- (1-((s_L[i,t-1]*lam_Lf)+((1-s_L[i,t-1])*r_Lf[i,t-1]*lam_Ln)))*
        (1-((s_R[i,t-1]*lam_Rf)+((1-s_R[i,t-1])*r_Rf[i,t-1]*r_Rn[i,t-1]*lam_Rn)))
      
      
      pt[3,1,i,t] <- ((s_L[i,t-1]*lam_Lf)+((1-s_L[i,t-1])*r_Lf[i,t-1]*r_Ln[i,t-1]*lam_Ln))*
        ((s_R[i,t-1]*lam_Rf)+((1-s_R[i,t-1])*r_Rf[i,t-1]*lam_Rn))
      
      pt[3,2,i,t] <- ((s_L[i,t-1]*lam_Lf)+((1-s_L[i,t-1])*r_Lf[i,t-1]*r_Ln[i,t-1]*lam_Ln))*
        (1-((s_R[i,t-1]*lam_Rf)+((1-s_R[i,t-1])*r_Rf[i,t-1]*lam_Rn)))
      
      pt[3,3,i,t] <- (1-((s_L[i,t-1]*lam_Lf)+((1-s_L[i,t-1])*r_Lf[i,t-1]*r_Ln[i,t-1]*lam_Ln)))*
        ((s_R[i,t-1]*lam_Rf)+((1-s_R[i,t-1])*r_Rf[i,t-1]*r_Rn[i,t-1]*lam_Rn))
      
      pt[3,4,i,t] <- (1-((s_L[i,t-1]*lam_Lf)+((1-s_L[i,t-1])*r_Lf[i,t-1]*r_Ln[i,t-1]*lam_Ln)))*
        (1-((s_R[i,t-1]*lam_Rf)+((1-s_R[i,t-1])*r_Rf[i,t-1]*lam_Rn)))
      
      
      pt[4,1,i,t] <- 0
      
      pt[4,2,i,t] <- 0
      
      pt[4,3,i,t] <- 0
      
      pt[4,4,i,t] <- 1
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
      k[i,t] <- ind_1(g[i,t])
      m[i,t] <- ind_2(a[i,t],k[i,t])
      h[i,t] ~ dcat(po[m[i, t], 1:2])
      
    }
  }
})

my.data <- list(h = cap + 1, a = alive, g = tag,
                s_L = s_indi_L, s_R = s_indi_R,
                r_Lf = retag_Lf, r_Ln = retag_Ln,
                r_Rf = retag_Rf, r_Rn = retag_Rn)

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
                     inprob_g=c(0.34, 0.33, 0.33, 0))

initial.values <- list(phi = runif(1,0,1), 
                       p = runif(1,0,1), 
                       lam_Lf = runif(1,0,1),
                       lam_Rf = runif(1,0,1),
                       lam_Ln = runif(1,0,1),
                       lam_Rn = runif(1,0,1),
                       a = a_inits,
                       g = g_inits)

parameters.to.save <- c("phi", "p", "lam_Lf", "lam_Rf",
                        "lam_Ln", "lam_Rn")

#parameters.to.save <- c("phi", "p", "lam_L", "lam_R", "g_1", "g_2", "h")

n.iter <- 50000
n.burnin <- 2000
n.chains <- 3

# m <- nimbleModel(code = code_m1,
#                  constants = my.constants,
#                  data = my.data,
#                  inits = initial.values)
# m$calculate()
# 
# m$g
# m$logProb_g

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

pdf(file = "ms_wopt_mp_lam2.pdf")
MCMCplot(samples, HPD = T)
dev.off()

s <- MCMCsummary(samples, round = 5)
MCMCtrace(samples,pdf = T,open_pdf = F,filename = "ms_wopt_lam2", ind = TRUE,
          Rhat = FALSE, n.eff = FALSE)
write.csv(s, file = "ms_wopt_lam2_sum.csv")