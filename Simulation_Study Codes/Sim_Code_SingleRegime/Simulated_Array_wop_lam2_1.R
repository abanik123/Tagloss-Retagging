library(nimble)
library(MCMCvis)
library(tidyverse)
library(dplyr)

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
    
    if (x < 1 | x > length(prob)) {
      print("I found a problem with x = ", x)
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

# x = 1 is the alive state and y < 4 means at least one tag exist

ind_11 <- nimbleFunction(
  run = function(x = double(0), y=double(0)) {
    out_2 <- 2
    if (x == 1 & y < 4) {
      out_2 <- 1
    }
    returnType(double(0))
    return(out_2)
  }
)



# Define the model code for SIMULATION

sim_function <- function(params){
  
  iter = params$iter
  
  ############################
  ## Functions
  ############################
  
  ## Alive & tag indicator for capture
  ind_1 <- function(a, g) {
    if (a == 1 && g < 4) return(1)
    return(2)
  }
  
  
  ############################
  ## Dimensions
  ############################
  N      <- 2000
  nyears <- 6
  
  ############################
  ## Parameters
  ############################
  
  ## Survival
  phi <- params$phi
  
  ## Detection
  p <- params$p
  
  ## Retention probabilities
  lam_Lf <- params$lam_Lf   # fresh left
  lam_Rf <- params$lam_Rf   # fresh right
  
  lam_Ln <- params$lam_Ln   # non-fresh left
  lam_Rn <- params$lam_Rn   # non-fresh right
  
  ############################
  ## Helper functions
  ############################
  
  # Latent state process probabilities
  
  px <- matrix(NA, nrow = 2, ncol = 2)
  
  px[1, 1] <- phi
  px[1, 2] <- 1 - phi
  
  px[2, 1] <- 0
  px[2, 2] <- 1
  
  # Observation process probabilities
  
  po <- matrix(NA, nrow = 2, ncol = 2)
  
  po[1, 1] <- 1 - p
  po[1, 2] <- p
  
  po[2, 1] <- 1
  po[2, 2] <- 0
  
  # Tag state process probabilities
  pt <- array(NA, dim = c(4, 4, N, nyears))
  
  ############################
  ## Storage
  ############################
  a <- matrix(NA, nrow = N, ncol = nyears)
  g <- matrix(NA, nrow = N, ncol = nyears)
  h <- matrix(NA, nrow = N, ncol = nyears)
  m <- matrix(NA, nrow = N, ncol = nyears)
  
  # Always 1
  r_Lf <- r_Rf <- r_Ln <- r_Rn <- matrix(0, N, nyears)
  
  # Freshness
  s_L  <- matrix(0, N, nyears)
  s_R  <- matrix(0, N, nyears)
  
  ############################
  ## First capture times
  ############################
  first <- rep(NA, N)
  for (i in 1:N) {
    first[i] <- sample(1:nyears,1)
  }
  
  ############################
  ## Initialization
  ############################
  for (i in 1:N) {
    
    # Initial tag state
    g[i, first[i]] <- 1  # at least one tag
    a[i, first[i]] <- 1
    h[i, first[i]] <- 2
    
    # Initial s_L, s_R
    if (g[i, first[i]] == 1) {       # both tags
      s_L[i, first[i]] <- 1
      s_R[i, first[i]] <- 1
    } 
    
    if (first[i] > 1) {
      g[i, 1:(first[i]-1)] <- NA
      a[i, 1:(first[i]-1)] <- 2
      h[i, 1:(first[i]-1)] <- NA
    }
  }
  
  
  # Create permanent zero flags for each individual
  s_L_zero_permanent <- rep(FALSE, N)
  s_R_zero_permanent <- rep(FALSE, N)
  
  r_Lf_permanent <- rep(FALSE, N)
  r_Rf_permanent <- rep(FALSE, N)
  
  ############################
  ## Main simulation
  ############################
  for (i in 1:N) {
    
    if (first[i] < nyears) {
      
      for (t in (first[i]+1):nyears) {
        
        if (h[i,t-1] == 2 && s_L[i,t-1] == 1) {
          
          # Only sample if not permanently zero
          if (!s_L_zero_permanent[i]) {
            s_L[i, t] <- sample(c(1,0), 1, prob = c(lam_Lf, 1 - lam_Lf))
            if (s_L[i, t] == 0) s_L_zero_permanent[i] <- TRUE
          } else {
            s_L[i, t] <- 0
          }
          
          if (!s_R_zero_permanent[i]) {
            s_R[i, t] <- sample(c(1,0), 1, prob = c(lam_Rf, 1 - lam_Rf))
            if (s_R[i, t] == 0) s_R_zero_permanent[i] <- TRUE
          } else {
            s_R[i, t] <- 0
          }
          
        } else {
          # If not captured, keep 1 only if it hasn't gone permanently zero
          s_L[i, t] <- ifelse(s_L_zero_permanent[i], 0, 1)
          s_R[i, t] <- ifelse(s_R_zero_permanent[i], 0, 1)
        }
        
        # Inside for (i in 1:N) and for (t in ...)
        # ---------------------------------------
        # Reset non-fresh indicators (event-based)
        r_Ln[i, t-1] <- 0
        r_Rn[i, t-1] <- 0
        
        # LEFT tag loss observed (state 3 = left missing)
        if (h[i, t-1] == 2 && g[i, t-1] == 3) {
          
          # Non-fresh indicator (event only)
          r_Ln[i, t-1] <- 1
          
          # Fresh indicator (absorbing)
          if (!r_Lf_permanent[i]) {
            r_Lf[i, t-1] <- 1
            r_Lf_permanent[i] <- TRUE
          }
        }
        
        # RIGHT tag loss observed (state 2 = right missing)
        if (h[i, t-1] == 2 && g[i, t-1] == 2) {
          
          # Non-fresh indicator (event only)
          r_Rn[i, t-1] <- 1
          
          # Fresh indicator (absorbing)
          if (!r_Rf_permanent[i]) {
            r_Rf[i, t-1] <- 1
            r_Rf_permanent[i] <- TRUE
          }
        }
        
        # Enforce absorbing structure explicitly
        if (r_Lf_permanent[i]) r_Lf[i, t-1] <- 1
        if (r_Rf_permanent[i]) r_Rf[i, t-1] <- 1
        
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
        
        a[i,t] = rcat(n=1, px[a[i, t-1], 1:2])
        g[i,t] = rcat(n=1, pt[g[i, t - 1], 1:4,i,t])
        m[i,t] = ind_1(a[i,t],g[i,t])
        h[i,t] = rcat(n=1, po[m[i, t], 1:2])
        
      }
    }
  }
  
  cap <- h
  
  cap[cap == 1] <- 0
  cap[cap == 2] <- 1
  
  cap[is.na(cap)] <- 0
  
  known.state.cjs <- function(ch){
    state <- ch
    for (i in 1:dim(ch)[1]){
      n1 <- min(which(ch[i,]==1))
      n2 <- max(which(ch[i,]==1))
      state[i,n1:n2] <- 1
      #state[i,1:n1-1] <- NA
    }
    state[state==0] <- NA
    return(state)
  }
  
  alive <- known.state.cjs(cap)
  
  g_1 <- g
  
  g_1[cap != 1] <- NA
  
  r_Ln[,] <- 0
  r_Ln[cap == 1 & g_1 == 3] <- 1
  
  r_Rn[,] <- 0
  r_Rn[cap == 1 & g_1 == 2] <- 1
  
  fill_g1_zero <- function(x) {
    idx <- which(!is.na(x))
    if (length(idx) == 0) return(x)
    
    # Before first known value
    if (idx[1] > 1) {
      x[1:(idx[1] - 1)] <- 0
    }
    
    # Between known values
    if (length(idx) > 1) {
      for (k in 1:(length(idx) - 1)) {
        if (idx[k] + 1 <= idx[k + 1] - 1) {
          x[(idx[k] + 1):(idx[k + 1] - 1)] <- 0
        }
      }
    }
    
    x
  }
  
  g_1_f <- g_1
  for (i in 1:nrow(g_1)) {
    g_1_f[i, ] <- fill_g1_zero(g_1[i, ])
  }
  
  process_row <- function(a, b, c) {
    
    original <- a
    nT <- length(a)
    
    for (i in 1:(nT - 1)) {
      
      if (is.na(a[i])) next
      
      next_idx <- which(!is.na(a[(i + 1):nT]))
      if (length(next_idx) == 0) next
      
      j <- i + next_idx[1]
      
      # ---------- RULE 1 ----------
      if (a[i] == 1 && a[j] == 1) {
        a[(i + 1):(j - 1)] <- ifelse(
          is.na(a[(i + 1):(j - 1)]), 1, a[(i + 1):(j - 1)]
        )
      }
      
      # ---------- RULE 2 ----------
      if (a[i] == 2) {
        if (!is.na(c[i]) && c[i] == 1 && a[j] == 1) {
          a[(i + 1):(j - 1)] <- ifelse(
            is.na(a[(i + 1):(j - 1)]), 1, a[(i + 1):(j - 1)]
          )
        }
        if ((is.na(c[i]) || c[i] == 0) && a[j] == 2) {
          a[(i + 1):(j - 1)] <- ifelse(
            is.na(a[(i + 1):(j - 1)]), 2, a[(i + 1):(j - 1)]
          )
        }
      }
      
      # ---------- RULE 3 ----------
      if (a[i] == 3) {
        if (!is.na(b[i]) && b[i] == 1 && a[j] == 1) {
          a[(i + 1):(j - 1)] <- ifelse(
            is.na(a[(i + 1):(j - 1)]), 1, a[(i + 1):(j - 1)]
          )
        }
        if ((is.na(b[i]) || b[i] == 0) && a[j] == 3) {
          a[(i + 1):(j - 1)] <- ifelse(
            is.na(a[(i + 1):(j - 1)]), 3, a[(i + 1):(j - 1)]
          )
        }
      }
      
      # ---------- RULE 4 ----------
      if (a[i] == 4) {
        
        bb <- ifelse(is.na(b[i]), 0, b[i])
        cc <- ifelse(is.na(c[i]), 0, c[i])
        
        if (bb == 1 && cc == 1 && a[j] == 1) {
          a[(i + 1):(j - 1)] <- ifelse(
            is.na(a[(i + 1):(j - 1)]), 1, a[(i + 1):(j - 1)]
          )
        }
        if (bb == 1 && cc == 0 && a[j] == 2) {
          a[(i + 1):(j - 1)] <- ifelse(
            is.na(a[(i + 1):(j - 1)]), 2, a[(i + 1):(j - 1)]
          )
        }
        if (bb == 0 && cc == 1 && a[j] == 3) {
          a[(i + 1):(j - 1)] <- ifelse(
            is.na(a[(i + 1):(j - 1)]), 3, a[(i + 1):(j - 1)]
          )
        }
        if (bb == 0 && cc == 0 && a[j] == 4) {
          a[(i + 1):(j - 1)] <- ifelse(
            is.na(a[(i + 1):(j - 1)]), 4, a[(i + 1):(j - 1)]
          )
        }
      }
    }
    
    # Preserve original observed values
    a[!is.na(original)] <- original[!is.na(original)]
    
    return(a)
  }
  
  # Convert 0 → NA
  A <- g_1_f; A[A == 0] <- NA
  B <- r_Ln;   B[B == 0] <- NA
  C <- r_Rn;   C[C == 0] <- NA
  
  out <- A
  
  for (r in 1:nrow(A)) {
    out[r, ] <- process_row(
      as.numeric(A[r, ]),
      as.numeric(B[r, ]),
      as.numeric(C[r, ])
    )
  }
  
  convert_first1_to_ones <- function(A, na_as_zero = TRUE) {
    mat <- as.matrix(A)
    mat_num <- apply(mat, 2, function(x) as.numeric(as.character(x)))
    if (na_as_zero) {
      mat_num[is.na(mat_num)] <- 0   # optional: treat NA as 0
    }
    B_mat <- t(apply(mat_num, 1, cummax))
    B <- as.data.frame(B_mat, stringsAsFactors = FALSE)
    colnames(B) <- colnames(A)
    rownames(B) <- rownames(A)
    B[] <- lapply(B, function(x) as.integer(x))
    return(B)
  }
  
  r_Lf <- as.matrix(convert_first1_to_ones(r_Ln))
  r_Rf <- as.matrix(convert_first1_to_ones(r_Rn))
  
  absorbing_zero_from_first_one <- function(A) {
    
    N      <- nrow(A)
    Tt     <- ncol(A)
    
    B <- matrix(1, N, Tt)
    
    for (i in 1:N) {
      
      hit <- FALSE
      
      for (t in 1:Tt) {
        
        if (!hit && A[i, t] == 1) {
          hit <- TRUE
        }
        
        if (hit) {
          B[i, t] <- 0
        }
      }
    }
    
    return(B)
  }
  
  s_indi_L <- absorbing_zero_from_first_one(r_Lf)
  s_indi_R <- absorbing_zero_from_first_one(r_Rf)
  
  ####################################################################
  ## ------------------Finalizing Datasets------------------------- ##
  ####################################################################
  
  cap <- cap
  alive <- alive
  tag <- out
  
  retag_Lf <- r_Lf
  retag_Rf <- r_Rf
  
  retag_Ln <- r_Ln
  retag_Rn <- r_Rn
  
  retag_Lf <- retag_Lf[, -ncol(retag_Lf)]
  retag_Rf <- retag_Rf[, -ncol(retag_Rf)]
  
  retag_Ln <- retag_Ln[, -ncol(retag_Ln)]
  retag_Rn <- retag_Rn[, -ncol(retag_Rn)]
  
  s_indi_L <- s_indi_L[, -ncol(s_indi_L)]
  s_indi_R <- s_indi_R[, -ncol(s_indi_R)]
  
  keep <- first < nyears
  
  cap     <- cap[keep, , drop = FALSE]
  alive    <- alive[keep, , drop = FALSE]
  tag    <- tag[keep, , drop = FALSE]
  
  retag_Lf    <- retag_Lf[keep, , drop = FALSE]
  retag_Rf    <- retag_Rf[keep, , drop = FALSE]
  
  retag_Ln    <- retag_Ln[keep, , drop = FALSE]
  retag_Rn    <- retag_Rn[keep, , drop = FALSE]
  
  s_indi_L    <- s_indi_L[keep, , drop = FALSE]
  s_indi_R    <- s_indi_R[keep, , drop = FALSE]
  
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
      for (t in (f_1[i]+1):nyears) {
        
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
      
      for (t in (f_1[i]+1):nyears) { 
        
        a[i,t] ~ dcat(px[a[i, t - 1], 1:2])
        g[i,t] ~ dcatt(pt[g[i, t - 1], 1:4,i,t])
        m[i,t] <- ind_11(a[i,t],g[i,t])
        h[i,t] ~ dcat(po[m[i,t], 1:2])
        
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
  
  initial.values <- list(phi = 0.5, 
                         p = 0.3, 
                         lam_Lf = 0.6,
                         lam_Rf = 0.6,
                         lam_Ln = 0.6,
                         lam_Rn = 0.6,
                         a = a_inits,
                         g = g_inits)
  
  parameters.to.save <- c("phi", "p", "lam_Lf", "lam_Rf",
                          "lam_Ln", "lam_Rn")
  
  #parameters.to.save <- c("phi", "p", "lam_L", "lam_R", "g_1", "g_2", "h")
  
  n.iter <- 30000
  n.burnin <- 3000
  n.chains <- 2
  
  mcmc.multistate <- nimbleMCMC(code = code_m1, 
                                constants = my.constants,
                                data = my.data,              
                                inits = initial.values,
                                monitors = parameters.to.save,
                                niter = n.iter,
                                nburnin = n.burnin, 
                                nchains = n.chains,
                                summary = TRUE, WAIC = TRUE)
  
  samples_sim <- mcmc.multistate$samples
  
  sim_lam_Lf <- MCMCsummary(samples_sim, params = "lam_Lf", round = 5)$mean
  sim_lam_Rf <- MCMCsummary(samples_sim, params = "lam_Rf", round = 5)$mean
  
  sim_lam_Ln <- MCMCsummary(samples_sim, params = "lam_Ln", round = 5)$mean
  sim_lam_Rn <- MCMCsummary(samples_sim, params = "lam_Rn", round = 5)$mean
  
  sim_phi <- MCMCsummary(samples_sim, params = "phi", round = 5)$mean
  sim_p <- MCMCsummary(samples_sim, params = "p", round = 5)$mean
  
  d_b = data.frame(sim_lam_Lf, sim_lam_Rf,
                   sim_lam_Ln, sim_lam_Rn,
                   sim_phi, sim_p)
  
  write.csv(d_b, paste0("lam_Lf", sim_lam_Lf, "_lam_Rf",sim_lam_Rf,
                        "_phi",sim_phi,"_p",sim_p,
                        "_Sim_",iter,".csv"), row.names = TRUE)
  
  return(d_b)
}

args <- commandArgs(trailingOnly = TRUE)
iter <- as.numeric(args[1])
params= list(p = 0.4, phi = 0.7,
             lam_Lf = 0.97, lam_Rf = 0.97,
             lam_Ln = 0.7, lam_Rn = 0.5,
             iter=iter)
sim_function(params)