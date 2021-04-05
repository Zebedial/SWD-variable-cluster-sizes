##########################################################################
# Key Functions Needed to Run the Main Simulation When the True Correlation 
# Structure is Exponential Decay (ED)
##########################################################################


##########################################################################
# Using cluster-period mean GEE to compute the variance of treatment effect 
# estimator, when the working correlation structure correctly specifies the
# true (ED) correlation structure
#
# Input:
# n: a (number of clusters) x (number of periods) size matrix of cluster 
#    period sizes
# scheme: a vector of randomization scheme whose elements are the number of 
#         clusters randomized in each step, e.g. c(5,5,5,5)
# delta: treatment effect in scale of log(OR)
# beta: a vector of period effect in scale of log(OR) with length (number of
#       periods)
# alpha0: WP-ICC alpha_0
# rho: decay rate from 0 to 1
#
# Output: a computed variance of treatment effect estimator
###########################################################################

sumVAR.exp <- function(n, scheme, delta, beta, alpha0, rho) {
  t <- ncol(n)
  I <- nrow(n)
  trtSeq <- matrix(0, t - 1, t)
  trtSeq[upper.tri(trtSeq)] <- 1
  trtSeq <- trtSeq[rep(1:nrow(trtSeq), scheme), ]
  if (sum(scheme) != I) {
    stop("invalid randomization scheme")
  }
  Omega <- matrix(0, t + 1, t + 1)
  
  # loop through each cluster
  for (i in 1:I) {
    # design matrix
    period <- rep(1:t)
    X <- trtSeq[i, ]
    for (j in 1:t) {
      X <- cbind(X, as.numeric(period == j))
    }
    
    # unique marginal means
    gmu <- c(X %*% c(delta, beta))
    mu <- plogis(gmu)
    
    
    B = matrix(0, t, t)
    v = mu * (1 - mu)
    ni = n[i, ]
    for (j in 1:(t - 1)) {
      B[j, j] <- (v[j]/ni[j]) * (1 + (ni[j] - 1) * alpha0)
      for (k in (j + 1):t) {
        B[j, k] <- sqrt(v[j] * v[k]) * alpha0 * rho^(abs(k - j))
      }
    }
    B[t, t] <- (v[t]/ni[t]) * (1 + (ni[t] - 1) * alpha0)
    B[lower.tri(B)] = t(B)[lower.tri(B)]
    invV <- solve(B)
    
    D <- X * (mu * (1 - mu))
    Omega <- Omega + t(D) %*% invV %*% D
  }
  vardelta <- solve(Omega)[1, 1]
  return(vardelta)
}

###########################################################################
# Using cluster-period mean GEE to compute the variance of treatment effect 
# estimator, when the working correlation structure misspecified the true  
# (ED) correlation structure by an independence (IND) correlation structure
#
# Input:
# n: a (number of clusters) x (number of periods) size matrix of cluster 
#    period sizes
# scheme: a vector of randomization scheme whose elements are the number of  
#         clusters randomized in each step, e.g. c(5,5,5,5)
# delta: treatment effect in scale of log(OR)
# beta: a vector of period effect in scale of log(OR) with length (number 
#       of periods)
# alpha0: WP-ICC alpha_0
# rho: decay rate from 0 to 1
#
# Output: a computed variance of treatment effect estimator
###########################################################################

sumVAR.exp.naive <- function(n, scheme, delta, beta, alpha0, rho) {
  t <- ncol(n)
  I <- nrow(n)
  trtSeq <- matrix(0, t - 1, t)
  trtSeq[upper.tri(trtSeq)] <- 1
  trtSeq <- trtSeq[rep(1:nrow(trtSeq), scheme),]
  if (sum(scheme) != I) {
    stop("invalid randomization scheme")
  }
  Model <- matrix(0, t + 1, t + 1)
  Omega <- matrix(0, t + 1, t + 1)
  
  # loop through each cluster
  for (i in 1:I) {
    # design matrix
    period <- rep(1:t)
    X <- trtSeq[i,]
    for (j in 1:t) {
      X <- cbind(X, as.numeric(period == j))
    }
    
    # unique marginal means
    gmu <- c(X %*% c(delta, beta))
    mu <- plogis(gmu)
    
    
    B = matrix(0, t, t)
    v = mu * (1 - mu)
    ni = n[i,]
    for (j in 1:(t - 1)) {
      B[j, j] <- (v[j] / ni[j]) * (1 + (ni[j] - 1) * alpha0)
      for (k in (j + 1):t) {
        B[j, k] <- sqrt(v[j] * v[k]) * alpha0 * rho ^ (abs(k - j))
      }
    }
    B[t, t] <- (v[t] / ni[t]) * (1 + (ni[t] - 1) * alpha0)
    B[lower.tri(B)] = t(B)[lower.tri(B)]
    # invV<-solve(B)
    
    D <- X * (mu * (1 - mu))
    invV1 <- diag(ni / (mu * (1 - mu)))
    Model <- Model + t(D) %*% invV1 %*% D
    Omega <- Omega + t(D) %*% invV1 %*% B %*% invV1 %*% D
  }
  invModel <- solve(Model)
  vardelta <- invModel %*% Omega %*% invModel
  return(vardelta[1, 1])
}


##########################################################################
# Construct the cluster-by-period sample size matrix of size (number of 
# clusters) x (number of periods) with between cluster imbalance
# 
# Input:
# I: number of clusters (I)
# t: number of periods (J)
# cv: coefficient of variation of the cluster-specific mean cluster-period 
#     size 
# cmean: mean cluster period size 
#
# Output: a matrix of cluster-by-period sample size with only between 
#         cluster imbalance parameterized by CV
##########################################################################

n.variable <- function(I, t, cv, cmean) {
  if (cv == 0) {
    return(matrix(rep(cmean, I * t), nrow = I, ncol = t))
  } else {
    N.raw <- rgamma(I, 1/cv^2, rate = 1/(cmean * cv^2))
    N <- as.integer(N.raw * (cmean/mean(N.raw)))
    N[N < 5] = 5
    return(matrix(rep(N, t), nrow = I, ncol = t))
  }
}

##########################################################################
# Construct the cluster-by-period sample size matrix of size (number of 
# clusters) x (number of periods) with no imbalance introduced
#
# Input:
# I: number of clusters (I)
# t: number of periods (J)
# cmean: mean cluster period size 
#
# Output: a matrix of cluster-by-period sample size with equal cluster-
#         period size
##########################################################################

n.equal <- function(I, t, cmean) {
  return(matrix(rep(cmean, I * t), nrow = I, ncol = t))
}

##########################################################################
# Construct the randomization scheme 
#
# Input:
# I: number of clusters (I)
# t: number of periods (J), assume t-1 to be a divisor of I in our main 
#    simulation 
#
# Output: a vector of number of clusters randomized at each step
##########################################################################


scheme <- function(I, t) {
  arm <- t - 1
  return(rep(I/arm, arm))
}

##########################################################################
# Truncated multinomial distribution 
#
# Input:
# n: sum of cluster-period sizes in a specific cluster
# p: probability vector of length (number of period) 
# a: lower bound for sampling
#
# Output: a vector of sampled cluster-period sizes for a specific cluster
##########################################################################

trmultinom <- function(n, p, a) {
  N <- as.vector(rmultinom(1, n, p))
  while (min(N) < a) {
    N <- as.vector(rmultinom(1, n, p))
  }
  return(N)
}

##########################################################################
# Construct the cluster-by-period sample size matrix of size (number of  
# clusters) x (number of periods) with both between cluster imbalance and 
# within-cluster imbalance (pattern 1 to pattern 4)
#
# Input:
# I: number of clusters (I)
# t: number of periods (J)
# cv: coefficient of variation of the cluster-specific mean cluster-period 
#     size 
# cmean: mean cluster period size 
# case: a integer ranges from 1 to 4 that specifies pattern 1 to pattern 4 
#       within-cluster imbalance
# p1: initial value of the probability vector (0.2 for J=3, 0.1 for J=5, 
#     0.05 for J=13, where J is the number of periods)
#
# Output: a matrix of cluster-by-period sample size with both between 
#         cluster imbalance parameterized by CV and the within cluster 
#         imbalance
##########################################################################

n.variable.twoway <- function(I, t, cv, cmean, case, p1) {
  if (cv == 0) {
    if (case == 1) {
      p <- rep(1/t, t)
      N <- rep(cmean, I)
      M <- NULL
      for (i in 1:I) {
        M <- c(M, trmultinom(n = t * N[i], p = p, a = 2))
      }
    }
    if (case == 2) {
      d = (2 - 2 * p1 * t)/(t * (t - 1))
      p <- seq(p1, p1 + (t - 1) * d, length.out = t)
      N <- rep(cmean, I)
      M <- NULL
      for (i in 1:I) {
        M <- c(M, trmultinom(n = t * N[i], p = p, a = 2))
      }
    }
    if (case == 3) {
      d = (2 - 2 * p1 * t)/(t * (t - 1))
      p <- seq(p1 + (t - 1) * d, p1, length.out = t)
      N <- rep(cmean, I)
      M <- NULL
      for (i in 1:I) {
        M <- c(M, trmultinom(n = t * N[i], p = p, a = 2))
      }
    }
    if (case == 4) {
      d = (2 - 2 * p1 * t)/(t * (t - 1))
      p <- seq(p1, p1 + (t - 1) * d, length.out = t)
      N <- rep(cmean, I)
      M <- NULL
      for (i in 1:I) {
        M <- c(M, trmultinom(n = t * N[i], p = sample(p), a = 2))
      }
    }
    return(matrix(M, nrow = I, ncol = t))
  } else {
    if (case == 1) {
      p <- rep(1/t, t)
      N.raw <- rgamma(I, 1/cv^2, rate = 1/(cmean * cv^2))
      N <- as.integer(N.raw * (cmean/mean(N.raw)))
      N[N < 5] = 5
      M <- NULL
      for (i in 1:I) {
        M <- c(M, trmultinom(n = t * N[i], p = p, a = 2))
      }
    }
    if (case == 2) {
      d = (2 - 2 * p1 * t)/(t * (t - 1))
      p <- seq(p1, p1 + (t - 1) * d, length.out = t)
      N.raw <- rgamma(I, 1/cv^2, rate = 1/(cmean * cv^2))
      N <- as.integer(N.raw * (cmean/mean(N.raw)))
      N[N < 5] = 5
      M <- NULL
      for (i in 1:I) {
        M <- c(M, trmultinom(n = t * N[i], p = p, a = 2))
      }
    }
    if (case == 3) {
      d = (2 - 2 * p1 * t)/(t * (t - 1))
      p <- seq(p1 + (t - 1) * d, p1, length.out = t)
      N.raw <- rgamma(I, 1/cv^2, rate = 1/(cmean * cv^2))
      N <- as.integer(N.raw * (cmean/mean(N.raw)))
      N[N < 5] = 5
      M <- NULL
      for (i in 1:I) {
        M <- c(M, trmultinom(n = t * N[i], p = p, a = 2))
      }
    }
    if (case == 4) {
      d = (2 - 2 * p1 * t)/(t * (t - 1))
      p <- seq(p1, p1 + (t - 1) * d, length.out = t)
      N.raw <- rgamma(I, 1/cv^2, rate = 1/(cmean * cv^2))
      N <- as.integer(N.raw * (cmean/mean(N.raw)))
      N[N < 5] = 5
      M <- NULL
      for (i in 1:I) {
        M <- c(M, trmultinom(n = t * N[i], p = p, a = 2))
      }
    }
    return(matrix(M, nrow = I, ncol = t))
  }
}