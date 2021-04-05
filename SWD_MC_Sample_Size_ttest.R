#############################################################
# Apply MC Sample size calculation algorithm to WEPT trial as
# described in section 7 and produce Table 3 in the main article
#############################################################

dir1 = "C:\\Users\\ASUS\\Desktop\\YCAS\\Simulation\\Nest Exchangeable Correlation"
dir2 = "C:\\Users\\ASUS\\Desktop\\YCAS\\Simulation\\Exponential Decay Correlation"
source(paste0(dir1,"/NEX_functions.R"))
source(paste0(dir2,"/ED_functions.R"))


##################################################################
# New function that is able to specify appropriate scheme when it  
# is impossible to set equal number of clusters in each treatment 
# sequence (see section 7)
#
# Input:
# I: number of clusters (I)
# t: number of periods (J), assume t-1 to be a divisor of I in our 
#    main simulation 
#
# Output: a vector of number of clusters randomized at each step
##################################################################

scheme <- function(I, t) {
  arm <- t - 1
  if (I %% arm == 0) {
    out <- rep(I / arm, arm)
  }
  else {
    I <- I - 2
    res <- I %% arm
    base <- (I - res) / arm
    add1 <- c(rep(1, res), rep(0, arm - res))
    add2 <- c(1, rep(0, arm - 2), 1)
    out <- rep(base, arm) + add1 + add2
  }
  return(out)
}


####################################################################
# Compute mean variance under designs of unequal cluster-period sizes
#
# Input:
# I: number of clusters (I)
# t: number of periods (J)
# cv: degree of between cluster imbalance, coefficient of variation
# cmean: mean cluster-period size
# beta: a vector of period effect in scale of log(OR) with length 
#       (number of periods)
# delta: treatment effect in scale of log(OR)
# Corr: specification of true correlation structure, "NEX" or "ED"
# ICC: a set of ICC parameter
# naive: 0 for true correctly specified working correlation structure,
#        1 for independence working correlation structure
# pattern: a integer ranges from 0 to 4 that specifies no within-
#          cluster imbalance orpattern 1 to pattern 4 within-cluster
#          imbalance
# p1: initial value of the probability vector (0.2 for J=3, 0.1 for 
#     J=5, 0.05 for J=13, where J is the number of periods)       
# seed: a seed to control reproducible output with default value 8888
# nsims: times for the simulation runs with default value 1000
#
# Output: mean variance under designs of unequal cluster-period sizes
#####################################################################

var_uneq <- function(I, t, cv, cmean, beta, delta, Corr, ICC, naive,
                     pattern, p1, nsims=1000, seed=8888) {
  
  var <- numeric(nsims)
  if (Corr == "NEX" & naive == 0 & pattern == 0) {
    for (i in 1:nsims) {
      set.seed(seed + i)
      var[i] <-
        sumVAR(n.variable(I, t, cv, cmean),
               scheme(I, t),
               delta,
               beta,
               c(ICC[1], ICC[2]))
    }
  }
  else if (Corr == "NEX" & naive == 1 & pattern == 0) {
    for (i in 1:nsims) {
      set.seed(seed + i)
      var[i] <-
        sumVAR.naive(n.variable(I, t, cv, cmean),
                     scheme(I, t),
                     delta,
                     beta,
                     c(ICC[1], ICC[2]))
    }
  }
  else if (Corr == "NEX" & naive == 0 & pattern != 0) {
    for (i in 1:nsims) {
      set.seed(seed + i)
      var[i] <-
        sumVAR(
          n.variable.twoway(I, t, cv, cmean, case = pattern, p1 = p1),
          scheme(I, t),
          delta,
          beta,
          c(ICC[1], ICC[2])
        )
    }
  }
  else if (Corr == "NEX" & naive == 1 & pattern != 0) {
    for (i in 1:nsims) {
      set.seed(seed + i)
      var[i] <-
        sumVAR.naive(
          n.variable.twoway(I, t, cv, cmean, case = pattern, p1 = p1),
          scheme(I, t),
          delta,
          beta,
          c(ICC[1], ICC[2])
        )
    }
  }
  else if (Corr == "ED" & naive == 0 & pattern == 0) {
    for (i in 1:nsims) {
      set.seed(seed + i)
      var[i] <-
        sumVAR.exp(n.variable(I, t, cv, cmean),
                   scheme(I, t),
                   delta,
                   beta,
                   ICC[1],
                   ICC[2])
    }
  }
  else if (Corr == "ED" & naive == 1 & pattern == 0) {
    for (i in 1:nsims) {
      set.seed(seed + i)
      var[i] <-
        sumVAR.exp.naive(n.variable(I, t, cv, cmean),
                         scheme(I, t),
                         delta,
                         beta,
                         ICC[1],
                         ICC[2])
    }
  }
  else if (Corr == "ED" & naive == 0 & pattern != 0) {
    for (i in 1:nsims) {
      set.seed(seed + i)
      var[i] <-
        sumVAR.exp(
          n.variable.twoway(I, t, cv, cmean, case = pattern, p1 = p1),
          scheme(I, t),
          delta,
          beta,
          ICC[1],
          ICC[2]
        )
    }
  }
  else if (Corr == "ED" & naive == 1 & pattern != 0) {
    for (i in 1:nsims) {
      set.seed(seed + i)
      var[i] <-
        sumVAR.exp.naive(
          n.variable.twoway(I, t, cv, cmean, case = pattern, p1 = p1),
          scheme(I, t),
          delta,
          beta,
          ICC[1],
          ICC[2]
        )
    }
  }
  var_mean <- mean(var)
  return(var_mean)
}

####################################################################
# Compute mean variance under designs of equal cluster-period sizes
#
# Input:
# I: number of clusters (I)
# t: number of periods (J)
# cmean: mean cluster-period size
# beta: a vector of period effect in scale of log(OR) with length 
#       (number of periods)
# delta: treatment effect in scale of log(OR)
# Corr: specification of true correlation structure, "NEX" or "ED"
# ICC: a set of ICC parameter
# naive: 0 for true correctly specified working correlation structure,  
#        1 for independence working correlation structure
#
# Output: mean variance under designs of equal cluster-period sizes
#####################################################################

var_eq <- function(I, t, cmean, beta, delta, Corr, ICC, naive) {
  if (Corr == "NEX" & naive == 0) {
    var <-
      sumVAR(n.equal(I, t, cmean), scheme(I, t), delta, beta, 
             c(ICC[1], ICC[2]))
  }
  else if (Corr == "NEX" & naive == 1) {
    var <-
      sumVAR.naive(n.equal(I, t, cmean), scheme(I, t), delta, beta, 
                   c(ICC[1], ICC[2]))
  }
  else if (Corr == "ED" & naive == 0) {
    var <-
      sumVAR.exp(n.equal(I, t, cmean), scheme(I, t), delta, beta, 
                 ICC[1], ICC[2])
  }
  else if (Corr == "ED" & naive == 1) {
    var <-
      sumVAR.exp.naive(n.equal(I, t, cmean), scheme(I, t), delta,
                       beta, ICC[1], ICC[2])
  }
  return(var)
}

######################################################################
# Apply MC sample size calculation algorithm to find optimal number of 
# cluster (I) when the cluster-period sizes are unequal
#
# Input:
# I: number of clusters (I)
# t: number of periods (J)
# cv: degree of between cluster imbalance, coefficient of variation
# cmean: mean cluster-period size
# beta: a vector of period effect in scale of log(OR) with length 
#       (number of periods)
# delta: treatment effect in scale of log(OR)
# Corr: specification of true correlation structure, "NEX" or "ED"
# ICC: a set of ICC parameter
# naive: 0 for true correctly specified working correlation structure,  
#        1 for independence working correlation structure
# pattern: a integer ranges from 0 to 4 that specifies no within-cluster 
#          imbalance or pattern 1 to pattern 4 within-cluster imbalance
# p1: initial value of the probability vector (0.2 for J=3, 0.1 for J=5, 
#     0.05 for J=13, where J is the number of periods)       
# e1: type I error rate
# e2: type II error rate
#
# Output: mean variance under designs of unequal cluster-period sizes
######################################################################

find_I <- function(I, t, cv, cmean, beta, delta, Corr, ICC, naive=0, 
                   pattern=0, p1=0.1, e1=0.05, e2=0.2) {
  left_data <- NULL
  right_data <- NULL
  left <- I
  var <-
    var_uneq(I = left, t, cv, cmean, beta, delta, Corr, ICC, naive, 
             pattern, p1)
  right <-
    ceiling((qt(e1 / 2, df = left - 2) + qt(e2, df = left - 2)) ^ 2 * 
              var * left / (delta ^ 2))
  print(c(left, right))
  left_data <- c(left_data, left)
  right_data <- c(right_data, right)
  while (left < right) {
    left <- right
    var <-
      var_uneq(I = left, t, cv, cmean, beta, delta, Corr, ICC, naive, 
               pattern, p1)
    right <-
      ceiling((qt(e1 / 2, df = left - 2) + qt(e2, df = left - 2)) ^ 2 *
                var * left / (delta ^ 2))
    print(c(left, right))
    left_data <- c(left_data, left)
    right_data <- c(right_data, right)
  }
  while (left >= right) {
    left <- left - 1
    var <-
      var_uneq(I = left, t, cv, cmean, beta, delta, Corr, ICC, naive, 
               pattern, p1)
    right <-
      ceiling((qt(e1 / 2, df = left - 2) + qt(e2, df = left - 2)) ^ 2 * 
                var * left / (delta ^ 2))
    print(c(left, right))
    left_data <- c(left_data, left)
    right_data <- c(right_data, right)
  }
  return(list(I = left_data[length(left_data) - 1], data = 
                cbind(left_data, right_data)))
}

#######################################################################
# Apply MC sample size calculation algorithm to find optimal number of 
# cluster (I) when the cluster-period sizes are equal
#
# Input:
# I: number of clusters (I)
# t: number of periods (J)
# cmean: mean cluster-period size
# beta: a vector of period effect in scale of log(OR) with length 
#       (number of periods)
# delta: treatment effect in scale of log(OR)
# Corr: specification of true correlation structure, "NEX" or "ED"
# ICC: a set of ICC parameter
# naive: 0 for true correctly specified working correlation structure,  
#        1 for independence working correlation structure
# e1: type I error rate
# e2: type II error rate
#
# Output: mean variance under designs of unequal cluster-period sizes
######################################################################


find_I_eq <- function(I, t, cmean, beta, delta, Corr, ICC, naive=0, 
                      e1=0.05, e2=0.2) {
  left_data <- NULL
  right_data <- NULL
  left <- I
  var <- var_eq(I = left, t, cmean, beta, delta, Corr, ICC, naive)
  right <-
    ceiling((qt(e1 / 2, df = left - 2) + qt(e2, df = left - 2)) ^ 2 * 
              var * left / (delta ^ 2))
  print(c(left, right))
  left_data <- c(left_data, left)
  right_data <- c(right_data, right)
  while (left < right) {
    left <- right
    var <- var_eq(I = left, t, cmean, beta, delta, Corr, ICC, naive)
    right <-
      ceiling((qt(e1 / 2, df = left - 2) + qt(e2, df = left - 2)) ^ 2 * 
                var * left / (delta ^ 2))
    print(c(left, right))
    left_data <- c(left_data, left)
    right_data <- c(right_data, right)
  }
  while (left >= right) {
    left <- left - 1
    var <- var_eq(I = left, t, cmean, beta, delta, Corr, ICC, naive)
    right <-
      ceiling((qt(e1 / 2, df = left - 2) + qt(e2, df = left - 2)) ^ 2 * 
                var * left / (delta ^ 2))
    print(c(left, right))
    left_data <- c(left_data, left)
    right_data <- c(right_data, right)
  }
  return(list(I = left_data[length(left_data) - 1], data = 
                cbind(left_data, right_data)))
}


#######################################
# Apply MC sample size algorithm and
# produce Table 3 in the main article
#
# See section 7 in the main article
#######################################


#mean cluster period size of WEPT trial
cmean <- 305

## Simple exchangeable

#beta <- c(-2.443, -2.454, -2.535, -2.609, -2.537)
beta <- rep(-2.5, 5)
#beta <- rep(-2.95,5)
#delta <- -0.14
delta <- log(0.7)


cv <- c(0, 0.25,0.75,1.25)
t <- 5

alpha0 = alpha1 =0.007
ICC <- c(alpha0, alpha1)

I <- 30

#var <- var_eq(I=12, t, cmean, beta, delta, Corr="NEX", ICC, naive=0)

#ceiling(((qnorm(0.05/2)+qnorm(0.2))^2*var*I)/(delta^2))

I1 <- find_I_eq(I, t, cmean, beta, delta, Corr="NEX", ICC)$I

I2 <- find_I_eq(I, t, cmean, beta, delta, Corr="NEX", ICC, naive=1)$I


#I1 <- 12

#I2 <- 31

EX_fixed_0 <- find_I(I1, t, cv[1], cmean, beta, delta, Corr="NEX", 
                     ICC)$I
EX_fixed_0_naive <- find_I(I2, t, cv[1], cmean, beta, delta, 
                           Corr="NEX", ICC, naive=1)$I
EX_p2_0 <- find_I(I1, t, cv[1], cmean, beta, delta, Corr="NEX", 
                  ICC, pattern=2)$I
EX_p2_0_naive <- find_I(I2, t, cv[1], cmean, beta, delta, Corr="NEX", 
                        ICC, naive=1, pattern=2)$I
EX_p4_0 <- find_I(I1, t, cv[1], cmean, beta, delta, Corr="NEX", ICC, 
                  pattern=4)$I
EX_p4_0_naive <- find_I(I2, t, cv[1], cmean, beta, delta, Corr="NEX", 
                        ICC, naive=1, pattern=4)$I

EX_fixed_1 <- find_I(I1, t, cv[2], cmean, beta, delta, Corr="NEX", 
                     ICC)$I
EX_fixed_1_naive <- find_I(I2, t, cv[2], cmean, beta, delta, 
                           Corr="NEX", ICC, naive=1)$I
EX_p2_1 <- find_I(I1, t, cv[2], cmean, beta, delta, Corr="NEX", ICC, 
                  pattern=2)$I
EX_p2_1_naive <- find_I(I2, t, cv[2], cmean, beta, delta, Corr="NEX", 
                        ICC, naive=1, pattern=2)$I
EX_p4_1 <- find_I(I1, t, cv[2], cmean, beta, delta, Corr="NEX", ICC, 
                  pattern=4)$I
EX_p4_1_naive <- find_I(I2, t, cv[2], cmean, beta, delta, Corr="NEX", 
                        ICC, naive=1, pattern=4)$I

EX_fixed_2 <- find_I(I1, t, cv[3], cmean, beta, delta, Corr="NEX",
                     ICC)$I
EX_fixed_2_naive <- find_I(I2, t, cv[3], cmean, beta, delta, 
                           Corr="NEX", ICC, naive=1)$I
EX_p2_2 <- find_I(I1, t, cv[3], cmean, beta, delta, Corr="NEX", ICC, 
                  pattern=2)$I
EX_p2_2_naive <- find_I(I2, t, cv[3], cmean, beta, delta, Corr="NEX", 
                        ICC, naive=1, pattern=2)$I
EX_p4_2 <- find_I(I1, t, cv[3], cmean, beta, delta, Corr="NEX", ICC, 
                  pattern=4)$I
EX_p4_2_naive <- find_I(I2, t, cv[3], cmean, beta, delta, Corr="NEX", 
                        ICC, naive=1, pattern=4)$I

EX_fixed_3 <- find_I(I1, t, cv[4], cmean, beta, delta, Corr="NEX", 
                     ICC)$I
EX_fixed_3_naive <- find_I(I2, t, cv[4], cmean, beta, delta, 
                           Corr="NEX", ICC, naive=1)$I
EX_p2_3 <- find_I(I1, t, cv[4], cmean, beta, delta, Corr="NEX", ICC, 
                  pattern=2)$I
EX_p2_3_naive <- find_I(I2, t, cv[4], cmean, beta, delta, Corr="NEX",
                        ICC, naive=1, pattern=2)$I
EX_p4_3 <- find_I(I1, t, cv[4], cmean, beta, delta, Corr="NEX", ICC, 
                  pattern=4)$I
EX_p4_3_naive <- find_I(I2, t, cv[4], cmean, beta, delta, Corr="NEX",
                        ICC, naive=1, pattern=4)$I




##NEX

#beta <- c(-2.446, -2.439, -2.495, -2.606, -2.535)
#delta <- -0.142
alpha0 <- 0.007
alpha1 <- 0.0035
ICC <- c(alpha0, alpha1)

I <- 30

I1 <- find_I_eq(I, t, cmean, beta, delta, Corr = "NEX", ICC)$I
I2 <- find_I_eq(I, t, cmean, beta, delta, Corr = "NEX", ICC, 
                naive = 1)$I

NEX_fixed_0 <- find_I(I1, t, cv[1], cmean, beta, delta, Corr="NEX", 
                      ICC)$I
NEX_fixed_0_naive <- find_I(I2, t, cv[1], cmean, beta, delta, Corr="NEX", 
                            ICC, naive=1)$I
NEX_p2_0 <- find_I(I1, t, cv[1], cmean, beta, delta, Corr="NEX", ICC, 
                   pattern=2)$I
NEX_p2_0_naive <- find_I(I2, t, cv[1], cmean, beta, delta, Corr="NEX", 
                         ICC, naive=1, pattern=2)$I
NEX_p4_0 <- find_I(I1, t, cv[1], cmean, beta, delta, Corr="NEX", ICC, 
                   pattern=4)$I
NEX_p4_0_naive <- find_I(I2, t, cv[1], cmean, beta, delta, Corr="NEX", 
                         ICC, naive=1, pattern=4)$I

NEX_fixed_1 <- find_I(I1, t, cv[2], cmean, beta, delta, Corr="NEX", 
                      ICC)$I
NEX_fixed_1_naive <- find_I(I2, t, cv[2], cmean, beta, delta, Corr="NEX", 
                            ICC, naive=1)$I
NEX_p2_1 <- find_I(I1, t, cv[2], cmean, beta, delta, Corr="NEX", ICC, 
                   pattern=2)$I
NEX_p2_1_naive <- find_I(I2, t, cv[2], cmean, beta, delta, Corr="NEX", 
                         ICC, naive=1, pattern=2)$I
NEX_p4_1 <- find_I(I1, t, cv[2], cmean, beta, delta, Corr="NEX", ICC, 
                   pattern=4)$I
NEX_p4_1_naive <- find_I(I2, t, cv[2], cmean, beta, delta, Corr="NEX", 
                         ICC, naive=1, pattern=4)$I

NEX_fixed_2 <- find_I(I1, t, cv[3], cmean, beta, delta, Corr="NEX", 
                      ICC)$I
NEX_fixed_2_naive <- find_I(I2, t, cv[3], cmean, beta, delta, 
                            Corr="NEX", ICC, naive=1)$I
NEX_p2_2 <- find_I(I1, t, cv[3], cmean, beta, delta, Corr="NEX", ICC, 
                   pattern=2)$I
NEX_p2_2_naive <- find_I(I2, t, cv[3], cmean, beta, delta, Corr="NEX", 
                         ICC, naive=1, pattern=2)$I
NEX_p4_2 <- find_I(I1, t, cv[3], cmean, beta, delta, Corr="NEX", ICC, 
                   pattern=4)$I
NEX_p4_2_naive <- find_I(I2, t, cv[3], cmean, beta, delta, Corr="NEX", 
                         ICC, naive=1, pattern=4)$I

NEX_fixed_3 <- find_I(I1, t, cv[4], cmean, beta, delta, Corr="NEX", 
                      ICC)$I
NEX_fixed_3_naive <- find_I(I2, t, cv[4], cmean, beta, delta, Corr="NEX", 
                            ICC, naive=1)$I
NEX_p2_3 <- find_I(I1, t, cv[4], cmean, beta, delta, Corr="NEX", ICC, 
                   pattern=2)$I
NEX_p2_3_naive <- find_I(I2, t, cv[4], cmean, beta, delta, Corr="NEX", 
                         ICC, naive=1, pattern=2)$I
NEX_p4_3 <- find_I(I1, t, cv[4], cmean, beta, delta, Corr="NEX", ICC,
                   pattern=4)$I
NEX_p4_3_naive <- find_I(I2, t, cv[4], cmean, beta, delta, Corr="NEX", 
                         ICC, naive=1, pattern=4)$I

##EXP

#beta <- c(-2.437, -2.444, -2.508, -2.613, -2.552)
#delta <- -0.124
alpha0 <- 0.007
rho <- 0.7
ICC <- c(alpha0, rho)

I <- 30

I1 <- find_I_eq(I, t, cmean, beta, delta, Corr="ED", ICC)$I
I2 <- find_I_eq(I, t, cmean, beta, delta, Corr="ED", ICC, naive=1)$I

ED_fixed_0 <- find_I(I1, t, cv[1], cmean, beta, delta, Corr="ED", 
                     ICC)$I
ED_fixed_0_naive <- find_I(I2, t, cv[1], cmean, beta, delta, Corr="ED", 
                           ICC, naive=1)$I
ED_p2_0 <- find_I(I1, t, cv[1], cmean, beta, delta, Corr="ED", ICC, 
                  pattern=2)$I
ED_p2_0_naive <- find_I(I2, t, cv[1], cmean, beta, delta, Corr="ED", ICC,
                        naive=1, pattern=2)$I
ED_p4_0 <- find_I(I1, t, cv[1], cmean, beta, delta, Corr="ED", ICC, 
                  pattern=4)$I
ED_p4_0_naive <- find_I(I2, t, cv[1], cmean, beta, delta, Corr="ED", ICC,
                        naive=1, pattern=4)$I

ED_fixed_1 <- find_I(I1, t, cv[2], cmean, beta, delta, Corr="ED", 
                     ICC)$I
ED_fixed_1_naive <- find_I(I2, t, cv[2], cmean, beta, delta, Corr="ED",
                           ICC, naive=1)$I
ED_p2_1 <- find_I(I1, t, cv[2], cmean, beta, delta, Corr="ED", ICC, 
                  pattern=2)$I
ED_p2_1_naive <- find_I(I2, t, cv[2], cmean, beta, delta, Corr="ED", 
                        ICC, naive=1, pattern=2)$I
ED_p4_1 <- find_I(I1, t, cv[2], cmean, beta, delta, Corr="ED", ICC,
                  pattern=4)$I
ED_p4_1_naive <- find_I(I2, t, cv[2], cmean, beta, delta, Corr="ED", 
                        ICC, naive=1, pattern=4)$I

ED_fixed_2 <- find_I(I1, t, cv[3], cmean, beta, delta, Corr="ED", 
                     ICC)$I
ED_fixed_2_naive <- find_I(I2, t, cv[3], cmean, beta, delta, Corr="ED", 
                           ICC, naive=1)$I
ED_p2_2 <- find_I(I1, t, cv[3], cmean, beta, delta, Corr="ED", ICC, 
                  pattern=2)$I
ED_p2_2_naive <- find_I(I2, t, cv[3], cmean, beta, delta, Corr="ED", 
                        ICC, naive=1, pattern=2)$I
ED_p4_2 <- find_I(I1, t, cv[3], cmean, beta, delta, Corr="ED", ICC, 
                  pattern=4)$I
ED_p4_2_naive <- find_I(I2, t, cv[3], cmean, beta, delta, Corr="ED",
                        ICC, naive=1, pattern=4)$I

ED_fixed_3 <- find_I(I1, t, cv[4], cmean, beta, delta, Corr="ED", 
                     ICC)$I
ED_fixed_3_naive <- find_I(I2, t, cv[4], cmean, beta, delta, Corr="ED",
                           ICC, naive=1)$I
ED_p2_3 <- find_I(I1, t, cv[4], cmean, beta, delta, Corr="ED", ICC,
                  pattern=2)$I
ED_p2_3_naive <- find_I(I2, t, cv[4], cmean, beta, delta, Corr="ED",
                        ICC, naive=1, pattern=2)$I
ED_p4_3 <- find_I(I1, t, cv[4], cmean, beta, delta, Corr="ED", ICC,
                  pattern=4)$I
ED_p4_3_naive <- find_I(I2, t, cv[4], cmean, beta, delta, Corr="ED", 
                        ICC, naive=1, pattern=4)$I


##Table

EX_data <- rbind(c(EX_fixed_0, EX_p2_0, EX_p4_0, EX_fixed_0_naive, 
                   EX_p2_0_naive, EX_p4_0_naive),
                 c(EX_fixed_1, EX_p2_1, EX_p4_1, EX_fixed_1_naive,
                   EX_p2_1_naive, EX_p4_1_naive),
                 c(EX_fixed_2, EX_p2_2, EX_p4_2, EX_fixed_2_naive, 
                   EX_p2_2_naive, EX_p4_2_naive),
                 c(EX_fixed_3, EX_p2_3, EX_p4_3, EX_fixed_3_naive,
                   EX_p2_3_naive, EX_p4_3_naive))

NEX_data <- rbind(c(NEX_fixed_0, NEX_p2_0, NEX_p4_0, NEX_fixed_0_naive, 
                    NEX_p2_0_naive, NEX_p4_0_naive),
                  c(NEX_fixed_1, NEX_p2_1, NEX_p4_1, NEX_fixed_1_naive,
                    NEX_p2_1_naive, NEX_p4_1_naive),
                  c(NEX_fixed_2, NEX_p2_2, NEX_p4_2, NEX_fixed_2_naive, 
                    NEX_p2_2_naive, NEX_p4_2_naive),
                  c(NEX_fixed_3, NEX_p2_3, NEX_p4_3, NEX_fixed_3_naive,
                    NEX_p2_3_naive, NEX_p4_3_naive))

ED_data <- rbind(c(ED_fixed_0, ED_p2_0, ED_p4_0, ED_fixed_0_naive, 
                   ED_p2_0_naive, ED_p4_0_naive),
                 c(ED_fixed_1, ED_p2_1, ED_p4_1, ED_fixed_1_naive, 
                   ED_p2_1_naive, ED_p4_1_naive),
                 c(ED_fixed_2, ED_p2_2, ED_p4_2, ED_fixed_2_naive, 
                   ED_p2_2_naive, ED_p4_2_naive),
                 c(ED_fixed_3, ED_p2_3, ED_p4_3, ED_fixed_3_naive, 
                   ED_p2_3_naive, ED_p4_3_naive))

out_data <- data.frame(rbind(EX_data, NEX_data, ED_data))

colnames(out_data) <- c("fixed", "pattern 2", "pattern 4", "fixed", 
                        "pattern 2", "pattern 4")

#Path for saving the table
write.csv(out_data, "C:\\Users\\ASUS\\Desktop\\YCAS\\Application\\Washington_EPT_table_ttest.csv")