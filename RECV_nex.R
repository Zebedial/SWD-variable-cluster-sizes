#########################################################################
# Obtain Original Data for Relative Efficiency vs. CV Plots When the True 
# Correlation Structure is Nested Exchangeable (NEX)
#########################################################################

#For Running in HPC
mainDir = '/ysm-gpfs/home/zt223/ycas'

#Local path
#setwd("C:\\Users\\ASUS\\Desktop\\YCAS\\Simulation\\Nest Exchangeable Correlation")
source("./NEX_functions.R")


#########################################################################
# Raw data of simulated Relative efficiency given a whole set of design 
# parameters and marginal model parameters (See input below), when no 
# within-cluster imbalance is introduced. "REsumtable1" is used under the 
# case that the working correlation structure correctly specifies the true 
# one (NEX); "REsumtable2" is used when the working correlation model 
# misspecifies the true correlation structure by an indepdendence 
# correlation model (IND)
#
# Input:
# I: number of clusters (I)
# t: number of periods (J)
# cv.num: a vector of a finite sequence of CVs range from 0 (no between 
#         period imbalance) to 1.5
# cmean: mean cluster-period size
# beta: a vector of period effect in scale of log(OR) with length (number 
#       of periods)
# alpha0: WP-ICC
# alpha1: BP-ICC
# delta: treatment effect in scale of log(OR)
# seed: a seed to control reproducible output with default value 8888
# nsims: times for the simulation runs with default value 1000
#
# Output: a table whose rows record 1000 REs computed based on a specifc 
#         CV and other required parameters
#########################################################################

REsumtable1 <- function(I, t, cv.num, cmean, beta, alpha0, alpha1, delta, 
                        seed = 8888, nsims = 1000) {
  
  RE <- matrix(NA, length(cv.num), nsims)
  for (i in 1:length(cv.num)) {
    cv <- cv.num[i]
    for (j in 1:nsims) {
      set.seed(seed + j)
      RE[i, j] <-
        sumVAR(n.equal(I, t, cmean),
               scheme(I, t),
               delta,
               beta,
               c(alpha0, alpha1)) / sumVAR(n.variable(I, t, cv, cmean),
                                           scheme(I, t),
                                           delta,
                                           beta,
                                           c(alpha0, alpha1))
    }
    #print(i)
  }
  return(RE)
}

REsumtable2 <- function(I, t, cv.num, cmean, beta, alpha0, alpha1, delta, 
                        seed = 8888, nsims = 1000) {
  
  RE <- matrix(NA, length(cv.num), nsims)
  for (i in 1:length(cv.num)) {
    cv <- cv.num[i]
    for (j in 1:nsims) {
      set.seed(seed + j)
      RE[i, j] <-
        sumVAR.naive(n.equal(I, t, cmean),
                     scheme(I, t),
                     delta,
                     beta,
                     c(alpha0, 
                       alpha1)) / sumVAR.naive(n.variable(I, t, cv, 
                                                          cmean),
                                               scheme(I, t),
                                               delta,
                                               beta,
                                               c(alpha0, alpha1))
                                                         
    }
    #print(i)
  }
  return(RE)
}

##########################################################################
# Construct datasets used to produce the RE vs CV plots, e.g. Figure 2 in 
# the main article, when no within-cluster imbalance is introduced.
# "REvsCV.table" is used under the case that the working correlation structure 
# correctly specifies the true one (NEX); "REvsCV.naive.table" is used when  
# the working correlation model misspecifies the true correlation structure 
# by an indepdendence correlation model (IND)
#
# Input:
# I: number of clusters (I)
# t: number of periods (J)
# cv.num: a vector of a finite sequence of CVs range from 0 (no between 
#         period imbalance) to 1.5
# cmean: mean cluster-period size
# base: baseline prevalence
# end: end prevalence, e.g. equal to baseline prevalence means constant 
#      prevalence
# alpha0: WP-ICC 
# alpha1: a vector of BP-ICCS given the WP-ICC (three in the main simulation)
# delta: treatment effect in scale of log(OR)
#
# Output: a table with minimum, lower quartile (25 percentile), median, 
#         mean, upper quartile (75 percentile), and maximum of RE as a 
#         function of CV and WP-ICC given a set of other required parameters
##########################################################################

REvsCV.table <- function(I, t, cv.num, cmean, base, end, alpha0, alpha1, 
                         delta) {
  p = seq(base, end, length.out = t)
  beta = log(p / (1 - p))
  REsummary <- NULL
  for (i in 1:length(alpha1)) {
    REbase <- REsumtable1(I, t, cv.num, cmean, beta, alpha0, alpha1[i], 
                          delta)
    RE_min <- apply(REbase, 1, quantile)[1, ]
    RE_25 <- apply(REbase, 1, quantile)[2, ]
    RE_median <- apply(REbase, 1, quantile)[3, ]
    RE_mean <- apply(REbase, 1, mean)
    RE_75 <- apply(REbase, 1, quantile)[4, ]
    RE_max <- apply(REbase, 1, quantile)[5, ]
    REsummary <-
      rbind(REsummary,
            cbind(cv.num, RE_min, RE_25, RE_median, RE_mean, RE_75, RE_max))
  }
  ICC <- rep(c(1, 2, 3), each = length(cv.num))
  REtable <- data.frame(cbind(REsummary, ICC))
  REtable$ICC <- as.factor(REtable$ICC)
  return(REtable)
}

REvsCV.naive.table <- function(I, t, cv.num, cmean, base, end, alpha0, 
                               alpha1, delta) {
  p = seq(base, end, length.out = t)
  beta = log(p / (1 - p))
  REsummary <- NULL
  for (i in 1:length(alpha1)) {
    REbase <- REsumtable2(I, t, cv.num, cmean, beta, alpha0, alpha1[i], 
                          delta)
    RE_min <- apply(REbase, 1, quantile)[1, ]
    RE_25 <- apply(REbase, 1, quantile)[2, ]
    RE_median <- apply(REbase, 1, quantile)[3, ]
    RE_mean <- apply(REbase, 1, mean)
    RE_75 <- apply(REbase, 1, quantile)[4, ]
    RE_max <- apply(REbase, 1, quantile)[5, ]
    REsummary <-
      rbind(REsummary,
            cbind(cv.num, RE_min, RE_25, RE_median, RE_mean, RE_75, RE_max))
  }
  ICC <- rep(c(1, 2, 3), each = length(cv.num))
  REtable <- data.frame(cbind(REsummary, ICC))
  REtable$ICC <- as.factor(REtable$ICC)
  return(REtable)
}

############################################################################
# Raw data of simulated Relative efficiency given a whole set of design 
# parameters and marginal model parameters (See input below), when within-cluster 
# imbalance is introduced. "REvsCV.twoway1" is used under the case that the 
# working correlation structure correctly specifies the true one (NEX); 
# "REvsCV.twoway2" is used when the working correlation model misspecifies 
# the true correlation structure by an indepdendence correlation model (IND)
#
# Input:
# I: number of clusters (I)
# t: number of periods (J)
# cv.num: a vector of a finite sequence of CVs range from 0 (no between period 
#         imbalance) to 1.5
# cmean: mean cluster-period size
# beta: a vector of period effect in scale of log(OR) with length (number of 
#       periods)
# alpha0: WP-ICC
# alpha1: BP-ICC
# delta: treatment effect in scale of log(OR)
# case: a integer ranges from 1 to 4 that specifies pattern 1 to pattern 4 
#       within-cluster imbalance
# p1: initial value of the probability vector (0.2 for J=3, 0.1 for J=5, 0.05
#     for J=13, where J is the number of periods)
# seed: a seed to control reproducible output with default value 8888
# nsims: times for the simulation runs with default value 1000
#
# Output: a table whose rows record 1000 REs computed based on a specifc CV 
#         and other required parameters
#############################################################################

REvsCV.twoway1 <- function(I, t, cv.num, cmean, beta, alpha0, alpha1, delta, 
                           case, p1, seed = 8888, nsims = 1000) {
  
  RE <- matrix(NA, length(cv.num), nsims)
  for (i in 1:length(cv.num)) {
    cv <- cv.num[i]
    for (j in 1:nsims) {
      set.seed(seed + j)
      RE[i, j] <-
        sumVAR(n.equal(I, t, cmean),
               scheme(I, t),
               delta,
               beta,
               c(alpha0, alpha1)) / sumVAR(
                 n.variable.twoway(I, t, cv, cmean, case, p1),
                 scheme(I, t),
                 delta,
                 beta,
                 c(alpha0, alpha1)
               )
    }
    #print(i)
  }
  return(RE)
}

REvsCV.twoway2 <- function(I, t, cv.num, cmean, beta, alpha0, alpha1, delta, 
                           case, p1, seed=8888, nsims=1000) {
  set.seed(seed)
  RE <- matrix(NA, length(cv.num), nsims)
  for (i in 1:length(cv.num)) {
    cv <- cv.num[i]
    for (j in 1:nsims) {
      set.seed(seed + j)
      RE[i, j] <-
        sumVAR.naive(n.equal(I, t, cmean),
                     scheme(I, t),
                     delta,
                     beta,
                     c(alpha0, alpha1)) / sumVAR.naive(
                       n.variable.twoway(I, t, cv, cmean, case, p1),
                       scheme(I, t),
                       delta,
                       beta,
                       c(alpha0, alpha1)
                     )
    }
    #print(i)
  }
  return(RE)
}


########################################################################
# Construct datasets used to produce the RE vs CV plots, e.g. Figure 2  
# in the main article, when a given pattern of  within-cluster imbalance 
# is introduced. "REvsCV.twoway.table" is used under the case that the  
# working correlation structure correctly specifies the true one (NEX); 
# "REvsCV.twoway.naive.table" is used when the working correlation model 
# misspecifies the true correlation structure by an indepdendence 
# correlation model (IND)
#
# Input:
# I: number of clusters (I)
# t: number of periods (J)
# cv.num: a vector of a finite sequence of CVs range from 0 (no between 
#         period imbalance) to 1.5
# cmean: mean cluster-period size
# base: baseline prevalence
# end: end prevalence, e.g. equal to baseline prevalence means constant 
#      prevalence
# alpha0: WP-ICC 
# alpha1: a vector of BP-ICCS given the WP-ICC (three in the main simulation)
# delta: treatment effect in scale of log(OR)
# case: a integer ranges from 1 to 4 that specifies pattern 1 to pattern 
#       4 within-cluster imbalance
# p1: initial value of the probability vector (0.2 for J=3, 0.1 for J=5, 
#     0.05 for J=13, where J is the number of periods)
#
# Output: a table with minimum, lower quartile (25 percentile), median, 
#         mean, upper quartile (75 percentile), and maximum of RE as a 
#         function of CV and WP-ICC given a set of other required parameters
#########################################################################
REvsCV.twoway.table <- function(I, t, cv.num, cmean, base, end, alpha0, 
                                alpha1, delta, case, p1=0.1) {
  p = seq(base, end, length.out = t)
  beta = log(p / (1 - p))
  REsummary <- NULL
  for (i in 1:length(alpha1)) {
    REbase <-
      REvsCV.twoway1(I, t, cv.num, cmean, beta, alpha0, alpha1[i], delta, 
                     case, p1)
    RE_min <- apply(REbase, 1, quantile)[1, ]
    RE_25 <- apply(REbase, 1, quantile)[2, ]
    RE_median <- apply(REbase, 1, quantile)[3, ]
    RE_mean <- apply(REbase, 1, mean)
    RE_75 <- apply(REbase, 1, quantile)[4, ]
    RE_max <- apply(REbase, 1, quantile)[5, ]
    REsummary <-
      rbind(REsummary,
            cbind(cv.num, RE_min, RE_25, RE_median, RE_mean, RE_75, RE_max))
  }
  ICC <- rep(c(1, 2, 3), each = length(cv.num))
  REtable <- data.frame(cbind(REsummary, ICC))
  return(REtable)
}

REvsCV.twoway.naive.table <- function(I, t, cv.num, cmean, base, end, 
                                      alpha0, alpha1, delta, case, p1=0.1) {
  p = seq(base, end, length.out = t)
  beta = log(p / (1 - p))
  REsummary <- NULL
  for (i in 1:length(alpha1)) {
    REbase <-
      REvsCV.twoway2(I, t, cv.num, cmean, beta, alpha0, alpha1[i], delta, 
                     case, p1)
    RE_min <- apply(REbase, 1, quantile)[1, ]
    RE_25 <- apply(REbase, 1, quantile)[2, ]
    RE_median <- apply(REbase, 1, quantile)[3, ]
    RE_mean <- apply(REbase, 1, mean)
    RE_75 <- apply(REbase, 1, quantile)[4, ]
    RE_max <- apply(REbase, 1, quantile)[5, ]
    REsummary <-
      rbind(REsummary,
            cbind(cv.num, RE_min, RE_25, RE_median, RE_mean, RE_75, RE_max))
  }
  ICC <- rep(c(1, 2, 3), each = length(cv.num))
  REtable <- data.frame(cbind(REsummary, ICC))
  return(REtable)
}

#Path for saving the data sets
#setwd("C:\\Users\\ASUS\\Desktop\\YCAS\\result tables\\NEX_plots")

I <- c(12, 24, 48, 96)
t <- c(3, 5, 13)
cv <- seq(0, 1.5, by = 0.05)
cmean <- c(50, 100, 300)
base = 0.3
end = 0.3
alpha0 <- c(0.05, 0.1, 0.2)
delta = log(0.35)

#################################
## Produce data used to construct
## Figure 2 in the main article
#################################

REvsCV_table_1.1 <-
  REvsCV.table(
    I = I[1],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta
  )
write.csv(REvsCV_table_1.1, file = paste0("./REvsCV_1.1.csv"))

REvsCV_table_1.4 <-
  REvsCV.table(
    I = I[4],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta
  )
write.csv(REvsCV_table_1.4, file = paste0("./REvsCV_1.4.csv"))

REvsCV.naive_table_1.1 <-
  REvsCV.naive.table(
    I = I[1],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta
  )
write.csv(REvsCV.naive_table_1.1,
          file = paste0("./REvsCV.naive_1.1.csv"))

REvsCV.naive_table_1.4 <-
  REvsCV.naive.table(
    I = I[4],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta
  )
write.csv(REvsCV.naive_table_1.4, 
          file = paste0("./REvsCV.naive_1.4.csv"))


#################################
## Produce data used to construct
## Figure 3 in the main article
#################################

#pattern4
REvsCV_table_1.4.1 <-
  REvsCV.twoway.table(
    I = I[1],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta,
    case = 4,
    p1 = 0.1
  )
write.csv(REvsCV_table_1.4.1,
          file = paste0("./REvsCV_1.4.1.csv"))

REvsCV_table_1.4.4 <-
  REvsCV.twoway.table(
    I = I[4],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta,
    case = 4,
    p1 = 0.1
  )
write.csv(REvsCV_table_1.4.4, 
          file = paste0("./REvsCV_1.4.4.csv"))

REvsCV.naive_table_1.4.1 <-
  REvsCV.twoway.naive.table(
    I = I[1],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta,
    case = 4,
    p1 = 0.1
  )
write.csv(REvsCV.naive_table_1.4.1,
          file = paste0("./REvsCV.naive_1.4.1.csv"))

REvsCV.naive_table_1.4.4 <-
  REvsCV.twoway.naive.table(
    I = I[4],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta,
    case = 4,
    p1 = 0.1
  )
write.csv(REvsCV.naive_table_1.4.4,
          file = paste0("./REvsCV.naive_1.4.4.csv"))

#################################
## Produce data used to construct
## Web Figure 1
#################################

#pattern1
REvsCV_table_1.1.1 <-
  REvsCV.twoway.table(
    I = I[1],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta,
    case = 1,
    p1 = 0.1
  )
write.csv(REvsCV_table_1.1.1,
          file = paste0("./REvsCV_1.1.1.csv"))

REvsCV_table_1.1.4 <-
  REvsCV.twoway.table(
    I = I[4],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta,
    case = 1,
    p1 = 0.1
  )
write.csv(REvsCV_table_1.1.4,
          file = paste0("./REvsCV_1.1.4.csv"))

REvsCV.naive_table_1.1.1 <-
  REvsCV.twoway.naive.table(
    I = I[1],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta,
    case = 1,
    p1 = 0.1
  )
write.csv(REvsCV.naive_table_1.1.1,
          file = paste0("./REvsCV.naive_1.1.1.csv"))

REvsCV.naive_table_1.1.4 <-
  REvsCV.twoway.naive.table(
    I = I[4],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta,
    case = 1,
    p1 = 0.1
  )
write.csv(REvsCV.naive_table_1.1.4,
          file = paste0("./REvsCV.naive_1.1.4.csv"))

#################################
## Produce data used to construct
## Web Figure 2
#################################

#pattern2
REvsCV_table_1.2.1 <-
  REvsCV.twoway.table(
    I = I[1],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta,
    case = 2,
    p1 = 0.1
  )
write.csv(REvsCV_table_1.2.1, 
          file = paste0("./REvsCV_1.2.1.csv"))

REvsCV_table_1.2.4 <-
  REvsCV.twoway.table(
    I = I[4],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta,
    case = 2,
    p1 = 0.1
  )
write.csv(REvsCV_table_1.2.4, 
          file = paste0("./REvsCV_1.2.4.csv"))

REvsCV.naive_table_1.2.1 <-
  REvsCV.twoway.naive.table(
    I = I[1],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta,
    case = 2,
    p1 = 0.1
  )
write.csv(REvsCV.naive_table_1.2.1,
          file = paste0("./REvsCV.naive_1.2.1.csv"))

REvsCV.naive_table_1.2.4 <-
  REvsCV.twoway.naive.table(
    I = I[4],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta,
    case = 2,
    p1 = 0.1
  )
write.csv(REvsCV.naive_table_1.2.4,
          file = paste0("./REvsCV.naive_1.2.4.csv"))

#################################
## Produce data used to construct
## Web Figure 3
#################################

#pattern3
REvsCV_table_1.3.1 <-
  REvsCV.twoway.table(
    I = I[1],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta,
    case = 3,
    p1 = 0.1
  )
write.csv(REvsCV_table_1.3.1, 
          file = paste0("./REvsCV_1.3.1.csv"))

REvsCV_table_1.3.4 <-
  REvsCV.twoway.table(
    I = I[4],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta,
    case = 3,
    p1 = 0.1
  )
write.csv(REvsCV_table_1.3.4, 
          file = paste0("./REvsCV_1.3.4.csv"))

REvsCV.naive_table_1.3.1 <-
  REvsCV.twoway.naive.table(
    I = I[1],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta,
    case = 3,
    p1 = 0.1
  )
write.csv(REvsCV.naive_table_1.3.1,
          file = paste0("./REvsCV.naive_1.3.1.csv"))

REvsCV.naive_table_1.3.4 <-
  REvsCV.twoway.naive.table(
    I = I[4],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta,
    case = 3,
    p1 = 0.1
  )
write.csv(REvsCV.naive_table_1.3.4, 
          file = paste0("./REvsCV.naive_1.3.4.csv"))

#################################
## Produce data used to construct
## Web Figure 12
#################################

REvsCV_table_1.2 <-
  REvsCV.table(
    I = I[2],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta
  )
write.csv(REvsCV_table_1.2, 
          file = paste0("./REvsCV_1.2.csv"))

REvsCV_table_1.3 <-
  REvsCV.table(
    I = I[3],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta
  )
write.csv(REvsCV_table_1.3, 
          file = paste0("./REvsCV_1.3.csv"))

REvsCV.naive_table_1.2 <-
  REvsCV.naive.table(
    I = I[2],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta
  )
write.csv(REvsCV.naive_table_1.2, 
          file = paste0("./REvsCV.naive_1.2.csv"))

REvsCV.naive_table_1.3 <-
  REvsCV.naive.table(
    I = I[3],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta
  )
write.csv(REvsCV.naive_table_1.3, 
          file = paste0("./REvsCV.naive_1.3.csv"))













#################################
## Produce data used to conduct
## experiments that do not showed
## in the article
#################################

###############################
## No Within cluster variability
###############################

REvsCV_table_1.5 <-
  REvsCV.table(
    I = I[3],
    t = t[1],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta
  )
REvsCV_table_1.6 <-
  REvsCV.table(
    I = I[3],
    t = t[3],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta
  )

REvsCV_table_2.1 <-
  REvsCV.table(
    I = I[1],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[2],
    alpha1 = c(0.001, alpha0[2] / 2, alpha0[2]),
    delta = delta
  )
REvsCV_table_2.2 <-
  REvsCV.table(
    I = I[2],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[2],
    alpha1 = c(0.001, alpha0[2] / 2, alpha0[2]),
    delta = delta
  )
REvsCV_table_2.3 <-
  REvsCV.table(
    I = I[3],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[2],
    alpha1 = c(0.001, alpha0[2] / 2, alpha0[2]),
    delta = delta
  )
REvsCV_table_2.4 <-
  REvsCV.table(
    I = I[4],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[2],
    alpha1 = c(0.001, alpha0[2] / 2, alpha0[2]),
    delta = delta
  )
REvsCV_table_2.5 <-
  REvsCV.table(
    I = I[3],
    t = t[1],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[2],
    alpha1 = c(0.001, alpha0[2] / 2, alpha0[2]),
    delta = delta
  )
REvsCV_table_2.6 <-
  REvsCV.table(
    I = I[3],
    t = t[3],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[2],
    alpha1 = c(0.001, alpha0[2] / 2, alpha0[2]),
    delta = delta
  )

REvsCV_table_3.1 <-
  REvsCV.table(
    I = I[1],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[3],
    alpha1 = c(0.001, alpha0[3] / 2, alpha0[3]),
    delta = delta
  )
REvsCV_table_3.2 <-
  REvsCV.table(
    I = I[2],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[3],
    alpha1 = c(0.001, alpha0[3] / 2, alpha0[3]),
    delta = delta
  )
REvsCV_table_3.3 <-
  REvsCV.table(
    I = I[3],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[3],
    alpha1 = c(0.001, alpha0[3] / 2, alpha0[3]),
    delta = delta
  )
REvsCV_table_3.4 <-
  REvsCV.table(
    I = I[4],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[3],
    alpha1 = c(0.001, alpha0[3] / 2, alpha0[3]),
    delta = delta
  )
REvsCV_table_3.5 <-
  REvsCV.table(
    I = I[3],
    t = t[1],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[3],
    alpha1 = c(0.001, alpha0[3] / 2, alpha0[3]),
    delta = delta
  )
REvsCV_table_3.6 <-
  REvsCV.table(
    I = I[3],
    t = t[3],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[3],
    alpha1 = c(0.001, alpha0[3] / 2, alpha0[3]),
    delta = delta
  )


write.csv(REvsCV_table_1.5, file = paste0("./REvsCV_1.5.csv"))
write.csv(REvsCV_table_1.6, file = paste0("./REvsCV_1.6.csv"))

write.csv(REvsCV_table_2.1, file = paste0("./REvsCV_2.1.csv"))
write.csv(REvsCV_table_2.2, file = paste0("./REvsCV_2.2.csv"))
write.csv(REvsCV_table_2.3, file = paste0("./REvsCV_2.3.csv"))
write.csv(REvsCV_table_2.4, file = paste0("./REvsCV_2.4.csv"))
write.csv(REvsCV_table_2.5, file = paste0("./REvsCV_2.5.csv"))
write.csv(REvsCV_table_2.6, file = paste0("./REvsCV_2.6.csv"))

write.csv(REvsCV_table_3.1, file = paste0("./REvsCV_3.1.csv"))
write.csv(REvsCV_table_3.2, file = paste0("./REvsCV_3.2.csv"))
write.csv(REvsCV_table_3.3, file = paste0("./REvsCV_3.3.csv"))
write.csv(REvsCV_table_3.4, file = paste0("./REvsCV_3.4.csv"))
write.csv(REvsCV_table_3.5, file = paste0("./REvsCV_3.5.csv"))
write.csv(REvsCV_table_3.6, file = paste0("./REvsCV_3.6.csv"))


REvsCV.naive_table_1.5 <-
  REvsCV.naive.table(
    I = I[3],
    t = t[1],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta
  )
REvsCV.naive_table_1.6 <-
  REvsCV.naive.table(
    I = I[3],
    t = t[3],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta
  )

REvsCV.naive_table_2.1 <-
  REvsCV.naive.table(
    I = I[1],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[2],
    alpha1 = c(0.001, alpha0[2] / 2, alpha0[2]),
    delta = delta
  )
REvsCV.naive_table_2.2 <-
  REvsCV.naive.table(
    I = I[2],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[2],
    alpha1 = c(0.001, alpha0[2] / 2, alpha0[2]),
    delta = delta
  )
REvsCV.naive_table_2.3 <-
  REvsCV.naive.table(
    I = I[3],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[2],
    alpha1 = c(0.001, alpha0[2] / 2, alpha0[2]),
    delta = delta
  )
REvsCV.naive_table_2.4 <-
  REvsCV.naive.table(
    I = I[4],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[2],
    alpha1 = c(0.001, alpha0[2] / 2, alpha0[2]),
    delta = delta
  )
REvsCV.naive_table_2.5 <-
  REvsCV.naive.table(
    I = I[3],
    t = t[1],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[2],
    alpha1 = c(0.001, alpha0[2] / 2, alpha0[2]),
    delta = delta
  )
REvsCV.naive_table_2.6 <-
  REvsCV.naive.table(
    I = I[3],
    t = t[3],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[2],
    alpha1 = c(0.001, alpha0[2] / 2, alpha0[2]),
    delta = delta
  )

REvsCV.naive_table_3.1 <-
  REvsCV.naive.table(
    I = I[1],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[3],
    alpha1 = c(0.001, alpha0[3] / 2, alpha0[3]),
    delta = delta
  )
REvsCV.naive_table_3.2 <-
  REvsCV.naive.table(
    I = I[2],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[3],
    alpha1 = c(0.001, alpha0[3] / 2, alpha0[3]),
    delta = delta
  )
REvsCV.naive_table_3.3 <-
  REvsCV.naive.table(
    I = I[3],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[3],
    alpha1 = c(0.001, alpha0[3] / 2, alpha0[3]),
    delta = delta
  )
REvsCV.naive_table_3.4 <-
  REvsCV.naive.table(
    I = I[4],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[3],
    alpha1 = c(0.001, alpha0[3] / 2, alpha0[3]),
    delta = delta
  )
REvsCV.naive_table_3.5 <-
  REvsCV.naive.table(
    I = I[3],
    t = t[1],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[3],
    alpha1 = c(0.001, alpha0[3] / 2, alpha0[3]),
    delta = delta
  )
REvsCV.naive_table_3.6 <-
  REvsCV.naive.table(
    I = I[3],
    t = t[3],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[3],
    alpha1 = c(0.001, alpha0[3] / 2, alpha0[3]),
    delta = delta
  )


write.csv(REvsCV.naive_table_1.5, file = paste0("./REvsCV.naive_1.5.csv"))
write.csv(REvsCV.naive_table_1.6, file = paste0("./REvsCV.naive_1.6.csv"))

write.csv(REvsCV.naive_table_2.1, file = paste0("./REvsCV.naive_2.1.csv"))
write.csv(REvsCV.naive_table_2.2, file = paste0("./REvsCV.naive_2.2.csv"))
write.csv(REvsCV.naive_table_2.3, file = paste0("./REvsCV.naive_2.3.csv"))
write.csv(REvsCV.naive_table_2.4, file = paste0("./REvsCV.naive_2.4.csv"))
write.csv(REvsCV.naive_table_2.5, file = paste0("./REvsCV.naive_2.5.csv"))
write.csv(REvsCV.naive_table_2.6, file = paste0("./REvsCV.naive_2.6.csv"))

write.csv(REvsCV.naive_table_3.1, file = paste0("./REvsCV.naive_3.1.csv"))
write.csv(REvsCV.naive_table_3.2, file = paste0("./REvsCV.naive_3.2.csv"))
write.csv(REvsCV.naive_table_3.3, file = paste0("./REvsCV.naive_3.3.csv"))
write.csv(REvsCV.naive_table_3.4, file = paste0("./REvsCV.naive_3.4.csv"))
write.csv(REvsCV.naive_table_3.5, file = paste0("./REvsCV.naive_3.5.csv"))
write.csv(REvsCV.naive_table_3.6, file = paste0("./REvsCV.naive_3.6.csv"))

###############################
## Within cluster variability
###############################

REvsCV_table_1.1.2 <-
  REvsCV.twoway.table(
    I = I[2],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta,
    case = 1,
    p1 = 0.1
  )
REvsCV_table_1.1.3 <-
  REvsCV.twoway.table(
    I = I[3],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta,
    case = 1,
    p1 = 0.1
  )

REvsCV_table_1.2.2 <-
  REvsCV.twoway.table(
    I = I[2],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta,
    case = 2,
    p1 = 0.1
  )
REvsCV_table_1.2.3 <-
  REvsCV.twoway.table(
    I = I[3],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta,
    case = 2,
    p1 = 0.1
  )

REvsCV_table_1.3.2 <-
  REvsCV.twoway.table(
    I = I[2],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta,
    case = 3,
    p1 = 0.1
  )
REvsCV_table_1.3.3 <-
  REvsCV.twoway.table(
    I = I[3],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta,
    case = 3,
    p1 = 0.1
  )


REvsCV_table_1.4.2 <-
  REvsCV.twoway.table(
    I = I[2],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta,
    case = 4,
    p1 = 0.1
  )
REvsCV_table_1.4.3 <-
  REvsCV.twoway.table(
    I = I[3],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta,
    case = 4,
    p1 = 0.1
  )

write.csv(REvsCV_table_1.1.2, file = paste0("./REvsCV_1.1.2.csv"))
write.csv(REvsCV_table_1.1.3, file = paste0("./REvsCV_1.1.3.csv"))

write.csv(REvsCV_table_1.2.2, file = paste0("./REvsCV_1.2.2.csv"))
write.csv(REvsCV_table_1.2.3, file = paste0("./REvsCV_1.2.3.csv"))

write.csv(REvsCV_table_1.3.2, file = paste0("./REvsCV_1.3.2.csv"))
write.csv(REvsCV_table_1.3.3, file = paste0("./REvsCV_1.3.3.csv"))

write.csv(REvsCV_table_1.4.2, file = paste0("./REvsCV_1.4.2.csv"))
write.csv(REvsCV_table_1.4.3, file = paste0("./REvsCV_1.4.3.csv"))

REvsCV.naive_table_1.1.2 <-
  REvsCV.twoway.naive.table(
    I = I[2],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta,
    case = 1,
    p1 = 0.1
  )
REvsCV.naive_table_1.1.3 <-
  REvsCV.twoway.naive.table(
    I = I[3],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta,
    case = 1,
    p1 = 0.1
  )

REvsCV.naive_table_1.2.2 <-
  REvsCV.twoway.naive.table(
    I = I[2],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta,
    case = 2,
    p1 = 0.1
  )
REvsCV.naive_table_1.2.3 <-
  REvsCV.twoway.naive.table(
    I = I[3],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta,
    case = 2,
    p1 = 0.1
  )

REvsCV.naive_table_1.3.2 <-
  REvsCV.twoway.naive.table(
    I = I[2],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta,
    case = 3,
    p1 = 0.1
  )
REvsCV.naive_table_1.3.3 <-
  REvsCV.twoway.naive.table(
    I = I[3],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta,
    case = 3,
    p1 = 0.1
  )

REvsCV.naive_table_1.4.2 <-
  REvsCV.twoway.naive.table(
    I = I[2],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta,
    case = 4,
    p1 = 0.1
  )
REvsCV.naive_table_1.4.3 <-
  REvsCV.twoway.naive.table(
    I = I[3],
    t = t[2],
    cv.num = cv,
    cmean = cmean[2],
    base = base,
    end = end,
    alpha0 = alpha0[1],
    alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]),
    delta = delta,
    case = 4,
    p1 = 0.1
  )

write.csv(REvsCV.naive_table_1.1.2, file = paste0("./REvsCV.naive_1.1.2.csv"))
write.csv(REvsCV.naive_table_1.1.3, file = paste0("./REvsCV.naive_1.1.3.csv"))

write.csv(REvsCV.naive_table_1.2.2, file = paste0("./REvsCV.naive_1.2.2.csv"))
write.csv(REvsCV.naive_table_1.2.3, file = paste0("./REvsCV.naive_1.2.3.csv"))

write.csv(REvsCV.naive_table_1.3.2, file = paste0("./REvsCV.naive_1.3.2.csv"))
write.csv(REvsCV.naive_table_1.3.3, file = paste0("./REvsCV.naive_1.3.3.csv"))

write.csv(REvsCV.naive_table_1.4.2, file = paste0("./REvsCV.naive_1.4.2.csv"))
write.csv(REvsCV.naive_table_1.4.3, file = paste0("./REvsCV.naive_1.4.3.csv"))
