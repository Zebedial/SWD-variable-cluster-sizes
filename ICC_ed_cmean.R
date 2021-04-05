#############################################
## Median RE under different alpha0, alpha1/alpha0 
## and mean cluster-period sizes {50, 300}
## when the true correlation structure is ED 
#############################################

mainDir = '/ysm-gpfs/home/zt223/ycas/cmean'

#Local path
#setwd("C:\\Users\\ASUS\\Desktop\\YCAS\\Simulation\\Exponential Decay Correlation")
source("./ED_functions.R")

################################################################################
# Tables used to produce median RE vs. rho plots. The true correlation structure 
# is ED, and no within-cluster imbalance is introduced. Four values of WP-ICC 
# alpha0 are used in each table (see section 4 in the main article)
#
# Input:
# I: number of clusters (I)
# t: number of periods (J)
# cv: coefficient of variation, the degree of between cluster imbalance
# cmean: mean cluster-period size
# delta: treatment effect in scale of log(OR)
# base: baseline prevalence
# end: end prevalence, e.g. equal to baseline prevalence means constant 
#      prevalence
# type: working correlation structure. 1 for correctly specified working 
#       correlation structure, 2 for independence working correlation structure
# seed: a seed to control reproducible output with default value 8888
# nsims: times for the simulation runs with default value 1000
#
# Output: a table with minimum, lower quartile (25 percentile), median, mean, 
#         upper quartile (75 percentile), and maximum of RE as a number of 
#         clusters and CV specified in section 4 of the main article given a set 
#         of other required parameters
################################################################################

icc_ed_fixed <- function(I, CV, t, cmean, delta, base, end, type, nsims = 1000, 
                         seed = 8888) {
  
  p <- seq(base, end, length.out = t)
  beta <- log(p / (1 - p))
  
  rho.num <- seq(0.025, 1, length.out = 40)
  RE.1 <- matrix(NA, length(rho.num), nsims)
  RE.2 <- matrix(NA, length(rho.num), nsims)
  RE.3 <- matrix(NA, length(rho.num), nsims)
  RE.4 <- matrix(NA, length(rho.num), nsims)
  
  alpha0 <- c(0.01, 0.05, 0.1, 0.2)
  if (type == 1) {
    for (i in 1:length(rho.num)) {
      rho <- rho.num[i]
      for (j in 1:nsims) {
        set.seed(seed + j)
        RE.1[i, j] <-
          sumVAR.exp(n.equal(I, t, cmean),
                     scheme(I, t),
                     delta,
                     beta,
                     alpha0[1],
                     rho) / sumVAR.exp(n.variable(I, t, CV, cmean),
                                       scheme(I, t),
                                       delta,
                                       beta,
                                       alpha0[1],
                                       rho)
        RE.2[i, j] <-
          sumVAR.exp(n.equal(I, t, cmean),
                     scheme(I, t),
                     delta,
                     beta,
                     alpha0[2],
                     rho) / sumVAR.exp(n.variable(I, t, CV, cmean),
                                       scheme(I, t),
                                       delta,
                                       beta,
                                       alpha0[2],
                                       rho)
        RE.3[i, j] <-
          sumVAR.exp(n.equal(I, t, cmean),
                     scheme(I, t),
                     delta,
                     beta,
                     alpha0[3],
                     rho) / sumVAR.exp(n.variable(I, t, CV, cmean),
                                       scheme(I, t),
                                       delta,
                                       beta,
                                       alpha0[3],
                                       rho)
        RE.4[i, j] <-
          sumVAR.exp(n.equal(I, t, cmean),
                     scheme(I, t),
                     delta,
                     beta,
                     alpha0[4],
                     rho) / sumVAR.exp(n.variable(I, t, CV, cmean),
                                       scheme(I, t),
                                       delta,
                                       beta,
                                       alpha0[4],
                                       rho)
      }
      print(i)
    }
  }
  if (type == 2) {
    for (i in 1:length(rho.num)) {
      rho <- rho.num[i]
      for (j in 1:nsims) {
        set.seed(seed + j)
        RE.1[i, j] <-
          sumVAR.exp.naive(n.equal(I, t, cmean),
                           scheme(I, t),
                           delta,
                           beta,
                           alpha0[1],
                           rho) / sumVAR.exp.naive(n.variable(I, t, CV, 
                                                              cmean),
                                                   scheme(I, t),
                                                   delta,
                                                   beta,
                                                   alpha0[1],
                                                   rho)
        RE.2[i, j] <-
          sumVAR.exp.naive(n.equal(I, t, cmean),
                           scheme(I, t),
                           delta,
                           beta,
                           alpha0[2],
                           rho) / sumVAR.exp.naive(n.variable(I, t, CV, 
                                                              cmean),
                                                   scheme(I, t),
                                                   delta,
                                                   beta,
                                                   alpha0[2],
                                                   rho)
        RE.3[i, j] <-
          sumVAR.exp.naive(n.equal(I, t, cmean),
                           scheme(I, t),
                           delta,
                           beta,
                           alpha0[3],
                           rho) / sumVAR.exp.naive(n.variable(I, t, CV, 
                                                              cmean),
                                                   scheme(I, t),
                                                   delta,
                                                   beta,
                                                   alpha0[3],
                                                   rho)
        RE.4[i, j] <-
          sumVAR.exp.naive(n.equal(I, t, cmean),
                           scheme(I, t),
                           delta,
                           beta,
                           alpha0[4],
                           rho) / sumVAR.exp.naive(n.variable(I, t, CV, 
                                                              cmean),
                                                   scheme(I, t),
                                                   delta,
                                                   beta,
                                                   alpha0[4],
                                                   rho)
      }
    }
  }
  RE.1_min <- apply(RE.1, 1, quantile)[1, ]
  RE.1_25 <- apply(RE.1, 1, quantile)[2, ]
  RE.1_median <- apply(RE.1, 1, quantile)[3, ]
  RE.1_75 <- apply(RE.1, 1, quantile)[4, ]
  RE.1_max <- apply(RE.1, 1, quantile)[5, ]
  RE.1_mean <- apply(RE.1, 1, mean)
  
  RE.2_min <- apply(RE.2, 1, quantile)[1, ]
  RE.2_25 <- apply(RE.2, 1, quantile)[2, ]
  RE.2_median <- apply(RE.2, 1, quantile)[3, ]
  RE.2_75 <- apply(RE.2, 1, quantile)[4, ]
  RE.2_max <- apply(RE.2, 1, quantile)[5, ]
  RE.2_mean <- apply(RE.2, 1, mean)
  
  RE.3_min <- apply(RE.3, 1, quantile)[1, ]
  RE.3_25 <- apply(RE.3, 1, quantile)[2, ]
  RE.3_median <- apply(RE.3, 1, quantile)[3, ]
  RE.3_75 <- apply(RE.3, 1, quantile)[4, ]
  RE.3_max <- apply(RE.3, 1, quantile)[5, ]
  RE.3_mean <- apply(RE.3, 1, mean)
  
  RE.4_min <- apply(RE.4, 1, quantile)[1, ]
  RE.4_25 <- apply(RE.4, 1, quantile)[2, ]
  RE.4_median <- apply(RE.4, 1, quantile)[3, ]
  RE.4_75 <- apply(RE.4, 1, quantile)[4, ]
  RE.4_max <- apply(RE.4, 1, quantile)[5, ]
  RE.4_mean <- apply(RE.4, 1, mean)
  
  Rho <- rep(rho.num, 4)
  Alpha0 <- rep(alpha0, each = length(rho.num))
  
  REdata <-
    rbind(
      cbind(RE.1_median, RE.1_25, RE.1_75, RE.1_mean, RE.1_min, RE.1_max),
      cbind(RE.2_median, RE.2_25, RE.2_75, RE.2_mean, RE.2_min, RE.2_max),
      cbind(RE.3_median, RE.3_25, RE.3_75, RE.3_mean, RE.3_min, RE.3_max),
      cbind(RE.4_median, RE.4_25, RE.4_75, RE.4_mean, RE.4_min, RE.4_max)
    )
  colnames(REdata) <-
    c("RE_Median", "RE_25", "RE_75", "RE_Mean", "RE_min", "RE_max")
  REdata <- data.frame(cbind(Alpha0, Rho, REdata))
  cvid <- ((CV - 0.25) / 0.5) + 1
  if (type == 1) {
    write.csv(REdata,
              file = paste0("./icc_table_ed_", I, "_", cmean, "_", cvid, ".csv"))
  }
  if (type == 2) {
    write.csv(REdata,
              file = paste0("./icc_table_ed_", I, "_", cmean, "_", cvid, 
                            "_naive.csv"))
  }
}

####################################################################################
# Tables used to produce median RE vs. rho plots. The true correlation structure is 
# ED, and patterns of within-cluster imbalance are introduced. Four values of WP-ICC 
# alpha0 are used in each table (see section 4 in the main article)
#
# Input:
# I: number of clusters (I)
# t: number of periods (J)
# cv: coefficient of variation, the degree of between cluster imbalance
# cmean: mean cluster-period size
# delta: treatment effect in scale of log(OR)
# base: baseline prevalence
# end: end prevalence, e.g. equal to baseline prevalence means constant prevalence
# type: working correlation structure. 1 for correctly specified working correlation 
#       structure, 2 for independence working correlation structure
# initial: p1 in the probability vector for the multinomial sampling
# case: pattern of within-cluster. 1 for constant, 2 for increasing, 3 for decreasing
#       4 for randomly permuted
# seed: a seed to control reproducible output with default value 8888
# nsims: times for the simulation runs with default value 1000
#
# Output: a table with minimum, lower quartile (25 percentile), median, mean, upper
#         quartile (75 percentile), and maximum of RE as a number of clusters and CV
#         specified in section 4 of the main article given a set of other required 
#         parameters
####################################################################################

icc_ed_var <- function(I, CV, t, cmean, delta, base, end, type, initial, case, 
                       nsims = 1000, seed = 8888) {
  
  p <- seq(base, end, length.out = t)
  beta <- log(p / (1 - p))
  
  rho.num <- seq(0.025, 1, length.out = 40)
  RE.1 <- matrix(NA, length(rho.num), nsims)
  RE.2 <- matrix(NA, length(rho.num), nsims)
  RE.3 <- matrix(NA, length(rho.num), nsims)
  RE.4 <- matrix(NA, length(rho.num), nsims)
  
  alpha0 <- c(0.01, 0.05, 0.1, 0.2)
  if (type == 1) {
    for (i in 1:length(rho.num)) {
      rho <- rho.num[i]
      for (j in 1:nsims) {
        set.seed(seed + j)
        RE.1[i, j] <-
          sumVAR.exp(n.equal(I, t, cmean),
                     scheme(I, t),
                     delta,
                     beta,
                     alpha0[1],
                     rho) / sumVAR.exp(
                       n.variable.twoway(I, t, CV, cmean, case = case, p1 = initial),
                       scheme(I, t),
                       delta,
                       beta,
                       alpha0[1],
                       rho
                     )
        RE.2[i, j] <-
          sumVAR.exp(n.equal(I, t, cmean),
                     scheme(I, t),
                     delta,
                     beta,
                     alpha0[2],
                     rho) / sumVAR.exp(
                       n.variable.twoway(I, t, CV, cmean, case = case, p1 = initial),
                       scheme(I, t),
                       delta,
                       beta,
                       alpha0[2],
                       rho
                     )
        RE.3[i, j] <-
          sumVAR.exp(n.equal(I, t, cmean),
                     scheme(I, t),
                     delta,
                     beta,
                     alpha0[3],
                     rho) / sumVAR.exp(
                       n.variable.twoway(I, t, CV, cmean, case = case, p1 = initial),
                       scheme(I, t),
                       delta,
                       beta,
                       alpha0[3],
                       rho
                     )
        RE.4[i, j] <-
          sumVAR.exp(n.equal(I, t, cmean),
                     scheme(I, t),
                     delta,
                     beta,
                     alpha0[4],
                     rho) / sumVAR.exp(
                       n.variable.twoway(I, t, CV, cmean, case = case, p1 = initial),
                       scheme(I, t),
                       delta,
                       beta,
                       alpha0[4],
                       rho
                     )
      }
      print(i)
    }
  }
  if (type == 2) {
    for (i in 1:length(rho.num)) {
      rho <- rho.num[i]
      for (j in 1:nsims) {
        set.seed(seed + j)
        RE.1[i, j] <-
          sumVAR.exp.naive(n.equal(I, t, cmean),
                           scheme(I, t),
                           delta,
                           beta,
                           alpha0[1],
                           rho) / sumVAR.exp.naive(
                             n.variable.twoway(I, t, CV, cmean, case = case, 
                                               p1 = initial),
                             scheme(I, t),
                             delta,
                             beta,
                             alpha0[1],
                             rho
                           )
        RE.2[i, j] <-
          sumVAR.exp.naive(n.equal(I, t, cmean),
                           scheme(I, t),
                           delta,
                           beta,
                           alpha0[2],
                           rho) / sumVAR.exp.naive(
                             n.variable.twoway(I, t, CV, cmean, case = case, 
                                               p1 = initial),
                             scheme(I, t),
                             delta,
                             beta,
                             alpha0[2],
                             rho
                           )
        RE.3[i, j] <-
          sumVAR.exp.naive(n.equal(I, t, cmean),
                           scheme(I, t),
                           delta,
                           beta,
                           alpha0[3],
                           rho) / sumVAR.exp.naive(
                             n.variable.twoway(I, t, CV, cmean, case = case, 
                                               p1 = initial),
                             scheme(I, t),
                             delta,
                             beta,
                             alpha0[3],
                             rho
                           )
        RE.4[i, j] <-
          sumVAR.exp.naive(n.equal(I, t, cmean),
                           scheme(I, t),
                           delta,
                           beta,
                           alpha0[4],
                           rho) / sumVAR.exp.naive(
                             n.variable.twoway(I, t, CV, cmean, case = case, 
                                               p1 = initial),
                             scheme(I, t),
                             delta,
                             beta,
                             alpha0[4],
                             rho
                           )
      }
    }
  }
  RE.1_min <- apply(RE.1, 1, quantile)[1, ]
  RE.1_25 <- apply(RE.1, 1, quantile)[2, ]
  RE.1_median <- apply(RE.1, 1, quantile)[3, ]
  RE.1_75 <- apply(RE.1, 1, quantile)[4, ]
  RE.1_max <- apply(RE.1, 1, quantile)[5, ]
  RE.1_mean <- apply(RE.1, 1, mean)
  
  RE.2_min <- apply(RE.2, 1, quantile)[1, ]
  RE.2_25 <- apply(RE.2, 1, quantile)[2, ]
  RE.2_median <- apply(RE.2, 1, quantile)[3, ]
  RE.2_75 <- apply(RE.2, 1, quantile)[4, ]
  RE.2_max <- apply(RE.2, 1, quantile)[5, ]
  RE.2_mean <- apply(RE.2, 1, mean)
  
  RE.3_min <- apply(RE.3, 1, quantile)[1, ]
  RE.3_25 <- apply(RE.3, 1, quantile)[2, ]
  RE.3_median <- apply(RE.3, 1, quantile)[3, ]
  RE.3_75 <- apply(RE.3, 1, quantile)[4, ]
  RE.3_max <- apply(RE.3, 1, quantile)[5, ]
  RE.3_mean <- apply(RE.3, 1, mean)
  
  RE.4_min <- apply(RE.4, 1, quantile)[1, ]
  RE.4_25 <- apply(RE.4, 1, quantile)[2, ]
  RE.4_median <- apply(RE.4, 1, quantile)[3, ]
  RE.4_75 <- apply(RE.4, 1, quantile)[4, ]
  RE.4_max <- apply(RE.4, 1, quantile)[5, ]
  RE.4_mean <- apply(RE.4, 1, mean)
  
  Rho <- rep(rho.num, 4)
  Alpha0 <- rep(alpha0, each = length(rho.num))
  
  REdata <-
    rbind(
      cbind(RE.1_median, RE.1_25, RE.1_75, RE.1_mean, RE.1_min, RE.1_max),
      cbind(RE.2_median, RE.2_25, RE.2_75, RE.2_mean, RE.2_min, RE.2_max),
      cbind(RE.3_median, RE.3_25, RE.3_75, RE.3_mean, RE.3_min, RE.3_max),
      cbind(RE.4_median, RE.4_25, RE.4_75, RE.4_mean, RE.4_min, RE.4_max)
    )
  colnames(REdata) <-
    c("RE_Median", "RE_25", "RE_75", "RE_Mean", "RE_min", "RE_max")
  REdata <- data.frame(cbind(Alpha0, Rho, REdata))
  cvid <- ((CV - 0.25) / 0.5) + 1
  if (type == 1) {
    write.csv(REdata,
              file = paste0(
                "./icc_table_ed_pattern",
                case,
                "_",
                I,
                "_",
                cmean,
                "_",
                cvid,
                ".csv"
              ))
  }
  if (type == 2) {
    write.csv(
      REdata,
      file = paste0(
        "./icc_table_ed_pattern",
        case,
        "_",
        I,
        "_",
        cmean,
        "_",
        cvid,
        "_naive.csv"
      )
    )
  }
}

################################
## Data used to produce Web table 37
## I=12, true working correlation
################################

icc_ed_fixed(I=12, CV=0.25, t=5, cmean=50, delta=log(0.35), base=0.3, 
             end=0.3, type=1)


icc_ed_fixed(I=12, CV=0.75, t=5, cmean=50, delta=log(0.35), base=0.3, 
             end=0.3, type=1)


icc_ed_fixed(I=12, CV=1.25, t=5, cmean=50, delta=log(0.35), base=0.3, 
             end=0.3, type=1)



icc_ed_fixed(I=12, CV=0.25, t=5, cmean=300, delta=log(0.35), base=0.3, 
             end=0.3, type=1)


icc_ed_fixed(I=12, CV=0.75, t=5, cmean=300, delta=log(0.35), base=0.3, 
             end=0.3, type=1)


icc_ed_fixed(I=12, CV=1.25, t=5, cmean=300, delta=log(0.35), base=0.3, 
             end=0.3, type=1)


################################
## Data used to produce Web table 38
## I=96, true working correlation
################################

icc_ed_fixed(I=96, CV=0.25, t=5, cmean=50, delta=log(0.35), base=0.3, 
             end=0.3, type=1)


icc_ed_fixed(I=96, CV=0.75, t=5, cmean=50, delta=log(0.35), base=0.3, 
             end=0.3, type=1)


icc_ed_fixed(I=96, CV=1.25, t=5, cmean=50, delta=log(0.35), base=0.3, 
             end=0.3, type=1)



icc_ed_fixed(I=96, CV=0.25, t=5, cmean=300, delta=log(0.35), base=0.3, 
             end=0.3, type=1)


icc_ed_fixed(I=96, CV=0.75, t=5, cmean=300, delta=log(0.35), base=0.3, 
             end=0.3, type=1)

icc_ed_fixed(I=96, CV=1.25, t=5, cmean=300, delta=log(0.35), base=0.3, 
             end=0.3, type=1)


################################
## Data used to produce Web table 41
## I=12, IND working correlation
################################

icc_ed_fixed(I=12, CV=0.25, t=5, cmean=50, delta=log(0.35), base=0.3, 
             end=0.3, type=2)


icc_ed_fixed(I=12, CV=0.75, t=5, cmean=50, delta=log(0.35), base=0.3, 
             end=0.3, type=2)


icc_ed_fixed(I=12, CV=1.25, t=5, cmean=50, delta=log(0.35), base=0.3, 
             end=0.3, type=2)



icc_ed_fixed(I=12, CV=0.25, t=5, cmean=300, delta=log(0.35), base=0.3, 
             end=0.3, type=2)


icc_ed_fixed(I=12, CV=0.75, t=5, cmean=300, delta=log(0.35), base=0.3, 
             end=0.3, type=2)


icc_ed_fixed(I=12, CV=1.25, t=5, cmean=300, delta=log(0.35), base=0.3, 
             end=0.3, type=2)


################################
## Data used to produce Web table 42
## I=96, IND working correlation
################################

icc_ed_fixed(I=96, CV=0.25, t=5, cmean=50, delta=log(0.35), base=0.3, 
             end=0.3, type=2)


icc_ed_fixed(I=96, CV=0.75, t=5, cmean=50, delta=log(0.35), base=0.3, 
             end=0.3, type=2)


icc_ed_fixed(I=96, CV=1.25, t=5, cmean=50, delta=log(0.35), base=0.3, 
             end=0.3, type=2)



icc_ed_fixed(I=96, CV=0.25, t=5, cmean=300, delta=log(0.35), base=0.3, 
             end=0.3, type=2)


icc_ed_fixed(I=96, CV=0.75, t=5, cmean=300, delta=log(0.35), base=0.3, 
             end=0.3, type=2)


icc_ed_fixed(I=96, CV=1.25, t=5, cmean=300, delta=log(0.35), base=0.3, 
             end=0.3, type=2)


################################
## Data used to produce Web table 39
## I=12, true working correlation
################################

icc_ed_var(I=12, CV=0.25, t=5, cmean=50, delta=log(0.35), base=0.3, end=0.3, 
           type=1, initial=0.1, case=4)


icc_ed_var(I=12, CV=0.75, t=5, cmean=50, delta=log(0.35), base=0.3, end=0.3, 
           type=1, initial=0.1, case=4)


icc_ed_var(I=12, CV=1.25, t=5, cmean=50, delta=log(0.35), base=0.3, end=0.3, 
           type=1, initial=0.1, case=4)



icc_ed_var(I=12, CV=0.25, t=5, cmean=300, delta=log(0.35), base=0.3, end=0.3,
           type=1, initial=0.1, case=4)


icc_ed_var(I=12, CV=0.75, t=5, cmean=300, delta=log(0.35), base=0.3, end=0.3, 
           type=1, initial=0.1, case=4)


icc_ed_var(I=12, CV=1.25, t=5, cmean=300, delta=log(0.35), base=0.3, end=0.3,
           type=1, initial=0.1, case=4)


################################
## Data used to produce Web table 40
## I=96, true working correlation
################################

icc_ed_var(I=96, CV=0.25, t=5, cmean=50, delta=log(0.35), base=0.3, end=0.3,
           type=1, initial=0.1, case=4)


icc_ed_var(I=96, CV=0.75, t=5, cmean=50, delta=log(0.35), base=0.3, end=0.3,
           type=1, initial=0.1, case=4)


icc_ed_var(I=96, CV=1.25, t=5, cmean=50, delta=log(0.35), base=0.3, end=0.3, 
           type=1, initial=0.1, case=4)



icc_ed_var(I=96, CV=0.25, t=5, cmean=300, delta=log(0.35), base=0.3, end=0.3, 
           type=1, initial=0.1, case=4)


icc_ed_var(I=96, CV=0.75, t=5, cmean=300, delta=log(0.35), base=0.3, end=0.3, 
           type=1, initial=0.1, case=4)

icc_ed_var(I=96, CV=1.25, t=5, cmean=300, delta=log(0.35), base=0.3, end=0.3, 
           type=1, initial=0.1, case=4)


################################
## Data used to produce Web table 43
## I=12, IND working correlation
################################

icc_ed_var(I=12, CV=0.25, t=5, cmean=50, delta=log(0.35), base=0.3, end=0.3, 
           type=2, initial=0.1, case=4)


icc_ed_var(I=12, CV=0.75, t=5, cmean=50, delta=log(0.35), base=0.3, end=0.3, 
           type=2, initial=0.1, case=4)


icc_ed_var(I=12, CV=1.25, t=5, cmean=50, delta=log(0.35), base=0.3, end=0.3, 
           type=2, initial=0.1, case=4)



icc_ed_var(I=12, CV=0.25, t=5, cmean=300, delta=log(0.35), base=0.3, end=0.3, 
           type=2, initial=0.1, case=4)


icc_ed_var(I=12, CV=0.75, t=5, cmean=300, delta=log(0.35), base=0.3, end=0.3, 
           type=2, initial=0.1, case=4)


icc_ed_var(I=12, CV=1.25, t=5, cmean=300, delta=log(0.35), base=0.3, end=0.3, 
           type=2, initial=0.1, case=4)


################################
## Data used to produce Web table 44
## I=96, IND working correlation
################################

icc_ed_var(I=96, CV=0.25, t=5, cmean=50, delta=log(0.35), base=0.3, end=0.3, 
           type=2, initial=0.1, case=4)


icc_ed_var(I=96, CV=0.75, t=5, cmean=50, delta=log(0.35), base=0.3, end=0.3, 
           type=2, initial=0.1, case=4)


icc_ed_var(I=96, CV=1.25, t=5, cmean=50, delta=log(0.35), base=0.3, end=0.3, 
           type=2, initial=0.1, case=4)



icc_ed_var(I=96, CV=0.25, t=5, cmean=300, delta=log(0.35), base=0.3, end=0.3,
           type=2, initial=0.1, case=4)


icc_ed_var(I=96, CV=0.75, t=5, cmean=300, delta=log(0.35), base=0.3, end=0.3, 
           type=2, initial=0.1, case=4)


icc_ed_var(I=96, CV=1.25, t=5, cmean=300, delta=log(0.35), base=0.3, end=0.3, 
           type=2, initial=0.1, case=4)

