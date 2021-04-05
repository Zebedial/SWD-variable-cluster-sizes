#############################################
## Sensitivity to treatment effect, baseline 
## prevalence and prevalence trend under patten
## 4 within-cluster imbalance, ED true correlation
## (Web Tables 20, 21, 22, 23)
#############################################

dir="C:\\Users\\ASUS\\Desktop\\YCAS\\result tables\\EXP_tables"
#setwd("C:\\Users\\ASUS\\Desktop\\YCAS\\Simulation\\Exponential Decay Correlation")
source("./ED_functions.R")

####################################################################################
# Construct datasets used to produce the summary table of the impact of treatment effect
# baseline prevalence and prevalence trend, when patten 4 within-cluster imbalance is 
# introduced. The output is a subset of the final table.
# 
# Input:
# I: number of clusters (I)
# t: number of periods (J)
# cmean: mean cluster-period size
# alpha0: WP-ICC 
# rho: decay rate
# delta: treatment effect in scale of log(OR)
# type: working correlation structure. 1 for correctly specified working correlation 
#       structure, 2 for independence working correlation structure
# case: pattern of within-cluster. 1 for constant, 2 for increasing, 3 for decreasing
#       4 for randomly permuted
# seed: a seed to control reproducible output with default value 8888
# nsims: times for the simulation runs with default value 1000
#
# Output: a table with minimum, lower quartile (25 percentile), median, mean, upper
#         quartile (75 percentile), and maximum of RE as a function of CV, baseline
#         prevalence, and prevalence trend specified in section 4 of the main article
#         given a set of other required parameters
####################################################################################


prev.table <- function(I, t, cmean, alpha0, rho, delta, type, case, nsims = 1000, 
                       seed = 8888) {
  
  CV <- rep(c(0.25, 0.75, 1.25), each = 2)
  Base <- rep(c(0.1, 0.3), 3)
  
  RE.constant <- matrix(NA, 6, nsims)
  RE.increase <- matrix(NA, 6, nsims)
  RE.decrease <- matrix(NA, 6, nsims)
  
  for (i in 1:6) {
    base <- Base[i]
    end.constant <- Base[i]
    end.increase <- Base[i] * 1.2
    end.decrease <- Base[i] * 0.8
    
    p1 <- seq(base, end.constant, length.out = t)
    beta.constant <- log(p1 / (1 - p1))
    p2 <- seq(base, end.increase, length.out = t)
    beta.increase <- log(p2 / (1 - p2))
    p3 <- seq(base, end.decrease, length.out = t)
    beta.decrease <- log(p3 / (1 - p3))
    
    if (type == 1) {
      for (j in 1:nsims) {
        set.seed(seed + j)
        RE.constant[i, j] <-
          sumVAR.exp(n.equal(I, t, cmean),
                     scheme(I, t),
                     delta,
                     beta.constant,
                     alpha0,
                     rho) / sumVAR.exp(
                       n.variable.twoway(I, t, CV[i], cmean, case, 
                                         p1 = 0.1),
                       scheme(I, t),
                       delta,
                       beta.constant,
                       alpha0,
                       rho
                     )
        RE.increase[i, j] <-
          sumVAR.exp(n.equal(I, t, cmean),
                     scheme(I, t),
                     delta,
                     beta.increase,
                     alpha0,
                     rho) / sumVAR.exp(
                       n.variable.twoway(I, t, CV[i], cmean, case, 
                                         p1 = 0.1),
                       scheme(I, t),
                       delta,
                       beta.increase,
                       alpha0,
                       rho
                     )
        RE.decrease[i, j] <-
          sumVAR.exp(n.equal(I, t, cmean),
                     scheme(I, t),
                     delta,
                     beta.decrease,
                     alpha0,
                     rho) / sumVAR.exp(
                       n.variable.twoway(I, t, CV[i], cmean, case, 
                                         p1 = 0.1),
                       scheme(I, t),
                       delta,
                       beta.decrease,
                       alpha0,
                       rho
                     )
      }
    }
    
    else if (type == 2) {
      for (j in 1:nsims) {
        set.seed(seed + j)
        RE.constant[i, j] <-
          sumVAR.exp.naive(n.equal(I, t, cmean),
                           scheme(I, t),
                           delta,
                           beta.constant,
                           alpha0,
                           rho) / sumVAR.exp.naive(
                             n.variable.twoway(I, t, CV[i], cmean, case, 
                                               p1 = 0.1),
                             scheme(I, t),
                             delta,
                             beta.constant,
                             alpha0,
                             rho
                           )
        RE.increase[i, j] <-
          sumVAR.exp.naive(n.equal(I, t, cmean),
                           scheme(I, t),
                           delta,
                           beta.increase,
                           alpha0,
                           rho) / sumVAR.exp.naive(
                             n.variable.twoway(I, t, CV[i], cmean, case, 
                                               p1 = 0.1),
                             scheme(I, t),
                             delta,
                             beta.increase,
                             alpha0,
                             rho
                           )
        RE.decrease[i, j] <-
          sumVAR.exp.naive(n.equal(I, t, cmean),
                           scheme(I, t),
                           delta,
                           beta.decrease,
                           alpha0,
                           rho) / sumVAR.exp.naive(
                             n.variable.twoway(I, t, CV[i], cmean, case, 
                                               p1 = 0.1),
                             scheme(I, t),
                             delta,
                             beta.decrease,
                             alpha0,
                             rho
                           )
      }
    }
  }
  
  RE.constant_min <- apply(RE.constant, 1, quantile)[1, ]
  RE.constant_25 <- apply(RE.constant, 1, quantile)[2, ]
  RE.constant_median <- apply(RE.constant, 1, quantile)[3, ]
  RE.constant_75 <- apply(RE.constant, 1, quantile)[4, ]
  RE.constant_max <- apply(RE.constant, 1, quantile)[5, ]
  RE.constant_mean <- apply(RE.constant, 1, mean)
  
  RE.increase_min <- apply(RE.increase, 1, quantile)[1, ]
  RE.increase_25 <- apply(RE.increase, 1, quantile)[2, ]
  RE.increase_median <- apply(RE.increase, 1, quantile)[3, ]
  RE.increase_75 <- apply(RE.increase, 1, quantile)[4, ]
  RE.increase_max <- apply(RE.increase, 1, quantile)[5, ]
  RE.increase_mean <- apply(RE.increase, 1, mean)
  
  RE.decrease_min <- apply(RE.decrease, 1, quantile)[1, ]
  RE.decrease_25 <- apply(RE.decrease, 1, quantile)[2, ]
  RE.decrease_median <- apply(RE.decrease, 1, quantile)[3, ]
  RE.decrease_75 <- apply(RE.decrease, 1, quantile)[4, ]
  RE.decrease_max <- apply(RE.decrease, 1, quantile)[5, ]
  RE.decrease_mean <- apply(RE.decrease, 1, mean)
  
  Prevalence_compare <-
    data.frame(
      cbind(
        CV,
        Base,
        RE.constant_median,
        RE.constant_25,
        RE.constant_75,
        RE.increase_median,
        RE.increase_25,
        RE.increase_75,
        RE.decrease_median,
        RE.decrease_25,
        RE.decrease_75,
        RE.constant_mean,
        RE.constant_min,
        RE.constant_max,
        RE.increase_mean,
        RE.increase_min,
        RE.increase_max,
        RE.decrease_mean,
        RE.decrease_min,
        RE.decrease_max
      )
    )
  return(Prevalence_compare)
}

############################################################
## Scenario 1: I=24, pattern 4 within-cluster imbalance, rho=0.9
############################################################

treatment <- rep(c("log(0.35)", "log(0.75)"), each = 6)
prevtable <-
  rbind(
    prev.table(
      I = 24,
      t = 5,
      cmean = 100,
      alpha0 = 0.05,
      rho = 0.9,
      delta = log(0.35),
      type = 1,
      case = 4
    ),
    prev.table(
      I = 24,
      t = 5,
      cmean = 100,
      alpha0 = 0.05,
      rho = 0.9,
      delta = log(0.75),
      type = 1,
      case = 4
    )
  )
prevtable <- cbind(treatment, prevtable)

prev_out1 <- prevtable[, 1:3]
prev_out2 <- matrix(0, ncol = 3, nrow = 12)

for (i in 1:12) {
  for (j in 1:3) {
    prev_out2[i, j] <-
      paste0(round(prevtable[i, 3 * j + 1], 3),
             " (",
             round(prevtable[i, j * 3 + 2], 3),
             ", ",
             round(prevtable[i, j * 3 + 3], 3),
             ")")
  }
}
prev_out <- cbind(prev_out1, prev_out2)
colnames(prev_out) <-
  c("Treatment",
    "CV",
    "Baseline",
    "Constant",
    "Increaisng",
    "Decreasing")

write.csv(prev_out, file = paste0(dir, "/nuisance_ed_pattern4_1.csv"))

############################################################
## Scenario 2: I=24, pattern 4 within-cluster imbalance, rho=0.5
############################################################

treatment <- rep(c("log(0.35)", "log(0.75)"), each = 6)
prevtable <-
  rbind(
    prev.table(
      I = 24,
      t = 5,
      cmean = 100,
      alpha0 = 0.05,
      rho = 0.5,
      delta = log(0.35),
      type = 1,
      case = 4
    ),
    prev.table(
      I = 24,
      t = 5,
      cmean = 100,
      alpha0 = 0.05,
      rho = 0.5,
      delta = log(0.75),
      type = 1,
      case = 4
    )
  )
prevtable <- cbind(treatment, prevtable)

prev_out1 <- prevtable[, 1:3]
prev_out2 <- matrix(0, ncol = 3, nrow = 12)

for (i in 1:12) {
  for (j in 1:3) {
    prev_out2[i, j] <-
      paste0(round(prevtable[i, 3 * j + 1], 3),
             " (",
             round(prevtable[i, j * 3 + 2], 3),
             ", ",
             round(prevtable[i, j * 3 + 3], 3),
             ")")
  }
}
prev_out <- cbind(prev_out1, prev_out2)
colnames(prev_out) <-
  c("Treatment",
    "CV",
    "Baseline",
    "Constant",
    "Increaisng",
    "Decreasing")

write.csv(prev_out, file = paste0(dir, "/nuisance_ed_pattern4_2.csv"))

#############################################################
## Scenario 3: I=24, pattern 4 within-cluster imbalance, IND working 
##            correlation, rho=0.9
#############################################################

treatment <- rep(c("log(0.35)", "log(0.75)"), each = 6)
prevtable <-
  rbind(
    prev.table(
      I = 24,
      t = 5,
      cmean = 100,
      alpha0 = 0.05,
      rho = 0.9,
      delta = log(0.35),
      type = 2,
      case = 4
    ),
    prev.table(
      I = 24,
      t = 5,
      cmean = 100,
      alpha0 = 0.05,
      rho = 0.9,
      delta = log(0.75),
      type = 2,
      case = 4
    )
  )
prevtable <- cbind(treatment, prevtable)

prev_out1 <- prevtable[, 1:3]
prev_out2 <- matrix(0, ncol = 3, nrow = 12)

for (i in 1:12) {
  for (j in 1:3) {
    prev_out2[i, j] <-
      paste0(round(prevtable[i, 3 * j + 1], 3),
             " (",
             round(prevtable[i, j * 3 + 2], 3),
             ", ",
             round(prevtable[i, j * 3 + 3], 3),
             ")")
  }
}
prev_out <- cbind(prev_out1, prev_out2)
colnames(prev_out) <-
  c("Treatment",
    "CV",
    "Baseline",
    "Constant",
    "Increaisng",
    "Decreasing")

write.csv(prev_out, file = paste0(dir, "/nuisance_ed_pattern4_1_naive.csv"))

#############################################################
## Scenario 4: I=24, pattern 4 within-cluster imbalance, IND working 
##            correlation, rho=0.5
#############################################################

treatment <- rep(c("log(0.35)", "log(0.75)"), each = 6)
prevtable <-
  rbind(
    prev.table(
      I = 24,
      t = 5,
      cmean = 100,
      alpha0 = 0.05,
      rho = 0.5,
      delta = log(0.35),
      type = 2,
      case = 4
    ),
    prev.table(
      I = 24,
      t = 5,
      cmean = 100,
      alpha0 = 0.05,
      rho = 0.5,
      delta = log(0.75),
      type = 2,
      case = 4
    )
  )
prevtable <- cbind(treatment, prevtable)

prev_out1 <- prevtable[, 1:3]
prev_out2 <- matrix(0, ncol = 3, nrow = 12)

for (i in 1:12) {
  for (j in 1:3) {
    prev_out2[i, j] <-
      paste0(round(prevtable[i, 3 * j + 1], 3),
             " (",
             round(prevtable[i, j * 3 + 2], 3),
             ", ",
             round(prevtable[i, j * 3 + 3], 3),
             ")")
  }
}
prev_out <- cbind(prev_out1, prev_out2)
colnames(prev_out) <-
  c("Treatment",
    "CV",
    "Baseline",
    "Constant",
    "Increaisng",
    "Decreasing")

write.csv(prev_out, file = paste0(dir, "/nuisance_ed_pattern4_2_naive.csv"))


library(xtable)
dir = "C:\\Users\\ASUS\\Desktop\\YCAS\\result tables\\EXP_tables"

xtable(read.csv(paste0(dir, "/nuisance_ed_pattern4_1.csv"), header = T))
xtable(read.csv(paste0(dir, "/nuisance_ed_pattern4_2.csv"), header = T))
xtable(read.csv(paste0(
  dir, "/nuisance_ed_pattern4_1_naive.csv"
), header = T))
xtable(read.csv(paste0(
  dir, "/nuisance_ed_pattern4_2_naive.csv"
), header = T))
