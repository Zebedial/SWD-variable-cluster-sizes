#######################################################################
# Illustration of Impact of number of period (J) under NEX correlation 
# structure (Table 2 in the main article and Web Tables 1, 2, 3)
#######################################################################

setwd("C:\\Users\\ASUS\\Desktop\\YCAS\\Simulation\\Nest Exchangeable Correlation")
source(".\\NEX_functions.R")

#Install xtable if needed
#library(xtable)

#######################################################################
# Minimum, lower quartile (25 percentile), median, mean, upper quartile 
# (75 percentile), and maximum of RE were computed based on different 
# values of CV, different patterns of within-cluster variability, given 
# other required parameters, when the true correlation model is nested 
# exchangeable (NEX)
#
# Input:
# I: number of clusters (I)
# t: number of periods (J)
# cmean: mean cluster-period size
# alpha0: WP-ICC
# alpha1: BP-ICC 
# delta: treatment effect in scale of log(OR)
# base: baseline prevalence
# end: end prevalence, e.g. equal to baseline prevalence means constant 
#      prevalence
# type: 1 for the case of correctly specified working correlation model, 
#       2 for the working independence scenario
# initial: initial value of the probability vector (0.2 for J=3, 0.1 for 
#          J=5, 0.05 for J=13, where J is the number of periods)
# nsims: times for the simulation runs with default value 1000
# seed: a seed to control reproducible output with default value 8888
#
# Output: 
# a table of minimum lower quartile (25 percentile), median, mean, upper
# quartile (75 percentile), and maximum of RE as a function of different 
# CV and patterns of within-cluster imbalance
#########################################################################

CV_table <- function(I, t, cmean, alpha0, alpha1, delta, base, end, type, 
                     initial, nsims=1000, seed=8888) {
  CV <- c(0.25, 0.5, 0.75, 1, 1.25, 1.5)
  p <- seq(base, end, length.out = t)
  beta <- log(p / (1 - p))
  
  RE.0 <- matrix(NA, 6, nsims)
  RE.1 <- matrix(NA, 6, nsims)
  RE.2 <- matrix(NA, 6, nsims)
  RE.3 <- matrix(NA, 6, nsims)
  RE.4 <- matrix(NA, 6, nsims)
  
  if (type == 1) {
    for (i in 1:6) {
      for (j in 1:nsims) {
        set.seed(seed + j)
        RE.0[i, j] <-
          sumVAR(n.equal(I, t, cmean),
                 scheme(I, t),
                 delta,
                 beta,
                 c(alpha0, alpha1)) / sumVAR(n.variable(I, t, CV[i],
                                                        cmean),
                                             scheme(I, t),
                                             delta,
                                             beta,
                                             c(alpha0, alpha1))
        RE.1[i, j] <-
          sumVAR(n.equal(I, t, cmean),
                 scheme(I, t),
                 delta,
                 beta,
                 c(alpha0, alpha1)) / sumVAR(
                   n.variable.twoway(I, t, CV[i], cmean, case = 1, 
                                     p1 = initial),
                   scheme(I, t),
                   delta,
                   beta,
                   c(alpha0, alpha1)
                 )
        RE.2[i, j] <-
          sumVAR(n.equal(I, t, cmean),
                 scheme(I, t),
                 delta,
                 beta,
                 c(alpha0, alpha1)) / sumVAR(
                   n.variable.twoway(I, t, CV[i], cmean, case = 2, 
                                     p1 = initial),
                   scheme(I, t),
                   delta,
                   beta,
                   c(alpha0, alpha1)
                 )
        RE.3[i, j] <-
          sumVAR(n.equal(I, t, cmean),
                 scheme(I, t),
                 delta,
                 beta,
                 c(alpha0, alpha1)) / sumVAR(
                   n.variable.twoway(I, t, CV[i], cmean, case = 3, 
                                     p1 = initial),
                   scheme(I, t),
                   delta,
                   beta,
                   c(alpha0, alpha1)
                 )
        RE.4[i, j] <-
          sumVAR(n.equal(I, t, cmean),
                 scheme(I, t),
                 delta,
                 beta,
                 c(alpha0, alpha1)) / sumVAR(
                   n.variable.twoway(I, t, CV[i], cmean, case = 4,
                                     p1 = initial),
                   scheme(I, t),
                   delta,
                   beta,
                   c(alpha0, alpha1)
                 )
      }
    }
  }
  
  if (type == 2) {
    for (i in 1:6) {
      for (j in 1:nsims) {
        set.seed(seed + j)
        RE.0[i, j] <-
          sumVAR.naive(n.equal(I, t, cmean),
                       scheme(I, t),
                       delta,
                       beta,
                       c(alpha0, alpha1)) / sumVAR.naive(
                         n.variable(I, t, CV[i], cmean),
                         scheme(I, t),
                         delta,
                         beta,
                         c(alpha0, alpha1)
                       ) 
                                                           
        RE.1[i, j] <-
          sumVAR.naive(n.equal(I, t, cmean),
                       scheme(I, t),
                       delta,
                       beta,
                       c(alpha0, alpha1)) / sumVAR.naive(
                         n.variable.twoway(I, t, CV[i], cmean, case = 1, 
                                           p1 = initial),
                         scheme(I, t),
                         delta,
                         beta,
                         c(alpha0, alpha1)
                       )
        RE.2[i, j] <-
          sumVAR.naive(n.equal(I, t, cmean),
                       scheme(I, t),
                       delta,
                       beta,
                       c(alpha0, alpha1)) / sumVAR.naive(
                         n.variable.twoway(I, t, CV[i], cmean, case = 2, 
                                           p1 = initial),
                         scheme(I, t),
                         delta,
                         beta,
                         c(alpha0, alpha1)
                       )
        RE.3[i, j] <-
          sumVAR.naive(n.equal(I, t, cmean),
                       scheme(I, t),
                       delta,
                       beta,
                       c(alpha0, alpha1)) / sumVAR.naive(
                         n.variable.twoway(I, t, CV[i], cmean, case = 3, 
                                           p1 = initial),
                         scheme(I, t),
                         delta,
                         beta,
                         c(alpha0, alpha1)
                       )
        RE.4[i, j] <-
          sumVAR.naive(n.equal(I, t, cmean),
                       scheme(I, t),
                       delta,
                       beta,
                       c(alpha0, alpha1)) / sumVAR.naive(
                         n.variable.twoway(I, t, CV[i], cmean, case = 4, 
                                           p1 = initial),
                         scheme(I, t),
                         delta,
                         beta,
                         c(alpha0, alpha1)
                       )
      }
    }
  }
  
  RE.0_min <- apply(RE.0, 1, quantile)[1, ]
  RE.0_25 <- apply(RE.0, 1, quantile)[2, ]
  RE.0_median <- apply(RE.0, 1, quantile)[3, ]
  RE.0_75 <- apply(RE.0, 1, quantile)[4, ]
  RE.0_max <- apply(RE.0, 1, quantile)[5, ]
  RE.0_mean <- apply(RE.0, 1, mean)
  
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
  
  Pattern_compare <-
    data.frame(
      cbind(
        CV,
        RE.0_median,
        RE.0_25,
        RE.0_75,
        RE.0_mean,
        RE.0_min,
        RE.0_max,
        RE.1_median,
        RE.1_25,
        RE.1_75,
        RE.1_mean,
        RE.1_min,
        RE.1_max,
        RE.2_median,
        RE.2_25,
        RE.2_75,
        RE.2_mean,
        RE.2_min,
        RE.2_max,
        RE.3_median,
        RE.3_25,
        RE.3_75,
        RE.3_mean,
        RE.3_min,
        RE.3_max,
        RE.4_median,
        RE.4_25,
        RE.4_75,
        RE.4_mean,
        RE.4_min,
        RE.4_max
      )
    )
  return(Pattern_compare)
}

#Path for output tables
setwd("C:/Users/ASUS/Desktop/YCAS/result tables/NEX_tables")


####################################################################################
# Combine tables produced by the function above to obtain the tables like Table 2
# showed in the main article. Four tables are produced below, corresponding to the 
# scenarios with J=12, 24, 48, and 96
####################################################################################


#############################################
## J and CV table (Scenario 1) Web Table 1
#############################################

alpha0 <- 0.05
alpha1 <- 0.025
I <- 12

J3 <-
  CV_table(
    I = I,
    t = 3,
    cmean = 100,
    alpha0 = alpha0,
    alpha1 = alpha1,
    delta = log(0.35),
    base = 0.3,
    end = 0.3,
    type = 1,
    initial = 0.2
  )
J5 <-
  CV_table(
    I = I,
    t = 5,
    cmean = 100,
    alpha0 = alpha0,
    alpha1 = alpha1,
    delta = log(0.35),
    base = 0.3,
    end = 0.3,
    type = 1,
    initial = 0.1
  )
J13 <-
  CV_table(
    I = I,
    t = 13,
    cmean = 100,
    alpha0 = alpha0,
    alpha1 = alpha1,
    delta = log(0.35),
    base = 0.3,
    end = 0.3,
    type = 1,
    initial = 0.05
  )

J3.naive <-
  CV_table(
    I = I,
    t = 3,
    cmean = 100,
    alpha0 = alpha0,
    alpha1 = alpha1,
    delta = log(0.35),
    base = 0.3,
    end = 0.3,
    type = 2,
    initial = 0.2
  )
J5.naive <-
  CV_table(
    I = I,
    t = 5,
    cmean = 100,
    alpha0 = alpha0,
    alpha1 = alpha1,
    delta = log(0.35),
    base = 0.3,
    end = 0.3,
    type = 2,
    initial = 0.1
  )
J13.naive <-
  CV_table(
    I = I,
    t = 13,
    cmean = 100,
    alpha0 = alpha0,
    alpha1 = alpha1,
    delta = log(0.35),
    base = 0.3,
    end = 0.3,
    type = 2,
    initial = 0.05
  )

WCor <- rep(c("TRUE", "IND"), each = 18)
J <- rep(c(3, 5, 13, 3, 5, 13), each = 6)
CV.data <-
  round(rbind(J3, J5, J13, J3.naive, J5.naive, J13.naive), digits = 3)

CV.data.new <- matrix(0, ncol = 6, nrow = 36)
CV.data.new[, 1] <- CV.data[, 1]
colnames(CV.data.new) <-
  c("CV",
    "No Within-cluster imbalance",
    "Pattern 1",
    "Pattern 2",
    "Pattern 3",
    "Pattern 4")

for (i in 1:36) {
  for (j in 1:5) {
    CV.data.new[i, j + 1] <-
      paste0(CV.data[i, j * 6 - 4], " (", CV.data[i, j * 6 - 3], ", ", 
             CV.data[i, j * 6 - 2], ")")
                                                                                
  }
}

J_CV_table1 <- data.frame(cbind(WCor, J, CV.data.new))

write.csv(J_CV_table1, "./J_CV_table_12.csv")

#xtable(J_CV_table1[seq(1,36,by=2),-7])



#############################################
## J and CV table (Scenario 2) Table 2 in
## the main article
#############################################

alpha0 <- 0.05
alpha1 <- 0.025
I <- 24

J3 <-
  CV_table(
    I = I,
    t = 3,
    cmean = 100,
    alpha0 = alpha0,
    alpha1 = alpha1,
    delta = log(0.35),
    base = 0.3,
    end = 0.3,
    type = 1,
    initial = 0.2
  )
J5 <-
  CV_table(
    I = I,
    t = 5,
    cmean = 100,
    alpha0 = alpha0,
    alpha1 = alpha1,
    delta = log(0.35),
    base = 0.3,
    end = 0.3,
    type = 1,
    initial = 0.1
  )
J13 <-
  CV_table(
    I = I,
    t = 13,
    cmean = 100,
    alpha0 = alpha0,
    alpha1 = alpha1,
    delta = log(0.35),
    base = 0.3,
    end = 0.3,
    type = 1,
    initial = 0.05
  )

J3.naive <-
  CV_table(
    I = I,
    t = 3,
    cmean = 100,
    alpha0 = alpha0,
    alpha1 = alpha1,
    delta = log(0.35),
    base = 0.3,
    end = 0.3,
    type = 2,
    initial = 0.2
  )
J5.naive <-
  CV_table(
    I = I,
    t = 5,
    cmean = 100,
    alpha0 = alpha0,
    alpha1 = alpha1,
    delta = log(0.35),
    base = 0.3,
    end = 0.3,
    type = 2,
    initial = 0.1
  )
J13.naive <-
  CV_table(
    I = I,
    t = 13,
    cmean = 100,
    alpha0 = alpha0,
    alpha1 = alpha1,
    delta = log(0.35),
    base = 0.3,
    end = 0.3,
    type = 2,
    initial = 0.05
  )

WCor <- rep(c("TRUE", "IND"), each = 18)
J <- rep(c(3, 5, 13, 3, 5, 13), each = 6)
CV.data <-
  round(rbind(J3, J5, J13, J3.naive, J5.naive, J13.naive), digits = 3)

CV.data.new <- matrix(0, ncol = 6, nrow = 36)
CV.data.new[, 1] <- CV.data[, 1]
colnames(CV.data.new) <-
  c("CV",
    "No Within-cluster imbalance",
    "Pattern 1",
    "Pattern 2",
    "Pattern 3",
    "Pattern 4")

for (i in 1:36) {
  for (j in 1:5) {
    CV.data.new[i, j + 1] <-
      paste0(CV.data[i, j * 6 - 4], " (", CV.data[i, j * 6 - 3], ", ", 
             CV.data[i, j * 6 - 2], ")")
                                                                                 
  }
}

J_CV_table2 <- data.frame(cbind(WCor, J, CV.data.new))

write.csv(J_CV_table2, "./J_CV_table_24.csv")

#xtable(J_CV_table2[seq(1,36,by=2),-7])


#############################################
## J and CV table (Scenario 3) Web Table 2
#############################################

alpha0 <- 0.05
alpha1 <- 0.025
I <- 48

J3 <-
  CV_table(
    I = I,
    t = 3,
    cmean = 100,
    alpha0 = alpha0,
    alpha1 = alpha1,
    delta = log(0.35),
    base = 0.3,
    end = 0.3,
    type = 1,
    initial = 0.2
  )
J5 <-
  CV_table(
    I = I,
    t = 5,
    cmean = 100,
    alpha0 = alpha0,
    alpha1 = alpha1,
    delta = log(0.35),
    base = 0.3,
    end = 0.3,
    type = 1,
    initial = 0.1
  )
J13 <-
  CV_table(
    I = I,
    t = 13,
    cmean = 100,
    alpha0 = alpha0,
    alpha1 = alpha1,
    delta = log(0.35),
    base = 0.3,
    end = 0.3,
    type = 1,
    initial = 0.05
  )

J3.naive <-
  CV_table(
    I = I,
    t = 3,
    cmean = 100,
    alpha0 = alpha0,
    alpha1 = alpha1,
    delta = log(0.35),
    base = 0.3,
    end = 0.3,
    type = 2,
    initial = 0.2
  )
J5.naive <-
  CV_table(
    I = I,
    t = 5,
    cmean = 100,
    alpha0 = alpha0,
    alpha1 = alpha1,
    delta = log(0.35),
    base = 0.3,
    end = 0.3,
    type = 2,
    initial = 0.1
  )
J13.naive <-
  CV_table(
    I = I,
    t = 13,
    cmean = 100,
    alpha0 = alpha0,
    alpha1 = alpha1,
    delta = log(0.35),
    base = 0.3,
    end = 0.3,
    type = 2,
    initial = 0.05
  )

WCor <- rep(c("TRUE", "IND"), each = 18)
J <- rep(c(3, 5, 13, 3, 5, 13), each = 6)
CV.data <-
  round(rbind(J3, J5, J13, J3.naive, J5.naive, J13.naive), digits = 3)

CV.data.new <- matrix(0, ncol = 6, nrow = 36)
CV.data.new[, 1] <- CV.data[, 1]
colnames(CV.data.new) <-
  c("CV",
    "No Within-cluster imbalance",
    "Pattern 1",
    "Pattern 2",
    "Pattern 3",
    "Pattern 4")

for (i in 1:36) {
  for (j in 1:5) {
    CV.data.new[i, j + 1] <-
      paste0(CV.data[i, j * 6 - 4], " (", CV.data[i, j * 6 - 3], ", ",
             CV.data[i, j * 6 - 2], ")")
                                                                                
  }
}

J_CV_table3 <- data.frame(cbind(WCor, J, CV.data.new))

write.csv(J_CV_table3, "./J_CV_table_48.csv")

#xtable(J_CV_table3[seq(1,36,by=2),-7])


#############################################
## J and CV table (Scenario 4) Web Table 3
#############################################

alpha0 <- 0.05
alpha1 <- 0.025
I <- 96

J3 <-
  CV_table(
    I = I,
    t = 3,
    cmean = 100,
    alpha0 = alpha0,
    alpha1 = alpha1,
    delta = log(0.35),
    base = 0.3,
    end = 0.3,
    type = 1,
    initial = 0.2
  )
J5 <-
  CV_table(
    I = I,
    t = 5,
    cmean = 100,
    alpha0 = alpha0,
    alpha1 = alpha1,
    delta = log(0.35),
    base = 0.3,
    end = 0.3,
    type = 1,
    initial = 0.1
  )
J13 <-
  CV_table(
    I = I,
    t = 13,
    cmean = 100,
    alpha0 = alpha0,
    alpha1 = alpha1,
    delta = log(0.35),
    base = 0.3,
    end = 0.3,
    type = 1,
    initial = 0.05
  )

J3.naive <-
  CV_table(
    I = I,
    t = 3,
    cmean = 100,
    alpha0 = alpha0,
    alpha1 = alpha1,
    delta = log(0.35),
    base = 0.3,
    end = 0.3,
    type = 2,
    initial = 0.2
  )
J5.naive <-
  CV_table(
    I = I,
    t = 5,
    cmean = 100,
    alpha0 = alpha0,
    alpha1 = alpha1,
    delta = log(0.35),
    base = 0.3,
    end = 0.3,
    type = 2,
    initial = 0.1
  )
J13.naive <-
  CV_table(
    I = I,
    t = 13,
    cmean = 100,
    alpha0 = alpha0,
    alpha1 = alpha1,
    delta = log(0.35),
    base = 0.3,
    end = 0.3,
    type = 2,
    initial = 0.05
  )

WCor <- rep(c("TRUE", "IND"), each = 18)
J <- rep(c(3, 5, 13, 3, 5, 13), each = 6)
CV.data <-
  round(rbind(J3, J5, J13, J3.naive, J5.naive, J13.naive), digits = 3)

CV.data.new <- matrix(0, ncol = 6, nrow = 36)
CV.data.new[, 1] <- CV.data[, 1]
colnames(CV.data.new) <-
  c("CV",
    "No Within-cluster imbalance",
    "Pattern 1",
    "Pattern 2",
    "Pattern 3",
    "Pattern 4")

for (i in 1:36) {
  for (j in 1:5) {
    CV.data.new[i, j + 1] <-
      paste0(CV.data[i, j * 6 - 4], " (", CV.data[i, j * 6 - 3], ", ", 
             CV.data[i, j * 6 - 2], ")")
                                                                                 
  }
}

J_CV_table4 <- data.frame(cbind(WCor, J, CV.data.new))

write.csv(J_CV_table4, "./J_CV_table_96.csv")

#xtable(J_CV_table4[seq(1,36,by=2),-7])


