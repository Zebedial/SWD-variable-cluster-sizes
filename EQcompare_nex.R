####################################################################################
# Obtain Original Data for Alternative Relative Efficiency (under true vs. 
# independence working correlaton structure) Defined in Section 8 When the True 
# Correlation Structure is Nested Exchangeable (NEX)
####################################################################################

#Local path
#setwd("C:\\Users\\ASUS\\Desktop\\YCAS\\Simulation\\Nest Exchangeable Correlation")
source("./NEX_functions.R")

dir="C:\\Users\\ASUS\\Desktop\\YCAS\\result tables\\NEX_plots"
#dir="C:\\Users\\ASUS\\Desktop\\YCAS\\result tables\\EXP_plots"

####################################################################################
# Tables used to produce RE vs. alpha1/alpha0 plots. The true correlation structure 
# is NEX. Four values of WP-ICC alpha0 are used in each table (see section 4 in the 
# main article).
#
# Input:
# I: number of clusters (I)
# t: number of periods (J)
# cmean: mean cluster-period size
# delta: treatment effect in scale of log(OR)
# base: baseline prevalence
# end: end prevalence, e.g. equal to baseline prevalence means constant prevalence
#
# Output: a table with RE (defined in section 8) as a number of clusters, number of 
#         periods, and cluster-period size (constant) given a set of other required 
#         parameters
####################################################################################

eq_nex <- function(I, t, cmean, delta, base, end) {
  p <- seq(base, end, length.out = t)
  beta <- log(p / (1 - p))
  
  ratio <- seq(0.025, 1, length.out = 40)
  RE.1 <- numeric(length(ratio))
  RE.2 <- numeric(length(ratio))
  RE.3 <- numeric(length(ratio))
  RE.4 <- numeric(length(ratio))
  
  alpha0 <- c(0.01, 0.05, 0.1, 0.2)
  for (i in 1:length(ratio)) {
    alpha1 <- ratio[i] * alpha0
    RE.1[i] <-
      sumVAR(n.equal(I, t, cmean),
             scheme(I, t),
             delta,
             beta,
             c(alpha0[1], alpha1[1])) / sumVAR.naive(n.equal(I, t, cmean),
                                                     scheme(I, t),
                                                     delta,
                                                     beta,
                                                     c(alpha0[1], alpha1[1]))
    RE.2[i] <-
      sumVAR(n.equal(I, t, cmean),
             scheme(I, t),
             delta,
             beta,
             c(alpha0[2], alpha1[2])) / sumVAR.naive(n.equal(I, t, cmean),
                                                     scheme(I, t),
                                                     delta,
                                                     beta,
                                                     c(alpha0[2], alpha1[2]))
    RE.3[i] <-
      sumVAR(n.equal(I, t, cmean),
             scheme(I, t),
             delta,
             beta,
             c(alpha0[3], alpha1[3])) / sumVAR.naive(n.equal(I, t, cmean),
                                                     scheme(I, t),
                                                     delta,
                                                     beta,
                                                     c(alpha0[3], alpha1[3]))
    RE.4[i] <-
      sumVAR(n.equal(I, t, cmean),
             scheme(I, t),
             delta,
             beta,
             c(alpha0[4], alpha1[4])) / sumVAR.naive(n.equal(I, t, cmean),
                                                     scheme(I, t),
                                                     delta,
                                                     beta,
                                                     c(alpha0[4], alpha1[4]))
    
  }
  
  
  Ratio <- rep(ratio, 4)
  Alpha0 <- rep(alpha0, each = length(ratio))
  
  RE <- c(RE.1, RE.2, RE.3, RE.4)
  REdata <- data.frame(cbind(Alpha0, Ratio, RE))
  write.csv(REdata,
            file = paste0(dir, "/EQcompare_nex_", I, "_", t, "_", cmean, ".csv"))
}

#################################
## Produce data used to construct
## Web Figure 45
#################################


eq_nex(I=12, t=5, cmean=100, delta=log(0.35), base=0.3, end=0.3)

eq_nex(I=96, t=5, cmean=100, delta=log(0.35), base=0.3, end=0.3)

eq_nex(I=12, t=13, cmean=100, delta=log(0.35), base=0.3, end=0.3)

eq_nex(I=96, t=13, cmean=100, delta=log(0.35), base=0.3, end=0.3)

eq_nex(I=12, t=5, cmean=300, delta=log(0.35), base=0.3, end=0.3)

eq_nex(I=96, t=5, cmean=300, delta=log(0.35), base=0.3, end=0.3)








#################################
## Produce data used to conduct 
## experiments that are not showed 
## in the article
#################################

eq_nex(I=12, t=3, cmean=100, delta=log(0.35), base=0.3, end=0.3)

eq_nex(I=96, t=3, cmean=100, delta=log(0.35), base=0.3, end=0.3)

eq_nex(I=12, t=5, cmean=50, delta=log(0.35), base=0.3, end=0.3)

eq_nex(I=96, t=5, cmean=50, delta=log(0.35), base=0.3, end=0.3)