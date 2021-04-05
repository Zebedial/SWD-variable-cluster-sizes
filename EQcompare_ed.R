####################################################################################
# Obtain Original Data for Alternative Relative Efficiency (under true vs. 
# independence working correlaton structure) Defined in Section 8 When the True 
# Correlation Structure is Exponential Decay (ED)
####################################################################################


#Local path
#setwd("C:\\Users\\ASUS\\Desktop\\YCAS\\Simulation\\Exponential Decay Correlation")
source("./ED_functions.R")

#dir="C:\\Users\\ASUS\\Desktop\\YCAS\\result tables\\NEX_plots"
dir="C:\\Users\\ASUS\\Desktop\\YCAS\\result tables\\EXP_plots"

####################################################################################
# Tables used to produce RE vs. alpha1/alpha0 plots. The true correlation structure 
# is ED. Four values of WP-ICC alpha0 are used in each table (see section 4 in the 
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

eq_ed <- function(I, t, cmean, delta, base, end) {
  p <- seq(base, end, length.out = t)
  beta <- log(p / (1 - p))
  
  rho.num <- seq(0.025, 1, length.out = 40)
  RE.1 <- numeric(length(rho.num))
  RE.2 <- numeric(length(rho.num))
  RE.3 <- numeric(length(rho.num))
  RE.4 <- numeric(length(rho.num))
  
  alpha0 <- c(0.01, 0.05, 0.1, 0.2)
  for (i in 1:length(rho.num)) {
    rho <- rho.num[i]
    RE.1[i] <-
      sumVAR.exp(n.equal(I, t, cmean), scheme(I, t), delta, beta, alpha0[1], rho) /
      sumVAR.exp.naive(n.equal(I, t, cmean), scheme(I, t), delta, beta, 
                       alpha0[1], rho)
    RE.2[i] <-
      sumVAR.exp(n.equal(I, t, cmean), scheme(I, t), delta, beta, alpha0[2], rho) /
      sumVAR.exp.naive(n.equal(I, t, cmean), scheme(I, t), delta, beta, 
                       alpha0[2], rho)
    RE.3[i] <-
      sumVAR.exp(n.equal(I, t, cmean), scheme(I, t), delta, beta, alpha0[3], rho) /
      sumVAR.exp.naive(n.equal(I, t, cmean), scheme(I, t), delta, beta, 
                       alpha0[3], rho)
    RE.4[i] <-
      sumVAR.exp(n.equal(I, t, cmean), scheme(I, t), delta, beta, alpha0[4], rho) /
      sumVAR.exp.naive(n.equal(I, t, cmean), scheme(I, t), delta, beta,
                       alpha0[4], rho)
    
  }
  
  
  Rho <- rep(rho.num, 4)
  Alpha0 <- rep(alpha0, each = length(rho.num))
  
  RE <- c(RE.1, RE.2, RE.3, RE.4)
  REdata <- data.frame(cbind(Alpha0, Rho, RE))
  write.csv(REdata,
            file = paste0(dir, "/EQcompare_ed_", I, "_", t, "_", cmean, ".csv"))
}

#################################
## Produce data used to construct
## Web Figure 46
#################################

eq_ed(I=12, t=5, cmean=100, delta=log(0.35), base=0.3, end=0.3)

eq_ed(I=96, t=5, cmean=100, delta=log(0.35), base=0.3, end=0.3)

eq_ed(I=12, t=13, cmean=100, delta=log(0.35), base=0.3, end=0.3)

eq_ed(I=96, t=13, cmean=100, delta=log(0.35), base=0.3, end=0.3)

eq_ed(I=12, t=5, cmean=300, delta=log(0.35), base=0.3, end=0.3)

eq_ed(I=96, t=5, cmean=300, delta=log(0.35), base=0.3, end=0.3)







#################################
## Produce data used to conduct 
## experiments that are not showed 
## in the article
#################################

eq_ed(I=12, t=3, cmean=100, delta=log(0.35), base=0.3, end=0.3)

eq_ed(I=96, t=3, cmean=100, delta=log(0.35), base=0.3, end=0.3)

eq_ed(I=12, t=5, cmean=50, delta=log(0.35), base=0.3, end=0.3)

eq_ed(I=96, t=5, cmean=50, delta=log(0.35), base=0.3, end=0.3)