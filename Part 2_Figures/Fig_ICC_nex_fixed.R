#######################################
# Produce Figure 4 in the main article
#######################################
library(ggplot2)
library(ggpubr)
#Path for storing data
dir="C:\\Users\\ASUS\\Desktop\\YCAS\\result tables\\NEX_plots"
#dir="C:\\Users\\ASUS\\Desktop\\YCAS\\result tables\\EXP_plots"


setwd("C:\\Users\\ASUS\\Desktop\\YCAS\\Figures")
pdf("ICC_nex_fixed_1.pdf", width=10, height=10.5, paper="special")

#Read in data

ICC_1 <- read.csv(paste0(dir, "/icc_table_nex_12_1.csv"), header = T)
ICC_2 <- read.csv(paste0(dir, "/icc_table_nex_96_1.csv"), header = T)
ICC_3 <- read.csv(paste0(dir, "/icc_table_nex_12_2.csv"), header = T)
ICC_4 <- read.csv(paste0(dir, "/icc_table_nex_96_2.csv"), header = T)
ICC_5 <- read.csv(paste0(dir, "/icc_table_nex_12_3.csv"), header = T)
ICC_6 <- read.csv(paste0(dir, "/icc_table_nex_96_3.csv"), header = T)

#Colors for legend
colcode <- c("#490092", "#b66dff", "#009292", "#db6d00")


ICC_1$ICC <- as.factor(ICC_1$Alpha0)
ICC_2$ICC <- as.factor(ICC_2$Alpha0)
ICC_3$ICC <- as.factor(ICC_3$Alpha0)
ICC_4$ICC <- as.factor(ICC_4$Alpha0)
ICC_5$ICC <- as.factor(ICC_5$Alpha0)
ICC_6$ICC <- as.factor(ICC_6$Alpha0)

#plot
iccplot <- function(table, I, CV, alpha0 = c(0.01, 0.05, 0.1, 0.2)) {
  my.labs <-
    list(
      bquote(alpha[0] ==  ~ .(alpha0[1])),
      bquote(alpha[0] ==  ~ .(alpha0[2])),
      bquote(alpha[0] ==  ~ .(alpha0[3])),
      bquote(alpha[0] ==  ~ .(alpha0[4]))
    )
  
  ggplot(table, aes(
    x = Ratio,
    y = RE_Median,
    group = ICC,
    color = ICC
  )) +
    geom_line(size = 1) + scale_color_manual(values = colcode, labels =
                                               my.labs) +
    scale_y_continuous(limits = c(0.27, 1.03),
                       breaks = seq(0.3, 1, 0.1)) +
    labs(
      title = paste0("I = ", I, ", CV =", CV, ", 
                     Working correlation = NEX"),
      y = "Median of Relative Efficiency",
      x = expression(alpha[1] / alpha[0])
    ) +
    guides(
      shape = guide_legend(override.aes = list(size = 1)),
      color = guide_legend(override.aes = list(size = 0.5)),
      fill = guide_legend(title = "ICC")
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      axis.text = element_text(size = 10),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 10),
      legend.key.size = unit(1, "lines"),
      legend.position = "none"
    )
}



I <- c(12, 24, 48, 96)
CV <- c(0.25, 0.75, 1.25)


p1 <- iccplot(ICC_1, I = I[1], CV = CV[1])
p2 <- iccplot(ICC_2, I = I[4], CV = CV[1])
p3 <- iccplot(ICC_3, I = I[1], CV = CV[2])
p4 <- iccplot(ICC_4, I = I[4], CV = CV[2])
p5 <- iccplot(ICC_5, I = I[1], CV = CV[3])
p6 <- iccplot(ICC_6, I = I[4], CV = CV[3])



#Legend
iccplot.null <- function(table, alpha0 = c(0.01, 0.05, 0.1, 0.2)) {
  my.labs <-
    list(
      bquote(alpha[0] ==  ~ .(alpha0[1])),
      bquote(alpha[0] ==  ~ .(alpha0[2])),
      bquote(alpha[0] ==  ~ .(alpha0[3])),
      bquote(alpha[0] ==  ~ .(alpha0[4]))
    )
  ggplot(table, aes(
    x = Ratio,
    y = RE_Median,
    group = ICC,
    color = ICC
  )) +
    geom_line(size = 1) + scale_color_manual(values = colcode,
                                             labels = my.labs,
                                             name = "WP-ICC") +
    theme_void() +
    lims(x = c(0, 0), y = c(0, 0)) +
    guides(shape = guide_legend(override.aes = list(size = 1)),
           color = guide_legend(override.aes = list(size = 0.5))) +
    theme(
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 10),
      legend.key.size = unit(1.2, "lines")
    )
}


p7 <- iccplot.null(ICC_1)

figure <-
  ggarrange(
    p1,
    p2,
    p7,
    p3,
    p4,
    p7,
    p5,
    p6,
    p7,
    labels = c("(A)", "(B)", "", "(C)", "(D)", "", "(E)", "(F)", ""),
    font.label = list(size = 15),
    ncol = 3,
    nrow = 3,
    widths = c(3, 3, 0.8)
  )
figure

dev.off()


#######################################
# Produce Figure 5 in the main article
#######################################

dir="C:\\Users\\ASUS\\Desktop\\YCAS\\result tables\\NEX_plots"
#dir="C:\\Users\\ASUS\\Desktop\\YCAS\\result tables\\EXP_plots"


setwd("C:\\Users\\ASUS\\Desktop\\YCAS\\Figures")
pdf("ICC_nex_fixed_1_naive.pdf", width=10, height=10.5, paper="special")

#Read in data

ICC_1 <-
  read.csv(paste0(dir, "/icc_table_nex_12_1_naive.csv"), header = T)
ICC_2 <-
  read.csv(paste0(dir, "/icc_table_nex_96_1_naive.csv"), header = T)
ICC_3 <-
  read.csv(paste0(dir, "/icc_table_nex_12_2_naive.csv"), header = T)
ICC_4 <-
  read.csv(paste0(dir, "/icc_table_nex_96_2_naive.csv"), header = T)
ICC_5 <-
  read.csv(paste0(dir, "/icc_table_nex_12_3_naive.csv"), header = T)
ICC_6 <-
  read.csv(paste0(dir, "/icc_table_nex_96_3_naive.csv"), header = T)



ICC_1$ICC <- as.factor(ICC_1$Alpha0)
ICC_2$ICC <- as.factor(ICC_2$Alpha0)
ICC_3$ICC <- as.factor(ICC_3$Alpha0)
ICC_4$ICC <- as.factor(ICC_4$Alpha0)
ICC_5$ICC <- as.factor(ICC_5$Alpha0)
ICC_6$ICC <- as.factor(ICC_6$Alpha0)

#plot
iccplot <- function(table, I, CV, alpha0 = c(0.01, 0.05, 0.1, 0.2)) {
  my.labs <-
    list(
      bquote(alpha[0] ==  ~ .(alpha0[1])),
      bquote(alpha[0] ==  ~ .(alpha0[2])),
      bquote(alpha[0] ==  ~ .(alpha0[3])),
      bquote(alpha[0] ==  ~ .(alpha0[4]))
    )
  
  ggplot(table, aes(
    x = Ratio,
    y = RE_Median,
    group = ICC,
    color = ICC
  )) +
    geom_line(size = 1) + scale_color_manual(values = colcode, labels =
                                               my.labs) +
    scale_y_continuous(limits = c(0.27, 1.03),
                       breaks = seq(0.3, 1, 0.1)) +
    labs(
      title = paste0("I = ", I, ", CV =", CV, ", 
                     Working correlation = IND"),
      y = "Median of Relative Efficiency",
      x = expression(alpha[1] / alpha[0])
    ) +
    guides(
      shape = guide_legend(override.aes = list(size = 1)),
      color = guide_legend(override.aes = list(size = 0.5)),
      fill = guide_legend(title = "ICC")
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      axis.text = element_text(size = 10),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 10),
      legend.key.size = unit(1, "lines"),
      legend.position = "none"
    )
}



I <- c(12, 24, 48, 96)
CV <- c(0.25, 0.75, 1.25)


p1 <- iccplot(ICC_1, I = I[1], CV = CV[1])
p2 <- iccplot(ICC_2, I = I[4], CV = CV[1])
p3 <- iccplot(ICC_3, I = I[1], CV = CV[2])
p4 <- iccplot(ICC_4, I = I[4], CV = CV[2])
p5 <- iccplot(ICC_5, I = I[1], CV = CV[3])
p6 <- iccplot(ICC_6, I = I[4], CV = CV[3])



#Legend
iccplot.null <- function(table, alpha0 = c(0.01, 0.05, 0.1, 0.2)) {
  my.labs <-
    list(
      bquote(alpha[0] ==  ~ .(alpha0[1])),
      bquote(alpha[0] ==  ~ .(alpha0[2])),
      bquote(alpha[0] ==  ~ .(alpha0[3])),
      bquote(alpha[0] ==  ~ .(alpha0[4]))
    )
  ggplot(table, aes(
    x = Ratio,
    y = RE_Median,
    group = ICC,
    color = ICC
  )) +
    geom_line(size = 1) + scale_color_manual(values = colcode,
                                             labels = my.labs,
                                             name = "WP-ICC") +
    theme_void() +
    lims(x = c(0, 0), y = c(0, 0)) +
    guides(shape = guide_legend(override.aes = list(size = 1)),
           color = guide_legend(override.aes = list(size = 0.5))) +
    theme(
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 10),
      legend.key.size = unit(1.2, "lines")
    )
}


p7 <- iccplot.null(ICC_1)

figure <-
  ggarrange(
    p1,
    p2,
    p7,
    p3,
    p4,
    p7,
    p5,
    p6,
    p7,
    labels = c("(A)", "(B)", "", "(C)", "(D)", "", "(E)", "(F)", ""),
    font.label = list(size = 15),
    ncol = 3,
    nrow = 3,
    widths = c(3, 3, 0.8)
  )
figure

dev.off()



