#######################################
# Produce Web Figure 46
#######################################
library(ggplot2)
library(ggpubr)
#dir="C:\\Users\\ASUS\\Desktop\\YCAS\\result tables\\NEX_plots"
dir="C:\\Users\\ASUS\\Desktop\\YCAS\\result tables\\EXP_plots"

#Colors for legend
colcode <- c("#490092", "#b66dff", "#009292", "#db6d00")


setwd("C:\\Users\\ASUS\\Desktop\\YCAS\\Figures")
pdf("EQcompare_ED.pdf", width=10, height=13.5, paper="special")

#Read in data

ICC_1 <- 
  read.csv(paste0(dir, "/EQcompare_ed_12_5_100.csv"), header = T)
ICC_2 <-
  read.csv(paste0(dir, "/EQcompare_ed_96_5_100.csv"), header = T)
ICC_3 <-
  read.csv(paste0(dir, "/EQcompare_ed_12_13_100.csv"), header = T)
ICC_4 <-
  read.csv(paste0(dir, "/EQcompare_ed_96_13_100.csv"), header = T)
ICC_5 <-
  read.csv(paste0(dir, "/EQcompare_ed_12_5_300.csv"), header = T)
ICC_6 <-
  read.csv(paste0(dir, "/EQcompare_ed_96_5_300.csv"), header = T)


ICC_1$ICC <- as.factor(ICC_1$Alpha0)
ICC_2$ICC <- as.factor(ICC_2$Alpha0)
ICC_3$ICC <- as.factor(ICC_3$Alpha0)
ICC_4$ICC <- as.factor(ICC_4$Alpha0)
ICC_5$ICC <- as.factor(ICC_5$Alpha0)
ICC_6$ICC <- as.factor(ICC_6$Alpha0)

#plot
iccplot <- function(table,
                    I,
                    t,
                    cmean,
                    alpha0 = c(0.01, 0.05, 0.1, 0.2)) {
  my.labs <-
    list(
      bquote(alpha[0] ==  ~ .(alpha0[1])),
      bquote(alpha[0] ==  ~ .(alpha0[2])),
      bquote(alpha[0] ==  ~ .(alpha0[3])),
      bquote(alpha[0] ==  ~ .(alpha0[4]))
    )
  
  ggplot(table, aes(
    x = Rho,
    y = RE,
    group = ICC,
    color = ICC
  )) +
    geom_line(size = 1) + scale_color_manual(values = colcode, labels =
                                               my.labs) +
    scale_y_continuous(limits = c(0, 1.03), breaks = seq(0, 1, 0.1)) +
    labs(
      title = bquote("I = " * .(I) * ", J =" * .(t) * ", " * 
                       bar(n) * " = " * .(cmean)),
      y = "Relative Efficiency",
      x = expression(rho)
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




p1 <- iccplot(ICC_1, I = 12, t = 5, cmean = 100)
p2 <- iccplot(ICC_2, I = 96, t = 5, cmean = 100)
p3 <- iccplot(ICC_3, I = 12, t = 13, cmean = 100)
p4 <- iccplot(ICC_4, I = 96, t = 13, cmean = 100)
p5 <- iccplot(ICC_5, I = 12, t = 5, cmean = 300)
p6 <- iccplot(ICC_6, I = 96, t = 5, cmean = 300)



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
    x = Rho,
    y = RE,
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


