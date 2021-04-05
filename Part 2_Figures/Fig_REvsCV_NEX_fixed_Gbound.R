#######################################
# Produce Web Figure 47
#######################################
library(ggplot2)
library(ggpubr)
dir="C:\\Users\\ASUS\\Desktop\\YCAS\\result tables\\NEX_plots"
#dir="C:\\Users\\ASUS\\Desktop\\YCAS\\result tables\\EXP_plots"


setwd("C:\\Users\\ASUS\\Desktop\\YCAS\\Figures")
pdf("REvsCV_NEX_fixed_1_Gbound.pdf", width=10, height=7, 
    paper="special")


#CV vs RE

REvsCV1 <- read.csv(paste0(dir, "/REvsCV_1.1.csv"), header = T)
REvsCV2 <- read.csv(paste0(dir, "/REvsCV_1.4.csv"), header = T)
REvsCV1.naive <-
  read.csv(paste0(dir, "/REvsCV.naive_1.1.csv"), header = T)
REvsCV2.naive <-
  read.csv(paste0(dir, "/REvsCV.naive_1.4.csv"), header = T)

#color for lines
colcode <- c("#009292", "#b66dff", "#db6d00")


#RE vs CV plot
RECV.plot1 <- function(table, I, alpha1) {
  table$ICC <- as.factor(table$ICC)
  my.labs <-
    list(bquote(alpha[1] ==  ~ .(alpha1[1])),
         bquote(alpha[1] ==  ~ .(alpha1[2])),
         bquote(alpha[1] ==  ~ .(alpha1[3])))
  
  ggplot(table, aes(x = cv.num, y = RE_median, group = ICC)) +
    geom_ribbon(aes(
      ymin = RE_25,
      ymax = RE_75,
      color = ICC
    ), alpha = 0.2)  +
    geom_line(aes(color = ICC)) + scale_y_continuous(limits = c(0.27, 1.03),
                                                     breaks = seq(0.3, 1, 0.1)) +
    labs(
      title = paste0("I = ", I, ", Working correlation = NEX"),
      y = "Relative Efficiency",
      x = "CV"
    ) +
    scale_colour_manual(values = colcode,
                        labels = my.labs,
                        name = "BP-ICC") +
    guides(shape = guide_legend(override.aes = list(size = 1)),
           color = guide_legend(override.aes = list(size = 0.5))) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      axis.text = element_text(size = 10),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 8),
      legend.key.size = unit(1.2, "lines"),
      legend.position = "none"
    ) + stat_function(
      fun = ~ 1 / (1 + (.x) ^ 2),
      linetype = 2,
      col = "#9B9B9B",
      size = 0.7
    )
}

#

#RE vs CV plot
RECV.plot2 <- function(table, I, alpha1) {
  table$ICC <- as.factor(table$ICC)
  my.labs <-
    list(bquote(alpha[1] ==  ~ .(alpha1[1])),
         bquote(alpha[1] ==  ~ .(alpha1[2])),
         bquote(alpha[1] ==  ~ .(alpha1[3])))
  
  ggplot(table, aes(x = cv.num, y = RE_median, group = ICC)) +
    geom_ribbon(aes(
      ymin = RE_25,
      ymax = RE_75,
      color = ICC
    ), alpha = 0.2) +
    geom_line(aes(color = ICC)) + 
    scale_y_continuous(limits = c(0.27, 1.03), breaks = seq(0.3, 1, 0.1)) +
                                                     
    labs(
      title = paste0("I = ", I, ", Working correlation = IND"),
      y = "Relative Efficiency",
      x = "CV"
    ) +
    scale_colour_manual(values = colcode,
                        labels = my.labs,
                        name = "BP-ICC") +
    guides(shape = guide_legend(override.aes = list(size = 1)),
           color = guide_legend(override.aes = list(size = 0.5))) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      axis.text = element_text(size = 10),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 8),
      legend.key.size = unit(1.2, "lines"),
      legend.position = "none"
    ) + stat_function(
      fun = ~ 1 / (1 + (.x) ^ 2),
      linetype = 2,
      col = "#9B9B9B",
      size = 0.7
    )
}



I <- c(12, 24, 48, 96)
t <- c(3, 5, 13)
cv <- seq(0, 1.5, by = 0.05)
cmean <- c(50, 100, 200, 500)
base = 0.3
end = 0.3
alpha0 <- c(0.05, 0.1, 0.2)
alpha1 <- c(0.001, alpha0 / 2, alpha0)
delta = log(0.35)

p1 <-
  RECV.plot1(REvsCV1,
             I = I[1],
             alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]))
p2 <-
  RECV.plot1(REvsCV2,
             I = I[4],
             alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]))

p3 <-
  RECV.plot2(REvsCV1.naive,
             I = I[1],
             alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]))
p4 <-
  RECV.plot2(REvsCV2.naive,
             I = I[4],
             alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]))

#Legend
RECV.null <- function(table) {
  table$ICC <- as.factor(table$ICC)
  my.labs <-
    list(bquote(alpha[1] ==  ~ .(alpha1[1])),
         bquote(alpha[1] ==  ~ .(alpha1[2])),
         bquote(alpha[1] ==  ~ .(alpha1[3])))
  
  ggplot(table, aes(x = cv.num, y = RE_median, group = ICC)) +
    geom_ribbon(aes(
      ymin = RE_25,
      ymax = RE_75,
      color = ICC
    ), alpha = 0.2) +
    geom_line(aes(color = ICC)) +
    theme_void() +
    lims(x = c(0, 0), y = c(0, 0)) +
    scale_colour_manual(values = colcode,
                        labels = my.labs,
                        name = "BP-ICC") +
    guides(shape = guide_legend(override.aes = list(size = 1)),
           color = guide_legend(override.aes = list(size = 0.5))) +
    theme(
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 10),
      legend.key.size = unit(1.2, "lines")
    )
  
}



p5 <- RECV.null(REvsCV1)

figure <-
  ggarrange(
    p1,
    p2,
    p5,
    p3,
    p4,
    p5,
    labels = c("(A)", "(B)", "", "(C)", "(D)", ""),
    font.label = list(size = 15),
    ncol = 3,
    nrow = 2,
    widths = c(3, 3, 0.9)
  )

figure

dev.off()

#######################################
# Produce Web Figure 48
#######################################

dir="C:\\Users\\ASUS\\Desktop\\YCAS\\result tables\\NEX_plots"
#dir="C:\\Users\\ASUS\\Desktop\\YCAS\\result tables\\EXP_plots"


setwd("C:\\Users\\ASUS\\Desktop\\YCAS\\Figures")
pdf("REvsCV_NEX_fixed_2_Gbound.pdf", width=10, height=7, 
    paper="special")

#CV vs RE

REvsCV1 <- read.csv(paste0(dir, "/REvsCV_1.2.csv"), header = T)
REvsCV2 <- read.csv(paste0(dir, "/REvsCV_1.3.csv"), header = T)
REvsCV1.naive <-
  read.csv(paste0(dir, "/REvsCV.naive_1.2.csv"), header = T)
REvsCV2.naive <-
  read.csv(paste0(dir, "/REvsCV.naive_1.3.csv"), header = T)


#RE vs CV plot
RECV.plot1 <- function(table, I, alpha1) {
  table$ICC <- as.factor(table$ICC)
  my.labs <-
    list(bquote(alpha[1] ==  ~ .(alpha1[1])),
         bquote(alpha[1] ==  ~ .(alpha1[2])),
         bquote(alpha[1] ==  ~ .(alpha1[3])))
  
  ggplot(table, aes(x = cv.num, y = RE_median, group = ICC)) +
    geom_ribbon(aes(
      ymin = RE_25,
      ymax = RE_75,
      color = ICC
    ), alpha = 0.2) +
    geom_line(aes(color = ICC)) + scale_y_continuous(limits = c(0.27, 1.03),
                                                     breaks = seq(0.3, 1, 0.1)) +
    labs(
      title = paste0("I = ", I, ", Working correlation = NEX"),
      y = "Relative Efficiency",
      x = "CV"
    ) +
    scale_colour_manual(values = colcode,
                        labels = my.labs,
                        name = "BP-ICC") +
    guides(shape = guide_legend(override.aes = list(size = 1)),
           color = guide_legend(override.aes = list(size = 0.5))) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      axis.text = element_text(size = 10),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 8),
      legend.key.size = unit(1.2, "lines"),
      legend.position = "none"
    ) + stat_function(
      fun = ~ 1 / (1 + (.x) ^ 2),
      linetype = 2,
      col = "#9B9B9B",
      size = 0.7
    )
}

#RE vs CV plot
RECV.plot2 <- function(table, I, alpha1) {
  table$ICC <- as.factor(table$ICC)
  my.labs <-
    list(bquote(alpha[1] ==  ~ .(alpha1[1])),
         bquote(alpha[1] ==  ~ .(alpha1[2])),
         bquote(alpha[1] ==  ~ .(alpha1[3])))
  
  ggplot(table, aes(x = cv.num, y = RE_median, group = ICC)) +
    geom_ribbon(aes(
      ymin = RE_25,
      ymax = RE_75,
      color = ICC
    ), alpha = 0.2) +
    geom_line(aes(color = ICC)) + 
    scale_y_continuous(limits = c(0.27, 1.03), breaks = seq(0.3, 1, 0.1)) +
                                                     
    labs(
      title = paste0("I = ", I, ", Working correlation = IND"),
      y = "Relative Efficiency",
      x = "CV"
    ) +
    scale_colour_manual(values = colcode,
                        labels = my.labs,
                        name = "BP-ICC") +
    guides(shape = guide_legend(override.aes = list(size = 1)),
           color = guide_legend(override.aes = list(size = 0.5))) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      axis.text = element_text(size = 10),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 8),
      legend.key.size = unit(1.2, "lines"),
      legend.position = "none"
    ) + stat_function(
      fun = ~ 1 / (1 + (.x) ^ 2),
      linetype = 2,
      col = "#9B9B9B",
      size = 0.7
    )
}



I <- c(12, 24, 48, 96)
t <- c(3, 5, 13)
cv <- seq(0, 1.5, by = 0.05)
cmean <- c(50, 100, 200, 500)
base = 0.3
end = 0.3
alpha0 <- c(0.05, 0.1, 0.2)
alpha1 <- c(0.001, alpha0 / 2, alpha0)
delta = log(0.35)

p1 <-
  RECV.plot1(REvsCV1,
             I = I[2],
             alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]))
p2 <-
  RECV.plot1(REvsCV2,
             I = I[3],
             alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]))

p3 <-
  RECV.plot2(REvsCV1.naive,
             I = I[2],
             alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]))
p4 <-
  RECV.plot2(REvsCV2.naive,
             I = I[3],
             alpha1 = c(0.001, alpha0[1] / 2, alpha0[1]))

#Legend
RECV.null <- function(table) {
  table$ICC <- as.factor(table$ICC)
  my.labs <-
    list(bquote(alpha[1] ==  ~ .(alpha1[1])),
         bquote(alpha[1] ==  ~ .(alpha1[2])),
         bquote(alpha[1] ==  ~ .(alpha1[3])))
  
  ggplot(table, aes(x = cv.num, y = RE_median, group = ICC)) +
    geom_ribbon(aes(
      ymin = RE_25,
      ymax = RE_75,
      color = ICC
    ), alpha = 0.2) +
    geom_line(aes(color = ICC)) +
    theme_void() +
    lims(x = c(0, 0), y = c(0, 0)) +
    scale_colour_manual(values = colcode,
                        labels = my.labs,
                        name = "BP-ICC") +
    guides(shape = guide_legend(override.aes = list(size = 1)),
           color = guide_legend(override.aes = list(size = 0.5))) +
    theme(
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 10),
      legend.key.size = unit(1.2, "lines")
    )
  
}



p5 <- RECV.null(REvsCV1)

figure <-
  ggarrange(
    p1,
    p2,
    p5,
    p3,
    p4,
    p5,
    labels = c("(A)", "(B)", "", "(C)", "(D)", ""),
    font.label = list(size = 15),
    ncol = 3,
    nrow = 2,
    widths = c(3, 3, 0.9)
  )

figure

dev.off()