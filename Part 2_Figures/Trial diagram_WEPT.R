####################################################
# Produce the cluster-by-period diagram of the WEPT
# trial (Figure 1 in the main article)
####################################################

library(ggplot2)

setwd("C:\\Users\\ASUS\\Desktop\\YCAS\\Figures")
pdf("Trial_diagram_EPT.pdf", width=10, height=10, paper="special")


## Washington EPT cluster-period sizes
LHJ1 <- c(608, 787, 1075, 1553, 1370)
LHJ2 <- c(374, 482, 591, 616, 557)
LHJ3 <- c(26, 67, 115, 174, 94)
LHJ4 <- c(493, 590, 532, 785, 961)
LHJ5 <- c(64, 261, 176, 133, 90)
LHJ6 <- c(45, 212, 199, 79, 131)
LHJ7 <- c(359, 396, 329, 397, 354)
LHJ8 <- c(168, 272, 396, 476, 423)
LHJ9 <- c(221, 270, 276, 356, 225)
LHJ10 <- c(26, 44, 65, 83, 96)
LHJ11 <- c(48, 67, 73, 120, 166)
LHJ12 <- c(231, 318, 293, 305, 352)
LHJ13 <- c(333, 848, 1189, 1325, 1050)
LHJ14 <- c(135, 223, 250, 382, 855)
LHJ15 <- c(117, 132, 129, 224, 316)
LHJ16 <- c(17, 58, 79, 99, 59)
LHJ17 <- c(49, 44, 51, 83, 54)
LHJ18 <- c(36, 57, 62, 79, 131)
LHJ19 <- c(421, 568, 550, 844, 1050)
LHJ20 <- c(78, 58, 41, 52, 48)
LHJ21 <- c(19, 59, 109, 122, 101)
LHJ22 <- c(77, 78, 94, 122, 133)

WEPT <- rbind(LHJ1, LHJ2, LHJ3, LHJ4, LHJ5, LHJ6,
              LHJ7, LHJ8, LHJ9, LHJ10, LHJ11, LHJ12,
              LHJ13, LHJ14, LHJ15, LHJ16, LHJ17, LHJ18,
              LHJ19, LHJ20, LHJ21, LHJ22)

#Treatment indicator of each cluster period
WEPT.trt <- matrix(c(0, 1, 1, 1, 1,
                     0, 1, 1, 1, 1,
                     0, 1, 1, 1, 1,
                     0, 1, 1, 1, 1,
                     0, 1, 1, 1, 1,
                     0, 1, 1, 1, 1,
                     0, 0, 1, 1, 1,
                     0, 0, 1, 1, 1,
                     0, 0, 1, 1, 1,
                     0, 0, 1, 1, 1,
                     0, 0, 1, 1, 1,
                     0, 0, 1, 1, 1,
                     0, 0, 0, 1, 1,
                     0, 0, 0, 1, 1,
                     0, 0, 0, 1, 1,
                     0, 0, 0, 1, 1,
                     0, 0, 0, 1, 1,
                     0, 0, 0, 1, 1,
                     0, 0, 0, 0, 1,
                     0, 0, 0, 0, 1,
                     0, 0, 0, 0, 1,
                     0, 0, 0, 0, 1), nrow=22, ncol=5, byrow = T)

WEPT.data <- NULL
for (i in 1:22) {
  for (j in 1:5) {
    WEPT.data <- rbind(WEPT.data, c(i, j, WEPT[i, j],
                                    WEPT.trt[i, j]))
  }
}
colnames(WEPT.data) <- c("LHJ", "Period", "N", "Treatment")

WEPT.data <- data.frame(WEPT.data)

WEPT.data$Treatment <- factor(WEPT.data$Treatment, levels=c(0,1), 
                              labels = c("Control", "Intervention"))

## cluster by period diagram
breaks <- levels(unique(WEPT.data$Treatment))
p1 <- ggplot(WEPT.data, aes(Period, LHJ)) +
  geom_tile(aes(fill = Treatment)) +
  geom_text(aes(label = N), size = 5) +
  scale_fill_manual(
    values = c("#5cb8ff", "#52ffa3"),
    name = "Treatment",
    breaks = breaks,
    labels = breaks
  ) +
  xlab("Period") + ylab("Cluster (LHJ)") + 
  scale_y_continuous(limits = c(0.5, 22.5), breaks = seq(1, 22, 1)) +
  #labs(fill = as.expression(bquote(n[ij]))) +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(
      size = 14,
      angle = 0,
      vjust = 0.5
    ),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    axis.text = element_text(size = 15),
    legend.position = "right"
  )


p1

dev.off()