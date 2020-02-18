if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")
library(dplyr)

temp <- filter(temp, temp$gene != "GAD1")

aabiosynthesis_handcur <- temp

x <- filter(x, x$adj.P.Val < 0.05)
y <- x[ which(x$Yorf1 %in% temp$name), ]
z <- x[ which(x$Yorf1 %in% translationcontrol_handcur$name), ]
a <- x[ which(x$Yorf1 %in% tRNAaasynthetase_handcur$name), ]
b <- x[ which(x$Yorf1 %in% polIIIsubunit_handcur$name), ]
c <- x[ which(x$Yorf1 %in% tRNAprocessing_handcur$name), ]

#cummulative distribution

par(pty='s')  # force the plot to be square before we start

plot(ecdf(-x$logFC),
     xlim=c(-2,2),
     xlab="Log 2 Fold Change",
     ylab="Cumulative Proportion",
     main="Cummulative Distribution of GO subsets",
     lwd = 2,
     col="gray52")
lines(ecdf(-y$logFC),
      do.points = FALSE,
      verticals = TRUE, 
      lwd = 2,
      col="springgreen")
lines(ecdf(-z$logFC),
      do.points = FALSE,
      verticals = TRUE,
      lwd = 2,
      col="darkgreen")
lines(ecdf(-a$logFC),
      do.points = FALSE,
      verticals = TRUE,
      lwd = 2,
      col="olivedrab3")
lines(ecdf(-b$logFC),
      do.points = FALSE,
      verticals = TRUE,
      lwd = 2,
      col="tomato")
lines(ecdf(-c$logFC),
      do.points = FALSE,
      verticals = TRUE,
      lwd = 2,
      col="tomato4")

legend('bottomright', 
       legend=c("DAX","SMI","CAC","FTSE"),  # text in the legend
       col=c("seashell2","khaki","purple","red"),  # point colors
       pch=15)  # specify the point type to be a square

write.table(temp, "~/Downloads/aabiosynthesis_handcur", sep="\t")

