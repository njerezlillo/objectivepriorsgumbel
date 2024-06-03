load("simu1.RData")

pdf("geweke1.pdf", width = 10, height = 6)
par(mar = c(4, 4.4, 1.5, 2), mfrow = c(2, 3))
x_axis <- paste(seq(5, 60, by = 5))

box1 <- data.frame(t(gematriz1))
colnames(box1) <- x_axis
boxplot(box1, xlab = "n", col = "dodgerblue3",
        ylab = expression(paste("Geweke test"," (", mu,")")))

box3 <- data.frame(t(gematriz3))
colnames(box3) <- x_axis
boxplot(box3, xlab = "n", col = "goldenrod3",
        ylab = expression(paste("Geweke test"," (", mu,")")))

box5 <- data.frame(t(gematriz5))
colnames(box5) <- x_axis
boxplot(box5, xlab = "n", col = "darkgreen",
        ylab = expression(paste("Geweke test"," (", mu,")")))

box2 <- data.frame(t(gematriz2))
colnames(box2) <- x_axis
boxplot(box2, xlab = "n", col = "dodgerblue3",
        ylab = expression(paste("Geweke test"," (", varphi,")")))

box4 <- data.frame(t(gematriz4))
colnames(box4) <- x_axis
boxplot(box4, xlab = "n", col = "goldenrod3",
        ylab = expression(paste("Geweke test"," (", varphi,")")))

box6 <- data.frame(t(gematriz6))
colnames(box6) <- x_axis
boxplot(box6, xlab = "n", col = "darkgreen",
        ylab = expression(paste("Geweke test"," (", varphi,")")))
dev.off()