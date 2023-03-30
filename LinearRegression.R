firstfn <- "/Volumes/Seagate Backup Plus Drive/WBCT Reconstruction/segmentation images/coronal_wbct/profileResampled.csv"
secondfn <- "/Volumes/Seagate Backup Plus Drive/WBCT Reconstruction/segmentation images/coronal_clinical/profile.csv"

firstData <- read.csv(firstfn)
secondData <- read.csv(secondfn)

firstBmd <- firstData$BMD800
secondBmd <- secondData$BMD800

firstLength <- length(firstBmd)
secondLength <- length(secondBmd)

if (firstLength <= secondLength) {
  secondBmd <- secondBmd[1:firstLength]
} else {
  firstBmd <- firstBmd[1:secondLength]
}

model <- lm(secondBmd~firstBmd)

title <- expression("BMD 800 mgHA/cm"^3)
xLabel <- "WBCT.Ph"
yLabel <- "cCT.Ph"

plot(secondBmd~firstBmd, main = title, xlab = xLabel, ylab = yLabel, pch=19)
abline(model, col='red', lwd=2)
legend("topleft", 4, legend=c(paste("y =", format(model$coefficients[1],digits=7), "+", format(model$coefficients[2],digits=5),"x"), as.expression(bquote(r^2 * " = " * .(summary(model)$r.squared)))), bty='n')
summary(model)
