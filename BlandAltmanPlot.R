firstfn <- "/Volumes/Seagate Backup Plus Drive/WBCT Reconstruction/segmentation images/coronal_clinical/profile.csv"
secondfn <- "/Volumes/Seagate Backup Plus Drive/WBCT Reconstruction/segmentation images/coronal_wbct/profileResampled.csv"

firstData <- read.csv(firstfn)
secondData <- read.csv(secondfn)

firstBmd <- firstData$BMD100
secondBmd <- secondData$BMD100

firstLength <- length(firstBmd)
secondLength <- length(secondBmd)

if (firstLength <= secondLength) {
  secondBmd <- secondBmd[1:firstLength]
} else {
  firstBmd <- firstBmd[1:secondLength]
}

df <- data.frame(bmd1 = firstBmd, bmd2 = secondBmd)

df$mean <- rowMeans(df)

df$difference <- df$bmd1 - df$bmd2

meanDiff <- mean(df$difference)
print(meanDiff)

model <- lm(difference~mean, data = df)

newX <- seq(min(df$mean), max(df$mean), length.out=length(df$difference))
prediction <- predict(model, newdata = data.frame(mean=newX), interval='confidence')

title <- expression("BMD 100 mgHA/cm"^3)
xLabel <- "Mean (cCT.Ph and WBCT.Ph)"
yLabel <- "Error (cCT.Ph - WBCT.Ph)"

plot(difference~mean, data = df, main = title, xlab = xLabel, ylab = yLabel, pch=19)
text(x=max(df$mean)-0.05, y=meanDiff, bquote(bar(x) * " = " * .(meanDiff)))
abline(model, col='red', lwd=2)
abline(h=0,col='grey',lty=1, lwd=1) 
lines(newX, prediction[,3], col='blue', lty = 2, lwd = 2)
lines(newX, prediction[,2], col='blue', lty = 2, lwd = 2)
summary(model)
