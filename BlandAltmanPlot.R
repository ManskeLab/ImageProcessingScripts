argv <- commandArgs(TRUE)

firstfn <- "/Volumes/Seagate Backup Plus Drive/WBCT Reconstruction/segmentation images/coronal_clinical/profile.csv"
secondfn <- "/Volumes/Seagate Backup Plus Drive/WBCT Reconstruction/segmentation images/PKOA_112_clinical/profile.csv"

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

df <- data.frame(bmd1 = firstBmd, bmd2 = secondBmd)

df$mean <- rowMeans(df)

df$difference <- df$bmd1 - df$bmd2

meanDiff <- mean(df$difference)

model <- lm(df$difference~df$mean)
model

# confidence interval of mean difference
zStar <- qnorm(0.975)
lowerLimit <- meanDiff - zStar * sd(df$difference)
upperLimit <- meanDiff + zStar * sd(df$difference)

plot(df$mean, df$difference, main = "Bland-Altman Plot", xlab = "Mean", ylab = "Difference", ylim = c(lowerLimit - 8, upperLimit + 5))
abline(h=0,col='grey',lty=1, lwd=1) 
abline(h=meanDiff,col='blue',lty=1, lwd=1.5)
abline(h=lowerLimit,col='red',lty=2, lwd=1.5)
abline(h=upperLimit,col='red',lty=2, lwd=1.5)
abline(model)

summary(model)