firstfn <- "~/Documents/wbctbmddata/segmentation images/coronal_wbct/profileResampled.csv"
secondfn <- "~/Documents/wbctbmddata/segmentation images/coronal_clinical/profile.csv"

firstData <- read.csv(firstfn)
secondData <- read.csv(secondfn)

firstBmd <- firstData$BMD100
secondBmd <- secondData$BMD100

firstLength <- length(firstBmd)
secondLength <- length(secondBmd)
minLength <- min(firstLength, secondLength)

firstBmd <- firstBmd[1:minLength]
secondBmd <- secondBmd[1:minLength]

df <- data.frame(bmd1 = firstBmd, bmd2 = secondBmd)

df$mean <- rowMeans(df)

df$difference <- df$bmd1 - df$bmd2

meanDiff <- mean(df$difference)
print(meanDiff)

model <- lm(difference~mean, data = df)

tStar <- qt(0.975, df = minLength - 1)

lowerCI <- meanDiff - tStar * sd(df$difference)
upperCI <- meanDiff + tStar * sd(df$difference)

newX <- seq(min(df$mean), max(df$mean), length.out=length(df$difference))
prediction <- predict(model, newdata = data.frame(mean=newX), interval='confidence')

title <- expression("BMD 100 mgHA/cm"^3)
xLabel <- "Mean (WBCT.Pa and cCT.Pa)"
yLabel <- "Error (WBCT.Pa - cCT.Pa)"

plot(difference~mean, data = df, main = title, xlab = xLabel, ylab = yLabel, pch=19, ylim=c(-70,70))
text(x=max(df$mean), y=meanDiff, bquote(bar(x) * " = " * .(format(meanDiff, digits=4))))
abline(model, col='red', lwd=2)
abline(h=0,col='grey',lty=1, lwd=1) 
abline(h=lowerCI, col='blue', lty = 2, lwd = 2)
abline(h=upperCI, col='blue', lty = 2, lwd = 2)
# lines(newX, prediction[,3], col='blue', lty = 2, lwd = 2)
# lines(newX, prediction[,2], col='blue', lty = 2, lwd = 2)
summary(model)
