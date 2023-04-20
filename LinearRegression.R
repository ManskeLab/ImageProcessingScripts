firstfn <- "~/Documents/wbctbmddata/segmentation images/coronal_wbct/profileResampled.csv"
secondfn <- "~/Documents/wbctbmddata/segmentation images/coronal_clinical/profile.csv"

firstData <- read.csv(firstfn)
secondData <- read.csv(secondfn)

firstBmd100 <- firstData$BMD100
secondBmd100 <- secondData$BMD100
firstBmd400 <- firstData$BMD400
secondBmd400 <- secondData$BMD400
firstBmd800 <- firstData$BMD800
secondBmd800 <- secondData$BMD800

firstLength <- length(firstBmd800)
secondLength <- length(secondBmd800)

minLength <- min(firstLength, secondLength)

firstBmd100 <- firstBmd100[1:minLength]
secondBmd100 <- secondBmd100[1:minLength]
firstBmd400 <- firstBmd400[1:minLength]
secondBmd400 <- secondBmd400[1:minLength]
firstBmd800 <- firstBmd800[1:minLength]
secondBmd800 <- secondBmd800[1:minLength]

firstBmd <- c(firstBmd100, firstBmd400, firstBmd800)
secondBmd <- c(secondBmd100, secondBmd400, secondBmd800)

model <- lm(secondBmd~firstBmd)

newX <- seq(min(firstBmd), max(firstBmd), length.out=length(firstBmd))
prediction <- predict(model, newdata = data.frame(mean=newX), interval='confidence')

xLabel <- "WBCT.Ph"
yLabel <- "cCT.Ph"
title <- paste(xLabel, "vs.", yLabel) # expression("BMD 800 mgHA/cm"^3)

plot(secondBmd~firstBmd, main = title, xlab = xLabel, ylab = yLabel, pch=19)
abline(model, col='red', lwd=2)
lines(newX, prediction[,3], col='blue', lty = 2, lwd = 2)
lines(newX, prediction[,2], col='blue', lty = 2, lwd = 2)
legend("topleft", 4, legend=c(paste("y =", format(model$coefficients[1],digits=7), "+", format(model$coefficients[2],digits=5),"x"), as.expression(bquote(r^2 * " = " * .(summary(model)$r.squared)))), bty='n')
summary(model)
