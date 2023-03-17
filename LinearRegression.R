firstfn <- "/Volumes/Seagate Backup Plus Drive/WBCT Reconstruction/segmentation images/coronal_wbct/profile.csv"
secondfn <- "/Volumes/Seagate Backup Plus Drive/WBCT Reconstruction/segmentation images/PKOA_112_wbct/profile.csv"

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
xLabel <- "Mean (cCT.Ph and cCT.Pa)"
yLabel <- "Error (cCT.Ph - cCT.Pa)"

legend()

plot(firstBmd, secondBmd, main = title, xlab = xLabel, ylab = yLabel)