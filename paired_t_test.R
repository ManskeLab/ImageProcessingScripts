argv <- commandArgs(TRUE)

firstfn <- argv[1]
secondfn <- argv[2]

firstData <- read.csv(firstfn)
secondData <- read.csv(secondfn)

bmd100First <- firstData$BMD100
bmd100Second <- secondData$BMD100
first100Length <- length(bmd100First)
second100Length <- length(bmd100Second)
if (first100Length <= second100Length) {
  bmd100Second <- bmd100Second[1:first100Length]
} else {
  bmd100First <- bmd100First[1:second100Length]
}

bmd400First <- firstData$BMD400
bmd400Second <- secondData$BMD400
first400Length <- length(bmd400First)
second400Length <- length(bmd400Second)
if (first400Length <= second400Length) {
  bmd400Second <- bmd400Second[1:first400Length]
} else {
  bmd400First <- bmd400First[1:second400Length]
}

bmd800First <- firstData$BMD800
bmd800Second <- secondData$BMD800
first800Length <- length(bmd800First)
second800Length <- length(bmd800Second)
if (first800Length <= second800Length) {
  bmd800Second <- bmd800Second[1:first800Length]
} else {
  bmd800First <- bmd800First[1:second800Length]
}
df <- data.frame(firstBmd100=bmd100First, secondBmd100=bmd100Second, firstBmd400=bmd400First, secondBmd400=bmd400Second, firstBmd800=bmd800First, secondBmd800=bmd800Second)

print(t.test(df$firstBmd100, df$secondBmd100, paired = TRUE))
print(t.test(df$firstBmd400, df$secondBmd400, paired = TRUE))
print(t.test(df$firstBmd800, df$secondBmd800, paired = TRUE))
