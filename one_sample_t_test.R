argv <- commandArgs(TRUE)

profileFilename <- argv[1]
profile <- read.csv(profileFilename)

print(paste("Mean of BMD 100:", mean(profile$BMD100)))
print(paste("SD of BMD 100: ", sd(profile$BMD100)))
print(paste("n =", length(profile$BMD100)))
print(paste("Mean of BMD 400: ", mean(profile$BMD400)))
print(paste("SD of BMD 400: ", sd(profile$BMD400)))
print(paste("n =", length(profile$BMD400)))
print(paste("Mean of BMD 800: ", mean(profile$BMD800)))
print(paste("SD of BMD 800: ", sd(profile$BMD800)))
print(paste("n =", length(profile$BMD800)))

print(t.test(profile$BMD100, mu = 100))
print(t.test(profile$BMD400, mu = 400))
print(t.test(profile$BMD800, mu = 800))