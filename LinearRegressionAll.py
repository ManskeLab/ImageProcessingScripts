import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as poly
import sys, os
from uniform_seg import abline

# [clinical, wbct, wbct resampled, clinical, wbct, wbct resampled]
allProfiles100 = [None] * 6
allProfiles400 = [None] * 6
allProfiles800 = [None] * 6
allPlots = {100: allProfiles100, 400: allProfiles400, 800: allProfiles800}

scanToName = {
  0: "cCT.Ph",
  1: "WBCT.Ph",
  2: "WBCT.Ph",
  3: "cCT.Pa",
  4: "WBCT.Pa",
  5: "WBCT.Pa"
}

# 0: phantom clinical, 1: phantom wbct, 2: phantom wbct resampled, 3: patient clinical, 4: patient wbct, 5: patient wbct resampled
combinations = [(2,0), (1,4), (0,3), (5,3)]

def abline(intercept, slope, color, xAxis):
  """Plot a line from slope and intercept"""
  y_vals = intercept + slope * xAxis
  plt.plot(xAxis, y_vals, '-', label=f"$y = {slope:.2f}x+{intercept:.2f}$", color=color)

def plot(comb: tuple, directory: str) -> None:
  index1 = comb[0]
  index2 = comb[1]

  firstProfile100 = allProfiles100[index1]
  firstProfile400 = allProfiles400[index1]
  firstProfile800 = allProfiles800[index1]
  secondProfile100 = allProfiles100[index2]
  secondProfile400 = allProfiles400[index2]
  secondProfile800 = allProfiles800[index2]

  minLength = min(len(firstProfile100), len(secondProfile100))

  firstProfile100 = firstProfile100[1:minLength]
  firstProfile400 = firstProfile400[1:minLength]
  firstProfile800 = firstProfile800[1:minLength]
  secondProfile100 = secondProfile100[1:minLength]
  secondProfile400 = secondProfile400[1:minLength]
  secondProfile800 = secondProfile800[1:minLength]

  xProfile = np.concatenate((firstProfile100, firstProfile400, firstProfile800))
  yProfile = np.concatenate((secondProfile100, secondProfile400, secondProfile800))
  model = poly.polyfit(xProfile, yProfile, 1)

  x_axis = np.arange(min(xProfile), max(xProfile), 0.5)

  plt.title("%s vs. %s" % (scanToName[index1], scanToName[index2]))
  plt.xlabel("%s" % scanToName[index1])
  plt.ylabel("%s" % scanToName[index2])
  plt.plot(firstProfile100, secondProfile100, 'o',
    color='blue', markersize=1
  )
  plt.plot(firstProfile400, secondProfile400, 'o',
    color='green', markersize=1
  )
  plt.plot(firstProfile800, secondProfile800, 'o',
    color='orange', markersize=1
  )

  plt.plot(x_axis, poly.polyval(x_axis, model), 'k--')

  # x_i = xProfile
  # y_i = yProfile
  # \hat{y}_i = a + bx
  # \bar{y} = \frac{\sum_{i=1}^numCalibration y_i}{numCalibration}
  # RSE = \sqrt{\frac{\sum_{i=1}^numCalibration (y_i-\hat{y}_i)}{df}}

  y_hat = poly.polyval(xProfile, model)
  y_bar = np.mean(yProfile)
  rss = np.sum((yProfile - y_hat)**2)
  tss = np.sum((yProfile - y_bar)**2)
  R2 = 1 - rss/tss

  plt.text(min(xProfile), max(yProfile), "$y = %.2f + %.2fx$" % (model[0], model[1]), fontsize=7)
  plt.text(min(xProfile), max(yProfile) - 50, "$R^2 = %.4f$" % R2, fontsize=7)

  abline(0, 1, "grey", x_axis)

  filename = "%s_%s_reg.png" % (scanToName[index1], scanToName[index2])
  plt.savefig(os.path.join(directory, filename))
  plt.clf()

def main():
  if (len(sys.argv) != 2):
    print("Usage: python BlandAltmanPlotAll.py <directory>")
    sys.exit(1)
  directory = sys.argv[1]

  newDirName = "LinearRegressionPlots"
  allPlotsDir = os.path.join(directory, newDirName)
  os.makedirs(allPlotsDir, exist_ok=True)

  firstIteration = True
  for subdir, dirs, files in os.walk(directory):
    if firstIteration:
      firstIteration = False
    else:
      if "notes.txt" not in files:
        continue

      isWBCT = False
      isPatient = False
      with open(os.path.join(subdir, "notes.txt")) as f:
        line = f.readline()
        line = line.rstrip().split()
        if "WBCT" == line[0]:
          isWBCT = True

        if "patient" in line:
          isPatient = True
      
      index = 0
      if isWBCT:
        if isPatient:
          index = 4
        else:
          index = 1
      else:
        if isPatient:
          index = 3
        else:
          index = 0
      profileFilename = "profile.csv"

      dfProfile = pd.read_csv(os.path.join(subdir, profileFilename))

      allProfiles100[index] = dfProfile["BMD100"].to_numpy()
      allProfiles400[index] = dfProfile["BMD400"].to_numpy()
      allProfiles800[index] = dfProfile["BMD800"].to_numpy()

      if isWBCT:
        profileFilename = "profileResampled.csv"
        dfProfile = pd.read_csv(os.path.join(subdir, profileFilename))
        allProfiles100[index + 1] = dfProfile["BMD100"].to_numpy()
        allProfiles400[index + 1] = dfProfile["BMD400"].to_numpy()
        allProfiles800[index + 1] = dfProfile["BMD800"].to_numpy()

  for c in combinations:
    plot(c, os.path.join(directory, allPlotsDir))
 
if __name__ == "__main__":
  main()