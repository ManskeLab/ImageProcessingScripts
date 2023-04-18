import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys

colors = ["r", "b", "m", "orange"] # color order
spacing = 0.625

allProfiles100 = []
allProfiles400 = []
allProfiles800 = []
allPlots = {100: allProfiles100, 400: allProfiles400, 800: allProfiles800}
allCalibratedLimit100 = []
allCalibratedLimit400 = []
allCalibratedLimit800 = []
calibratedLimit = {100: allCalibratedLimit100, 400: allCalibratedLimit400, 800: allCalibratedLimit800}
trueLimits = {}
x_axes = []

nameToProperName = {"clinical": "Clinical", "wbct": "WBCT"}

def plot(bmd: int, x_axis: np.array, directory: str) -> None:
  allProfiles = allPlots[bmd]
  calibratedMeans = calibratedLimit[bmd]

  for i in range(len(allProfiles)):
    plt.plot(x_axis, allProfiles[i][len(allProfiles[i]) - len(x_axis):], "-", color=colors[i])
    plt.plot(x_axis, calibratedMeans[i][: len(x_axis)], "--", color=colors[i])
  plt.plot(x_axis, trueLimits[bmd], "k--")

  if (bmd == 100):
    plt.ylim(20, 180)
  elif (bmd == 400):
    plt.ylim(320, 480)
  else:
    plt.ylim(700, 860)
  plt.title("Uniformity Profile: BMD%d" % bmd)
  plt.xlabel("Distance (mm)")
  plt.ylabel("BMD (mgHA/cm3)")
  filename = "all_uniformity_profile_BMD%d.png" % bmd
  plt.savefig(os.path.join(directory, filename))
  plt.clf()

def main():
  if (len(sys.argv) != 2):
    print("Usage: python plotAll.py <directory>")
    sys.exit(1)
  directory = sys.argv[1]

  newDirName = "allPlots"
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
      with open(os.path.join(subdir, "notes.txt")) as f:
        line = f.readline()
        line = line.rstrip().split()
        if "WBCT" == line[0]:
          isWBCT = True

      if isWBCT:
        profileFilename = "profileResampled.csv"
        cubesFilename = "cubesResampled.csv"
      else:
        profileFilename = "profile.csv"
        cubesFilename = "cubes.csv"

      dfProfile = pd.read_csv(os.path.join(subdir, profileFilename))
      dfCubes = pd.read_csv(os.path.join(subdir, cubesFilename))

      current_x_axis = np.zeros_like(dfProfile["BMD100"].to_numpy())
      for i in range(len(current_x_axis)):
        current_x_axis[i] = i * spacing
      x_axes.append(current_x_axis)

      allCalibratedLimit100.append(np.full(len(current_x_axis), dfCubes.iloc[3]["BMD100"]))
      allCalibratedLimit400.append(np.full(len(current_x_axis), dfCubes.iloc[3]["BMD400"]))
      allCalibratedLimit800.append(np.full(len(current_x_axis), dfCubes.iloc[3]["BMD800"]))

      allProfiles100.append(dfProfile["BMD100"].to_numpy())
      allProfiles400.append(dfProfile["BMD400"].to_numpy())
      allProfiles800.append(dfProfile["BMD800"].to_numpy())

  x_axis = sorted(x_axes, key=len)[0]
  trueLimits[100] = np.full(len(x_axis), 100)
  trueLimits[400] = np.full(len(x_axis), 400)
  trueLimits[800] = np.full(len(x_axis), 800)
  plot(100, x_axis, os.path.join(directory, allPlotsDir))
  plot(400, x_axis, os.path.join(directory, allPlotsDir))
  plot(800, x_axis, os.path.join(directory, allPlotsDir))
 
if __name__ == "__main__":
  main()
