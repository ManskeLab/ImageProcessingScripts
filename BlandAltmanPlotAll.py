import numpy as np
import pandas as pd
import sys, os

# [clinical, wbct, wbct resampled, clinical, wbct, wbct resampled]
allProfiles100 = []
allProfiles400 = []
allProfiles800 = []
allPlots = {100: allProfiles100, 400: allProfiles400, 800: allProfiles800}

# 0: phantom clinical, 1: phantom wbct, 2: phantom wbct resampled, 3: patient clinical, 4: patient wbct, 5: patient wbct resampled
combinations = [(2,0), (1,4), (0,3), (5,3)]

def plot():
  for c in combinations:
    print(c)

def main():
  if (len(sys.argv) != 2):
    print("Usage: python BlandAltmanPlotAll.py <directory>")
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

      profileFilename = "profile.csv"
      # cubesFilename = "cubes.csv"

      dfProfile = pd.read_csv(os.path.join(subdir, profileFilename))
      # dfCubes = pd.read_csv(os.path.join(subdir, cubesFilename))

      allProfiles100.append(dfProfile["BMD100"].to_numpy())
      allProfiles400.append(dfProfile["BMD400"].to_numpy())
      allProfiles800.append(dfProfile["BMD800"].to_numpy())

      if isWBCT:
        profileFilename = "profileResampled.csv"
        dfProfile = pd.read_csv(os.path.join(subdir, profileFilename))
        allProfiles100.append(dfProfile["BMD100"].to_numpy())
        allProfiles400.append(dfProfile["BMD400"].to_numpy())
        allProfiles800.append(dfProfile["BMD800"].to_numpy())

  plot(100, x_axis, os.path.join(directory, allPlotsDir))
  plot(400, x_axis, os.path.join(directory, allPlotsDir))
  plot(800, x_axis, os.path.join(directory, allPlotsDir))
 
if __name__ == "__main__":
  main()