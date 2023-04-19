import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as poly
import sys, os

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

  data = list(zip(firstProfile100, firstProfile400, firstProfile800,
    secondProfile100, secondProfile400, secondProfile800))
  df = pd.DataFrame(data)

  df["mean100"] = (firstProfile100 + secondProfile100) / 2
  df["mean400"] = (firstProfile400 + secondProfile400) / 2
  df["mean800"] = (firstProfile800 + secondProfile800) / 2

  df["difference100"] = firstProfile100 - secondProfile100
  df["difference400"] = firstProfile400 - secondProfile400
  df["difference800"] = firstProfile800 - secondProfile800

  meanDiff100 = np.mean(df["difference100"])
  meanDiff400 = np.mean(df["difference400"])
  meanDiff800 = np.mean(df["difference800"])

  model100 = poly.polyfit(df["mean100"], df["difference100"], 1)
  model400 = poly.polyfit(df["mean400"], df["difference400"], 1)
  model800 = poly.polyfit(df["mean800"], df["difference800"], 1)

  x_axis = np.linspace(min(df["mean100"]), max(df["mean800"]), len(df["mean100"]) + len(df["mean400"]) + len(df["mean800"]))

  plt.title("%s vs. %s" % (scanToName[index1], scanToName[index2]))
  plt.xlabel("Mean (%s and %s)" % (scanToName[index1], scanToName[index2]))
  plt.ylabel("Error (%s - %s)" % (scanToName[index1], scanToName[index2]))
  plt.ylim(-70,70)
  plt.plot(df["mean100"], df["difference100"], 'o')
  plt.plot(df["mean400"], df["difference400"], 'o')
  plt.plot(df["mean800"], df["difference800"], 'o')

  # filename = "%s_%s_BA.png" % (scanToName[index1], scanToName[index2])
  # plt.savefig(os.path.join(directory, filename))
  plt.show()

def main():
  if (len(sys.argv) != 2):
    print("Usage: python BlandAltmanPlotAll.py <directory>")
    sys.exit(1)
  directory = sys.argv[1]

  newDirName = "BlandAltmanPlots"
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