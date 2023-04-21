import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as poly
import sys, os
import scipy.stats

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

differencesToWrite = []

def plot(comb: tuple, directory: str) -> None:
  index1 = comb[0]
  index2 = comb[1]

  firstProfile100 = allProfiles100[index1]
  firstProfile400 = allProfiles400[index1]
  firstProfile800 = allProfiles800[index1]
  secondProfile100 = allProfiles100[index2]
  secondProfile400 = allProfiles400[index2]
  secondProfile800 = allProfiles800[index2]

  n = min(len(firstProfile100), len(secondProfile100))

  firstProfile100 = firstProfile100[0:n]
  firstProfile400 = firstProfile400[0:n]
  firstProfile800 = firstProfile800[0:n]
  secondProfile100 = secondProfile100[0:n]
  secondProfile400 = secondProfile400[0:n]
  secondProfile800 = secondProfile800[0:n]

  data = list(zip(firstProfile100, firstProfile400, firstProfile800,
    secondProfile100, secondProfile400, secondProfile800))
  df = pd.DataFrame(data)

  df["mean100"] = (firstProfile100 + secondProfile100) / 2
  df["mean400"] = (firstProfile400 + secondProfile400) / 2
  df["mean800"] = (firstProfile800 + secondProfile800) / 2

  df["difference100"] = firstProfile100 - secondProfile100
  df["difference400"] = firstProfile400 - secondProfile400
  df["difference800"] = firstProfile800 - secondProfile800

  allDifferences = np.concatenate((df["difference100"].to_numpy(),df["difference400"].to_numpy(),df["difference800"].to_numpy()))

  differencesToWrite.append((allDifferences, "%s_%s" % (scanToName[index1], scanToName[index2])))

  meanDiff100 = np.mean(df["difference100"])
  meanDiff400 = np.mean(df["difference400"])
  meanDiff800 = np.mean(df["difference800"])

  sd100 = np.std(df["difference100"], ddof=1)
  sd400 = np.std(df["difference400"], ddof=1)
  sd800 = np.std(df["difference800"], ddof=1)

  tStar = scipy.stats.t.ppf(0.975, df=n-1)

  individualN = 100
  x_axis100 = np.linspace(min(df["mean100"]), max(df["mean100"]), individualN)
  lowerCI100 = meanDiff100 - tStar * sd100
  upperCI100 = meanDiff100 + tStar * sd100

  x_axis400 = np.linspace(min(df["mean400"]), max(df["mean400"]), individualN)
  lowerCI400 = meanDiff400 - tStar * sd400
  upperCI400 = meanDiff400 + tStar * sd400

  x_axis800 = np.linspace(min(df["mean800"]), max(df["mean800"]), individualN)
  lowerCI800 = meanDiff800 - tStar * sd800
  upperCI800 = meanDiff800 + tStar * sd800

  plt.title("%s vs. %s - BMD %d" % (scanToName[index1], scanToName[index2], 100))
  plt.xlabel("Mean (%s and %s)" % (scanToName[index1], scanToName[index2]))
  plt.ylabel("Error (%s - %s)" % (scanToName[index1], scanToName[index2]))
  plt.ylim(-70,70)
  plt.plot(df["mean100"], df["difference100"], 'o',
    x_axis100, np.full(individualN, meanDiff100), '-',
    x_axis100, np.full(individualN, lowerCI100), '--',
    x_axis100, np.full(individualN, upperCI100), '--',
    color='cyan', markersize=4
  )
  filename = "%s_%s_BMD%d_BA.png" % (scanToName[index1], scanToName[index2], 100)
  plt.savefig(os.path.join(directory, filename))
  plt.clf()

  plt.title("%s vs. %s - BMD %d" % (scanToName[index1], scanToName[index2], 400))
  plt.xlabel("Mean (%s and %s)" % (scanToName[index1], scanToName[index2]))
  plt.ylabel("Error (%s - %s)" % (scanToName[index1], scanToName[index2]))
  plt.ylim(-70,70)
  plt.plot(df["mean400"], df["difference400"], 'o',
    x_axis400, np.full(individualN, meanDiff400), '-',
    x_axis400, np.full(individualN, lowerCI400), '--',
    x_axis400, np.full(individualN, upperCI400), '--',
    color='green', markersize=4
  )
  filename = "%s_%s_BMD%d_BA.png" % (scanToName[index1], scanToName[index2], 400)
  plt.savefig(os.path.join(directory, filename))
  plt.clf()

  plt.title("%s vs. %s - BMD %d" % (scanToName[index1], scanToName[index2], 800))
  plt.xlabel("Mean (%s and %s)" % (scanToName[index1], scanToName[index2]))
  plt.ylabel("Error (%s - %s)" % (scanToName[index1], scanToName[index2]))
  plt.ylim(-70,70)
  plt.plot(df["mean800"], df["difference800"], 'o',
    x_axis800, np.full(individualN, meanDiff800), '-',
    x_axis800, np.full(individualN, lowerCI800), '--',
    x_axis800, np.full(individualN, upperCI800), '--',
    color='orange', markersize=4
  )
  filename = "%s_%s_BMD%d_BA.png" % (scanToName[index1], scanToName[index2], 800)
  plt.savefig(os.path.join(directory, filename))
  plt.clf()

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

  d = {}
  for diff, name in differencesToWrite:
    d[name] = pd.Series(diff)

  df = pd.DataFrame(d)
  df.to_csv(os.path.join(directory, allPlotsDir, "differences.csv"))
 
if __name__ == "__main__":
  main()