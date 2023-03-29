import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from enum import Enum
from uniform_seg import PlotKey

plotTogether = {}

def dualPlotNoResampleSameXEnd(key1: PlotKey, key2: PlotKey, directory: str, bmd: int) -> None:
  (x_axis1, bmd_prof1, bmdTrueLimit1, bmd_limit1) = plotTogether[key1]
  (x_axis2, bmd_prof2, bmdTrueLimit2, bmd_limit2) = plotTogether[key2]

  if (x_axis1[0] <= x_axis2[0]):
    diff = x_axis2[0] - x_axis1[0]
    x_axis2 -= diff
  else:
    diff = x_axis1[0] - x_axis2[0]
    x_axis1 -= diff

  if (x_axis1[-1] <= x_axis2[-1]):
    xEnd = x_axis1[-1]
    smallestIndex = 0
    smallestDifference = abs(x_axis2[0] - xEnd)
    for i in range(1, len(x_axis2)):
      if (abs(x_axis2[i] - xEnd) < smallestDifference):
        smallestIndex = i
        smallestDifference = abs(x_axis2[i] - xEnd)
    x_axis2 = x_axis2[0:smallestIndex+1]
    bmd_prof2 = bmd_prof2[0:smallestIndex+1]
    bmdTrueLimit2 = bmdTrueLimit2[0:smallestIndex+1]
  else:
    xEnd = x_axis2[-1]
    smallestIndex = 0
    smallestDifference = abs(x_axis1[0] - xEnd)
    for i in range(1, len(x_axis1)):
      if (abs(x_axis1[i] - xEnd) < smallestDifference):
        smallestIndex = i
        smallestDifference = abs(x_axis1[i] - xEnd)
    x_axis1 = x_axis1[0:smallestIndex+1]
    bmd_prof1 = bmd_prof1[0:smallestIndex+1]
    bmdTrueLimit1 = bmdTrueLimit1[0:smallestIndex+1]
    bmd_limit1 = bmd_limit1[0:smallestIndex+1]

  plt.plot(x_axis1, bmd_prof1, 'r-', x_axis2, bmd_prof2, 'b-')
  plt.plot(x_axis1, bmdTrueLimit1, 'k--', x_axis2, bmdTrueLimit2, 'y--', x_axis1, bmd_limit1, 'm--')
  if (bmd == 100):
    plt.ylim(0, 160)
  elif (bmd == 400):
    plt.ylim(300, 460)
  else:
    plt.ylim(700, 860)
  plt.title("Uniformity Profile %d: %s vs. %s" % (bmd, key1.name, key2.name))
  plt.xlabel("Distance (mm)")
  plt.ylabel("BMD (mgHA/cm3)")
  firstName = " - ".join(key1.name.split("Coronal"))
  secondName = " - ".join(key2.name.split("Coronal"))
  plt.legend([firstName, secondName, firstName + " $\it{calib. avg}$", secondName + " $\it{calib. avg}$", "True BMD"], loc="best")
  filename = "%s_%s_BMD%d.png" % (key1.name, key2.name, bmd)
  plt.savefig(os.path.join(directory, filename))
  plt.clf()

def dualPlotResampled(key1: PlotKey, key2: PlotKey, directory: str, bmd: int) -> None:
  (x_axis1, bmd_prof1, bmdTrueLimit1, bmd_limit1) = plotTogether[key1]
  (x_axis2, bmd_prof2, bmdTrueLimit2, bmd_limit2) = plotTogether[key2]

  x_axis, bmd_800_limit = (x_axis1, bmd_limit1) if (len(x_axis1) <= len(x_axis2)) else (x_axis2, bmd_limit2)
  # bmd_prof1 = bmd_prof1[:len(x_axis)] # beginning of profile
  # bmd_prof2 = bmd_prof2[:len(x_axis)] # beginning of profile
  bmd_prof1 = bmd_prof1[len(bmd_prof1) - len(x_axis):len(bmd_prof1)] # end of profile
  bmd_prof2 = bmd_prof2[len(bmd_prof2) - len(x_axis):len(bmd_prof2)] # end of profile
  bmdTrueLimit1 = bmdTrueLimit1[:len(x_axis)]
  bmdTrueLimit2 = bmdTrueLimit2[:len(x_axis)]

  plt.plot(x_axis, bmd_prof1, 'r-', x_axis, bmd_prof2, 'b-')
  plt.plot(x_axis, bmdTrueLimit1, 'k--', x_axis, bmdTrueLimit2, 'y--', x_axis, bmd_800_limit, 'm--')
  if (bmd == 100):
    plt.ylim(0, 160)
  elif (bmd == 400):
    plt.ylim(300, 460)
  else:
    plt.ylim(700, 860)
  plt.title("Uniformity Profile %d: %s vs. %s" % (bmd, key1.name, key2.name))
  plt.xlabel("Distance (mm)")
  plt.ylabel("BMD (mgHA/cm3)")
  firstName = " - ".join(key1.name.split("Coronal"))
  secondName = " - ".join(key2.name.split("Coronal"))
  plt.legend([firstName, secondName, firstName + " $\it{calib. avg}$", secondName + " $\it{calib. avg}$", "True BMD"], loc="best")
  filename = "%s_%s_BMD%d.png" % (key1.name, key2.name, bmd)
  plt.savefig(os.path.join(directory, filename))
  plt.clf()
  
def plot(bmd: int, isResampled: bool, directory: str, saveDirectory: str):
  filename = "dualPlots%d.txt" % (bmd) if not isResampled else "dualPlots%dResampled.txt" % (bmd)
  with open(os.path.join(directory, filename), "r") as f:
    currentKey = None
    data = []
    for line in f:
      line = line.rstrip()
      if (line == ""):
        break
      if (len(line) == 1):
        if (currentKey is not None):
          plotTogether[currentKey] = data
          data = []
        for x in PlotKey:
          if (int(line) == x.value):
            currentKey = x
      else:
        arr = np.array(list(map(np.float64, line.split(','))))
        data.append(arr)
    if (len(data) > 0):
      plotTogether[currentKey] = data

  combinations = [
    (PlotKey.wbctCoronalPhantom, PlotKey.clinicalCoronalPhantom),
    (PlotKey.wbctCoronalPhantom, PlotKey.wbctCoronalPatient),
    (PlotKey.wbctCoronalPatient, PlotKey.clinicalCoronalPatient),
    (PlotKey.clinicalCoronalPhantom, PlotKey.clinicalCoronalPatient)
  ]
  for i in range(len(combinations)):
    if (isResampled):
      dualPlotResampled(combinations[i][0], combinations[i][1], saveDirectory, bmd)
    else:
      dualPlotNoResampleSameXEnd(combinations[i][0], combinations[i][1], saveDirectory, bmd)

def main():
  if (len(sys.argv) != 2):
    print("Usage: python dualPlot.py <directory>")
    sys.exit(1)
  directory = sys.argv[1]

  isResampled = 'y'
  print("Using resampled data? [y]/n")
  isResampled = input()
  if (isResampled == 'y'):
    isResampled = True
  elif (isResampled == 'n' or isResampled == ""):
    isResampled = False
  else:
    sys.exit(1)
  
  newDirName = "dual_plots" if not isResampled else "dual_plots_resampled"
  dualPlotsDirectory = os.path.join(directory, newDirName)
  os.makedirs(dualPlotsDirectory, exist_ok=True)
  plot(100, isResampled, directory, dualPlotsDirectory)
  plot(400, isResampled, directory, dualPlotsDirectory)
  plot(800, isResampled, directory, dualPlotsDirectory)
 
if __name__ == "__main__":
  main()
