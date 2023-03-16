import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from enum import Enum

plotTogether = {}

class PlotKey(Enum):
  clinicalAxialPhantom = 0
  wbctAxialPhantom = 1
  clinicalCoronalPhantom = 2
  wbctCoronalPhantom = 3
  wbctCoronalPatient = 4
  clinicalCoronalPatient = 5

# no resample, same xend
# if (x_axis1[0] <= x_axis2[0]):
  # diff = x_axis2[0] - x_axis1[0]
  # x_axis2 -= diff
# else:
  # diff = x_axis1[0] - x_axis2[0]
  # x_axis1 -= diff
# 
# if (x_axis1[-1] <= x_axis2[-1]):
  # xEnd = x_axis1[-1]
  # smallestIndex = 0
  # smallestDifference = abs(x_axis2[0] - xEnd)
  # for i in range(1, len(x_axis2)):
    # if (abs(x_axis2[i] - xEnd) < smallestDifference):
      # smallestIndex = i
      # smallestDifference = abs(x_axis2[i] - xEnd)
  # x_axis2 = x_axis2[0:smallestIndex+1]
  # bmd_800_prof2 = bmd_800_prof2[0:smallestIndex+1]
  # bmd800TrueLimit2 = bmd800TrueLimit2[0:smallestIndex+1]
# else:
  # xEnd = x_axis2[-1]
  # smallestIndex = 0
  # smallestDifference = abs(x_axis1[0] - xEnd)
  # for i in range(1, len(x_axis1)):
    # if (abs(x_axis1[i] - xEnd) < smallestDifference):
      # smallestIndex = i
      # smallestDifference = abs(x_axis1[i] - xEnd)
  # x_axis1 = x_axis1[0:smallestIndex+1]
  # bmd_800_prof1 = bmd_800_prof1[0:smallestIndex+1]
  # bmd800TrueLimit1 = bmd800TrueLimit1[0:smallestIndex+1]
  # bmd_800_limit1 = bmd_800_limit1[0:smallestIndex+1]

def dualPlot(key1: PlotKey, key2: PlotKey, directory: str) -> None:
  (x_axis1, bmd_800_prof1, bmd800TrueLimit1, bmd_800_limit1) = plotTogether[key1]
  (x_axis2, bmd_800_prof2, bmd800TrueLimit2, bmd_800_limit2) = plotTogether[key2]

  x_axis, bmd_800_limit = (x_axis1, bmd_800_limit1) if (len(x_axis1) <= len(x_axis2)) else (x_axis2, bmd_800_limit2)
  bmd_800_prof1 = bmd_800_prof1[:len(x_axis)]
  bmd_800_prof2 = bmd_800_prof2[:len(x_axis)]
  bmd800TrueLimit1 = bmd800TrueLimit1[:len(x_axis)]
  bmd800TrueLimit2 = bmd800TrueLimit2[:len(x_axis)]

  plt.plot(x_axis, bmd_800_prof1, 'r-', x_axis, bmd_800_prof2, 'b-')
  plt.plot(x_axis, bmd800TrueLimit1, 'k--', x_axis, bmd800TrueLimit2, 'y--', x_axis, bmd_800_limit, 'm-')
  # plt.plot(x_axis1, bmd_800_prof1, 'r-', x_axis2, bmd_800_prof2, 'b-')
  # plt.plot(x_axis1, bmd800TrueLimit1, 'k--', x_axis2, bmd800TrueLimit2, 'y--', x_axis1, bmd_800_limit1, 'm-')
  plt.ylim(700, 860)
  plt.title("Uniformity Profile: %s vs. %s" % (key1.name, key2.name))
  plt.xlabel("Distance (mm)")
  plt.ylabel("BMD (mgHA/cm3)")
  filename = "%s_%s.png" % (key1.name, key2.name)
  plt.savefig(os.path.join(directory, filename))
  plt.clf()

def main():
  if (len(sys.argv) != 2):
    print("Usage: python dualPlot.py <directory>")
    sys.exit(1)
  directory = sys.argv[1]
  
  with open(os.path.join(directory, "dualPlots.txt"), "r") as f:
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

  otherPlotsDirectory = os.path.join(directory, "other_plots")
  os.makedirs(otherPlotsDirectory, exist_ok=True)
  combinations = [
    (PlotKey.clinicalCoronalPhantom, PlotKey.wbctCoronalPhantom),
    (PlotKey.wbctCoronalPatient, PlotKey.wbctCoronalPhantom),
    (PlotKey.clinicalCoronalPatient, PlotKey.wbctCoronalPatient),
    (PlotKey.clinicalCoronalPatient, PlotKey.clinicalCoronalPhantom)
  ]
  for i in range(len(combinations)):
    dualPlot(combinations[i][0], combinations[i][1], otherPlotsDirectory)

 
if __name__ == "__main__":
  main()
