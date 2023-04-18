import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from enum import Enum
from uniform_seg import PlotKey

colors = ["blue", "red", "purple", "orange"]

plotTogether = {}

nameToProperName = {"clinical": "Clinical", "wbct": "WBCT"}

def plotAllFour(bmd: int, isResampled: bool, directory: str) -> None:
  pass
  
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
  print("Using resampled data? y/[n]")
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
