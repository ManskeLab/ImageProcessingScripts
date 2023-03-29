import time
start_time = time.time()

import vtk
import numpy as np
import numpy.polynomial.polynomial as poly
import pandas as pd
from vtk.util.numpy_support import vtk_to_numpy
import sys, os
from ccl3D import findConnectedComponents, findMultipleLabelConnectedComponents
from uniform_seg import Image, ImageBuilder
from uniform_seg import writeMapToDataFrameToCSVWithIndex, writeMapToDataFrameToCSV, mask_calc

def getCubesAndWriteValues(image: Image, isResampled: bool=False):
  imageData = vtk_to_numpy(image.getImage().GetOutput().GetPointData().GetScalars()).reshape(
    image.getImage().GetOutput().GetDimensions(), order='F'
  )

  segmentData = vtk_to_numpy(image.getSegmentation().GetOutput().GetPointData().GetScalars()).reshape(
    image.getSegmentation().GetOutput().GetDimensions(), order='F'
  )

  ((prof_100, prof_400, prof_800), attenuation_x, x_axis) = mask_calc(image, True)
  bmd_y = np.array([100, 400, 800])
  model = poly.polyfit(attenuation_x, bmd_y, 1)

  cubeDataNames = ["cube1_mean", "cube2_mean", "cube3_mean", "mean", "SD"]
  data = {"BMD100": [], "BMD400": [], "BMD800": []}

  cubeLabels = {4: "BMD800", 5: "BMD400", 6: "BMD100"}
  labels = [4, 5, 6]

  allLabelsToLabelToIndex = findMultipleLabelConnectedComponents(segmentData, labels)
  for label in labels:
    print("Label: %d" % label)
    labelToIndex = allLabelsToLabelToIndex[label]
    print(len(labelToIndex))

    allCubes = np.array([])
    for i in range(1, len(labelToIndex) + 1): # iterate 3 times
      cubeValues = np.zeros(len(labelToIndex[i]))
      for j in range(len(labelToIndex[i])): # iterate over each cube
        voxelValue = imageData[labelToIndex[i][j]]
        calibratedVoxelValue = poly.polyval(voxelValue, model)
        cubeValues[j] = calibratedVoxelValue
      data[cubeLabels[label]].append(np.mean(cubeValues))
      allCubes = np.append(allCubes, cubeValues)

    data[cubeLabels[label]].append(np.mean(allCubes))
    data[cubeLabels[label]].append(np.std(allCubes, ddof=1)) 
  
  if (isResampled):
    writeMapToDataFrameToCSVWithIndex(data, cubeDataNames, image.directory, "cubesResampled.csv")
  else:
    writeMapToDataFrameToCSVWithIndex(data, cubeDataNames, image.directory, "cubes.csv")

def getProfileAndWriteValues(image: Image, isResampled: bool=False):
  ((prof_100, prof_400, prof_800), attenuation_x, x_axis) = mask_calc(image, True)
  bmd_y = np.array([100, 400, 800])
  model = poly.polyfit(attenuation_x, bmd_y, 1)
  
  bmd_100_prof = model[0] + (prof_100 * model[1])
  bmd_400_prof = model[0] + (prof_400 * model[1])
  bmd_800_prof = model[0] + (prof_800 * model[1])

  data = {"BMD100": bmd_100_prof, "BMD400": bmd_400_prof, "BMD800": bmd_800_prof}
  if (isResampled):
    writeMapToDataFrameToCSV(data, image.directory, "profileResampled.csv")
  else:
    writeMapToDataFrameToCSV(data, image.directory, "profile.csv")
 
def main():
  directory = sys.argv[1] # of image

  if (len(sys.argv) != 2):
    print("Usage: python getComponents.py <directory>")
    sys.exit(1)
  directory = sys.argv[1]

  files = os.listdir(directory)
  if ("notes.txt" not in files):
    sys.exit(1)

  isResampled = 'y'
  print("Using resampled data? y/[n]")
  isResampled = input()
  if (isResampled == 'y'):
    isResampled = True
  elif (isResampled == 'n' or isResampled == ""):
    isResampled = False
  else:
    sys.exit(1)

  if (isResampled):
    imageFiles = list(filter(lambda x: x.endswith(".nii") and "resampled" in x, files))
  else:
    imageFiles = list(filter(lambda x: x.endswith(".nii") and "resampled" not in x, files))
  imageFilename, segmentationFilename = (imageFiles[0], imageFiles[1]) if "seg" in imageFiles[1] else (imageFiles[1], imageFiles[0])

  image = ImageBuilder.createImage(imageFilename, segmentationFilename, directory)
  getCubesAndWriteValues(image, isResampled=isResampled)
  getProfileAndWriteValues(image, isResampled=isResampled)

if __name__ == "__main__":
  main()
  end_time = time.time()
  print("--- Execution Time: %.2f seconds ---" % (end_time - start_time))
