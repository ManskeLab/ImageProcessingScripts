import time
start_time = time.time()

import numpy as np
from vtk.util.numpy_support import vtk_to_numpy
from uniform_seg import Image, ImageBuilder, mask_calc
import ccl3D, cclBFS3D
import sys, os
import numpy.polynomial.polynomial as poly

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
  labels = [4,5,6]
  allLabelsToLabelToIndex = ccl3D.findMultipleLabelConnectedComponents(segmentData, labels)
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

def main():
  directory = sys.argv[1] # of image

  if (len(sys.argv) != 2):
    print("Usage: python getComponents.py <directory>")
    sys.exit(1)
  directory = sys.argv[1]

  files = os.listdir(directory)
  if ("notes.txt" not in files):
    sys.exit(1)

  isResampled = False

  if (isResampled):
    imageFiles = list(filter(lambda x: x.endswith(".nii") and "resampled" in x, files))
  else:
    imageFiles = list(filter(lambda x: x.endswith(".nii") and "resampled" not in x, files))
  imageFilename, segmentationFilename = (imageFiles[0], imageFiles[1]) if "seg" in imageFiles[1] else (imageFiles[1], imageFiles[0])

  image = ImageBuilder.createImage(imageFilename, segmentationFilename, directory)
  getCubesAndWriteValues(image, isResampled=isResampled)

if __name__=='__main__':
  main()
  end_time = time.time()
  print("--- Exec time: %.2f second ---" % (end_time - start_time))
