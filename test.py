import vtk
import numpy as np
import numpy.polynomial.polynomial as poly
import pandas as pd
from vtk.util.numpy_support import vtk_to_numpy
import sys, os
from ConnectedComponentLabeling3D import findConnectedComponents
from uniform_seg import Image, ImageBuilder, Scanner, Plane, Alg
from uniform_seg import writeMapToDataFrameToCSVWithIndex, writeMapToDataFrameToCSV, mask_calc

imageFilename = sys.argv[1]
segmentFilename = sys.argv[2]
reader = vtk.vtkNIFTIImageReader()
reader.SetFileName(imageFilename)
reader.Update()
image = reader.GetOutput()
reader.SetFileName(segmentFilename)
reader.Update()
segment = reader.GetOutput()

imageData = vtk_to_numpy(image.GetPointData().GetScalars()).reshape(
  image.GetDimensions(), order='F'
)

segmentData = vtk_to_numpy(segment.GetPointData().GetScalars()).reshape(
  segment.GetDimensions(), order='F'
)

m = findConnectedComponents(imageData, 4)
print(len(m))

