import vtk
import sys
import os 

filename = sys.argv[1]
spacing = list(map(float, sys.argv[2:5]))

reader = vtk.vtkNIFTIImageReader()
reader.SetFileName(filename)
reader.Update()

resampledImage = vtk.vtkImageResample()
resampledImage.SetInputData(reader.GetOutput())
resampledImage.SetOutputSpacing(spacing)
resampledImage.Update()

newFilename = filename[:-4] + "_resampled.nii"

temp =  vtk.vtkNIFTIImageWriter()
temp.SetInputData(resampledImage.GetOutput())
temp.SetFileName(newFilename)
temp.SetNIFTIHeader(reader.GetNIFTIHeader())
temp.Write()

# writer = vtk.vtkImageWriter()
# writer.SetInputData(resampledImage.GetOutput())
# writer.SetFileName(newFilename)
# writer.Update()
