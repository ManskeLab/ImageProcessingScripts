"""
Cylinder profile analysis
-------------------------------------------------------------------------------
  This script will create a cylinder profile analysis from
  1. the image 2. the segmentation of a cylinder
  Will run for multiple images
-------------------------------------------------------------------------------
  Invoke script by:
  python uniform_seg.py <directory of directories with image and segmentation>
"""
import time
start_time = time.time()
import os
import vtk
import tkinter as tk
from tkinter import filedialog as fd
import numpy as np
import numpy.polynomial.polynomial as poly
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib.pyplot as plt
import sys # pass in arguments
from enum import Enum
import pandas as pd

numCalibration = 3
numEstimatedParameters = 2
df = numCalibration - numEstimatedParameters

attributeFile = "notes.txt"

calibrationCoefficients = []
calibrationPoints = []

isFirstTimeOpeningDualPlots = True

# keys
class PlotKey(Enum):
  clinicalAxialPhantom = 0
  wbctAxialPhantom = 1
  clinicalCoronalPhantom = 2
  wbctCoronalPhantom = 3
  wbctCoronalPatient = 4
  clinicalCoronalPatient = 5

class Scanner(Enum):
  CLINICAL = "clinical"
  WBCT = "WBCT"
  XCT2 = "XCT2"

class Plane(Enum):
  CORONAL = "coronal"
  SAGITTAL = "sagittal"
  AXIAL = "axial"

class Alg(Enum):
  NONE = "none"
  OLD = "old"
  NEW = "new"

class Image:
  def __init__(self, image, segmentation, directory: str, scanner: Scanner, plane: Plane, patient: bool=False, algUsed: Alg=Alg.NONE) -> None:
    self.directory = directory
    self.image = Image.create_reader(os.path.join(self.directory, image))
    self.segmentation = Image.create_reader(os.path.join(self.directory, segmentation))
    self.imageFilename = image
    self.segmentationFilename = segmentation
    self.scanner = scanner
    self.plane = plane
    self.spacing = Image.getSpacingFromScannerAndPlane(self.scanner, self.plane)
    self.patient = patient
    self.alg = Alg.NONE if self.scanner != Scanner.WBCT else algUsed

    print("  Scanner: %s" % self.scanner.value)
    print("  Plane: %s" % self.plane.value)
    print("  Algorithm used (WBCT): %s" % self.alg.value)
    print("  Patient in image: %s\n" % str(self.patient))
  
  def getImage(self):
    return self.image
  
  def getSegmentation(self):
    return self.segmentation

  def getImageFilename(self):
    return self.imageFilename
  
  def getSegmentationFilename(self):
    return self.segmentationFilename

  def getScanner(self) -> Scanner:
    return self.scanner

  def getPlane(self) -> Plane:
    return self.plane

  def getSpacing(self) -> float:
    return self.spacing

  def getPatient(self) -> bool:
    return self.patient

  def getAlg(self) -> Alg:
    return self.alg

  # Creater reader for chosen file type
  @staticmethod
  def create_reader(fn):
    if fn.endswith('.nii'): # Selects the NIFTI image reader if filename ends with '.nii'
      reader = vtk.vtkNIFTIImageReader()
      reader.SetFileName(fn)
      print("", fn)
    elif fn.endswith('.obj'):   # Selects the OBJreader if filename ends with '.obj'
      reader = vtk.vtkOBJReader()
      reader.SetFileName(fn)
      print("", fn)
    elif fn.endswith('.dcm'):    # Select the DICOM image reader for DCM folders
      reader = vtk.vtkDICOMImageReader()
      reader.SetDirectoryName(os.path.dirname(fn))
      print("", reader.GetDirectoryName())
    else:
      raise ValueError("Please select a dicom/nifti/obj file")
    reader.Update()
    return reader

  @staticmethod
  def getSpacingFromScannerAndPlane(typeOfScanner: Scanner, plane: Plane) -> float:
    spacing = 0
    if (typeOfScanner == Scanner.XCT2):
      spacing = 0.0607
    elif (typeOfScanner == Scanner.CLINICAL):
      if (plane == Plane.CORONAL or plane == Plane.SAGITTAL):
        spacing = 0.625
      else:
        spacing = 0.5352
    elif (typeOfScanner == Scanner.WBCT):
      spacing = 0.3
    else:
      spacing = 0

    return spacing

class ImageBuilder:
  @staticmethod
  def createImage(image, segmentation, directory: str) -> Image:
    scanner = None
    plane = None
    alg = Alg.NONE
    havePatient = False
    with open(os.path.join(directory, attributeFile), "r") as f:
      line = f.readline()
      line = line.rstrip().split()
      if ("WBCT" == line[0]):
        scanner = Scanner.WBCT
      elif ("clinical" == line[0]):
        scanner = Scanner.CLINICAL
      else:
        scanner = Scanner.XCT2

      if ("coronal" == line[1]):
        plane = Plane.CORONAL
      elif ("sagittal" == line[1]):
        plane = Plane.SAGITTAL
      else:
        plane = Plane.AXIAL

      if ("old" in line):
        alg = Alg.OLD
      elif ("new" in line):
        alg = Alg.NEW
      else:
        alg = Alg.NONE
      
      if ("patient" in line):
        havePatient = True
      else:
        havePatient = False
    return Image(image, segmentation, directory, scanner, plane, havePatient, alg)

def plotCalibration(x_axis, y_points, y, title, xlabel, ylabel, text, filename):
  plt.plot(x_axis, y_points, 'ks', x_axis, y,'k--')
  plt.title(title)
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  plt.text(1000, 200, text)
  plt.savefig(filename)
  plt.clf()

def plot(x_axis, y, limit, calibratedBMD, title, xlabel, ylabel, colour, filename):
  plt.plot(x_axis, y, colour + '-', x_axis, limit, 'm-', x_axis, calibratedBMD, 'k--')
  plt.title(title)
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  plt.savefig(filename)
  plt.clf()

# Obtain the line profiles and perform image calibration using segmentation image
def mask_calc(image: Image, calibrate = False):
  vtk_image_data = image.getImage().GetOutput()
  vtk_segment_data = image.getSegmentation().GetOutput()

  numpy_image_data = vtk_to_numpy(vtk_image_data.GetPointData().GetScalars()).reshape(vtk_image_data.GetDimensions(), order = 'F')

  numpy_segment_data = vtk_to_numpy(vtk_segment_data.GetPointData().GetScalars()).reshape(vtk_segment_data.GetDimensions(), order = 'F')
  extent = numpy_segment_data.shape

  prof_800 = np.zeros(extent[2]) if image.getPlane() == Plane.CORONAL or image.getPlane() == Plane.SAGITTAL else np.zeros(extent[0]) if image.getScanner() == Scanner.WBCT else np.zeros(extent[1])
  prof_400 = np.zeros(extent[2]) if image.getPlane() == Plane.CORONAL or image.getPlane() == Plane.SAGITTAL else np.zeros(extent[0]) if image.getScanner() == Scanner.WBCT else np.zeros(extent[1])
  prof_100 = np.zeros(extent[2]) if image.getPlane() == Plane.CORONAL or image.getPlane() == Plane.SAGITTAL else np.zeros(extent[0]) if image.getScanner() == Scanner.WBCT else np.zeros(extent[1])
  x_axis = np.zeros_like(prof_800)

  loopEnd = extent[2] if image.getPlane() == Plane.CORONAL or image.getPlane() == Plane.SAGITTAL else extent[0] if image.getScanner() == Scanner.WBCT else extent[1]

  for i in range(loopEnd):
    try:
      index_800 = np.where(numpy_segment_data[:,:,i] == 1) if image.getPlane() == Plane.CORONAL or image.getPlane() == Plane.SAGITTAL else np.where(numpy_segment_data[i,:,:] == 1) if image.getScanner() == Scanner.WBCT else np.where(numpy_segment_data[:,i,:] == 1)
      if np.size(index_800[0]) > 0:
        im_800_data = numpy_image_data[:,:,i][index_800] if image.getPlane() == Plane.CORONAL or image.getPlane() == Plane.SAGITTAL else numpy_image_data[i,:,:][index_800] if image.getScanner() == Scanner.WBCT else numpy_image_data[:,i,:][index_800]
        prof_800[i] = np.mean(im_800_data)
        x_axis[i] = i * image.getSpacing()

      index_400 = np.where(numpy_segment_data[:,:,i] == 2) if image.getPlane() == Plane.CORONAL or image.getPlane() == Plane.SAGITTAL else np.where(numpy_segment_data[i,:,:] == 2) if image.getScanner() == Scanner.WBCT else np.where(numpy_segment_data[:,i,:] == 2)
      if np.size(index_400[0]) > 0:
        im_400_data = numpy_image_data[:,:,i][index_400] if image.getPlane() == Plane.CORONAL or image.getPlane() == Plane.SAGITTAL else numpy_image_data[i,:,:][index_400] if image.getScanner() == Scanner.WBCT else numpy_image_data[:,i,:][index_400]
        prof_400[i] = np.mean(im_400_data)

      index_100 = np.where(numpy_segment_data[:,:,i] == 3) if image.getPlane() == Plane.CORONAL or image.getPlane() == Plane.SAGITTAL else np.where(numpy_segment_data[i,:,:] == 3) if image.getScanner() == Scanner.WBCT else np.where(numpy_segment_data[:,i,:] == 3)
      if np.size(index_100[0]) > 0:
        im_100_data = numpy_image_data[:,:,i][index_100] if image.getPlane() == Plane.CORONAL or image.getPlane() == Plane.SAGITTAL else numpy_image_data[i,:,:][index_100] if image.getScanner() == Scanner.WBCT else numpy_image_data[:,i,:][index_100]
        prof_100[i] = np.mean(im_100_data)
    except:
      continue
  
  start = np.nonzero(x_axis)[0][0]
  end = np.nonzero(x_axis)[0][-1]

  prof_800 = prof_800[ start : end + 1 ]
  prof_400 = prof_400[ start : end + 1 ]
  prof_100 = prof_100[ start : end + 1 ]
  x_axis = x_axis[ start : end + 1 ]

  if calibrate:
    bmd_800 = np.where(numpy_segment_data == 4)
    bmd_800_cal = np.mean(numpy_image_data[bmd_800])

    bmd_400 = np.where(numpy_segment_data == 5)
    bmd_400_cal = np.mean(numpy_image_data[bmd_400])

    bmd_100 = np.where(numpy_segment_data == 6)
    bmd_100_cal = np.mean(numpy_image_data[bmd_100])

    attenuation_x = np.array([bmd_100_cal, bmd_400_cal, bmd_800_cal])

  return ((prof_100, prof_400, prof_800), attenuation_x, x_axis)

def lineProfile(image: Image, directory: str) -> tuple:
  global isFirstTimeOpeningDualPlots
  # Allow user to import image file, creates reader and requests rendering type
  image_original_copy = image.getImage()

  # info to obtain the center of the original image for image reslicing
  (x_min, x_max, y_min, y_max, z_min, z_max) = image_original_copy.GetDataExtent()
  (space_x, space_y, space_z) = image_original_copy.GetDataSpacing()
  (x_orig, y_orig, z_orig) = image_original_copy.GetDataOrigin()

  calibrate = True
  ((prof_100, prof_400, prof_800), attenuation_x, x_axis) = mask_calc(image=image, calibrate=calibrate)
  # BMD calibration step
  bmd_y = np.array([100, 400, 800])
  model = poly.polyfit(attenuation_x, bmd_y, 1)
  
  xlabelCalibration = "Attenuation Units (AU)"
  xlabelBMD = "Distance (mm)"
  ylabel = "BMD (mgHA/cm3)"
  filenameCalibration = image.getImageFilename()[:-4] + '_calibration_plot.png'
  filename100 = image.getImageFilename()[:-4] + '_100mgHA.png'
  filename400 = image.getImageFilename()[:-4] + '_400mgHA.png'
  filename800 = image.getImageFilename()[:-4] + '_800mgHA.png'
  filenameCombined = image.getImageFilename()[:-4] + '_combined.png'

  plotCalibration(attenuation_x, bmd_y, poly.polyval(attenuation_x, model), 
    "BMD Calibration Plot", xlabelCalibration, ylabel,
    'y = ' + '{:.2f}'.format(model[0]) + ' + {:.2f}'.format(model[1]) + 'x',
    os.path.join(image.directory, filenameCalibration))

  calibratedBMD = poly.polyval(attenuation_x, model)

  calibrationCoefficients.append(model)
  calibrationPoints.append((attenuation_x, bmd_y))

  bmd_100_prof = model[0] + (prof_100 * model[1])
  bmd_400_prof = model[0] + (prof_400 * model[1])
  bmd_800_prof = model[0] + (prof_800 * model[1])

  stats100Prof = (np.mean(bmd_100_prof), np.std(bmd_100_prof))
  stats400Prof = (np.mean(bmd_400_prof), np.std(bmd_400_prof))
  stats800Prof = (np.mean(bmd_800_prof), np.std(bmd_800_prof))

  bmd_100_limit = np.zeros_like(bmd_100_prof)
  bmd_100_limit[:] = 100
  bmd_400_limit = np.zeros_like(bmd_400_prof)
  bmd_400_limit[:] = 400
  bmd_800_limit = np.zeros_like(bmd_800_prof)
  bmd_800_limit[:] = 800
  bmd100TrueLimit = np.zeros_like(bmd_100_prof)
  bmd100TrueLimit[:] = calibratedBMD[0]
  bmd400TrueLimit = np.zeros_like(bmd_400_prof)
  bmd400TrueLimit[:] = calibratedBMD[1]
  bmd800TrueLimit = np.zeros_like(bmd_800_prof)
  bmd800TrueLimit[:] = calibratedBMD[2]

  plot(x_axis, bmd_100_prof, bmd_100_limit, bmd100TrueLimit, "Uniformity profile (100 mgHA/cm3)", xlabelBMD, ylabel, 'b', os.path.join(image.directory, filename100))
  plot(x_axis, bmd_400_prof, bmd_400_limit, bmd400TrueLimit, "Uniformity profile (400 mgHA/cm3)", xlabelBMD, ylabel, 'g', os.path.join(image.directory, filename400))
  plot(x_axis, bmd_800_prof, bmd_800_limit, bmd800TrueLimit, "Uniformity profile (800 mgHA/cm3)", xlabelBMD, ylabel, 'r', os.path.join(image.directory, filename800))

  fig, (ax1, ax2, ax3) = plt.subplots(3,1)
  fig.subplots_adjust(hspace = 0.5)
  ax1.plot(x_axis, bmd_800_prof, 'r-', x_axis, bmd_800_limit, 'm', x_axis, bmd800TrueLimit, 'k--')
  ax2.plot(x_axis, bmd_400_prof, 'g-', x_axis, bmd_400_limit, 'm', x_axis, bmd400TrueLimit, 'k--')
  ax2.set_ylabel(ylabel)
  ax3.plot(x_axis, bmd_100_prof, 'b-', x_axis, bmd_100_limit, 'm', x_axis, bmd100TrueLimit, 'k--')
  ax3.set_xlabel(xlabelBMD)
  plt.savefig(os.path.join(image.directory, filenameCombined))
  plt.clf()

  # y_i = [100, 400, 800] -> bmd_y
  # \hat{y}_i = a + bx -> calibratedBMD
  # \bar{y} = \frac{\sum_{i=1}^numCalibration y_i}{numCalibration}
  # RSE = \sqrt{\frac{\sum_{i=1}^numCalibration (y_i-\hat{y}_i)}{df}}

  y_bar = np.mean(bmd_y)
  rss = np.sum((bmd_y - calibratedBMD)**2)
  tss = np.sum((bmd_y - y_bar)**2)
  R2 = 1 - rss/tss
  rse = np.sqrt(rss / df)

  result = (model, R2, rse, calibratedBMD, attenuation_x, (stats100Prof, stats400Prof, stats800Prof))

  # Store data for plots of double line profiles
  # line profile, true limit
  key = None
  if (image.getScanner() == Scanner.CLINICAL and image.getPlane() == Plane.AXIAL and image.getPatient() == False):
    key = PlotKey.clinicalAxialPhantom
  elif (image.getScanner() == Scanner.WBCT and image.getPlane() == Plane.AXIAL and image.getPatient() == False and image.getAlg() == Alg.NEW):
    key = PlotKey.wbctAxialPhantom
  elif (image.getScanner() == Scanner.CLINICAL and image.getPlane() == Plane.CORONAL and image.getPatient() == False):
    key = PlotKey.clinicalCoronalPhantom
  elif (image.getScanner() == Scanner.WBCT and image.getPlane() == Plane.CORONAL and image.getPatient() == False and image.getAlg() == Alg.NEW):
    key = PlotKey.wbctCoronalPhantom
  elif (image.getScanner() == Scanner.WBCT and image.getPlane() == Plane.CORONAL and image.getPatient() == True and image.getAlg() == Alg.NEW):
    key = PlotKey.wbctCoronalPatient
  elif (image.getScanner() == Scanner.CLINICAL and image.getPlane() == Plane.CORONAL and image.getPatient() == True):
    key = PlotKey.clinicalCoronalPatient

  resampledString = "resampled" # Make empty string if not resampled
  if (resampledString == "resampled" and image.getScanner() == Scanner.WBCT):
    resampledImageFilename = image.getImageFilename()[:-4] + "_resampled.nii"
    resampledSegFilename = image.getImageFilename()[:-4] + "_resampled_seg.nii"
    resampledImage = ImageBuilder.createImage(resampledImageFilename, resampledSegFilename, image.directory)
    ((prof_100, prof_400, prof_800), attenuation_x, x_axis) = mask_calc(resampledImage, calibrate)
    model = poly.polyfit(attenuation_x, bmd_y, 1)
    calibratedBMD = poly.polyval(attenuation_x, model)

    bmd_100_prof = model[0] + (prof_100 * model[1])
    bmd_100_limit = np.zeros_like(bmd_100_prof)
    bmd_100_limit[:] = 100
    bmd100TrueLimit = np.zeros_like(bmd_100_prof)
    bmd100TrueLimit[:] = calibratedBMD[0]
       
    bmd_400_prof = model[0] + (prof_400 * model[1])
    bmd_400_limit = np.zeros_like(bmd_400_prof)
    bmd_400_limit[:] = 400
    bmd400TrueLimit = np.zeros_like(bmd_400_prof)
    bmd400TrueLimit[:] = calibratedBMD[1]

    bmd_800_prof = model[0] + (prof_800 * model[1])
    bmd_800_limit = np.zeros_like(bmd_800_prof)
    bmd_800_limit[:] = 800
    bmd800TrueLimit = np.zeros_like(bmd_800_prof)
    bmd800TrueLimit[:] = calibratedBMD[2]

  if (key is not None):
    dual100 = os.path.join(directory, "dualPlots100%s.txt" % resampledString)
    dual400 = os.path.join(directory, "dualPlots400%s.txt" % resampledString)
    dual800 = os.path.join(directory, "dualPlots800%s.txt" % resampledString)

    if (isFirstTimeOpeningDualPlots):
      # Create new files
      with open(dual100, "w") as f1, open(dual400, "w") as f2, open(dual800, "w") as f3:
        isFirstTimeOpeningDualPlots = False
    
    with open(dual100, "a") as f1, open(dual400, "a") as f2, open(dual800, "a") as f3:
      line0 = str(key.value)
      x_axis_str = [str(x) for x in x_axis]
      line1 = ",".join(x_axis_str)
      bmd_100_prof_str = [str(x) for x in bmd_100_prof]
      line2 = ",".join(bmd_100_prof_str)
      bmd100TrueLimitStr = [str(x) for x in bmd100TrueLimit]
      line3 = ",".join(bmd100TrueLimitStr)
      bmd_100_limit_str = [str(x) for x in bmd_100_limit]
      line4 = ",".join(bmd_100_limit_str)
      f1.write(line0 + '\n')
      f1.write(line1 + '\n')
      f1.write(line2 + '\n')
      f1.write(line3 + '\n')
      f1.write(line4 + '\n')

      line0 = str(key.value)
      x_axis_str = [str(x) for x in x_axis]
      line1 = ",".join(x_axis_str)
      bmd_400_prof_str = [str(x) for x in bmd_400_prof]
      line2 = ",".join(bmd_400_prof_str)
      bmd400TrueLimitStr = [str(x) for x in bmd400TrueLimit]
      line3 = ",".join(bmd400TrueLimitStr)
      bmd_400_limit_str = [str(x) for x in bmd_400_limit]
      line4 = ",".join(bmd_400_limit_str)
      f2.write(line0 + '\n')
      f2.write(line1 + '\n')
      f2.write(line2 + '\n')
      f2.write(line3 + '\n')
      f2.write(line4 + '\n')

      line0 = str(key.value)
      x_axis_str = [str(x) for x in x_axis]
      line1 = ",".join(x_axis_str)
      bmd_800_prof_str = [str(x) for x in bmd_800_prof]
      line2 = ",".join(bmd_800_prof_str)
      bmd800TrueLimitStr = [str(x) for x in bmd800TrueLimit]
      line3 = ",".join(bmd800TrueLimitStr)
      bmd_800_limit_str = [str(x) for x in bmd_800_limit]
      line4 = ",".join(bmd_800_limit_str)
      f3.write(line0 + '\n')
      f3.write(line1 + '\n')
      f3.write(line2 + '\n')
      f3.write(line3 + '\n')
      f3.write(line4 + '\n')

  return result

def abline(intercept, slope, color):
  """Plot a line from slope and intercept"""
  x_vals = np.arange(-50, 2300, 0.5)
  y_vals = intercept + slope * x_vals
  plt.plot(x_vals, y_vals, '-', label=f"$y = {slope:.2f}x+{intercept:.2f}$", color=color)

def plotAllCalibration(directory: str) -> None:
  cmap = plt.get_cmap('gnuplot')
  colors = [cmap(i) for i in np.linspace(0, 1, len(calibrationCoefficients))]

  plt.title('BMD Calibration Plot')
  plt.xlabel('Attenuation Units (AU)')
  plt.ylabel('BMD (mgHA/cm3)')
  plt.xlim([0, 2300])
  plt.ylim([50, 850])

  for i in range(len(calibrationCoefficients)):
    plt.plot(calibrationPoints[i][0], calibrationPoints[i][1], 'o', color=colors[i])
    abline(calibrationCoefficients[i][0], calibrationCoefficients[i][1], colors[i])

  plt.legend(loc="best")

  plt.gcf().set_size_inches(8, 6)

  plt.savefig(os.path.join(directory, "calibration_plot_combined.png"), dpi=100)
  plt.clf()

def writeMapToDataFrameToCSVWithIndex(data: dict, index: list, directory: str, filename: str="data.csv") -> None:
  df = pd.DataFrame(data=data, index=index)
  df.to_csv(os.path.join(directory, filename))

def writeMapToDataFrameToCSV(data: dict, directory: str, filename: str="data.csv") -> None:
  df = pd.DataFrame(data=data)
  df.to_csv(os.path.join(directory, filename))

def main():
  print('\n************************************************************************')
  print(f'\nRunning BMD cylinder profile analysis...\n')

  if (len(sys.argv) != 2):
    print("Usage: python uniform_seg.py <directory>")
    sys.exit(1)
  directory = sys.argv[1]

  data = {
    "scanner": [], "plane": [],
    "intercept": [], "slope": [], "R2": [], "RSE": [],
    "calibratedBMD100": [], "calibratedBMD400": [], "calibratedBMD800": [],
    "meanBMD100CalibrationCubes": [], "meanBMD400CalibrationCubes": [], "meanBMD800CalibrationCubes": [],
    "mean100Profile": [], "stDev100Profile": [], "mean400Profile": [], "stDev400Profile": [], "mean800Profile": [], "stDev800Profile": []
  }
  index = []

  firstIteration = True
  for subdir, dirs, files in os.walk(directory):
    if (firstIteration):
      index = list(dirs)
      firstIteration = False
    else:
      imageFiles = list(filter(lambda x: x.endswith(".nii") and "resampled" not in x, files))

      if ("notes.txt" not in files):
        index.remove(os.path.basename(subdir))
        continue

      imageFilename, segmentationFilename = (imageFiles[0], imageFiles[1]) if "seg" in imageFiles[1] else (imageFiles[1], imageFiles[0])

      image = ImageBuilder.createImage(imageFilename, segmentationFilename, subdir)

      (model, R2, rse, calibratedBMD, attenuation_x, (stats100Prof, stats400Prof, stats800Prof)) = lineProfile(image, directory)

      data["scanner"].append(image.getScanner().name)
      data["plane"].append(image.getPlane().name)
      data["intercept"].append(model[0])
      data["slope"].append(model[1])
      data["R2"].append(R2)
      data["RSE"].append(rse)
      data["calibratedBMD100"].append(calibratedBMD[0])
      data["calibratedBMD400"].append(calibratedBMD[1])
      data["calibratedBMD800"].append(calibratedBMD[2])
      data["meanBMD100CalibrationCubes"].append(attenuation_x[0])
      data["meanBMD400CalibrationCubes"].append(attenuation_x[1])
      data["meanBMD800CalibrationCubes"].append(attenuation_x[2])
      data["mean100Profile"].append(stats100Prof[0])
      data["stDev100Profile"].append(stats100Prof[1])
      data["mean400Profile"].append(stats400Prof[0])
      data["stDev400Profile"].append(stats400Prof[1])
      data["mean800Profile"].append(stats800Prof[0])
      data["stDev800Profile"].append(stats800Prof[1])

  # writeMapToDataFrameToCSVWithIndex(data, index, directory)

  # plotAllCalibration(directory)

  print('\n Complete')
  print('\n************************************************************************')

if __name__ == '__main__':
  main()
  end_time = time.time()
  print("--- Execution Time: %.2f seconds ---" % (end_time - start_time))
