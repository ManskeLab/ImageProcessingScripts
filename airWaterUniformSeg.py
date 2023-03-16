import os
import vtk
import tkinter as tk
from tkinter import filedialog as fd
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib.pyplot as plt
import sys # pass in arguments
from enum import Enum
import pandas as pd

calibrationCoefficients = []
calibrationPoints = []

class Scanner(Enum):
    CLINICAL = 0
    WBCT = 1
    XCT2 = 2
    NONE = 3

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

def getSpacing(typeOfScanner: Scanner, isCoronal: bool) -> float:
    spacing = 0
    if (typeOfScanner == Scanner.XCT2):
        spacing = 0.0607
    elif (typeOfScanner == Scanner.CLINICAL):
        if (isCoronal):
            spacing = 0.625
        else:
            spacing = 0.5352
    elif (typeOfScanner == Scanner.WBCT):
        spacing = 0.3
    else:
        spacing = 0

    return spacing

def mask_calc(image, segment, scanner, calibrate):
    vtk_image_data = image.GetOutput()
    vtk_segment_data = segment.GetOutput()

    numpy_image_data = vtk_to_numpy(vtk_image_data.GetPointData().GetScalars()).reshape(vtk_image_data.GetDimensions(),order = 'F')

    numpy_segment_data = vtk_to_numpy(vtk_segment_data.GetPointData().GetScalars()).reshape(vtk_segment_data.GetDimensions(),order = 'F')
    extent = numpy_segment_data.shape

    spacingCoronal = getSpacing(scanner, isCoronal=True)
    spacingAxial = getSpacing(scanner, isCoronal=False)

    profCoronal = np.zeros(extent[2])
    profAxial = np.zeros(extent[0]) if scanner == Scanner.WBCT else np.zeros(extent[1])
    x_axisCoronal = np.zeros_like(profCoronal)
    x_axisAxial = np.zeros_like(profAxial)

    slicesCoronal = extent[2]
    slicesAxial = extent[0] if scanner == Scanner.WBCT else extent[1]
    
    for i in range(slicesCoronal):
        try:
            indexCoronal = np.where(numpy_segment_data[:,:,i] == 1)
            if (np.size(indexCoronal[0]) > 0):
                imageDataCoronal = numpy_image_data[:,:,i][indexCoronal]
                profCoronal[i] = np.mean(imageDataCoronal)
                x_axisCoronal[i] = i * spacingCoronal
        except:
            continue

    for i in range(slicesAxial):
        try:
            indexAxial = np.where(numpy_segment_data[i,:,:] == 2) if scanner == Scanner.WBCT else np.where(numpy_segment_data[i,:,:] == 2)
            if (np.size(indexAxial[0]) > 0):
                imageDataAxial = numpy_image_data[i,:,:][indexAxial] if scanner == Scanner.WBCT else numpy_image_data[i,:,:][indexAxial]
                profAxial[i] = np.mean(imageDataAxial)
                x_axisAxial[i] = i * spacingAxial
        except:
            continue

    startCoronal = np.nonzero(x_axisCoronal)[0][0]
    endCoronal = np.nonzero(x_axisCoronal)[0][-1]
    startAxial = np.nonzero(x_axisAxial)[0][0]
    endAxial = np.nonzero(x_axisAxial)[0][-1]

    profCoronal = profCoronal[startCoronal: endCoronal+1]
    x_axisCoronal = x_axisCoronal[startCoronal: endCoronal+1]
    profAxial = profAxial[startAxial: endAxial+1]
    x_axisAxial = x_axisAxial[startAxial: endAxial+1]

    if (calibrate):
        # water
        waterSeg = np.where(numpy_segment_data == 4)
        waterCal = np.mean(numpy_image_data[waterSeg])

        # air
        airSeg = np.where(numpy_segment_data == 5)
        airCal = np.mean(numpy_image_data[airSeg])

        attenuation_x = [waterCal, airCal]
        hounsfield_y = [0, -1000]
        coef_calib = np.polyfit(attenuation_x, hounsfield_y, 1)

        calibrationCoefficients.append(coef_calib)
        calibrationPoints.append((attenuation_x, hounsfield_y))

        coronal_prof_line = (profCoronal * coef_calib[0]) + coef_calib[1]
        axial_prof_line = (profAxial * coef_calib[0]) + coef_calib[1]

    avgDevList = [None] * 2
    avgDevList[0] = np.mean(np.abs(coronal_prof_line-hounsfield_y[0]))
    avgDevList[1] = np.mean(np.abs(axial_prof_line-hounsfield_y[1]))

    return ((coronal_prof_line, axial_prof_line), (x_axisCoronal, x_axisAxial), (waterCal, airCal), avgDevList)
    
def lineProfile(image: str, segment: str, typeOfScanner: Scanner, directory: str) -> tuple:
    # Allow user to import image file, creates reader and requests rendering type
    # image = get_filename_GUI()

    image_original = create_reader(os.path.join(directory, image))
    image_original_copy = image_original

    # segment = get_filename_GUI()
    segment_im = create_reader(os.path.join(directory, segment))

    # info to obtain the center of the original image for image reslicing
    (x_min, x_max, y_min, y_max, z_min, z_max) = image_original_copy.GetDataExtent()  # GetDimensions
    (space_x, space_y, space_z) = image_original_copy.GetDataSpacing()    # GetDataSpacing
    (x_orig, y_orig, z_orig) = image_original_copy.GetDataOrigin()    # GetDataOrigin

    calibrate = True
    (prof_line, x_axes, cal, avgDevList) = mask_calc(image_original, segment_im, typeOfScanner, calibrate)
    # BMD calibration step
    attenuation_x, hounsfield_y = calibrationPoints[-1] # last
    coef_calib = calibrationCoefficients[-1] # last
    poly1d_calib_fn = np.poly1d(coef_calib)

    plt.plot(attenuation_x, hounsfield_y, 'ks', attenuation_x, poly1d_calib_fn(attenuation_x),'k--')
    plt.title('Hounsfield Calibration Plot')
    plt.xlabel('Attenuation Units (AU)')
    plt.ylabel('HU')
    plt.text(1000, 200, 'y = ' + ' {:.2f}'.format(coef_calib[0]) + 'x' + ' + {:.2f}'.format(coef_calib[1]))
    filenameCalibration = image[:-4] + '_calibration_plot.png'
    plt.savefig(os.path.join(directory, filenameCalibration))
    plt.clf()
    # plt.show()

    # bmd_100_prof = (prof_100 * coef_calib[0]) + coef_calib[1]
    # bmd_400_prof = (prof_400 * coef_calib[0]) + coef_calib[1]
    # bmd_800_prof = (prof_800 * coef_calib[0]) + coef_calib[1]
    coronal_prof_line, axial_prof_line = prof_line
    x_axis_coronal, x_axis_axial = x_axes

    coronal_limit = np.zeros_like(coronal_prof_line)
    coronal_limit[:] = 0
    axial_limit = np.zeros_like(axial_prof_line)
    axial_limit[:] = 0

    # bmd_100_limit = np.zeros_like(bmd_100_prof)
    # bmd_100_limit[:] = 100
    # bmd_400_limit = np.zeros_like(bmd_400_prof)
    # bmd_400_limit[:] = 400
    # bmd_800_limit = np.zeros_like(bmd_800_prof)
    # bmd_800_limit[:] = 800

    plt.plot(x_axis_coronal, coronal_prof_line, 'b-', x_axis_coronal, coronal_limit, 'm-')
    plt.xlabel('Distance (mm)')
    # plt.ylabel('BMD (mgHA/cm3)')
    # plt.title('Uniformity profile (100 mgHA/cm3)')
    plt.ylabel('HU')
    plt.title('Uniformity Profile Coronal Cylinder in Water')
    coronalWaterFilename = image[:-4] + '_coronal_water.png'
    plt.savefig(os.path.join(directory, coronalWaterFilename))
    plt.clf()
    # plt.show()

    plt.plot(x_axis_axial, axial_prof_line, 'g-', x_axis_axial, axial_limit, 'm-')
    plt.xlabel('Distance (mm)')
    # plt.ylabel('BMD (mgHA/cm3)')
    # plt.title('Uniformity profile (400 mgHA/cm3)')
    plt.ylabel('HU')
    plt.title('Uniformity Profile Axial Cylinder in Water')
    axialWaterFilename = image[:-4] + '_axial_water.png'
    plt.savefig(os.path.join(directory, axialWaterFilename))
    plt.clf()
    # plt.show()

    # plt.plot(x_axis, bmd_800_prof, 'r-', x_axis, bmd_800_limit, 'm-')
    # plt.xlabel('Distance (mm)')
    # plt.ylabel('BMD (mgHA/cm3)')
    # plt.title('Uniformity profile (800 mgHA/cm3)')
    # filename800 = image[:-4] + '_800mgHA.png'
    # plt.savefig(os.path.join(directory, filename800))
    # plt.clf()
    # plt.show()

    """
    fig, (ax1, ax2, ax3) = plt.subplots(3,1)
    fig.subplots_adjust(hspace = 0.5)

    ax1.plot(x_axis, bmd_800_prof,'r-',x_axis,bmd_800_limit,'m') #, x_axis, water_profile_limit_1, 'k--', x_axis, water_profile_limit_2, 'k--')

    ax2.plot(x_axis, bmd_400_prof,'g-',x_axis,bmd_400_limit,'m') #, x_axis, water_profile_limit_1, 'k--', x_axis, water_profile_limit_2, 'k--')
    ax2.set_ylabel('BMD (mgHA/cm3)')

    ax3.plot(x_axis, bmd_100_prof,'b-',x_axis,bmd_100_limit,'m') #, x_axis, water_profile_limit_1, 'k--', x_axis, water_profile_limit_2, 'k--')
    ax3.set_xlabel('Distance (mm)')

    filenameCombined = image[:-4] + '_combined.png'
    plt.savefig(os.path.join(directory, filenameCombined))
    plt.clf()
    # plt.show()
    """

    calibrationList = [None] * 2
    calibrationList[0] = cal[0]
    calibrationList[1] = cal[1]

    return (avgDevList, calibrationList, coef_calib)

def abline(slope, intercept, x_points, y_points, color):
    """Plot a line from slope and intercept"""
    x_vals = np.arange(-50, 2300, 0.5)
    y_vals = intercept + slope * x_vals
    print(x_points, y_points)
    plt.plot(x_vals, y_vals, '-', label=f"$y = {slope:.2f}x+{intercept:.2f}$", color=color)

def plotAllCalibration(directory: str) -> None:
    cmap = plt.get_cmap('gnuplot')
    colors = [cmap(i) for i in np.linspace(0, 1, len(calibrationCoefficients))]

    plt.title('Hounsfield Calibration Plot')
    plt.xlabel('Attenuation Units (AU)')
    plt.ylabel('HU')
    plt.xlim([0, 2300])
    plt.ylim([50, 850])

    for i in range(len(calibrationCoefficients)):
        plt.plot(calibrationPoints[i][0], calibrationPoints[i][1], 'o', color=colors[i])
        abline(calibrationCoefficients[i][0], calibrationCoefficients[i][1], calibrationPoints[i][0], calibrationPoints[i][1], colors[i])

    plt.legend(loc="best")

    plt.gcf().set_size_inches(8, 6)

    plt.savefig(os.path.join(directory, "calibration_plot_combined.png"), dpi=100)
    plt.clf()

def writeToFile(data: dict, index: list, directory: str) -> None:
    df = pd.DataFrame(data=data, index=index)
    df.to_csv(os.path.join(directory, "data.csv"))

def main():
    print('\n************************************************************************')
    print(f'\nRunning HU cylinder profile analysis...\n')

    directory = sys.argv[1]

    data = {"avgDevCoronalWater": [], "avgDevAxialWater": [], "waterCal": [], "airCal": [], "slope": [], "intercept": []}
    index = []

    i = -1
    for subdir, dirs, files in os.walk(directory):
        if (i == -1):
            index = list(dirs)
        else:
            files = list(filter(lambda x: x.endswith(".nii"), files))

            image, segmentation = (files[0], files[1]) if "seg" in files[1] else (files[1], files[0])
            typeOfScanner = Scanner.CLINICAL if "clinical" in index[i] else Scanner.WBCT if "wbct" in index[i] else Scanner.XCT2 if "xct2" in index[i] else Scanner.NONE

            avgDevList, calibrationList, model = lineProfile(image, segmentation, typeOfScanner, subdir)
            
            data["avgDevCoronalWater"].append(avgDevList[0])
            data["avgDevAxialWater"].append(avgDevList[1])
            data["waterCal"].append(calibrationList[0])
            data["airCal"].append(calibrationList[1])
            data["slope"].append(model[0])
            data["intercept"].append(model[1])

        i += 1

    writeToFile(data, index, directory)

    plotAllCalibration(directory)

    print('\n Complete')
    print('\n************************************************************************')
    sys.exit(0)

if __name__ == '__main__':
    main()
