# Scripts to do image processing and statistical analysis for the WBCT BMD study and for general use

resample.py: resample an image to certain spacing.

uniform\_seg.py: contains many classes to hold image, segmentation, and information.  
Calibration equation, scanner, plane, R2, RSE, mean/sd of profile written to data.csv.  
Calibration plots, individual profile plots, all calibration plot are created.  
Writes information about plotting two profiles together to txt file for dual\_plot.py

dual\_plot.py: plots two profiles together for both resampled and unresampled images.  
Relies on txt file generated from uniform\_seg.py.

ccl3D/cclDFS3d.py: connected component labeling algorithms for both union-find and depth-first search methods.  
ccl3D can handle searching for multiple labels.

getComponents.py: utilizes ccl3D to search for individual calibration ROIs and writes values and each slice of the cylinder profile to csv.  

The below scripts use the result of getComponents.

one\_sample\_t\_test.R: performs a 1 sample t-test.

paired\_t\_test.R: performs a paired t-test.

plotAll.py: plots all profiles on the same plot (similar to dualPlot but for all and uses the getComponents profiles instead of txt file)

LinearRegression.R: plot linear regression. Can only handle one at a time.

BlandAltmanPlot.R: plot Bland-Altman. Can only handle one at a time.

BlandAltmanPlotAll.py: plots all Bland-Altman's individually in one go.

LinearRegressionAll.py: plots all linear regression on one plot in one go.

airWaterUniformSeg.py: for water phantom to calibrate to HU.

metadata.py: to read metadata from DICOM header.
