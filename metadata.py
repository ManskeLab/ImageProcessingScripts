import SimpleITK as sitk
import sys

imageFn = sys.argv[1]
image = sitk.ReadImage(imageFn)

for k in image.GetMetaDataKeys():
  v = image.GetMetaData(k)
  print(f'{k} = {v}')
