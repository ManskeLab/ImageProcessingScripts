import SimpleITK as sitk
import sys

imageFn = sys.argv[1]
image = sitk.ImageFileReader()
image.SetFileName(imageFn)
image.LoadPrivateTagsOn()
image.Execute()

for k in image.GetMetaDataKeys():
  v = image.GetMetaData(k)
  print(f'{k} = {v}')
