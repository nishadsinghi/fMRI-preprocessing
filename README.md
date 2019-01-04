# fMRI_PreProcessing
Perform slice time correction (using linear interpolation) and spatial smoothing of functional-MRI images in R.

## Slice Time Correction
Input: 
  - Image File name
  - TR
  - Target Time
  - Slice time acquisition file
  - Output File name (say output)

Output: 
  - Creates a file output.nii.gz which contains slice time corrected image corrected to the target time.

## Spatial Smoothing
Input:
  - Image File name
  - FWHM(in mm)(say k)
  - Output File name (say output.nii.gz)

Output:
  - Creates a file output.nii.gz spatially smoothed to k mm
  
## Required Packages
`oro.nifti`
