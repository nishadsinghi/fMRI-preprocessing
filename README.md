# fMRI_PreProcessing
Perform slice time correction (using linear interpolation) and spatial smoothing of functional-MRI images in R.
Input file must be in NIfTI format.

## Slice Time Correction
Input: 
  - Image File name
  - TR
  - Target Time
  - Slice time acquisition file
  - Output File name (say output)

Run: 
`Rscript sliceTimeCorrect.R inputFile <FWHM> outputFile`

Output: 
  - Creates a file output.nii.gz which contains slice time corrected image corrected to the target time.

## Spatial Smoothing
Input:
  - Image File name
  - FWHM(in mm)(say k)
  - Output File name (say output.nii.gz)

Run:
`Rscript spatialSmooth.R inputFile <TR> <Target Time> sliceTimeFile outputFileName`

Output:
  - Creates a file output.nii.gz spatially smoothed to k mm
  
## Required Packages
`oro.nifti`
