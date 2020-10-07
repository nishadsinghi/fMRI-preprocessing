# fMRI_PreProcessing
Perform slice time correction (using linear interpolation), temporal filtering and spatial smoothing of functional-MRI images using python.
Input file must be in NIfTI format.

## Usage:
`python preproc.py --inputFile --outputFile [optional: --xCoord --yCoord --zCoord --targetTime --lowCutoff --highCutoff --fwhm]`

  - xCoord, yCoord, zCoord: used to print the timeseries corresponding to this voxel.
  - targetTime: returns the brain scan for the given time instant using slice time correction.
  - lowCutoff, highCutoff: used to specify the band-pass filter for temporal filtering of fMRI data
  - fwhm: FWHM for spatial smoothing using Gaussian kernel.
