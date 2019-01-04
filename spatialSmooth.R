# Spatial Smoothing of a SINGLE VOLUME only
# FWHM is in mm

# Libraries
#install.packages("oro.nifti")
library(oro.nifti)

# Take input from arguments
args = commandArgs(trailingOnly = TRUE)
inputFile = args[1]
FWHM = as.numeric(args[2])
outputFile = args[3]

sigma = FWHM/2.355

# Load the image to be processed
nii = readNIfTI(inputFile, reorient = FALSE)
img_ = slot(nii, ".Data")
img = img_[, , , 1]
pd = pixdim(nii)
pd = pd[2:4]
dim = dim(img)

print("Everything loaded successfully!")

# Kernel Function
kernel = function(x, y, z){
  coeff = (2 * 3.14 * sigma * sigma)**(-3/2)
  return(coeff * exp(-1 * (1/(2 * sigma)**2)) * (x**2 + y**2 + z**2))
}

# Function for convolution
convolution = function(img, x, y, z){
  out = 0
  for(i in 1:x){
    for(j in 1:y){
      for(k in 1:z){
        out = out + img[i, j, k] * kernel(x-i*pd[1], y-j*pd[2], z-k*pd[3])
      }
    }
  }
  return(out)
}

output_img = array(0.0, dim)


for(i in 1:dim[1]){
  for(j in 1:dim[2]){
    for(k in 1:dim[3]){
      output_img[i, j, k] = convolution(img, i, j, k)
      print(paste(i, j, k, sep = " "))
    }
  }
}

nii_out = niftiarr(nii, output_img)
writeNIfTI(nii_out, outputFile)