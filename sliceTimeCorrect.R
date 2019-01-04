# Install libraries
# install.packages("oro.nifti")

# load libraries
library(oro.nifti)

# take inputs from the terminal
args = commandArgs(trailingOnly = TRUE)
inputFile = args[1]
TR = as.numeric(args[2])
target_time = as.numeric(args[3])
sliceTimeFile = args[4]
outputFile = args[5]

# Read data from the NIfTI image file
nim = readNIfTI(inputFile, reorient = FALSE)

# read the slice timings
sliceTime = read.delim2(sliceTimeFile)

# find the two instants just before and after target time
tar_time_lower = as.integer(target_time%/%TR)
tar_time_upper = tar_time_lower + 1

# take the data out of NIfTI file
data = slot(nim, ".Data")
img_lower = data[, , , tar_time_lower]
img_upper = data[, , , tar_time_upper]
d = dim(img_lower)
output = array(0, dim = d)

# linear interpolation
for(k in 1:d[3]){
  for(i in 1:d[1]){
    for(j in 1:d[2]){
      slope = as.double((img_upper[i, j, k] - img_lower[i, j, k])/TR)
      output[i, j, k] = img_upper[i, j, k] - (slope*(k - 1)*TR/d[3])
    }
  }
}

# write the output to '.nii.gz' as well as to '.txt'
nii = nifti(output)
writeNIfTI(nii, filename = outputFile)
write.table(output, file = paste(outputFile, ".txt"))
