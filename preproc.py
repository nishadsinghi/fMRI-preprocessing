import nibabel as nib
import numpy as np
from scipy.fftpack import fft, ifft
import argparse
import sys
import nibabel as nb


def generateTimeSeries(fMRIData, xCoordinate, yCoordinate, zCoordinate, outputTextFileName):
    outputTextFile = open(outputTextFileName, "a")
    timeSeries = fMRIData[xCoordinate, yCoordinate, zCoordinate, :]
    timeSeriesAsStr = [str(_) for _ in timeSeries]
    joinedString = " ".join(timeSeriesAsStr)
    outputTextFile.write(joinedString + "\n")


########################################################################################################################

def sliceTimeCorrectionForOneTimeSeries(y, sliceAcquisitionTime, TR, targetTime):
    if targetTime - sliceAcquisitionTime > 0:
        y1 = y[:-1]
        y2 = y[1:]
        coefficient = (targetTime - sliceAcquisitionTime) / TR
    elif targetTime - sliceAcquisitionTime < 0:
        initialValue = [2*y[0] - y[1]]
        y1 = np.concatenate((initialValue, y[:-2]))
        y2 = y[:-1]
        coefficient = (targetTime - sliceAcquisitionTime + TR) / TR

    yCorrectedWithoutLastVolume = y1 + coefficient * (y2 - y1)
    yCorrected = np.concatenate((yCorrectedWithoutLastVolume, [y[-1]]))

    return yCorrected


def sliceTimeCorrection(fMRIData, TR, targetTime, sliceTimeAcquisitionFile, outputTextFileName):
    outputTextFile = open(outputTextFileName, "a")

    TR = TR*1000
    if targetTime > TR or targetTime < 0:
        outputTextFile.write("SLICE TIME CORRECTION FAILURE")

        return fMRIData

    else:
        sliceAcquisitionTimes = np.loadtxt(sliceTimeAcquisitionFile)

        sliceTimeCorrectedImage = np.zeros(np.shape(fMRIData))
        for sliceIndex in range(np.shape(fMRIData)[2]):
            if sliceAcquisitionTimes[sliceIndex] > TR or sliceAcquisitionTimes[sliceIndex] < 0:
                outputTextFile.write("SLICE TIME CORRECTION FAILURE")

                return fMRIData
            for xCoordinate in range(np.shape(fMRIData)[0]):
                for yCoordinate in range(np.shape(fMRIData)[1]):
                    sliceTimeCorrectedImage[xCoordinate][yCoordinate][sliceIndex] = sliceTimeCorrectionForOneTimeSeries(fMRIData[xCoordinate][yCoordinate][sliceIndex], sliceAcquisitionTimes[sliceIndex], TR, targetTime)

        outputTextFile.write("SLICE TIME CORRECTION SUCCESS")
        outputTextFile.write("\n")
        return sliceTimeCorrectedImage


########################################################################################################################

# temporal filtering
class TemporalFilterForVoxel:
    def __init__(self, TR, getVoxelTimeSeries):
        self.TR = TR
        self.getVoxelTimeSeries = getVoxelTimeSeries

    def __call__(self, voxelCoordinates, highFreqCutoff, lowFreqCutoff):
        voxelTimeSeries = self.getVoxelTimeSeries(*voxelCoordinates)
        numSamplePoints = np.shape(voxelTimeSeries)[0]
        frequencies = np.concatenate((np.linspace(0.0, 1.0/(2.0*self.TR), numSamplePoints//2),
                                     np.linspace(-1.0/(2.0*self.TR), 0, numSamplePoints//2)))
        YFourier = fft(voxelTimeSeries)
        YFourier[np.abs(frequencies) > highFreqCutoff] = 0
        YFourier[np.abs(frequencies) < lowFreqCutoff] = 0
        filteredTimeSeries = ifft(YFourier)

        return filteredTimeSeries


def temporalFiltering(fMRIData, TR, high, low):
    highFreqCutoff = 1/low
    lowFreqCutoff = 1/high
    TR = TR

    if highFreqCutoff < lowFreqCutoff:
        temp = highFreqCutoff
        highFreqCutoff = lowFreqCutoff
        lowFreqCutoff = temp

    getVoxelTimeSeries = lambda xCoordinate, yCoordinate, zCoordinate: fMRIData[xCoordinate][yCoordinate][zCoordinate]
    temporalFilterForVoxel = TemporalFilterForVoxel(TR, getVoxelTimeSeries)

    temporallyFilteredData = np.zeros(np.shape(fMRIData))
    for xCoordinate in range(np.shape(fMRIData)[0]):
        for yCoordinate in range(np.shape(fMRIData)[1]):
            for zCoordinate in range(np.shape(fMRIData)[2]):
                temporallyFilteredData[xCoordinate][yCoordinate][zCoordinate][:] = \
                    temporalFilterForVoxel((xCoordinate, yCoordinate, zCoordinate), highFreqCutoff, lowFreqCutoff)

    return temporallyFilteredData


########################################################################################################################

kernelSize = 5

def applyKernel(sequence, kernel):
    kernelLength = np.size(kernel)
    sequenceLength = np.size(sequence)
    paddingLength = kernelLength // 2
    paddedSequence = np.concatenate((np.zeros(paddingLength), sequence, np.zeros(paddingLength)))

    output = np.zeros(sequenceLength)
    for startingPointer in range(0, sequenceLength):
        output[startingPointer] = np.sum(np.multiply(kernel, paddedSequence[startingPointer:startingPointer+kernelLength]))

    return output


def prepareGaussianKernel(kernelSize, fwhm, voxelSide):
    x = np.arange(int(-1*kernelSize/2), int(kernelSize/2)+1)
    sigma = fwhm / (np.sqrt(8 * np.log(2)) * voxelSide)
    y = np.exp(-(x) ** 2 / (2 * (sigma ** 2)))
    kernel = y/sum(y)

    return kernel


def spatialSmoothing(fMRIData, fwhm, voxelDimensions):
    outputData = np.zeros(np.shape(fMRIData))
    for timeStep in range(np.shape(fMRIData)[-1]):
        # apply kernel along x
        kernelX = prepareGaussianKernel(kernelSize, fwhm, voxelDimensions[0])
        outputX = np.zeros(np.shape(fMRIData)[:-1])
        for yCoordinate in range(np.shape(fMRIData)[1]):
            for zCoordinate in range(np.shape(fMRIData)[2]):
                outputX[:, yCoordinate, zCoordinate] = applyKernel(fMRIData[:, yCoordinate, zCoordinate, timeStep], kernelX)

        # apply kernel along y
        outputY = np.zeros(np.shape(fMRIData)[:-1])
        kernelY = prepareGaussianKernel(kernelSize, fwhm, voxelDimensions[1])
        for xCoordinate in range(np.shape(fMRIData)[0]):
            for zCoordinate in range(np.shape(fMRIData)[2]):
                outputX[xCoordinate, :, zCoordinate] = applyKernel(fMRIData[xCoordinate, :, zCoordinate, timeStep], kernelY)

        # apply kernel along z
        outputZ = np.zeros(np.shape(fMRIData)[:-1])
        kernelZ = prepareGaussianKernel(kernelSize, fwhm, voxelDimensions[2])
        for xCoordinate in range(np.shape(fMRIData)[0]):
            for yCoordinate in range(np.shape(fMRIData)[1]):
                outputX[xCoordinate, yCoordinate, :] = applyKernel(fMRIData[xCoordinate, yCoordinate, :, timeStep], kernelZ)

        meanOutputForTime = (outputX + outputY + outputZ)/3
        outputData[...,timeStep] = meanOutputForTime

    return fMRIData


########################################################################################################################

parser = argparse.ArgumentParser(description="preprocessing tool for fMRI images")
parser.add_argument("-i", "--inputFile", required=True)
parser.add_argument('-o', '--outputFile', required=True)
parser.add_argument('-ts', '--listTimeSeriesGeneration', nargs='+', default=[0, 0, 0])
parser.add_argument('-tc', '--listSliceTimeCorrection', nargs='+', default=[0, 0])
parser.add_argument('-tf', '--listTemporalFiltering', nargs='+', default=[0, 0])
parser.add_argument('-sm', '--fwhm', default=0)
args = parser.parse_args()

inputFile = str(args.inputFile)
outputFile = str(args.outputFile)

xCoordinate = int(args.listTimeSeriesGeneration[0])
yCoordinate = int(args.listTimeSeriesGeneration[1])
zCoordinate = int(args.listTimeSeriesGeneration[2])

targetTime = float(args.listSliceTimeCorrection[0])
sliceTimeAcquisitionFile = args.listSliceTimeCorrection[1]

high = float(args.listTemporalFiltering[0])
low = float(args.listTemporalFiltering[1])

fwhm = float(args.fwhm)

def main():
    img = nib.load(inputFile)
    fMRIData = img.get_data()
    header = img.header
    TR = header.get_zooms()[3]
    voxelDimensions = header.get_zooms()[:3]

    outputTextFileName = outputFile + '.txt'
    saveNII = False
    numSysArgs = len(sys.argv)
    for argIndex in range(numSysArgs):
        arg = sys.argv[argIndex]

        if arg == '-ts' or arg == '--listTimeSeriesGeneration':
            print("Time series generation")
            generateTimeSeries(fMRIData, xCoordinate, yCoordinate, zCoordinate, outputTextFileName)

        elif arg == '-tc' or arg == '--listSliceTimeCorrection':
            print("Slice time correction")
            fMRIData = sliceTimeCorrection(fMRIData, TR, targetTime, sliceTimeAcquisitionFile, outputTextFileName)
            saveNII = True

        elif arg == '-tf' or arg == '--listTemporalFiltering':
            print("temporal filtering")
            fMRIData = temporalFiltering(fMRIData, TR, high, low)
            saveNII = True

        elif arg == '-sm' or arg == '--fwhm':
            print("spatial smoothing")
            fMRIData = spatialSmoothing(fMRIData, fwhm, voxelDimensions)
            saveNII = True


    if saveNII:
        outputImg = nb.Nifti1Image(fMRIData, np.eye(4), header=header)
        outputNIIFileName = outputFile + '.nii.gz'
        nb.save(outputImg, outputNIIFileName)

if __name__ == '__main__':
    main()