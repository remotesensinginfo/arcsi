"""
Module that contains the ARCSIAbstractSensor class.
"""
############################################################################
#  arcsisensor.py
#
#  Copyright 2013 ARCSI.
#
#  ARCSI: 'Atmospheric and Radiometric Correction of Satellite Imagery'
#
#  ARCSI is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ARCSI is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with ARCSI.  If not, see <http://www.gnu.org/licenses/>.
#
#
# Purpose:  An abstract class which provides the base class for each 
#           sensor supported by the ARCSI system.
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 05/07/2013
# Version: 1.0
#
# History:
# Version 1.0 - Created.
#
############################################################################

# import abstract base class stuff
from abc import ABCMeta, abstractmethod
# Import the datetime module
import datetime
# Import the ARCSI exception class
from .arcsiexception import ARCSIException
# Import the python collections library
import collections
# Import python math module
import math
# Import the rsgislib module
import rsgislib
# Import the rsgislib imagecalc module
import rsgislib.imagecalc
# Import the RSGISLib segmentation Module
import rsgislib.segmentation
# Import the RSGISLib Raster GIS Module
import rsgislib.rastergis
# Import the RSGISLib Image Utils Module
import rsgislib.imageutils
#Import the OSGEO GDAL module
import osgeo.gdal as gdal
# Import the ARCSI utilities class
from .arcsiutils import ARCSIUtils
# Import OS path module for manipulating the file system 
import os.path
# Import the numpy module
import numpy
# Import the RIOS RAT library
from rios import rat
# Import the RIOS image reader
from rios.imagereader import ImageReader
# Import the RIOS image writer
from rios.imagewriter import ImageWriter
# Import the RBF interpolator from scipy
import scipy.interpolate.rbf

class ARCSIAbstractSensor (object):
    """
    An abstract class which represents a sensor and allows
    the various opperations required to be applied and standard
    variables (e.g., acqusiation date) stored and retrieved.
    """
    __metaclass__ = ABCMeta
        
    def __init__(self):
        self.sensor = "NA"
        self.acquisitionTime = datetime.datetime.today()
        self.latTL = 0.0
        self.lonTL = 0.0
        self.latTR = 0.0
        self.lonTR = 0.0
        self.latBL = 0.0
        self.lonBL = 0.0
        self.latBR = 0.0
        self.lonBR = 0.0
        self.latCentre = 0.0
        self.lonCentre = 0.0
        self.xTL = 0.0
        self.yTL = 0.0
        self.xTR = 0.0
        self.yTR = 0.0
        self.xBL = 0.0
        self.yBL = 0.0
        self.xBR = 0.0
        self.yBR = 0.0
        self.xCentre = 0.0
        self.yCentre = 0.0
        self.inWKT = ""
        self.solarZenith = 0.0
        self.solarAzimuth = 0.0
        self.senorZenith = 0.0
        self.senorAzimuth = 0.0
        self.epsgCodes = dict()
        self.epsgCodes["WGS84UTM1N"] = 32601
        self.epsgCodes["WGS84UTM2N"] = 32602
        self.epsgCodes["WGS84UTM3N"] = 32603
        self.epsgCodes["WGS84UTM4N"] = 32604
        self.epsgCodes["WGS84UTM5N"] = 32605
        self.epsgCodes["WGS84UTM6N"] = 32606
        self.epsgCodes["WGS84UTM7N"] = 32607
        self.epsgCodes["WGS84UTM8N"] = 32608
        self.epsgCodes["WGS84UTM9N"] = 32609
        self.epsgCodes["WGS84UTM10N"] = 32610
        self.epsgCodes["WGS84UTM11N"] = 32611
        self.epsgCodes["WGS84UTM12N"] = 32612
        self.epsgCodes["WGS84UTM13N"] = 32613
        self.epsgCodes["WGS84UTM14N"] = 32614
        self.epsgCodes["WGS84UTM15N"] = 32615
        self.epsgCodes["WGS84UTM16N"] = 32616
        self.epsgCodes["WGS84UTM17N"] = 32617
        self.epsgCodes["WGS84UTM18N"] = 32618
        self.epsgCodes["WGS84UTM19N"] = 32619
        self.epsgCodes["WGS84UTM20N"] = 32620
        self.epsgCodes["WGS84UTM21N"] = 32621
        self.epsgCodes["WGS84UTM22N"] = 32622
        self.epsgCodes["WGS84UTM23N"] = 32623
        self.epsgCodes["WGS84UTM24N"] = 32624
        self.epsgCodes["WGS84UTM25N"] = 32625
        self.epsgCodes["WGS84UTM26N"] = 32626
        self.epsgCodes["WGS84UTM27N"] = 32627
        self.epsgCodes["WGS84UTM28N"] = 32628
        self.epsgCodes["WGS84UTM29N"] = 32629
        self.epsgCodes["WGS84UTM30N"] = 32630
        self.epsgCodes["WGS84UTM31N"] = 32631
        self.epsgCodes["WGS84UTM32N"] = 32632
        self.epsgCodes["WGS84UTM33N"] = 32633
        self.epsgCodes["WGS84UTM34N"] = 32634
        self.epsgCodes["WGS84UTM35N"] = 32635
        self.epsgCodes["WGS84UTM36N"] = 32636
        self.epsgCodes["WGS84UTM37N"] = 32637
        self.epsgCodes["WGS84UTM38N"] = 32638
        self.epsgCodes["WGS84UTM39N"] = 32639
        self.epsgCodes["WGS84UTM40N"] = 32640
        self.epsgCodes["WGS84UTM41N"] = 32641
        self.epsgCodes["WGS84UTM42N"] = 32642
        self.epsgCodes["WGS84UTM43N"] = 32643
        self.epsgCodes["WGS84UTM44N"] = 32644
        self.epsgCodes["WGS84UTM45N"] = 32645
        self.epsgCodes["WGS84UTM46N"] = 32646
        self.epsgCodes["WGS84UTM47N"] = 32647
        self.epsgCodes["WGS84UTM48N"] = 32648
        self.epsgCodes["WGS84UTM49N"] = 32649
        self.epsgCodes["WGS84UTM50N"] = 32650
        self.epsgCodes["WGS84UTM51N"] = 32651
        self.epsgCodes["WGS84UTM52N"] = 32652
        self.epsgCodes["WGS84UTM53N"] = 32653
        self.epsgCodes["WGS84UTM54N"] = 32654
        self.epsgCodes["WGS84UTM55N"] = 32655
        self.epsgCodes["WGS84UTM56N"] = 32656
        self.epsgCodes["WGS84UTM57N"] = 32657
        self.epsgCodes["WGS84UTM58N"] = 32658
        self.epsgCodes["WGS84UTM59N"] = 32659
        self.epsgCodes["WGS84UTM60N"] = 32660
        self.epsgCodes["WGS84UTM1S"] = 32701
        self.epsgCodes["WGS84UTM2S"] = 32702
        self.epsgCodes["WGS84UTM3S"] = 32703
        self.epsgCodes["WGS84UTM4S"] = 32704
        self.epsgCodes["WGS84UTM5S"] = 32705
        self.epsgCodes["WGS84UTM6S"] = 32706
        self.epsgCodes["WGS84UTM7S"] = 32707
        self.epsgCodes["WGS84UTM8S"] = 32708
        self.epsgCodes["WGS84UTM9S"] = 32709
        self.epsgCodes["WGS84UTM10S"] = 32710
        self.epsgCodes["WGS84UTM11S"] = 32711
        self.epsgCodes["WGS84UTM12S"] = 32712
        self.epsgCodes["WGS84UTM13S"] = 32713
        self.epsgCodes["WGS84UTM14S"] = 32714
        self.epsgCodes["WGS84UTM15S"] = 32715
        self.epsgCodes["WGS84UTM16S"] = 32716
        self.epsgCodes["WGS84UTM17S"] = 32717
        self.epsgCodes["WGS84UTM18S"] = 32718
        self.epsgCodes["WGS84UTM19S"] = 32719
        self.epsgCodes["WGS84UTM20S"] = 32720
        self.epsgCodes["WGS84UTM21S"] = 32721
        self.epsgCodes["WGS84UTM22S"] = 32722
        self.epsgCodes["WGS84UTM23S"] = 32723
        self.epsgCodes["WGS84UTM24S"] = 32724
        self.epsgCodes["WGS84UTM25S"] = 32725
        self.epsgCodes["WGS84UTM26S"] = 32726
        self.epsgCodes["WGS84UTM27S"] = 32727
        self.epsgCodes["WGS84UTM28S"] = 32728
        self.epsgCodes["WGS84UTM29S"] = 32729
        self.epsgCodes["WGS84UTM30S"] = 32730
        self.epsgCodes["WGS84UTM31S"] = 32731
        self.epsgCodes["WGS84UTM32S"] = 32732
        self.epsgCodes["WGS84UTM33S"] = 32733
        self.epsgCodes["WGS84UTM34S"] = 32734
        self.epsgCodes["WGS84UTM35S"] = 32735
        self.epsgCodes["WGS84UTM36S"] = 32736
        self.epsgCodes["WGS84UTM37S"] = 32737
        self.epsgCodes["WGS84UTM38S"] = 32738
        self.epsgCodes["WGS84UTM39S"] = 32739
        self.epsgCodes["WGS84UTM40S"] = 32740
        self.epsgCodes["WGS84UTM41S"] = 32741
        self.epsgCodes["WGS84UTM42S"] = 32742
        self.epsgCodes["WGS84UTM43S"] = 32743
        self.epsgCodes["WGS84UTM44S"] = 32744
        self.epsgCodes["WGS84UTM45S"] = 32745
        self.epsgCodes["WGS84UTM46S"] = 32746
        self.epsgCodes["WGS84UTM47S"] = 32747
        self.epsgCodes["WGS84UTM48S"] = 32748
        self.epsgCodes["WGS84UTM49S"] = 32749
        self.epsgCodes["WGS84UTM50S"] = 32750
        self.epsgCodes["WGS84UTM51S"] = 32751
        self.epsgCodes["WGS84UTM52S"] = 32752
        self.epsgCodes["WGS84UTM53S"] = 32753
        self.epsgCodes["WGS84UTM54S"] = 32754
        self.epsgCodes["WGS84UTM55S"] = 32755
        self.epsgCodes["WGS84UTM56S"] = 32756
        self.epsgCodes["WGS84UTM57S"] = 32757
        self.epsgCodes["WGS84UTM58S"] = 32758
        self.epsgCodes["WGS84UTM59S"] = 32759
        self.epsgCodes["WGS84UTM60S"] = 32760
        
        
    
    @abstractmethod
    def extractHeaderParameters(self, inputHeader, wktStr): pass
    
    def defaultGenBaseOutFileName(self):
        """
        A function to generate a generic standard file
        base name which will be sensible.
        
        It is expected that individual sensors may override this function.
        """
        date = self.acquisitionTime.strftime("%Y%m%d")
        pos = "lat" + str(round(self.latCentre,)).replace('.', '').replace('-', '') + "lon" + str(round(self.lonCentre,2)).replace('.', '').replace('-', '')
        outname = self.sensor + "_" + date + "_" + pos
        return outname
        
    def generateOutputBaseName(self):
        """
        Provides a default implementation.
        """
        return self.defaultGenBaseOutFileName()
    
    def maskInputImages(self):
        return False
        
    def hasThermal(self):
        return False
    
    @abstractmethod
    def applyImageDataMask(self, inputHeader, outputPath, outputMaskName, outputImgName, outFormat): pass
    
    @abstractmethod
    def convertImageToRadiance(self, outputPath, outputReflName, outputThermalName, outFormat): pass
    
    @abstractmethod
    def generateImageSaturationMask(self, outputPath, outputName, outFormat): pass
    
    @abstractmethod
    def convertThermalToBrightness(self, inputRadImage, outputPath, outputName, outFormat): pass
    
    @abstractmethod
    def convertImageToTOARefl(self, inputRadImage, outputPath, outputName, outFormat): pass
    
    @abstractmethod
    def generateCloudMask(self, inputReflImage, inputSatImage, inputThermalImage, outputPath, outputName, outFormat, tmpPath): pass
    
    @abstractmethod
    def convertImageToSurfaceReflSglParam(self, inputRadImage, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF): pass

    @abstractmethod
    def calc6SCoefficients(self, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal): pass
    
    def buildElevation6SCoeffLUT(self, aeroProfile, atmosProfile, grdRefl, aotVal, useBRDF, surfaceAltitudeMin, surfaceAltitudeMax):
        elevLUTFeat = collections.namedtuple('ElevLUTFeat', ['Elev', 'Coeffs'])
        lut = list()
        elevRange = (surfaceAltitudeMax - surfaceAltitudeMin) / 100
        numElevSteps = int(math.ceil(elevRange) + 1)
        elevVal = surfaceAltitudeMin
        for i in range(numElevSteps):
            print("Building LUT Elevation ", elevVal)
            lut.append(elevLUTFeat(Elev=elevVal, Coeffs=self.calc6SCoefficients(aeroProfile, atmosProfile, grdRefl, (float(elevVal)/1000), aotVal, useBRDF)))
            elevVal = elevVal + 100
        return lut
    
    @abstractmethod
    def convertImageToSurfaceReflDEMElevLUT(self, inputRadImage, inputDEMFile, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, aotVal, useBRDF, surfaceAltitudeMin, surfaceAltitudeMax): pass
    
    def buildElevationAOT6SCoeffLUT(self, aeroProfile, atmosProfile, grdRefl, useBRDF, surfaceAltitudeMin, surfaceAltitudeMax, aotMin, aotMax):
        elevLUTFeat = collections.namedtuple('ElevLUTFeat', ['Elev', 'Coeffs'])
        aotLUTFeat = collections.namedtuple('AOTLUTFeat', ['AOT', 'Coeffs'])
        lut = list()
        elevRange = (surfaceAltitudeMax - surfaceAltitudeMin) / 100
        numElevSteps = int(math.ceil(elevRange) + 1)
        elevVal = surfaceAltitudeMin
        
        aotRange = (aotMax - aotMin) / 0.05
        numAOTSteps = int(math.ceil(aotRange) + 1) + 1
        aotVal = aotMin
        
        for i in range(numElevSteps):
            print("Building LUT Elevation ", elevVal)
            aotVal = aotMin
            aotCoeffLUT = list()
            for j in range(numAOTSteps):
                aotCoeffLUT.append(aotLUTFeat(AOT=aotVal,  Coeffs=self.calc6SCoefficients(aeroProfile, atmosProfile, grdRefl, (float(elevVal)/1000), aotVal, useBRDF)))
                aotVal = aotVal + 0.05
            lut.append(elevLUTFeat(Elev=elevVal, Coeffs=aotCoeffLUT))
            elevVal = elevVal + 100
        return lut

    @abstractmethod
    def convertImageToSurfaceReflAOTDEMElevLUT(self, inputRadImage, inputDEMFile, inputAOTImage, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, useBRDF, surfaceAltitudeMin, surfaceAltitudeMax, aotMin, aotMax): pass
    
    def calcDarkTargetOffsetsForBand(self, inputTOAImage, offsetImage, band, outFormat, histBinWidth, minObjSize, darkPxlPercentile, tmpDarkPxlsImg, tmpDarkPxlsClumpsImg, tmpDarkPxlsClumpsRMSmallImg, tmpDarkObjsImg):
        print("Band: ", band)
        bandHist = rsgislib.imagecalc.getHistogram(inputTOAImage, band, histBinWidth, False, 1, 10000)
        #print(bandHist)
        sumPxls = numpy.sum(bandHist[0])
        #print("Total Num Pixels: ", sumPxls)
        numPxlThreshold = sumPxls * darkPxlPercentile #0.001 # 1th percentile
        #print("Number of pixels = ", numPxlThreshold)
        
        pxlThreshold = 0
        pxlCount = 0
        for bin in bandHist[0]:
            pxlCount = pxlCount + bin
            if pxlCount < numPxlThreshold:
                pxlThreshold = pxlThreshold + histBinWidth
            else:
                break
        print("Image Band Threshold (For Dark Pixels) = ", pxlThreshold)
        
        dataType = rsgislib.TYPE_8UINT
        expression = str('(b1!=0)&&(b1<') + str(pxlThreshold) + str(')?1:0')
        bandDefns = []
        bandDefns.append(rsgislib.imagecalc.BandDefn('b1', inputTOAImage, band))
        rsgislib.imagecalc.bandMath(tmpDarkPxlsImg, expression, outFormat, dataType, bandDefns)
        rsgislib.segmentation.clump(tmpDarkPxlsImg, tmpDarkPxlsClumpsImg, outFormat, False, 0.0)
        rsgislib.rastergis.populateStats(tmpDarkPxlsClumpsImg, True, False)
        rsgislib.segmentation.rmSmallClumps(tmpDarkPxlsClumpsImg, tmpDarkPxlsClumpsRMSmallImg, minObjSize, outFormat)
        rsgislib.segmentation.relabelClumps(tmpDarkPxlsClumpsRMSmallImg, tmpDarkObjsImg, outFormat, False)
        rsgislib.rastergis.populateStats(tmpDarkObjsImg, True, False)
        stats2CalcTOA = list()
        stats2CalcTOA.append(rsgislib.rastergis.BandAttStats(band=(band), minField="MinTOARefl", meanField="MeanTOARefl"))
        rsgislib.rastergis.populateRATWithStats(inputTOAImage, tmpDarkObjsImg, stats2CalcTOA)
        
        ratDS = gdal.Open(tmpDarkObjsImg, gdal.GA_Update)
        Histogram = rat.readColumn(ratDS, "Histogram")
        selected = Histogram * 2
        selected[...] = 1
        selected[0] = 0
        rat.writeColumn(ratDS, "Selected", selected)
        ratDS = None
        
        rsgislib.rastergis.spatialLocation(tmpDarkObjsImg, "Eastings", "Northings")
        rsgislib.rastergis.selectClumpsOnGrid(tmpDarkObjsImg, "Selected", "SelectedGrid", "Eastings", "Northings", "MeanTOARefl", "min", 20, 20)
        
        ratDS = gdal.Open(tmpDarkObjsImg, gdal.GA_Update)
        Eastings = rat.readColumn(ratDS, "Eastings")
        Northings = rat.readColumn(ratDS, "Northings")
        MinTOARefl = rat.readColumn(ratDS, "MinTOARefl")
        SelectedGrid = rat.readColumn(ratDS, "SelectedGrid")
        ratDS = None
        
        Eastings = Eastings[SelectedGrid!=0]
        Northings = Northings[SelectedGrid!=0]
        MinTOARefl = MinTOARefl[SelectedGrid!=0]
        
        interpSmoothing = 10.0
        self.interpolateImageFromPointData(inputTOAImage, Eastings, Northings, MinTOARefl, offsetImage, outFormat, interpSmoothing)
        
        #rsgislib.rastergis.interpolateClumpValues2Image(tmpDarkObjsImg, "SelectedGrid", "Eastings", "Northings", "idwall", "MinTOARefl",  offsetImage, outFormat, rsgislib.TYPE_32FLOAT)
    
    def findPerBandDarkTargetsOffsets(self, inputTOAImage, numBands, outputPath, outputName, outFormat, tmpPath, minObjSize, darkPxlPercentile):
        try:
            arcsiUtils = ARCSIUtils()
            tmpBaseName = os.path.splitext(outputName)[0]
            binWidth = 1
                        
            bandDarkTargetOffsetImages = list()
            imgExtension = arcsiUtils.getFileExtension(outFormat)
            tmpDarkPxlsImg = os.path.join(tmpPath, tmpBaseName + "_darkpxls" + imgExtension)
            tmpDarkPxlsClumpsImg = os.path.join(tmpPath, tmpBaseName + "_darkclumps" + imgExtension)
            tmpDarkPxlsClumpsRMSmallImg = os.path.join(tmpPath, tmpBaseName + "_darkclumpsrmsmall" + imgExtension)
            tmpDarkObjsImg = os.path.join(tmpPath, tmpBaseName+"_darkobjs"+imgExtension)
            
            for band in range(numBands):
                offsetImage = os.path.join(tmpPath, tmpBaseName+"_darktargetoffs_b"+str(band+1)+imgExtension)
                self.calcDarkTargetOffsetsForBand(inputTOAImage, offsetImage, band+1, outFormat, binWidth, minObjSize, darkPxlPercentile, tmpDarkPxlsImg, tmpDarkPxlsClumpsImg, tmpDarkPxlsClumpsRMSmallImg, tmpDarkObjsImg)
                bandDarkTargetOffsetImages.append(offsetImage)
            
            outputImage = os.path.join(outputPath, tmpBaseName + "_dosuboffs" + imgExtension)
            print(outputImage)
            rsgislib.imageutils.stackImageBands(bandDarkTargetOffsetImages, None, outputImage, None, 0, outFormat, rsgislib.TYPE_32FLOAT)
            
            gdalDriver = gdal.GetDriverByName(outFormat)
            gdalDriver.Delete(tmpDarkPxlsImg)
            gdalDriver.Delete(tmpDarkPxlsClumpsImg)
            gdalDriver.Delete(tmpDarkPxlsClumpsRMSmallImg)
            gdalDriver.Delete(tmpDarkObjsImg)
            for image in bandDarkTargetOffsetImages:
                gdalDriver.Delete(image)
            
            return outputImage
            
        except Exception as e:
            raise e
    
    def performDOSOnSingleBand(self, inputTOAImage, band, outputPath, tmpBaseName, bandName, outFormat, tmpPath, minObjSize, darkPxlPercentile):
        try:
            arcsiUtils = ARCSIUtils()
            binWidth = 1
                        
            imgExtension = arcsiUtils.getFileExtension(outFormat)
            tmpDarkPxlsImg = os.path.join(tmpPath, tmpBaseName + "_darkpxls" + imgExtension)
            tmpDarkPxlsClumpsImg = os.path.join(tmpPath, tmpBaseName + "_darkclumps" + imgExtension)
            tmpDarkPxlsClumpsRMSmallImg = os.path.join(tmpPath, tmpBaseName + "_darkclumpsrmsmall" + imgExtension)
            tmpDarkObjsImg = os.path.join(tmpPath, tmpBaseName+"_darkobjs"+imgExtension)
            
            offsetImage = os.path.join(outputPath, tmpBaseName+"_darktargetoffs_b"+str(band)+imgExtension)
            self.calcDarkTargetOffsetsForBand(inputTOAImage, offsetImage, band, outFormat, binWidth, minObjSize, darkPxlPercentile, tmpDarkPxlsImg, tmpDarkPxlsClumpsImg, tmpDarkPxlsClumpsRMSmallImg, tmpDarkObjsImg)
            
            outputName = tmpBaseName + "DOS" + bandName + imgExtension
            outputImage = os.path.join(outputPath, outputName)
            print(outputImage)

            bandDefns = []
            bandDefns.append(rsgislib.imagecalc.BandDefn('Off', offsetImage, 1))
            bandDefns.append(rsgislib.imagecalc.BandDefn('TOA', inputTOAImage, band))
            expression = '(TOA==0)?0:(TOA-Off)<=0?1:TOA-Off'
            rsgislib.imagecalc.bandMath(outputImage, expression, outFormat, rsgislib.TYPE_16UINT, bandDefns)
                        
            gdalDriver = gdal.GetDriverByName(outFormat)
            gdalDriver.Delete(tmpDarkPxlsImg)
            gdalDriver.Delete(tmpDarkPxlsClumpsImg)
            gdalDriver.Delete(tmpDarkPxlsClumpsRMSmallImg)
            gdalDriver.Delete(tmpDarkObjsImg)
            
            return outputImage
            
        except Exception as e:
            raise e
            
    def performLocalDOSOnSingleBand(self, inputTOAImage, band, outputPath, tmpBaseName, bandName, outFormat, tmpPath, minObjSize, darkPxlPercentile, blockSize):
        try:
            arcsiUtils = ARCSIUtils()
            bandDarkTargetOffsetImages = list()
            imgExtension = arcsiUtils.getFileExtension(outFormat)
            tmpDarkTargetAllImage = os.path.join(tmpPath, tmpBaseName + "_darkpxls_allbands" + imgExtension)
            tmpDarkPxlsImg = os.path.join(tmpPath, tmpBaseName + "_darkpxls" + imgExtension)
            tmpDarkPxlsClumpsImg = os.path.join(tmpPath, tmpBaseName + "_darkclumps" + imgExtension)
            tmpDarkPxlsClumpsRMSmallImg = os.path.join(tmpPath, tmpBaseName + "_darkclumpsrmsmall" + imgExtension)
            tmpDarkObjsImg = os.path.join(tmpPath, tmpBaseName+"_darkobjs"+imgExtension)
            
            binWidth = 1
            self.findDOSLocalDarkTargets(inputTOAImage, tmpDarkTargetAllImage, blockSize, outFormat, binWidth, darkPxlPercentile)
                        
            print("Band ", band)
            rsgislib.imageutils.selectImageBands(tmpDarkTargetAllImage, tmpDarkPxlsImg, outFormat, rsgislib.TYPE_8UINT, [band]) 
            rsgislib.segmentation.clump(tmpDarkPxlsImg, tmpDarkPxlsClumpsImg, outFormat, False, 0.0)                
            rsgislib.rastergis.populateStats(tmpDarkPxlsClumpsImg, True, False)
            rsgislib.segmentation.rmSmallClumps(tmpDarkPxlsClumpsImg, tmpDarkPxlsClumpsRMSmallImg, minObjSize, outFormat)
            rsgislib.segmentation.relabelClumps(tmpDarkPxlsClumpsRMSmallImg, tmpDarkObjsImg, outFormat, False)
            rsgislib.rastergis.populateStats(tmpDarkObjsImg, True, False)
            stats2CalcTOA = list()
            stats2CalcTOA.append(rsgislib.rastergis.BandAttStats(band=(band+1), minField="MinTOARefl", meanField="MeanTOARefl"))
            rsgislib.rastergis.populateRATWithStats(inputTOAImage, tmpDarkObjsImg, stats2CalcTOA)

            ratDS = gdal.Open(tmpDarkObjsImg, gdal.GA_Update)
            Histogram = rat.readColumn(ratDS, "Histogram")
            selected = Histogram * 2
            selected[...] = 1
            selected[0] = 0
            rat.writeColumn(ratDS, "Selected", selected)
            ratDS = None

            rsgislib.rastergis.spatialLocation(tmpDarkObjsImg, "Eastings", "Northings")
            rsgislib.rastergis.selectClumpsOnGrid(tmpDarkObjsImg, "Selected", "SelectedGrid", "Eastings", "Northings", "MeanTOARefl", "min", 20, 20)
        
            print("Interpolating the offset image...")
        
            offsetImage = os.path.join(tmpPath, tmpBaseName+"_darktargetoffs_b"+str(band)+imgExtension)
        
            ratDS = gdal.Open(tmpDarkObjsImg, gdal.GA_Update)
            Eastings = rat.readColumn(ratDS, "Eastings")
            Northings = rat.readColumn(ratDS, "Northings")
            MinTOARefl = rat.readColumn(ratDS, "MinTOARefl")
            SelectedGrid = rat.readColumn(ratDS, "SelectedGrid")
            ratDS = None
        
            Eastings = Eastings[SelectedGrid!=0]
            Northings = Northings[SelectedGrid!=0]
            MinTOARefl = MinTOARefl[SelectedGrid!=0]
        
            interpSmoothing = 10.0
            self.interpolateImageFromPointData(inputTOAImage, Eastings, Northings, MinTOARefl, offsetImage, outFormat, interpSmoothing)
                                                
            outputName = tmpBaseName + "DOS" + bandName + imgExtension
            outputImage = os.path.join(outputPath, outputName)
            print(outputImage)

            bandDefns = []
            bandDefns.append(rsgislib.imagecalc.BandDefn('Off', offsetImage, 1))
            bandDefns.append(rsgislib.imagecalc.BandDefn('TOA', inputTOAImage, band))
            expression = '(TOA==0)?0:(TOA-Off)<=0?1:TOA-Off'
            rsgislib.imagecalc.bandMath(outputImage, expression, outFormat, rsgislib.TYPE_16UINT, bandDefns)
                        
            gdalDriver = gdal.GetDriverByName(outFormat)
            gdalDriver.Delete(tmpDarkPxlsImg)
            gdalDriver.Delete(tmpDarkPxlsClumpsImg)
            gdalDriver.Delete(tmpDarkPxlsClumpsRMSmallImg)
            gdalDriver.Delete(tmpDarkObjsImg)
            gdalDriver.Delete(tmpDarkTargetAllImage)
            gdalDriver.Delete(offsetImage)
            
            return outputImage
            
        except Exception as e:
            raise e
    
    def findDOSLocalDarkTargets(self, inputTOAImage, darkTargetImage, blockSize, outFormat, histBinWidth, darkPxlPercentile):
        reader = ImageReader(inputTOAImage, windowxsize=blockSize, windowysize=blockSize)
        writer = None
        for (info, block) in reader:
            out = numpy.zeros_like(block)
                
            # Iterate through the image bands
            for i in range(len(out)):
                minVal = numpy.min(block[i])
                maxVal = numpy.max(block[i])
                #print("Band: ", i+1, " Min = ", minVal, " Max = ", maxVal)
                if ((maxVal - minVal) > 5):
                    data = block[i].flatten()
                    data = data[data != 0]
                
                    minVal = numpy.min(data)
                    maxVal = numpy.max(data)
                
                    numBins = (math.ceil(maxVal - minVal) / histBinWidth) + 1
                    if data.shape[0] > ((blockSize*blockSize)*0.1):
                        #out[i,...] = 1 
                        histo, histoBins = numpy.histogram(data, bins=numBins, range=(float(minVal), float(maxVal)))
                        numValues = numpy.sum(histo)
                        numValsPercentile = math.floor(numValues * darkPxlPercentile)
                        #print("Num Values: ", numValues, " percentile: ", numValsPercentile)
                        #print("Histo: ", histo)
                        #print("Bins: ", histoBins)
                        binValCount = 0
                        threshold = 0.0
                        for n in range(histo.shape[0]):
                            if (binValCount + histo[n]) > numValsPercentile:
                                break
                            else:
                                binValCount =  binValCount + histo[n]
                                threshold = histoBins[n+1]
                        #print("Threshold = ", threshold)
                        out[i, ((block[i] < threshold) & (block[i] > 0))] = 1
                    else:
                        out[i,...] = 0 
            
            if writer is None:
                writer = ImageWriter(darkTargetImage, 
                                     info=info, 
                                     firstblock=out, 
                                     drivername=outFormat,
                                     creationoptions=[])
            else:
                writer.write(out)
        writer.close(calcStats=False)
    
    def findPerBandLocalDarkTargetsOffsets(self, inputTOAImage, numBands, outputPath, outputName, outFormat, tmpPath, blockSize, minObjSize, darkPxlPercentile):
        try:
            arcsiUtils = ARCSIUtils()
            tmpBaseName = os.path.splitext(outputName)[0]
            binWidth = 1
            
            bandDarkTargetOffsetImages = list()
            imgExtension = arcsiUtils.getFileExtension(outFormat)
            tmpDarkTargetAllImage = os.path.join(tmpPath, tmpBaseName + "_darkpxls_allbands" + imgExtension)
            tmpDarkPxlsImg = os.path.join(tmpPath, tmpBaseName + "_darkpxls" + imgExtension)
            tmpDarkPxlsClumpsImg = os.path.join(tmpPath, tmpBaseName + "_darkclumps" + imgExtension)
            tmpDarkPxlsClumpsRMSmallImg = os.path.join(tmpPath, tmpBaseName + "_darkclumpsrmsmall" + imgExtension)
            tmpDarkObjsImg = os.path.join(tmpPath, tmpBaseName+"_darkobjs"+imgExtension)
            
            self.findDOSLocalDarkTargets(inputTOAImage, tmpDarkTargetAllImage, blockSize, outFormat, binWidth, darkPxlPercentile)
                        
            for band in range(numBands):
                print("Band ", band+1)
                rsgislib.imageutils.selectImageBands(tmpDarkTargetAllImage, tmpDarkPxlsImg, outFormat, rsgislib.TYPE_8UINT, [band+1]) 
                rsgislib.segmentation.clump(tmpDarkPxlsImg, tmpDarkPxlsClumpsImg, outFormat, False, 0.0)                
                rsgislib.rastergis.populateStats(tmpDarkPxlsClumpsImg, True, False)
                rsgislib.segmentation.rmSmallClumps(tmpDarkPxlsClumpsImg, tmpDarkPxlsClumpsRMSmallImg, minObjSize, outFormat)
                rsgislib.segmentation.relabelClumps(tmpDarkPxlsClumpsRMSmallImg, tmpDarkObjsImg, outFormat, False)
                rsgislib.rastergis.populateStats(tmpDarkObjsImg, True, False)
                stats2CalcTOA = list()
                stats2CalcTOA.append(rsgislib.rastergis.BandAttStats(band=(band+1), minField="MinTOARefl", meanField="MeanTOARefl"))
                rsgislib.rastergis.populateRATWithStats(inputTOAImage, tmpDarkObjsImg, stats2CalcTOA)
        
                ratDS = gdal.Open(tmpDarkObjsImg, gdal.GA_Update)
                Histogram = rat.readColumn(ratDS, "Histogram")
                selected = Histogram * 2
                selected[...] = 1
                selected[0] = 0
                rat.writeColumn(ratDS, "Selected", selected)
                ratDS = None
        
                rsgislib.rastergis.spatialLocation(tmpDarkObjsImg, "Eastings", "Northings")
                rsgislib.rastergis.selectClumpsOnGrid(tmpDarkObjsImg, "Selected", "SelectedGrid", "Eastings", "Northings", "MeanTOARefl", "min", 20, 20)
                
                print("Interpolating the offset image...")
                
                offsetImage = os.path.join(tmpPath, tmpBaseName+"_darktargetoffs_b"+str(band+1)+imgExtension)
                
                ratDS = gdal.Open(tmpDarkObjsImg, gdal.GA_Update)
                Eastings = rat.readColumn(ratDS, "Eastings")
                Northings = rat.readColumn(ratDS, "Northings")
                MinTOARefl = rat.readColumn(ratDS, "MinTOARefl")
                SelectedGrid = rat.readColumn(ratDS, "SelectedGrid")
                ratDS = None
                
                Eastings = Eastings[SelectedGrid!=0]
                Northings = Northings[SelectedGrid!=0]
                MinTOARefl = MinTOARefl[SelectedGrid!=0]
                
                interpSmoothing = 10.0
                self.interpolateImageFromPointData(inputTOAImage, Eastings, Northings, MinTOARefl, offsetImage, outFormat, interpSmoothing)
                                                
                bandDarkTargetOffsetImages.append(offsetImage)
            
            outputImage = os.path.join(outputPath, tmpBaseName + "_dosuboffs" + imgExtension)
            print(outputImage)
            rsgislib.imageutils.stackImageBands(bandDarkTargetOffsetImages, None, outputImage, None, 0, outFormat, rsgislib.TYPE_32FLOAT)
            
            gdalDriver = gdal.GetDriverByName(outFormat)
            gdalDriver.Delete(tmpDarkPxlsImg)
            gdalDriver.Delete(tmpDarkPxlsClumpsImg)
            gdalDriver.Delete(tmpDarkPxlsClumpsRMSmallImg)
            gdalDriver.Delete(tmpDarkObjsImg)
            gdalDriver.Delete(tmpDarkTargetAllImage)
            for image in bandDarkTargetOffsetImages:
                gdalDriver.Delete(image)
            
            return outputImage
            
        except Exception as e:
            raise e
    
    @abstractmethod
    def convertImageToReflectanceDarkSubstract(self, inputTOAImage, outputPath, outputName, outFormat, tmpPath, globalDOS, dosOutRefl): pass

    @abstractmethod
    def findDDVTargets(self, inputTOAImage, outputPath, outputName, outFormat, tmpPath): pass

    @abstractmethod
    def estimateImageToAOD(self, inputRADImage, inputTOAImage, inputDEMFile, outputPath, outputName, outFormat, tmpPath, aeroProfile, atmosProfile, grdRefl, aotValMin, aotValMax): pass

    @abstractmethod
    def estimateImageToAODUsingDOS(self, inputRADImage, inputTOAImage, inputDEMFile, outputPath, outputName, outFormat, tmpPath, aeroProfile, atmosProfile, grdRefl, aotValMin, aotValMax, globalDOS, dosOutRefl): pass
    
    @abstractmethod
    def setBandNames(self, imageFile): pass
    
    def interpolateImageFromPointData(self, templateInImage, xVals, yVals, zVals, outputImage, outFormat, smoothingParam):
        print("Interpolating Image: Number of Features = ", xVals.shape[0])
        rbfi = scipy.interpolate.rbf.Rbf(xVals, yVals, zVals, function='linear', smooth=smoothingParam)
        
        reader = ImageReader(templateInImage, windowxsize=200, windowysize=200)
        writer = None
        for (info, block) in reader:
            pxlCoords = info.getBlockCoordArrays()
            interZ = rbfi(pxlCoords[0].flatten(), pxlCoords[1].flatten())
            out = numpy.reshape(interZ, block[0].shape)
            out = numpy.expand_dims(out, axis=0)
    
            if writer is None:
                writer = ImageWriter(outputImage, 
                                     info=info, 
                                     firstblock=out, 
                                     drivername=outFormat,
                                     creationoptions=[])
            else:
                writer.write(out)
        writer.close(calcStats=True)
        print("Interpolating Image - Complete")

    

