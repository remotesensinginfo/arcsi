"""
Module that contains the ARCSILandsat3MSSSensor class.
"""
############################################################################
#  arcsisensorlandsat.py
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
# Purpose:  A class for read the landsat sensor header file and applying
#           the pre-processing operations within ARCSI to the landsat 3 MSS
#           datasets.
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 13/07/2013
# Version: 1.0
#
# History:
# Version 1.0 - Created.
#
############################################################################

# Import updated print function into python 2.7
from __future__ import print_function
# Import updated division operator into python 2.7
from __future__ import division
# import abstract base class stuff
from .arcsisensor import ARCSIAbstractSensor
# Import the ARCSI exception class
from .arcsiexception import ARCSIException
# Import the ARCSI utilities class
from .arcsiutils import ARCSIUtils
# Import the datetime module
import datetime
# Import the GDAL/OGR spatial reference library
from osgeo import osr
from osgeo import ogr
# Import OS path module for manipulating the file system
import os.path
# Import the sys module
import sys
# Import the RSGISLib Image Calibration Module.
import rsgislib.imagecalibration
# Import the RSGISLib Image Utilities Module.
import rsgislib.imageutils
# Import the collections module
import collections
# Import the py6s module for running 6S from python.
import Py6S
# Import the python maths library
import math
# Import the numpy module
import numpy
# Import the GDAL python module
import osgeo.gdal as gdal
# Import the RIOS RAT module
from rios import rat
# Import JSON module
import json
# Import the solar angle tools from RSGISLib
import rsgislib.imagecalibration.solarangles

class ARCSILandsat3MSSSensor (ARCSIAbstractSensor):
    """
    A class which represents the landsat 3 MSS sensor to read
    header parameters and apply data processing operations.
    """
    def __init__(self, debugMode, inputImage):
        ARCSIAbstractSensor.__init__(self, debugMode, inputImage)
        self.sensor = "LS3MSS"
        self.band4File = ""
        self.band5File = ""
        self.band6File = ""
        self.band7File = ""
        self.row = 0
        self.path = 0

        self.b4CalMin = 0
        self.b4CalMax = 0
        self.b5CalMin = 0
        self.b5CalMax = 0
        self.b6CalMin = 0
        self.b6CalMax = 0
        self.b7CalMin = 0
        self.b7CalMax = 0

        self.b4MinRad = 0.0
        self.b4MaxRad = 0.0
        self.b5MinRad = 0.0
        self.b5MaxRad = 0.0
        self.b6MinRad = 0.0
        self.b6MaxRad = 0.0
        self.b7MinRad = 0.0
        self.b7MaxRad = 0.0

        self.sensorID = ""
        self.spacecraftID = ""
        self.cloudCover = 0.0
        self.cloudCoverLand = 0.0
        self.earthSunDistance = 0.0
        self.gridCellSizeRefl = 0.0

    def extractHeaderParameters(self, inputHeader, wktStr):
        """
        Understands and parses the Landsat MTL header files
        """
        try:
            if not self.userSpInputImage is None:
                raise ARCSIException("Landsat sensor cannot accept a user specified image file - only the images in the header file will be used.")
            self.headerFileName = os.path.split(inputHeader)[1]
            
            arcsiUtils = ARCSIUtils()

            print("Reading header file")
            hFile = open(inputHeader, 'r')
            headerParams = dict()
            for line in hFile:
                line = line.strip()
                if line:
                    lineVals = line.split('=')
                    if len(lineVals) == 2:
                        if (lineVals[0].strip() != "GROUP") or (lineVals[0].strip() != "END_GROUP"):
                            headerParams[lineVals[0].strip()] = lineVals[1].strip().replace('"','')
            hFile.close()
            print("Extracting Header Values")
            # Get the sensor info.
            if ((headerParams["SPACECRAFT_ID"].upper() == "LANDSAT_3") or (headerParams["SPACECRAFT_ID"].upper() == "LANDSAT3")) and (headerParams["SENSOR_ID"].upper() == "MSS"):
                self.sensor = "LS3MSS"
            else:
                raise ARCSIException("Do no recognise the spacecraft and sensor or combination.")

            self.sensorID = headerParams["SENSOR_ID"]
            self.spacecraftID = headerParams["SPACECRAFT_ID"]

            # Get row/path
            self.row = int(headerParams["WRS_ROW"])
            self.path = int(headerParams["WRS_PATH"])

            # Get date and time of the acquisition
            acData = headerParams["DATE_ACQUIRED"].split('-')
            acTime = headerParams["SCENE_CENTER_TIME"].split(':')
            secsTime = acTime[2].split('.')
            self.acquisitionTime = datetime.datetime(int(acData[0]), int(acData[1]), int(acData[2]), int(acTime[0]), int(acTime[1]), int(secsTime[0]))

            self.solarZenith = 90-arcsiUtils.str2Float(headerParams["SUN_ELEVATION"])
            self.solarAzimuth = arcsiUtils.str2Float(headerParams["SUN_AZIMUTH"])

            # Get the geographic lat/long corners of the image.
            self.latTL = arcsiUtils.str2Float(headerParams["CORNER_UL_LAT_PRODUCT"])
            self.lonTL = arcsiUtils.str2Float(headerParams["CORNER_UL_LON_PRODUCT"])
            self.latTR = arcsiUtils.str2Float(headerParams["CORNER_UR_LAT_PRODUCT"])
            self.lonTR = arcsiUtils.str2Float(headerParams["CORNER_UR_LON_PRODUCT"])
            self.latBL = arcsiUtils.str2Float(headerParams["CORNER_LL_LAT_PRODUCT"])
            self.lonBL = arcsiUtils.str2Float(headerParams["CORNER_LL_LON_PRODUCT"])
            self.latBR = arcsiUtils.str2Float(headerParams["CORNER_LR_LAT_PRODUCT"])
            self.lonBR = arcsiUtils.str2Float(headerParams["CORNER_LR_LON_PRODUCT"])

            # Get the projected X/Y corners of the image
            self.xTL = arcsiUtils.str2Float(headerParams["CORNER_UL_PROJECTION_X_PRODUCT"])
            self.yTL = arcsiUtils.str2Float(headerParams["CORNER_UL_PROJECTION_Y_PRODUCT"])
            self.xTR = arcsiUtils.str2Float(headerParams["CORNER_UR_PROJECTION_X_PRODUCT"])
            self.yTR = arcsiUtils.str2Float(headerParams["CORNER_UR_PROJECTION_Y_PRODUCT"])
            self.xBL = arcsiUtils.str2Float(headerParams["CORNER_LL_PROJECTION_X_PRODUCT"])
            self.yBL = arcsiUtils.str2Float(headerParams["CORNER_LL_PROJECTION_Y_PRODUCT"])
            self.xBR = arcsiUtils.str2Float(headerParams["CORNER_LR_PROJECTION_X_PRODUCT"])
            self.yBR = arcsiUtils.str2Float(headerParams["CORNER_LR_PROJECTION_Y_PRODUCT"])

            # Get projection
            inProj = osr.SpatialReference()
            if (headerParams["MAP_PROJECTION"] == "UTM") and (headerParams["DATUM"] == "WGS84") and (headerParams["ELLIPSOID"] == "WGS84"):
                utmZone = int(headerParams["UTM_ZONE"])
                utmCode = "WGS84UTM" + str(utmZone) + str("N")
                #print("UTM: ", utmCode)
                inProj.ImportFromEPSG(self.epsgCodes[utmCode])
            elif (headerParams["MAP_PROJECTION"] == "PS") and (headerParams["DATUM"] == "WGS84") and (headerParams["ELLIPSOID"] == "WGS84"):
                inProj.ImportFromWkt("PROJCS[\"PS WGS84\", GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563, AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433],AUTHORITY[\"EPSG\",\"4326\"]],PROJECTION[\"Polar_Stereographic\"],PARAMETER[\"latitude_of_origin\",-71],PARAMETER[\"central_meridian\",0],PARAMETER[\"scale_factor\",1],PARAMETER[\"false_easting\",0],PARAMETER[\"false_northing\",0],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]]]")
            else:
                raise ARCSIException("Expecting Landsat to be projected in UTM or PolarStereographic (PS) with datum=WGS84 and ellipsoid=WGS84.")

            if self.inWKT is "":
                self.inWKT = inProj.ExportToWkt()

            # Check image is square!
            if not ((self.xTL == self.xBL) and (self.yTL == self.yTR) and (self.xTR == self.xBR) and (self.yBL == self.yBR)):
                raise ARCSIException("Image is not square in projected coordinates.")

            self.xCentre = self.xTL + ((self.xTR - self.xTL)/2)
            self.yCentre = self.yBR + ((self.yTL - self.yBR)/2)

            self.lonCentre, self.latCentre = arcsiUtils.getLongLat(inProj, self.xCentre, self.yCentre)

            #print("Lat: " + str(self.latCentre) + " Long: " + str(self.lonCentre))

            filesDIR = os.path.dirname(inputHeader)

            self.band4File = os.path.join(filesDIR, headerParams["FILE_NAME_BAND_4"])
            self.band5File = os.path.join(filesDIR, headerParams["FILE_NAME_BAND_5"])
            self.band6File = os.path.join(filesDIR, headerParams["FILE_NAME_BAND_6"])
            self.band7File = os.path.join(filesDIR, headerParams["FILE_NAME_BAND_7"])

            self.b4CalMin = arcsiUtils.str2Float(headerParams["QUANTIZE_CAL_MIN_BAND_4"], 1.0)
            self.b4CalMax = arcsiUtils.str2Float(headerParams["QUANTIZE_CAL_MAX_BAND_4"], 255.0)
            self.b5CalMin = arcsiUtils.str2Float(headerParams["QUANTIZE_CAL_MIN_BAND_5"], 1.0)
            self.b5CalMax = arcsiUtils.str2Float(headerParams["QUANTIZE_CAL_MAX_BAND_5"], 255.0)
            self.b6CalMin = arcsiUtils.str2Float(headerParams["QUANTIZE_CAL_MIN_BAND_6"], 1.0)
            self.b6CalMax = arcsiUtils.str2Float(headerParams["QUANTIZE_CAL_MAX_BAND_6"], 255.0)
            self.b7CalMin = arcsiUtils.str2Float(headerParams["QUANTIZE_CAL_MIN_BAND_7"], 1.0)
            self.b7CalMax = arcsiUtils.str2Float(headerParams["QUANTIZE_CAL_MAX_BAND_7"], 255.0)

            self.b4MinRad = arcsiUtils.str2Float(headerParams["RADIANCE_MINIMUM_BAND_4"], 3.600)
            self.b4MaxRad = arcsiUtils.str2Float(headerParams["RADIANCE_MAXIMUM_BAND_4"], 234.600)
            self.b5MinRad = arcsiUtils.str2Float(headerParams["RADIANCE_MINIMUM_BAND_5"], 2.800)
            self.b5MaxRad = arcsiUtils.str2Float(headerParams["RADIANCE_MAXIMUM_BAND_5"], 164.200)
            self.b6MinRad = arcsiUtils.str2Float(headerParams["RADIANCE_MINIMUM_BAND_6"], 2.900)
            self.b6MaxRad = arcsiUtils.str2Float(headerParams["RADIANCE_MAXIMUM_BAND_6"], 146.200)
            self.b7MinRad = arcsiUtils.str2Float(headerParams["RADIANCE_MINIMUM_BAND_7"], 1.000)
            self.b7MaxRad = arcsiUtils.str2Float(headerParams["RADIANCE_MAXIMUM_BAND_7"], 121.700)

            if "CLOUD_COVER" in headerParams:
                self.cloudCover = arcsiUtils.str2Float(headerParams["CLOUD_COVER"], 0.0)
            if "CLOUD_COVER_LAND" in headerParams:
                self.cloudCoverLand = arcsiUtils.str2Float(headerParams["CLOUD_COVER_LAND"], 0.0)
            if "EARTH_SUN_DISTANCE" in headerParams:
                self.earthSunDistance = arcsiUtils.str2Float(headerParams["EARTH_SUN_DISTANCE"], 0.0)
            if "GRID_CELL_SIZE_REFLECTIVE" in headerParams:
                self.gridCellSizeRefl = arcsiUtils.str2Float(headerParams["GRID_CELL_SIZE_REFLECTIVE"], 60.0)

            fileDateStr = headerParams["FILE_DATE"].strip()
            fileDateStr = fileDateStr.replace('Z', '')
            self.fileDateObj = datetime.datetime.strptime(fileDateStr, "%Y-%m-%dT%H:%M:%S")

        except Exception as e:
            raise e

    def getSolarIrrStdSolarGeom(self):
        """
        Get Solar Azimuth and Zenith as standard geometry.
        Azimuth: N=0, E=90, S=180, W=270.
        """
        solarAz = rsgislib.imagecalibration.solarangles.getSolarIrrConventionSolarAzimuthFromUSGS(self.solarAzimuth)
        return (solarAz, self.solarZenith)

    def getSensorViewGeom(self):
        """
        Get sensor viewing angles
        returns (viewAzimuth, viewZenith)
        """
        return (0.0, 0.0)

    def generateOutputBaseName(self):
        """
        Provides an implementation for the landsat sensor
        """
        rowpath = "r" + str(self.row) + "p" + str(self.path)
        outname = self.defaultGenBaseOutFileName()
        outname = outname + str("_") + rowpath
        return outname

    def generateMetaDataFile(self, outputPath, outputFileName, productsStr, validMaskImage="", footprintCalc=False, calcdValuesDict=dict(), outFilesDict=dict()):
        """
        Generate file metadata.
        """
        outJSONFilePath = os.path.join(outputPath, outputFileName)
        jsonData = self.getJSONDictDefaultMetaData(productsStr, validMaskImage, footprintCalc, calcdValuesDict, outFilesDict)
        sensorInfo = jsonData['SensorInfo']
        sensorInfo['Row'] = self.row
        sensorInfo['Path'] = self.path
        sensorInfo['SensorID'] = self.sensorID
        sensorInfo['SpacecraftID'] = self.spacecraftID
        acqDict = jsonData['AcquasitionInfo']
        acqDict['EarthSunDistance'] = self.earthSunDistance
        imgInfo = dict()
        imgInfo['CloudCover'] = self.cloudCover
        imgInfo['CloudCoverLand'] = self.cloudCoverLand
        imgInfo['CellSizeRefl'] = self.gridCellSizeRefl
        jsonData['ImageInfo'] = imgInfo

        with open(outJSONFilePath, 'w') as outfile:
            json.dump(jsonData, outfile, sort_keys=True,indent=4, separators=(',', ': '), ensure_ascii=False)

    def expectedImageDataPresent(self):
        imageDataPresent = True

        if not os.path.exists(self.band4File):
            imageDataPresent = False
        if not os.path.exists(self.band5File):
            imageDataPresent = False
        if not os.path.exists(self.band6File):
            imageDataPresent = False
        if not os.path.exists(self.band7File):
            imageDataPresent = False

        return imageDataPresent

    def applyImageDataMask(self, inputHeader, outputPath, outputMaskName, outputImgName, outFormat, outWKTFile):
        raise ARCSIException("Landsat 3 does not provide any image masks, do not use the MASK option.")

    def mosaicImageTiles(self, outputPath):
        raise ARCSIException("Image data does not need mosaicking")

    def resampleImgRes(self, outputPath, resampleToLowResImg, resampleMethod='cubic', multicore=False):
        raise ARCSIException("Image data does not need resampling")

    def sharpenLowResRadImgBands(self, inputImg, outputImage, outFormat):
        raise ARCSIException("Image sharpening is not available for this sensor.")

    def generateValidImageDataMask(self, outputPath, outputMaskName, viewAngleImg, outFormat):
        print("Create the valid data mask")
        tmpBaseName = os.path.splitext(outputMaskName)[0]
        tmpValidPxlMsk = os.path.join(outputPath, tmpBaseName+'vldpxlmsk.kea')
        outputImage = os.path.join(outputPath, outputMaskName)
        inImages = [self.band4File, self.band5File, self.band6File, self.band7File]
        rsgislib.imageutils.genValidMask(inimages=inImages, outimage=tmpValidPxlMsk, gdalformat='KEA', nodata=0.0)
        rsgislib.rastergis.populateStats(tmpValidPxlMsk, True, False, True)
        # Check there is valid data
        ratDS = gdal.Open(tmpValidPxlMsk, gdal.GA_ReadOnly)
        Histogram = rat.readColumn(ratDS, "Histogram")
        ratDS = None
        if Histogram.shape[0] < 2:
            raise ARCSIException("There is no valid data in this image.")
        if not os.path.exists(viewAngleImg):
            rsgislib.rastergis.spatialExtent(clumps=tmpValidPxlMsk, minXX='MinXX', minXY='MinXY', maxXX='MaxXX', maxXY='MaxXY', minYX='MinYX', minYY='MinYY', maxYX='MaxYX', maxYY='MaxYY', ratband=1)
            rsgislib.imagecalibration.calcNadirImgViewAngle(tmpValidPxlMsk, viewAngleImg, 'KEA', 705000.0, 'MinXX', 'MinXY', 'MaxXX', 'MaxXY', 'MinYX', 'MinYY', 'MaxYX', 'MaxYY')
        rsgislib.imagecalc.imageMath(viewAngleImg, outputImage, 'b1<7.65?1:0', outFormat, rsgislib.TYPE_8UINT)
        rsgisUtils = rsgislib.RSGISPyUtils()
        rsgisUtils.deleteFileWithBasename(tmpValidPxlMsk)
        return outputImage

    def convertImageToRadiance(self, outputPath, outputReflName, outputThermalName, outFormat):
        print("Converting to Radiance")
        outputImage = os.path.join(outputPath, outputReflName)
        bandDefnSeq = list()
        lsBand = collections.namedtuple('LSBand', ['bandName', 'fileName', 'bandIndex', 'lMin', 'lMax', 'qCalMin', 'qCalMax'])
        bandDefnSeq.append(lsBand(bandName="Green", fileName=self.band4File, bandIndex=1, lMin=self.b4MinRad, lMax=self.b4MaxRad, qCalMin=self.b4CalMin, qCalMax=self.b4CalMax))
        bandDefnSeq.append(lsBand(bandName="Red", fileName=self.band5File, bandIndex=1, lMin=self.b5MinRad, lMax=self.b5MaxRad, qCalMin=self.b5CalMin, qCalMax=self.b5CalMax))
        bandDefnSeq.append(lsBand(bandName="NIR1", fileName=self.band6File, bandIndex=1, lMin=self.b6MinRad, lMax=self.b6MaxRad, qCalMin=self.b6CalMin, qCalMax=self.b6CalMax))
        bandDefnSeq.append(lsBand(bandName="NIR2", fileName=self.band7File, bandIndex=1, lMin=self.b7MinRad, lMax=self.b7MaxRad, qCalMin=self.b7CalMin, qCalMax=self.b7CalMax))
        rsgislib.imagecalibration.landsat2Radiance(outputImage, outFormat, bandDefnSeq)
        return outputImage, None

    def convertThermalToBrightness(self, inputRadImage, outputPath, outputName, outFormat, scaleFactor):
        raise ARCSIException("There are no thermal bands...")

    def convertImageToTOARefl(self, inputRadImage, outputPath, outputName, outFormat, scaleFactor):
        print("Converting to TOA")
        outputImage = os.path.join(outputPath, outputName)
        solarIrradianceVals = list()
        IrrVal = collections.namedtuple('SolarIrradiance', ['irradiance'])
        solarIrradianceVals.append(IrrVal(irradiance=1839.0))
        solarIrradianceVals.append(IrrVal(irradiance=1555.0))
        solarIrradianceVals.append(IrrVal(irradiance=1291.0))
        solarIrradianceVals.append(IrrVal(irradiance=887.9))
        rsgislib.imagecalibration.radiance2TOARefl(inputRadImage, outputImage, outFormat, rsgislib.TYPE_16UINT, scaleFactor, self.acquisitionTime.year, self.acquisitionTime.month, self.acquisitionTime.day, self.solarZenith, solarIrradianceVals)
        return outputImage

    def generateImageSaturationMask(self, outputPath, outputName, outFormat):
        print("Generate Saturation Image")
        outputImage = os.path.join(outputPath, outputName)

        lsBand = collections.namedtuple('LSBand', ['bandName', 'fileName', 'bandIndex', 'satVal'])
        bandDefnSeq = list()
        bandDefnSeq.append(lsBand(bandName="Green", fileName=self.band4File, bandIndex=1, satVal=self.b4CalMax))
        bandDefnSeq.append(lsBand(bandName="Red", fileName=self.band5File, bandIndex=1, satVal=self.b5CalMax))
        bandDefnSeq.append(lsBand(bandName="NIR1", fileName=self.band6File, bandIndex=1, satVal=self.b6CalMax))
        bandDefnSeq.append(lsBand(bandName="NIR2", fileName=self.band7File, bandIndex=1, satVal=self.b7CalMax))

        rsgislib.imagecalibration.saturatedPixelsMask(outputImage, outFormat, bandDefnSeq)

        return outputImage

    def generateCloudMask(self, inputReflImage, inputSatImage, inputThermalImage, inputViewAngleImg, inputValidImg, outputPath, outputName, outFormat, tmpPath, scaleFactor, cloud_msk_methods=None):
        raise ARCSIException("Cloud Masking Not Implemented for LS3 MSS.")

    def createCloudMaskDataArray(self, inImgDataArr):
        # Calc Whiteness
        meanArr = numpy.mean(inImgDataArr, axis=1)
        whitenessArr = numpy.absolute((inImgDataArr[...,0] - meanArr)/meanArr) + numpy.absolute((inImgDataArr[...,1] - meanArr)/meanArr) + numpy.absolute((inImgDataArr[...,2] - meanArr)/meanArr) + numpy.absolute((inImgDataArr[...,3] - meanArr)/meanArr)
        # Calc NDVI
        ndvi = (inImgDataArr[...,3] - inImgDataArr[...,1]) / (inImgDataArr[...,3] + inImgDataArr[...,1])
        
        # Create and populate the output array.
        inShape = inImgDataArr.shape
        outShape = [inShape[0], inShape[1]+3]    
        outArr = numpy.zeros(outShape, dtype=float)
        
        for i in range(inShape[1]):
            outArr[...,i] = inImgDataArr[...,i]
        
        idx = inShape[1]
        outArr[...,idx] = meanArr
        outArr[...,idx+1] = whitenessArr
        outArr[...,idx+2] = ndvi
        
        return outArr

    def defineDarkShadowImageBand(self):
        return 4

    def calc6SCoefficients(self, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF):
        sixsCoeffs = numpy.zeros((4, 6), dtype=numpy.float32)
        # Set up 6S model
        s = Py6S.SixS()
        s.atmos_profile = atmosProfile
        s.aero_profile = aeroProfile
        #s.ground_reflectance = Py6S.GroundReflectance.HomogeneousHapke(0.101, -0.263, 0.589, 0.046)
        s.ground_reflectance = grdRefl
        s.geometry = Py6S.Geometry.Landsat_TM()
        s.geometry.month = self.acquisitionTime.month
        s.geometry.day = self.acquisitionTime.day
        s.geometry.gmt_decimal_hour = float(self.acquisitionTime.hour) + float(self.acquisitionTime.minute)/60.0
        s.geometry.latitude = self.latCentre
        s.geometry.longitude = self.lonCentre
        s.altitudes = Py6S.Altitudes()
        s.altitudes.set_target_custom_altitude(surfaceAltitude)
        s.altitudes.set_sensor_satellite_level()
        if useBRDF:
            s.atmos_corr = Py6S.AtmosCorr.AtmosCorrBRDFFromRadiance(200)
        else:
            s.atmos_corr = Py6S.AtmosCorr.AtmosCorrLambertianFromRadiance(200)
        s.aot550 = aotVal

        # Band 1
        s.wavelength = Py6S.Wavelength(Py6S.SixSHelpers.PredefinedWavelengths.LANDSAT_MSS_B1)
        s.run()
        sixsCoeffs[0,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[0,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[0,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[0,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[0,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[0,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 2
        s.wavelength = Py6S.Wavelength(Py6S.SixSHelpers.PredefinedWavelengths.LANDSAT_MSS_B2)
        s.run()
        sixsCoeffs[1,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[1,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[1,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[1,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[1,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[1,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 3
        s.wavelength = Py6S.Wavelength(Py6S.SixSHelpers.PredefinedWavelengths.LANDSAT_MSS_B3)
        s.run()
        sixsCoeffs[2,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[2,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[2,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[2,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[2,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[2,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 4
        s.wavelength = Py6S.Wavelength(Py6S.SixSHelpers.PredefinedWavelengths.LANDSAT_MSS_B4)
        s.run()
        sixsCoeffs[3,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[3,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[3,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[3,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[3,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[3,5] = float(s.outputs.values['environmental_irradiance'])

        return sixsCoeffs

    def convertImageToSurfaceReflSglParam(self, inputRadImage, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF, scaleFactor):
        print("Converting to Surface Reflectance")
        outputImage = os.path.join(outputPath, outputName)

        imgBandCoeffs = list()

        sixsCoeffs = self.calc6SCoefficients(aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF)

        imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=1, aX=float(sixsCoeffs[0,0]), bX=float(sixsCoeffs[0,1]), cX=float(sixsCoeffs[0,2]), DirIrr=float(sixsCoeffs[0,3]), DifIrr=float(sixsCoeffs[0,4]), EnvIrr=float(sixsCoeffs[0,5])))
        imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=2, aX=float(sixsCoeffs[1,0]), bX=float(sixsCoeffs[1,1]), cX=float(sixsCoeffs[1,2]), DirIrr=float(sixsCoeffs[1,3]), DifIrr=float(sixsCoeffs[1,4]), EnvIrr=float(sixsCoeffs[1,5])))
        imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=3, aX=float(sixsCoeffs[2,0]), bX=float(sixsCoeffs[2,1]), cX=float(sixsCoeffs[2,2]), DirIrr=float(sixsCoeffs[2,3]), DifIrr=float(sixsCoeffs[2,4]), EnvIrr=float(sixsCoeffs[2,5])))
        imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=4, aX=float(sixsCoeffs[3,0]), bX=float(sixsCoeffs[3,1]), cX=float(sixsCoeffs[3,2]), DirIrr=float(sixsCoeffs[3,3]), DifIrr=float(sixsCoeffs[3,4]), EnvIrr=float(sixsCoeffs[3,5])))

        rsgislib.imagecalibration.apply6SCoeffSingleParam(inputRadImage, outputImage, outFormat, rsgislib.TYPE_16UINT, scaleFactor, 0, True, imgBandCoeffs)
        return outputImage

    def convertImageToSurfaceReflDEMElevLUT(self, inputRadImage, inputDEMFile, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, aotVal, useBRDF, surfaceAltitudeMin, surfaceAltitudeMax, scaleFactor, elevCoeffs=None):
        print("Converting to Surface Reflectance")
        outputImage = os.path.join(outputPath, outputName)

        if elevCoeffs is None:
            print("Build an LUT for elevation values.")
            elev6SCoeffsLUT = self.buildElevation6SCoeffLUT(aeroProfile, atmosProfile, grdRefl, aotVal, useBRDF, surfaceAltitudeMin, surfaceAltitudeMax)
            print("LUT has been built.")

            elevCoeffs = list()
            for elevLUT in elev6SCoeffsLUT:
                imgBandCoeffs = list()
                sixsCoeffs = elevLUT.Coeffs
                elevVal = elevLUT.Elev
                imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=1, aX=float(sixsCoeffs[0,0]), bX=float(sixsCoeffs[0,1]), cX=float(sixsCoeffs[0,2]), DirIrr=float(sixsCoeffs[0,3]), DifIrr=float(sixsCoeffs[0,4]), EnvIrr=float(sixsCoeffs[0,5])))
                imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=2, aX=float(sixsCoeffs[1,0]), bX=float(sixsCoeffs[1,1]), cX=float(sixsCoeffs[1,2]), DirIrr=float(sixsCoeffs[1,3]), DifIrr=float(sixsCoeffs[1,4]), EnvIrr=float(sixsCoeffs[1,5])))
                imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=3, aX=float(sixsCoeffs[2,0]), bX=float(sixsCoeffs[2,1]), cX=float(sixsCoeffs[2,2]), DirIrr=float(sixsCoeffs[2,3]), DifIrr=float(sixsCoeffs[2,4]), EnvIrr=float(sixsCoeffs[2,5])))
                imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=4, aX=float(sixsCoeffs[3,0]), bX=float(sixsCoeffs[3,1]), cX=float(sixsCoeffs[3,2]), DirIrr=float(sixsCoeffs[3,3]), DifIrr=float(sixsCoeffs[3,4]), EnvIrr=float(sixsCoeffs[3,5])))

                elevCoeffs.append(rsgislib.imagecalibration.ElevLUTFeat(Elev=float(elevVal), Coeffs=imgBandCoeffs))

        rsgislib.imagecalibration.apply6SCoeffElevLUTParam(inputRadImage, inputDEMFile, outputImage, outFormat, rsgislib.TYPE_16UINT, scaleFactor, 0, True, elevCoeffs)
        return outputImage, elevCoeffs


    def convertImageToSurfaceReflAOTDEMElevLUT(self, inputRadImage, inputDEMFile, inputAOTImage, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, useBRDF, surfaceAltitudeMin, surfaceAltitudeMax, aotMin, aotMax, scaleFactor, elevAOTCoeffs=None):
        print("Converting to Surface Reflectance")
        outputImage = os.path.join(outputPath, outputName)

        if elevAOTCoeffs is None:
            print("Build an LUT for elevation and AOT values.")
            elevAOT6SCoeffsLUT = self.buildElevationAOT6SCoeffLUT(aeroProfile, atmosProfile, grdRefl, useBRDF, surfaceAltitudeMin, surfaceAltitudeMax, aotMin, aotMax)

            elevAOTCoeffs = list()
            for elevLUT in elevAOT6SCoeffsLUT:
                elevVal = elevLUT.Elev
                aotLUT = elevLUT.Coeffs
                aot6SCoeffsOut = list()
                for aotFeat in aotLUT:
                    sixsCoeffs = aotFeat.Coeffs
                    aotVal = aotFeat.AOT
                    imgBandCoeffs = list()
                    imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=1, aX=float(sixsCoeffs[0,0]), bX=float(sixsCoeffs[0,1]), cX=float(sixsCoeffs[0,2]), DirIrr=float(sixsCoeffs[0,3]), DifIrr=float(sixsCoeffs[0,4]), EnvIrr=float(sixsCoeffs[0,5])))
                    imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=2, aX=float(sixsCoeffs[1,0]), bX=float(sixsCoeffs[1,1]), cX=float(sixsCoeffs[1,2]), DirIrr=float(sixsCoeffs[1,3]), DifIrr=float(sixsCoeffs[1,4]), EnvIrr=float(sixsCoeffs[1,5])))
                    imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=3, aX=float(sixsCoeffs[2,0]), bX=float(sixsCoeffs[2,1]), cX=float(sixsCoeffs[2,2]), DirIrr=float(sixsCoeffs[2,3]), DifIrr=float(sixsCoeffs[2,4]), EnvIrr=float(sixsCoeffs[2,5])))
                    imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=4, aX=float(sixsCoeffs[3,0]), bX=float(sixsCoeffs[3,1]), cX=float(sixsCoeffs[3,2]), DirIrr=float(sixsCoeffs[3,3]), DifIrr=float(sixsCoeffs[3,4]), EnvIrr=float(sixsCoeffs[3,5])))
                    aot6SCoeffsOut.append(rsgislib.imagecalibration.AOTLUTFeat(AOT=float(aotVal), Coeffs=imgBandCoeffs))
                elevAOTCoeffs.append(rsgislib.imagecalibration.ElevLUTFeat(Elev=float(elevVal), Coeffs=aot6SCoeffsOut))

        rsgislib.imagecalibration.apply6SCoeffElevAOTLUTParam(inputRadImage, inputDEMFile, inputAOTImage, outputImage, outFormat, rsgislib.TYPE_16UINT, scaleFactor, 0, True, elevAOTCoeffs)

        return outputImage, elevAOTCoeffs

    def run6SToOptimiseAODValue(self, aotVal, radBlueVal, predBlueVal, aeroProfile, atmosProfile, grdRefl, surfaceAltitude):
        """Used as part of the optimastion for identifying values of AOD"""
        print("Testing AOD Val: ", aotVal,)
        s = Py6S.SixS()

        s.atmos_profile = atmosProfile
        s.aero_profile = aeroProfile
        s.ground_reflectance = grdRefl
        s.geometry = Py6S.Geometry.Landsat_TM()
        s.geometry.month = self.acquisitionTime.month
        s.geometry.day = self.acquisitionTime.day
        s.geometry.gmt_decimal_hour = float(self.acquisitionTime.hour) + float(self.acquisitionTime.minute)/60.0
        s.geometry.latitude = self.latCentre
        s.geometry.longitude = self.lonCentre
        s.altitudes = Py6S.Altitudes()
        s.altitudes.set_target_custom_altitude(surfaceAltitude)
        s.altitudes.set_sensor_satellite_level()
        s.atmos_corr = Py6S.AtmosCorr.AtmosCorrLambertianFromRadiance(200)
        s.aot550 = aotVal

        # Band 1 (Blue!)
        s.wavelength = Py6S.Wavelength(Py6S.SixSHelpers.PredefinedWavelengths.LANDSAT_MSS_B1)
        s.run()
        aX = float(s.outputs.values['coef_xa'])
        bX = float(s.outputs.values['coef_xb'])
        cX = float(s.outputs.values['coef_xc'])
        tmpVal = (aX*radBlueVal)-bX;
        reflBlueVal = tmpVal/(1.0+cX*tmpVal)

        outDist = math.sqrt(math.pow((reflBlueVal - predBlueVal),2))
        print("\taX: ", aX, " bX: ", bX, " cX: ", cX, "     Dist = ", outDist)
        return outDist

    def findDDVTargets(self, inputTOAImage, outputPath, outputName, outFormat, tmpPath):
        print("Not implemented\n")
        sys.exit()

    def estimateImageToAODUsingDDV(self, inputTOAImage, outputPath, outputName, outFormat, tmpPath, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotValMin, aotValMax):
        print("Not implemented\n")
        sys.exit()

    def estimateImageToAODUsingDOS(self, inputRADImage, inputTOAImage, inputDEMFile, shadowMask, outputPath, outputName, outFormat, tmpPath, aeroProfile, atmosProfile, grdRefl, aotValMin, aotValMax, globalDOS, simpleDOS, dosOutRefl):
        try:
            print("Estimating AOD Using DOS")
            arcsiUtils = ARCSIUtils()
            outputAOTImage = os.path.join(outputPath, outputName)
            tmpBaseName = os.path.splitext(outputName)[0]
            imgExtension = arcsiUtils.getFileExtension(outFormat)

            dosBlueImage = ""
            minObjSize = 3
            darkPxlPercentile = 0.01
            blockSize = 1000
            if simpleDOS:
                outputDOSBlueName = tmpBaseName + "DOSBlue" + imgExtension
                dosBlueImage, bandOff = self.convertImageBandToReflectanceSimpleDarkSubtract(inputTOAImage, outputPath, outputDOSBlueName, outFormat, dosOutRefl, 1)
            elif globalDOS:
                dosBlueImage = self.performDOSOnSingleBand(inputTOAImage, 1, outputPath, tmpBaseName, "Blue", "KEA", tmpPath, minObjSize, darkPxlPercentile, dosOutRefl)
            else:
                dosBlueImage = self.performLocalDOSOnSingleBand(inputTOAImage, 1, outputPath, tmpBaseName, "Blue", "KEA", tmpPath, minObjSize, darkPxlPercentile, blockSize, dosOutRefl)

            thresImageClumpsFinal = os.path.join(tmpPath, tmpBaseName + "_clumps" + imgExtension)
            rsgislib.segmentation.segutils.runShepherdSegmentation(inputTOAImage, thresImageClumpsFinal, tmpath=tmpPath, gdalformat="KEA", numClusters=40, minPxls=10, bands=[1,2,3,4], processInMem=True)

            stats2CalcTOA = list()
            stats2CalcTOA.append(rsgislib.rastergis.BandAttStats(band=1, meanField="MeanElev"))
            rsgislib.rastergis.populateRATWithStats(inputDEMFile, thresImageClumpsFinal, stats2CalcTOA)

            stats2CalcTOA = list()
            stats2CalcTOA.append(rsgislib.rastergis.BandAttStats(band=1, meanField="MeanB1DOS"))
            rsgislib.rastergis.populateRATWithStats(dosBlueImage, thresImageClumpsFinal, stats2CalcTOA)

            stats2CalcRad = list()
            stats2CalcRad.append(rsgislib.rastergis.BandAttStats(band=1, meanField="MeanB1RAD"))
            stats2CalcRad.append(rsgislib.rastergis.BandAttStats(band=4, meanField="MeanB4RAD"))
            stats2CalcRad.append(rsgislib.rastergis.BandAttStats(band=3, meanField="MeanB3RAD"))
            rsgislib.rastergis.populateRATWithStats(inputRADImage, thresImageClumpsFinal, stats2CalcRad)

            ratDS = gdal.Open(thresImageClumpsFinal, gdal.GA_Update)
            Histogram = rat.readColumn(ratDS, "Histogram")
            MeanElev = rat.readColumn(ratDS, "MeanElev")

            MeanB4RAD = rat.readColumn(ratDS, "MeanB4RAD")
            MeanB3RAD = rat.readColumn(ratDS, "MeanB3RAD")

            radNDVI = (MeanB4RAD - MeanB3RAD)/(MeanB4RAD + MeanB3RAD)

            selected = Histogram * 2
            selected[...] = 0
            selected[radNDVI>0.2] = 1
            rat.writeColumn(ratDS, "Selected", selected)
            ratDS = None

            rsgislib.rastergis.spatialLocation(thresImageClumpsFinal, "Eastings", "Northings")
            rsgislib.rastergis.selectClumpsOnGrid(thresImageClumpsFinal, "Selected", "PredictAOTFor", "Eastings", "Northings", "MeanB1DOS", "min", 20, 20)

            ratDS = gdal.Open(thresImageClumpsFinal, gdal.GA_Update)
            MeanB1DOS = rat.readColumn(ratDS, "MeanB1DOS")
            MeanB1DOS = MeanB1DOS / 1000
            MeanB1RAD = rat.readColumn(ratDS, "MeanB1RAD")
            PredictAOTFor = rat.readColumn(ratDS, "PredictAOTFor")

            numAOTValTests = int(math.ceil((aotValMax - aotValMin)/0.05))+1

            if not numAOTValTests >= 1:
                raise ARCSIException("min and max AOT range are too close together, they need to be at least 0.05 apart.")

            cAOT = aotValMin
            cDist = 0.0
            minAOT = 0.0
            minDist = 0.0

            aotVals = numpy.zeros_like(MeanB1RAD, dtype=numpy.float)

            for i in range(len(MeanB1RAD)):
                if PredictAOTFor[i] == 1:
                    print("Predicting AOD for Segment ", i)
                    for j in range(numAOTValTests):
                        cAOT = aotValMin + (0.05 * j)
                        cDist = self.run6SToOptimiseAODValue(cAOT, MeanB1RAD[i], MeanB1DOS[i], aeroProfile, atmosProfile, grdRefl, MeanElev[i]/1000)
                        if j == 0:
                            minAOT = cAOT
                            minDist = cDist
                        elif cDist < minDist:
                            minAOT = cAOT
                            minDist = cDist
                    #predAOTArgs = (MinB1RAD[i], MeanB1DOS[i], aeroProfile, atmosProfile, grdRefl, MeanElev[i]/1000)
                    #res = minimize(self.run6SToOptimiseAODValue, minAOT, method='nelder-mead', options={'maxiter': 20, 'xtol': 0.001, 'disp': True}, args=predAOTArgs)
                    #aotVals[i] = res.x[0]
                    aotVals[i] = minAOT
                    print("IDENTIFIED AOT: ", aotVals[i])
                else:
                    aotVals[i] = 0
            rat.writeColumn(ratDS, "AOT", aotVals)

            Eastings = rat.readColumn(ratDS, "Eastings")
            Northings = rat.readColumn(ratDS, "Northings")
            ratDS = None

            Eastings = Eastings[PredictAOTFor!=0]
            Northings = Northings[PredictAOTFor!=0]
            aotVals = aotVals[PredictAOTFor!=0]

            interpSmoothing = 10.0
            self.interpolateImageFromPointData(inputTOAImage, Eastings, Northings, aotVals, outputAOTImage, outFormat, interpSmoothing, True, 0.05)

            if not self.debugMode:
                gdalDriver = gdal.GetDriverByName(outFormat)
                gdalDriver.Delete(thresImageClumpsFinal)
                gdalDriver.Delete(dosBlueImage)

            return outputAOTImage
        except Exception as e:
            raise e

    def estimateSingleAOTFromDOS(self, radianceImage, toaImage, inputDEMFile, tmpPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, minAOT, maxAOT, dosOutRefl):
        try:
            return self.estimateSingleAOTFromDOSBandImpl(radianceImage, toaImage, inputDEMFile, tmpPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, minAOT, maxAOT, dosOutRefl, 1)
        except Exception as e:
            raise

    def setBandNames(self, imageFile):
        """
        Set band names for Landsat 3 MSS

        Green: 500 - 600 nm
        Red: 600 - 700 nm
        NIR1: 700 - 800 nm
        NIR2: 800 - 1100 nm

        http://landsat.usgs.gov/about_landsat3.php

        """
        dataset = gdal.Open(imageFile, gdal.GA_Update)
        dataset.GetRasterBand(1).SetDescription("Green")
        dataset.GetRasterBand(2).SetDescription("Red")
        dataset.GetRasterBand(3).SetDescription("NIR1")
        dataset.GetRasterBand(4).SetDescription("NIR2")
        dataset = None

    def cleanLocalFollowProcessing(self):
        print("")


