"""
Module that contains the ARCSISPOT6Sensor class.
"""
############################################################################
#  arcsisensorspot6.py
#
#  Copyright 2017 ARCSI.
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
# Purpose:  A class for read the SPOT6 sensor header file and applying
#           the pre-processing operations within ARCSI to the SPOT6
#           datasets.
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 03/01/2017
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
# Import the OS module
import os
# Import OS path module for manipulating the file system
import os.path
# Import the RSGISLib Image Calibration Module.
import rsgislib.imagecalibration
# Import the RSGISLib Image Calculation Module
import rsgislib.imagecalc
# Import the RSGISLib Image Utilities Module
import rsgislib.imageutils
# Import the RSGISLib Vector Utilities Module
import rsgislib.vectorutils
# Import the collections module
import collections
# Import the py6s module for running 6S from python.
import Py6S
# Import the python maths library
import math
# Import python XML Parser
import xml.etree.ElementTree as ET
# Import the numpy module
import numpy
# Import the GDAL python module
import osgeo.gdal as gdal
# Import the python subprocess module - used to call commands line tools.
import subprocess
# Import the RIOS RAT module
from rios import rat

class ARCSISPOT6Sensor (ARCSIAbstractSensor):
    """
    A class which represents the SPOT6 sensor to read
    header parameters and apply data processing operations.
    """
    def __init__(self, debugMode, inputImage):
        ARCSIAbstractSensor.__init__(self, debugMode, inputImage)
        self.sensor = "SPOT6"
        self.fileName = ""
        self.inputImgType = ""
        self.inputImgFiles = []
        self.tiledInputImg = False
        self.inputImgProjed = True
        self.b0SolarIrr = 0.0
        self.b1SolarIrr = 0.0
        self.b2SolarIrr = 0.0
        self.b3SolarIrr = 0.0
        self.b0RadGain = 0.0
        self.b0RadBias = 0.0
        self.b1RadGain = 0.0
        self.b1RadBias = 0.0
        self.b2RadGain = 0.0
        self.b2RadBias = 0.0
        self.b3RadGain = 0.0
        self.b3RadBias = 0.0
        self.acqBitRange = 12
        self.prodBitRange = 12
        self.copiedKEADNImg = ''
        self.origDNImg = ''
        self.createdCopyKEADNImg = False
        self.cloudsMask = None
        self.roiMask = None
        

    def extractHeaderParameters(self, inputHeader, wktStr):
        """
        Understands and parses the SPOT 6 header file
        """
        try:
            self.headerFileName = os.path.split(inputHeader)[1]
            
            tree = ET.parse(inputHeader)
            root = tree.getroot()
            
            topLevelMetaIdent = root.find('Metadata_Identification')
            if topLevelMetaIdent == None:
                raise ARCSIException("Cannot open top level section \'Metadata_Identification\'")

            dimapVersion = topLevelMetaIdent.find('METADATA_FORMAT').attrib['version'].strip()
            if dimapVersion.split('.')[0] != '2':
                raise ARCSIException("Only DIMAP Version 2.X is supported by this reader; provided with: \'" + dimapVersion + "\'")

            metaSensorProfile = topLevelMetaIdent.find('METADATA_PROFILE').text.strip()
            if not ((metaSensorProfile == 'S6_SENSOR') or (metaSensorProfile == 'S6_ORTHO')):
                raise ARCSIException("Input file is not for SPOT 6, \'METADATA_PROFILE\' should be \'S6_SENSOR\' or \'S6_ORTHO\'; provided with: \'" + metaSensorProfile + "\'")
            self.inputImgType = metaSensorProfile

            topLevelDatasetIdent = root.find('Dataset_Identification')
            if topLevelDatasetIdent == None:
                raise ARCSIException("Cannot open top level section \'Dataset_Identification\'")

            topLevelDatasetContent = root.find('Dataset_Content')
            if topLevelDatasetContent == None:
                raise ARCSIException("Cannot open top level section \'Dataset_Content\'")

            topLevelProductInfo = root.find('Product_Information')
            if topLevelProductInfo == None:
                raise ARCSIException("Cannot open top level section \'Product_Information\'")

            topLevelCoordSystem = root.find('Coordinate_Reference_System')
            if topLevelCoordSystem == None:
                raise ARCSIException("Cannot open top level section \'Coordinate_Reference_System\'")

            topLevelGeoposition= root.find('Geoposition')
            if topLevelGeoposition == None:
                raise ARCSIException("Cannot open top level section \'Geoposition\'")

            topLevelProcessInfo = root.find('Processing_Information')
            if topLevelProcessInfo == None:
                raise ARCSIException("Cannot open top level section \'Processing_Information\'")

            topLevelRasterData = root.find('Raster_Data')
            if topLevelRasterData == None:
                raise ARCSIException("Cannot open top level section \'Raster_Data\'")

            topLevelRadiometricInfo = root.find('Radiometric_Data')
            if topLevelRadiometricInfo == None:
                raise ARCSIException("Cannot open top level section \'Radiometric_Data\'")

            topLevelGeometricInfo = root.find('Geometric_Data')
            if topLevelGeometricInfo == None:
                raise ARCSIException("Cannot open top level section \'Geometric_Data\'")

            topLevelQualityAssess = root.find('Quality_Assessment')
            if topLevelQualityAssess == None:
                raise ARCSIException("Cannot open top level section \'Quality_Assessment\'")

            topLevelQualityAssess = root.find('Quality_Assessment')
            if topLevelQualityAssess == None:
                raise ARCSIException("Cannot open top level section \'Quality_Assessment\'")

            topLevelDataSources = root.find('Dataset_Sources')
            if topLevelDataSources == None:
                raise ARCSIException("Cannot open top level section \'Dataset_Sources\'")
            
            # Get image file.
            filesDIR = os.path.dirname(inputHeader)
            if not self.userSpInputImage is None:
                self.fileName = os.path.abspath(self.userSpInputImage)
            else:
                if self.inputImgType == 'S6_SENSOR':
                    raise ARCSIException("You need to input an ortho-rectified product. Orthorectific / register your data and then use --imagefile to process this image file.")
                imgFileTags = topLevelRasterData.find('Data_Access').find('Data_Files')
                imgFiles = list()
                for child in imgFileTags:
                    if child.tag == 'Data_File':
                        imgFiles.append(os.path.join(filesDIR, imgFileTags.find('Data_File').find('DATA_FILE_PATH').attrib['href'].strip()))
                if len(imgFiles) == 0:
                    raise ARCSIException("Cannot find the input file.")
                elif len(imgFiles) == 1:
                    self.fileName = imgFiles[0]
                    self.tiledInputImg = False
                else:
                    self.inputImgFiles = imgFiles
                    self.tiledInputImg = True

            # Get acquasition time
            locGeoValsTag = None
            for child in topLevelGeometricInfo.find('Use_Area'):
                if child.tag == 'Located_Geometric_Values':
                    if child.find('LOCATION_TYPE').text.strip() == 'Center':
                        locGeoValsTag = child
                        break
            if locGeoValsTag == None:
                raise ARCSIException("Could not find \'Located_Geometric_Values\' for the centre of the image.")

            tmpAcquasitionTime = locGeoValsTag.find('TIME').text.strip()
            tmpAcquasitionTime = tmpAcquasitionTime.replace('Z', '')
            self.acquisitionTime = datetime.datetime.strptime(tmpAcquasitionTime, "%Y-%m-%dT%H:%M:%S.%f")

            # Get acquasition parameters - solar/view angles
            self.solarZenith = 90-float(locGeoValsTag.find('Solar_Incidences').find('SUN_ELEVATION').text.strip())
            self.solarAzimuth = float(locGeoValsTag.find('Solar_Incidences').find('SUN_AZIMUTH').text.strip())
            self.sensorZenith = 0.0 ## TODO: Not sure what to set this too!!
            self.sensorAzimuth = float(locGeoValsTag.find('Acquisition_Angles').find('AZIMUTH_ANGLE').text.strip())

            self.inputImgProjed = True
            if topLevelCoordSystem.find('Projected_CRS') == None:
                self.inputImgProjed = False
                if topLevelCoordSystem.find('Geodetic_CRS') == None:
                    raise ARCSIException("Could not find \'Projected_CRS\' or \'Geodetic_CRS\' therefore can't setup projection information.")
            else:
                self.inputImgProjed = True

            # Get coordinate system information
            epsgCode = -1
            if self.inputImgProjed:
                epsgCode = int(topLevelCoordSystem.find('Projected_CRS').find('PROJECTED_CRS_CODE').text.strip().split('::')[1])
            else:
                epsgCode = int(topLevelCoordSystem.find('Geodetic_CRS').find('GEODETIC_CRS_CODE').text.strip().split('::')[1])

            if epsgCode > 0:
                inProj = osr.SpatialReference()
                inProj.ImportFromEPSG(epsgCode)
                if self.inWKT == "":
                    self.inWKT = inProj.ExportToWkt()
            else:
                raise ARCSIException("EPSG code was not defined and therefore cannot defined the coordinate system and projection.")

            # Get BBOX
            tlTag = None
            trTag = None
            blTag = None
            brTag = None
            for child in topLevelDatasetContent.find('Dataset_Extent'):
                if child.tag == 'Vertex':
                    row = int(child.find('ROW').text.strip())
                    col = int(child.find('COL').text.strip())
                    if (row == 1) and (col == 1):
                        tlTag = child
                    elif (row == 1) and (col > 1):
                        trTag = child
                    elif (row > 1) and (col == 1):
                        blTag = child
                    elif (row > 1) and (col > 1):
                        brTag = child
                    else:
                        raise ARCSIException("Something strange has happened when parsing the header file, cannot find where vertex is located.")

            if tlTag == None:
                raise ARCSIException("Could not find the TL vertex tag.")
            if trTag == None:
                raise ARCSIException("Could not find the TR vertex tag.")
            if blTag == None:
                raise ARCSIException("Could not find the BL vertex tag.")
            if brTag == None:
                raise ARCSIException("Could not find the BR vertex tag.")

            self.latTL = float(tlTag.find('LAT').text.strip())
            self.lonTL = float(tlTag.find('LON').text.strip())
            self.latTR = float(trTag.find('LAT').text.strip())
            self.lonTR = float(trTag.find('LON').text.strip())
            self.latBL = float(blTag.find('LAT').text.strip())
            self.lonBL = float(blTag.find('LON').text.strip())
            self.latBR = float(brTag.find('LAT').text.strip())
            self.lonBR = float(brTag.find('LON').text.strip())
            self.latCentre = float(topLevelDatasetContent.find('Dataset_Extent').find('Center').find('LAT').text.strip())
            self.lonCentre = float(topLevelDatasetContent.find('Dataset_Extent').find('Center').find('LON').text.strip())

            if self.inputImgProjed:
                self.xTL = float(tlTag.find('X').text.strip())
                self.yTL = float(tlTag.find('Y').text.strip())
                self.xTR = float(trTag.find('X').text.strip())
                self.yTR = float(trTag.find('Y').text.strip())
                self.xBL = float(blTag.find('X').text.strip())
                self.yBL = float(blTag.find('Y').text.strip())
                self.xBR = float(brTag.find('X').text.strip())
                self.yBR = float(brTag.find('Y').text.strip())
                self.xCentre = float(topLevelDatasetContent.find('Dataset_Extent').find('Center').find('X').text.strip())
                self.yCentre = float(topLevelDatasetContent.find('Dataset_Extent').find('Center').find('Y').text.strip())
            
            # Find Radiometric Tags
            b0BandRadTag = None
            b1BandRadTag = None
            b2BandRadTag = None
            b3BandRadTag = None
            b0SolarIrrTag = None
            b1SolarIrrTag = None
            b2SolarIrrTag = None
            b3SolarIrrTag = None

            bandMeasTag = topLevelRadiometricInfo.find('Radiometric_Calibration').find('Instrument_Calibration').find('Band_Measurement_List')
            for child in bandMeasTag:
                if child.tag == 'Band_Radiance':
                    bandID = child.find('BAND_ID').text.strip()
                    if bandID == 'B0':
                        b0BandRadTag = child
                    elif bandID == 'B1':
                        b1BandRadTag = child
                    elif bandID == 'B2':
                        b2BandRadTag = child
                    elif bandID == 'B3':
                        b3BandRadTag = child

                elif child.tag == 'Band_Solar_Irradiance':
                    bandID = child.find('BAND_ID').text.strip()
                    if bandID == 'B0':
                        b0SolarIrrTag = child
                    elif bandID == 'B1':
                        b1SolarIrrTag = child
                    elif bandID == 'B2':
                        b2SolarIrrTag = child
                    elif bandID == 'B3':
                        b3SolarIrrTag = child

            # Get Solar Irradiance Values
            if b0SolarIrrTag == None:
                raise ARCSIException("Did not find B0 Solar Irradiance value")
            if b1SolarIrrTag == None:
                raise ARCSIException("Did not find B1 Solar Irradiance value")
            if b2SolarIrrTag == None:
                raise ARCSIException("Did not find B2 Solar Irradiance value")
            if b3SolarIrrTag == None:
                raise ARCSIException("Did not find B3 Solar Irradiance value")

            self.b0SolarIrr = float(b0SolarIrrTag.find('VALUE').text.strip())
            self.b1SolarIrr = float(b1SolarIrrTag.find('VALUE').text.strip())
            self.b2SolarIrr = float(b2SolarIrrTag.find('VALUE').text.strip())
            self.b3SolarIrr = float(b3SolarIrrTag.find('VALUE').text.strip())

            # Get band gain / bias values 
            if b0BandRadTag == None:
                raise ARCSIException("Did not find B0 radiance gain / bias")
            if b1BandRadTag == None:
                raise ARCSIException("Did not find B1 radiance gain / bias")
            if b2BandRadTag == None:
                raise ARCSIException("Did not find B2 radiance gain / bias")
            if b3BandRadTag == None:
                raise ARCSIException("Did not find B3 radiance gain / bias")

            self.b0RadGain = float(b0BandRadTag.find('GAIN').text.strip())
            self.b0RadBias = float(b0BandRadTag.find('BIAS').text.strip())
            self.b1RadGain = float(b1BandRadTag.find('GAIN').text.strip())
            self.b1RadBias = float(b1BandRadTag.find('BIAS').text.strip())
            self.b2RadGain = float(b2BandRadTag.find('GAIN').text.strip())
            self.b2RadBias = float(b2BandRadTag.find('BIAS').text.strip())
            self.b3RadGain = float(b3BandRadTag.find('GAIN').text.strip())
            self.b3RadBias = float(b3BandRadTag.find('BIAS').text.strip())

            # Find the bit ranges of the products and acquasition
            self.acqBitRange = int(topLevelRadiometricInfo.find('Dynamic_Range').find('ACQUISITION_RANGE').text.strip())
            self.prodBitRange = int(topLevelRadiometricInfo.find('Dynamic_Range').find('PRODUCT_RANGE').text.strip())

            rastDisTags = topLevelRasterData.find('Raster_Display')
            for child in rastDisTags:
                if child.tag == 'Special_Value':
                    if child.find('SPECIAL_VALUE_TEXT').text.strip() == 'NODATA':
                        self.inImgNoData = int(child.find('SPECIAL_VALUE_COUNT').text.strip())
                    elif child.find('SPECIAL_VALUE_TEXT').text.strip() == 'SATURATED':
                        self.inImgSatVal = int(child.find('SPECIAL_VALUE_COUNT').text.strip())

            for child in topLevelQualityAssess:
                if 'Cloud_Cotation' in child.find('MEASURE_NAME').text.strip():
                    self.cloudsMask = os.path.join(filesDIR, child.find('Quality_Mask').find('Component').find('COMPONENT_PATH').attrib['href'].strip())
                elif 'Area_Of_Interest' in child.find('MEASURE_NAME').text.strip():
                    self.roiMask = os.path.join(filesDIR, child.find('Quality_Mask').find('Component').find('COMPONENT_PATH').attrib['href'].strip())
            
            if not os.path.exists(self.cloudsMask):
                self.cloudsMask = None
            if not os.path.exists(self.roiMask):
                self.roiMask = None

        except Exception as e:
            raise e

    def getSolarIrrStdSolarGeom(self):
        """
        Get Solar Azimuth and Zenith as standard geometry.
        Azimuth: N=0, E=90, S=180, W=270.
        """
        return (self.solarAzimuth, self.solarZenith)

    def getSensorViewGeom(self):
        """
        Get sensor viewing angles
        returns (viewAzimuth, viewZenith)
        """
        return (self.sensorAzimuth, self.sensorZenith)

    def generateOutputBaseName(self):
        """
        Customises the generic name for the SPOT6 sensor
        """
        outname = self.defaultGenBaseOutFileName()
        if self.prodBitRange == 8:
            outname = outname + '_8bitDN'
        return outname

    def expectedImageDataPresent(self):
        imageDataPresent = True

        if not os.path.exists(self.fileName):
            imageDataPresent = False

        return imageDataPresent

    def applyImageDataMask(self, inputHeader, outputPath, outputMaskName, outputImgName, outFormat, outWKTFile):
        raise ARCSIException("SPOT6 does not provide any image masks, do not use the MASK option.")

    def imgNeedMosaicking(self):
        return self.tiledInputImg

    def mosaicImageTiles(self, outputPath):
        baseName = self.generateOutputBaseName()
        self.fileName = os.path.join(outputPath, baseName+'_mosic.kea')
        imageutils.createImageMosaic(self.inputImgFiles, self.fileName, self.inImgNoData, 0, 1, 0, 'KEA', rsgislib.imageutils.getRSGISLibDataType(self.inputImgFiles[0]))

    def resampleImgRes(self, outputPath, resampleToLowResImg, resampleMethod='cubic', multicore=False):
        raise ARCSIException("Image data does not need resampling")

    def sharpenLowResRadImgBands(self, inputImg, outputImage, outFormat):
        raise ARCSIException("Image sharpening is not available for this sensor.")

    def convertImageToRadiance(self, outputPath, outputReflName, outputThermalName, outFormat):
        print("Converting to Radiance")
        outputImage = os.path.join(outputPath, outputReflName)

        bandDefnSeq = list()
        spotBand = collections.namedtuple('Band', ['bandName', 'bandIndex', 'bias', 'gain'])
        bandDefnSeq.append(spotBand(bandName="Blue", bandIndex=3, bias=self.b2RadBias, gain=self.b2RadGain))
        bandDefnSeq.append(spotBand(bandName="Green", bandIndex=2, bias=self.b1RadBias, gain=self.b1RadGain))
        bandDefnSeq.append(spotBand(bandName="Red", bandIndex=1, bias=self.b0RadBias, gain=self.b0RadGain))
        bandDefnSeq.append(spotBand(bandName="NIR", bandIndex=4, bias=self.b3RadBias, gain=self.b3RadGain))
        rsgislib.imagecalibration.spot5ToRadiance(self.fileName, outputImage, outFormat, bandDefnSeq)

        return outputImage, None

    def generateImageSaturationMask(self, outputPath, outputName, outFormat):
        print("Generate Saturation Image")
        outputImage = os.path.join(outputPath, outputName)

        spotBand = collections.namedtuple('Band', ['bandName', 'fileName', 'bandIndex', 'satVal'])
        bandDefnSeq = list()
        self.inImgSatVal = float(self.inImgSatVal)
        bandDefnSeq.append(spotBand(bandName="Blue", fileName=self.fileName, bandIndex=3, satVal=self.inImgSatVal))
        bandDefnSeq.append(spotBand(bandName="Green", fileName=self.fileName, bandIndex=2, satVal=self.inImgSatVal))
        bandDefnSeq.append(spotBand(bandName="Red", fileName=self.fileName, bandIndex=1, satVal=self.inImgSatVal))
        bandDefnSeq.append(spotBand(bandName="NIR", fileName=self.fileName, bandIndex=4, satVal=self.inImgSatVal))

        rsgislib.imagecalibration.saturatedPixelsMask(outputImage, outFormat, bandDefnSeq)

        return outputImage

    def generateValidImageDataMask(self, outputPath, outputMaskName, viewAngleImg, outFormat):
        print("Generate valid image mask")
        # Test input image, if it is .jp2 then warp to KEA so processing can proceed.
        # Also set whether input file has projection correctly defined. If not, copy to KEA to proceed.
        arcsiUtils = ARCSIUtils()
        lclWktStr = arcsiUtils.getWKTProjFromImage(self.fileName)
        inFileExt = os.path.splitext(self.fileName)[1]
        if inFileExt.lower() == '.jp2':
            inFileName = os.path.splitext(os.path.basename(self.fileName))[0]
            self.copiedKEADNImg = os.path.join(outputPath, inFileName+'.kea')
            # Check input projection is defined if not define from header information.
            projCmd = ''
            createdWKTFile = False
            wktFileName = ''
            if (lclWktStr == '') or (lclWktStr == None):
                uidStr = arcsiUtils.uidGenerator()
                wktFileName = os.path.join(outputPath, uidStr+'_wkt.wkt')
                arcsiUtils.writeText2File(self.inWKT, wktFileName)
                createdWKTFile = True
                projCmd = ' -s_srs ' + wktFileName + ' -t_srs ' + wktFileName
            cmd = 'gdalwarp -of KEA -overwrite -r cubic ' + projCmd + ' -srcnodata ' + str(self.inImgNoData) + ' -dstnodata ' + str(self.inImgNoData) + ' "' + self.fileName + '" "' + self.copiedKEADNImg + '"'
            try:
                subprocess.call(cmd, shell=True)
            except OSError as e:
                raise ARCSIException('Could not warp image: ' + cmd)
            self.origDNImg = self.fileName
            self.fileName = self.copiedKEADNImg
            self.createdCopyKEADNImg = True
            if createdWKTFile:
                os.remove(wktFileName)
        elif (lclWktStr == '') or (lclWktStr == None):
            # Copy input file and defined projection if not define from header information.
            inFileName = os.path.splitext(os.path.basename(self.fileName))[0]
            self.copiedKEADNImg = os.path.join(outputPath, inFileName+'.kea')
            uidStr = arcsiUtils.uidGenerator()
            wktFileName = os.path.join(outputPath, uidStr+'_wkt.wkt')
            arcsiUtils.writeText2File(self.inWKT, wktFileName)
            cmd = 'gdal_translate -of KEA -a_nodata ' + str(self.inImgNoData) + ' -a_srs ' + wktFileName + ' "' + self.fileName + '" "' + self.copiedKEADNImg + '"'
            try:
                subprocess.call(cmd, shell=True)
            except OSError as e:
                raise ARCSIException('Could not translate image: ' + cmd)
            os.remove(wktFileName)
            self.origDNImg = self.fileName
            self.fileName = self.copiedKEADNImg
            self.createdCopyKEADNImg = True

        outputImage = os.path.join(outputPath, outputMaskName)
        rsgislib.imageutils.genValidMask(inimages=[self.fileName], outimage=outputImage, gdalformat=outFormat, nodata=self.inImgNoData)
        return outputImage

    def convertThermalToBrightness(self, inputRadImage, outputPath, outputName, outFormat, scaleFactor):
        raise ARCSIException("There are no thermal bands...")

    def convertImageToTOARefl(self, inputRadImage, outputPath, outputName, outFormat, scaleFactor):
        print("Converting to TOA")
        outputImage = os.path.join(outputPath, outputName)
        solarIrradianceVals = list()
        IrrVal = collections.namedtuple('SolarIrradiance', ['irradiance'])
        solarIrradianceVals.append(IrrVal(irradiance=self.b2SolarIrr))
        solarIrradianceVals.append(IrrVal(irradiance=self.b1SolarIrr))
        solarIrradianceVals.append(IrrVal(irradiance=self.b0SolarIrr))
        solarIrradianceVals.append(IrrVal(irradiance=self.b3SolarIrr))

        rsgislib.imagecalibration.radiance2TOARefl(inputRadImage, outputImage, outFormat, rsgislib.TYPE_16UINT, scaleFactor, self.acquisitionTime.year, self.acquisitionTime.month, self.acquisitionTime.day, self.solarZenith, solarIrradianceVals)
        return outputImage

    def generateCloudMask(self, inputReflImage, inputSatImage, inputThermalImage, inputViewAngleImg, inputValidImg, outputPath, outputName, outFormat, tmpPath, scaleFactor, cloud_msk_methods=None):
        outputCloudMaskImg = os.path.join(outputPath, outputName)
        if self.cloudsMask == None:
            rsgislib.imageutils.createCopyImage(inputReflImage, outputCloudMaskImg, 1, 0, outFormat, rsgislib.TYPE_8UINT)
        else:
            haveVecClouds = False
            tmpVectorDS = ogr.Open(self.cloudsMask)
            if tmpVectorDS != None:
                inVectorLayer = tmpVectorDS.GetLayer()
                if inVectorLayer != None:
                    if inVectorLayer.GetFeatureCount() > 0:
                        haveVecClouds = True
                    else:
                        haveVecClouds = False
                else:
                    haveVecClouds = False
            else:
                haveVecClouds = False
            if haveVecClouds:
                rsgislib.vectorutils.rasterise2Image(self.cloudsMask, inputReflImage, outputCloudMaskImg, gdalformat=outFormat, burnVal=1)
                rsgislib.rastergis.populateStats(clumps=outputCloudMaskImg, addclrtab=True, calcpyramids=True, ignorezero=True, ratband=1)
            else:
                rsgislib.imageutils.createCopyImage(inputReflImage, outputCloudMaskImg, 1, 0, outFormat, rsgislib.TYPE_8UINT)
        return outputCloudMaskImg

    def createCloudMaskDataArray(self, inImgDataArr):
        return inImgDataArr

    def defineDarkShadowImageBand(self):
        return 4

    def calc6SCoefficients(self, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF):
        sixsCoeffs = numpy.zeros((4, 6), dtype=numpy.float32)
        # Set up 6S model
        s = Py6S.SixS()
        s.atmos_profile = atmosProfile
        s.aero_profile = aeroProfile
        s.ground_reflectance = grdRefl
        s.geometry = Py6S.Geometry.User()
        s.geometry.solar_z = self.solarZenith
        s.geometry.solar_a = self.solarAzimuth
        s.geometry.view_z = self.sensorZenith
        s.geometry.view_a = self.sensorAzimuth
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

        # Band 1 - Blue
        s.wavelength = Py6S.Wavelength(0.4, 1, [0.000066,0.000060,0.000073,0.000070,0.000089,0.000064,0.000083,0.000055,0.000091,0.000139,0.000194,0.000330,0.001204,0.001908,0.003939,0.006938,0.018163,0.028061,0.067956,0.131090,0.264885,0.415586,0.614661,0.752286,0.857545,0.854075,0.853217,0.853871,0.875585,0.909005,0.933068,0.936660,0.933796,0.924426,0.920904,0.934558,0.963645,0.974364,0.966360,0.953720,0.958189,0.980703,1.000000,0.998082,0.984045,0.938674,0.833872,0.698591,0.431257,0.273925,0.104595,0.056546,0.027936,0.017517,0.008147,0.005309,0.003332,0.002576,0.001723,0.001304,0.000869,0.000679,0.000519,0.000474,0.000425,0.000388,0.000297,0.000245,0.000213,0.000211,0.000213,0.000210,0.000191,0.000172,0.000154,0.000150,0.000152,0.000153,0.000148,0.000141,0.000127,0.000116,0.000108,0.000102,0.000099,0.000099,0.000097,0.000098,0.000096,0.000096,0.000095,0.000095,0.000095,0.000096,0.000096,0.000095,0.000096,0.000095,0.000096,0.000098,0.000099,0.000102,0.000107,0.000111,0.000114,0.000118,0.000127,0.000135,0.000151,0.000157,0.000152,0.000131,0.000107,0.000099,0.000099,0.000108,0.000128,0.000141,0.000154,0.000152,0.000146,0.000132,0.000123,0.000108,0.000103,0.000096,0.000095,0.000092,0.000087,0.000087,0.000084,0.000085,0.000091,0.000091,0.000094,0.000085,0.000077,0.000076,0.000074,0.000073,0.000074,0.000071,0.000072,0.000071,0.000071,0.000071,0.000070,0.000069,0.000071,0.000072,0.000068,0.000067,0.000067,0.000066,0.000066,0.000066,0.000064,0.000063,0.000064,0.000061,0.000063,0.000062,0.000060,0.000060,0.000061,0.000060,0.000056,0.000056,0.000056,0.000055,0.000054,0.000055,0.000053,0.000055,0.000051,0.000050,0.000052,0.000049,0.000064,0.000069,0.000073,0.000072,0.000066,0.000071,0.000065,0.000069,0.000069,0.000074,0.000081,0.000088,0.000110,0.000128,0.000123,0.000186,0.000180,0.000226,0.000229,0.000219,0.000202,0.000190,0.000182,0.000144,0.000122,0.000105,0.000098,0.000089,0.000073,0.000063,0.000054,0.000052,0.000049,0.000046,0.000041,0.000041,0.000044,0.000041,0.000042,0.000041,0.000036,0.000032,0.000031,0.000030,0.000029,0.000028,0.000027,0.000026,0.000021,0.000019,0.000019,0.000016,0.000015,0.000015,0.000013,0.000011,0.000011,0.000009,0.000009,0.000008,0.000009,0.000009,0.000007])
        s.run()
        sixsCoeffs[0,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[0,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[0,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[0,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[0,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[0,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 2 - Green
        s.wavelength = Py6S.Wavelength(0.4, 1, [0.000010,0.000019,0.000029,0.000022,0.000067,0.000012,0.000044,0.000007,0.000005,0.000018,0.000000,0.000013,0.000020,0.000019,0.000001,0.000018,0.000019,0.000018,0.000013,0.000019,0.000019,0.000009,0.000012,0.000013,0.000009,0.000010,0.000011,0.000012,0.000014,0.000015,0.000018,0.000022,0.000032,0.000042,0.000061,0.000076,0.000105,0.000134,0.000202,0.000263,0.000380,0.000510,0.000979,0.001634,0.003279,0.005264,0.012676,0.025850,0.070359,0.130948,0.307512,0.474368,0.746747,0.897130,1.005725,1.022418,1.021159,1.017318,1.009553,1.004075,0.999629,1.001409,1.009106,1.010996,0.999011,0.984690,0.970225,0.970041,0.980672,0.992191,0.999125,0.990362,0.963916,0.910970,0.723618,0.576231,0.298406,0.147397,0.051017,0.028039,0.011787,0.006346,0.002507,0.001475,0.000808,0.000602,0.000419,0.000332,0.000230,0.000179,0.000129,0.000111,0.000099,0.000097,0.000096,0.000095,0.000092,0.000086,0.000083,0.000083,0.000092,0.000101,0.000113,0.000118,0.000112,0.000104,0.000095,0.000087,0.000075,0.000069,0.000064,0.000059,0.000053,0.000051,0.000051,0.000051,0.000050,0.000048,0.000047,0.000046,0.000045,0.000046,0.000047,0.000049,0.000050,0.000050,0.000049,0.000047,0.000053,0.000052,0.000071,0.000056,0.000040,0.000026,0.000006,0.000003,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000002,0.000001,0.000002,0.000005,0.000004,0.000010,0.000009,0.000006,0.000001,0.000001,0.000002,0.000002,0.000001,0.000001,0.000002,0.000001,0.000001,0.000001,0.000002,0.000002,0.000003,0.000003,0.000002,0.000002,0.000001,0.000001,0.000002,0.000002,0.000001,0.000002,0.000001,0.000002,0.000003,0.000003,0.000004,0.000006,0.000011,0.000010,0.000015,0.000014,0.000009,0.000014,0.000011,0.000010,0.000005,0.000006,0.000011,0.000011,0.000005,0.000005,0.000003,0.000013,0.000012,0.000012,0.000005,0.000009,0.000009,0.000016,0.000023,0.000024,0.000023,0.000020,0.000016,0.000014,0.000007,0.000008,0.000008,0.000012,0.000012,0.000014,0.000022,0.000024,0.000020,0.000015,0.000008,0.000005,0.000005,0.000005,0.000004,0.000003,0.000002,0.000002,0.000002,0.000002,0.000001,0.000001,0.000002,0.000002,0.000001,0.000001,0.000002,0.000001,0.000001,0.000001,0.000001,0.000001,0.000000,0.000001,0.000000])
        s.run()
        sixsCoeffs[1,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[1,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[1,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[1,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[1,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[1,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 3 - Red
        s.wavelength = Py6S.Wavelength(0.4, 1, [0.000012,0.000013,0.000014,0.000017,0.000024,0.000018,0.000013,0.000011,0.000004,0.000003,0.000005,0.000006,0.000007,0.000016,0.000017,0.000018,0.000031,0.000035,0.000048,0.000045,0.000058,0.000058,0.000036,0.000038,0.000040,0.000037,0.000034,0.000033,0.000045,0.000057,0.000057,0.000050,0.000060,0.000058,0.000069,0.000051,0.000020,0.000014,0.000008,0.000006,0.000004,0.000004,0.000004,0.000004,0.000004,0.000004,0.000004,0.000004,0.000004,0.000005,0.000006,0.000007,0.000008,0.000007,0.000008,0.000007,0.000008,0.000008,0.000008,0.000009,0.000010,0.000010,0.000010,0.000011,0.000013,0.000015,0.000018,0.000020,0.000028,0.000035,0.000043,0.000047,0.000051,0.000061,0.000104,0.000172,0.000297,0.000381,0.000491,0.000572,0.000884,0.001384,0.003100,0.005399,0.012507,0.019719,0.038710,0.064051,0.163370,0.307780,0.546333,0.717654,0.926693,0.972462,0.995125,1.000713,0.966180,0.929119,0.898687,0.900483,0.925482,0.951093,0.988735,1.007842,1.013651,1.000945,0.967365,0.944030,0.919796,0.913455,0.921124,0.933433,0.946550,0.936955,0.876512,0.816261,0.698532,0.609571,0.471508,0.377067,0.235139,0.155846,0.077590,0.047085,0.020647,0.011314,0.004546,0.002585,0.001250,0.000850,0.000550,0.000472,0.000403,0.000357,0.000298,0.000212,0.000125,0.000064,0.000024,0.000017,0.000013,0.000013,0.000096,0.000115,0.000128,0.000243,0.000127,0.000010,0.000005,0.000005,0.000004,0.000004,0.000005,0.000005,0.000004,0.000004,0.000003,0.000003,0.000003,0.000002,0.000003,0.000003,0.000003,0.000003,0.000004,0.000004,0.000004,0.000003,0.000004,0.000003,0.000005,0.000006,0.000005,0.000005,0.000004,0.000005,0.000005,0.000005,0.000004,0.000005,0.000006,0.000004,0.000006,0.000009,0.000012,0.000010,0.000012,0.000014,0.000011,0.000009,0.000013,0.000013,0.000012,0.000007,0.000014,0.000009,0.000015,0.000010,0.000021,0.000032,0.000057,0.000071,0.000074,0.000089,0.000097,0.000087,0.000067,0.000053,0.000045,0.000038,0.000025,0.000020,0.000015,0.000009,0.000009,0.000007,0.000005,0.000007,0.000007,0.000007,0.000006,0.000007,0.000008,0.000008,0.000005,0.000004,0.000002,0.000002,0.000005,0.000004,0.000002,0.000004,0.000003,0.000003,0.000003,0.000002,0.000001,0.000002,0.000001,0.000001,0.000002])
        s.run()
        sixsCoeffs[2,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[2,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[2,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[2,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[2,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[2,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 4 - NIR
        s.wavelength = Py6S.Wavelength(0.4, 1, [0.000002,0.000004,0.000005,0.000005,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000002,0.000002,0.000002,0.000005,0.000626,0.000848,0.000988,0.000035,0.000024,0.000037,0.000078,0.000108,0.000110,0.000050,0.000031,0.000021,0.000020,0.000031,0.000048,0.000080,0.000088,0.000070,0.000047,0.000023,0.000016,0.000017,0.000025,0.000036,0.000036,0.000028,0.000024,0.000031,0.000046,0.000103,0.000143,0.000177,0.000171,0.000166,0.000182,0.000194,0.000189,0.000164,0.000156,0.000161,0.000168,0.000172,0.000171,0.000168,0.000174,0.000204,0.000232,0.000262,0.000263,0.000275,0.000288,0.000339,0.000396,0.000538,0.000603,0.000590,0.000478,0.000252,0.000154,0.000087,0.000067,0.000050,0.000043,0.000033,0.000029,0.000023,0.000020,0.000017,0.000016,0.000014,0.000013,0.000013,0.000012,0.000012,0.000011,0.000012,0.000012,0.000011,0.000011,0.000011,0.000010,0.000011,0.000012,0.000012,0.000012,0.000012,0.000012,0.000012,0.000012,0.000013,0.000014,0.000016,0.000017,0.000020,0.000023,0.000028,0.000033,0.000041,0.000046,0.000058,0.000067,0.000085,0.000100,0.000130,0.000158,0.000217,0.000271,0.000397,0.000502,0.000724,0.000928,0.001371,0.001803,0.002803,0.003848,0.006455,0.009364,0.017006,0.025869,0.049567,0.076803,0.146158,0.231787,0.387710,0.513703,0.708353,0.822525,0.943603,0.984475,1.000000,0.995349,0.983444,0.976833,0.969379,0.964947,0.955247,0.947252,0.932992,0.922617,0.908670,0.900317,0.890406,0.884468,0.874879,0.866655,0.852033,0.840679,0.822377,0.810672,0.795384,0.785933,0.775054,0.768346,0.759152,0.752888,0.743113,0.735597,0.723345,0.714248,0.698655,0.687070,0.669891,0.659230,0.645375,0.638203,0.629921,0.624456,0.615019,0.607448,0.594257,0.583815,0.556993,0.522029,0.442705,0.378439,0.267302,0.189747,0.094132,0.054574,0.024255,0.014703,0.007443,0.005014,0.002936,0.002158,0.001430,0.001129,0.000835,0.000699,0.000550,0.000478,0.000395,0.000337,0.000277,0.000231,0.000184,0.000145,0.000110,0.000092,0.000064,0.000052,0.000024,0.000017,0.000013,0.000011,0.000013,0.000008,0.000005,0.000004,0.000003,0.000002,0.000003,0.000003,0.000002,0.000002,0.000003,0.000002,0.000001,0.000001,0.000001,0.000001,0.000001])
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
        # Set up 6S model
        s = Py6S.SixS()
        s.atmos_profile = atmosProfile
        s.aero_profile = aeroProfile
        s.ground_reflectance = grdRefl
        s.geometry = Py6S.Geometry.User()
        s.geometry.solar_z = self.solarZenith
        s.geometry.solar_a = self.solarAzimuth
        s.geometry.view_z = self.sensorZenith
        s.geometry.view_a = self.sensorAzimuth
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

        # Band Blue
        s.wavelength = Py6S.Wavelength(0.4, 1, [0.000066,0.000060,0.000073,0.000070,0.000089,0.000064,0.000083,0.000055,0.000091,0.000139,0.000194,0.000330,0.001204,0.001908,0.003939,0.006938,0.018163,0.028061,0.067956,0.131090,0.264885,0.415586,0.614661,0.752286,0.857545,0.854075,0.853217,0.853871,0.875585,0.909005,0.933068,0.936660,0.933796,0.924426,0.920904,0.934558,0.963645,0.974364,0.966360,0.953720,0.958189,0.980703,1.000000,0.998082,0.984045,0.938674,0.833872,0.698591,0.431257,0.273925,0.104595,0.056546,0.027936,0.017517,0.008147,0.005309,0.003332,0.002576,0.001723,0.001304,0.000869,0.000679,0.000519,0.000474,0.000425,0.000388,0.000297,0.000245,0.000213,0.000211,0.000213,0.000210,0.000191,0.000172,0.000154,0.000150,0.000152,0.000153,0.000148,0.000141,0.000127,0.000116,0.000108,0.000102,0.000099,0.000099,0.000097,0.000098,0.000096,0.000096,0.000095,0.000095,0.000095,0.000096,0.000096,0.000095,0.000096,0.000095,0.000096,0.000098,0.000099,0.000102,0.000107,0.000111,0.000114,0.000118,0.000127,0.000135,0.000151,0.000157,0.000152,0.000131,0.000107,0.000099,0.000099,0.000108,0.000128,0.000141,0.000154,0.000152,0.000146,0.000132,0.000123,0.000108,0.000103,0.000096,0.000095,0.000092,0.000087,0.000087,0.000084,0.000085,0.000091,0.000091,0.000094,0.000085,0.000077,0.000076,0.000074,0.000073,0.000074,0.000071,0.000072,0.000071,0.000071,0.000071,0.000070,0.000069,0.000071,0.000072,0.000068,0.000067,0.000067,0.000066,0.000066,0.000066,0.000064,0.000063,0.000064,0.000061,0.000063,0.000062,0.000060,0.000060,0.000061,0.000060,0.000056,0.000056,0.000056,0.000055,0.000054,0.000055,0.000053,0.000055,0.000051,0.000050,0.000052,0.000049,0.000064,0.000069,0.000073,0.000072,0.000066,0.000071,0.000065,0.000069,0.000069,0.000074,0.000081,0.000088,0.000110,0.000128,0.000123,0.000186,0.000180,0.000226,0.000229,0.000219,0.000202,0.000190,0.000182,0.000144,0.000122,0.000105,0.000098,0.000089,0.000073,0.000063,0.000054,0.000052,0.000049,0.000046,0.000041,0.000041,0.000044,0.000041,0.000042,0.000041,0.000036,0.000032,0.000031,0.000030,0.000029,0.000028,0.000027,0.000026,0.000021,0.000019,0.000019,0.000016,0.000015,0.000015,0.000013,0.000011,0.000011,0.000009,0.000009,0.000008,0.000009,0.000009,0.000007])
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

    def estimateImageToAODUsingDDV(self, inputRADImage, inputTOAImage, inputDEMFile, shadowMask, outputPath, outputName, outFormat, tmpPath, aeroProfile, atmosProfile, grdRefl, aotValMin, aotValMax):
        print("Not implemented\n")
        sys.exit()

    def estimateImageToAODUsingDOS(self, inputRADImage, inputTOAImage, inputDEMFile, shadowMask, outputPath, outputName, outFormat, tmpPath, aeroProfile, atmosProfile, grdRefl, aotValMin, aotValMax, globalDOS, simpleDOS, dosOutRefl):
        try:
            print("Estimating AOD Using DOS")
            arcsiUtils = ARCSIUtils()
            outputAOTImage = os.path.join(outputPath, outputName)
            tmpBaseName = os.path.splitext(outputName)[0]
            imgExtension = arcsiUtils.getFileExtension(outFormat)

            # Define Band Numbers
            blueBand = 1
            redBand = 3
            nirNamd = 4
            segBands = [1,2,3,4]

            dosBlueImage = ""
            minObjSize = 3
            darkPxlPercentile = 0.01
            blockSize = 1000
            if simpleDOS:
                outputDOSBlueName = tmpBaseName + "DOSBlue" + imgExtension
                dosBlueImage, bandOff = self.convertImageBandToReflectanceSimpleDarkSubtract(inputTOAImage, outputPath, outputDOSBlueName, outFormat, dosOutRefl, blueBand)
            elif globalDOS:
                dosBlueImage = self.performDOSOnSingleBand(inputTOAImage, blueBand, outputPath, tmpBaseName, "Blue", "KEA", tmpPath, minObjSize, darkPxlPercentile, dosOutRefl)
            else:
                dosBlueImage = self.performLocalDOSOnSingleBand(inputTOAImage, blueBand, outputPath, tmpBaseName, "Blue", "KEA", tmpPath, minObjSize, darkPxlPercentile, blockSize, dosOutRefl)

                
            thresImageClumpsFinal = os.path.join(tmpPath, tmpBaseName + "_clumps" + imgExtension)
            rsgislib.segmentation.segutils.runShepherdSegmentation(inputTOAImage, thresImageClumpsFinal, tmpath=tmpPath, gdalformat="KEA", numClusters=40, minPxls=10, bands=segBands, processInMem=True)

            stats2CalcTOA = list()
            stats2CalcTOA.append(rsgislib.rastergis.BandAttStats(band=1, meanField="MeanElev"))
            rsgislib.rastergis.populateRATWithStats(inputDEMFile, thresImageClumpsFinal, stats2CalcTOA)

            stats2CalcTOA = list()
            stats2CalcTOA.append(rsgislib.rastergis.BandAttStats(band=blueBand, meanField="MeanBlueDOS"))
            rsgislib.rastergis.populateRATWithStats(dosBlueImage, thresImageClumpsFinal, stats2CalcTOA)

            stats2CalcRad = list()
            stats2CalcRad.append(rsgislib.rastergis.BandAttStats(band=blueBand, meanField="MeanBlueRAD"))
            stats2CalcRad.append(rsgislib.rastergis.BandAttStats(band=nirNamd, meanField="MeanNIRRAD"))
            stats2CalcRad.append(rsgislib.rastergis.BandAttStats(band=redBand, meanField="MeanRedRAD"))
            rsgislib.rastergis.populateRATWithStats(inputRADImage, thresImageClumpsFinal, stats2CalcRad)

            ratDS = gdal.Open(thresImageClumpsFinal, gdal.GA_Update)
            Histogram = rat.readColumn(ratDS, "Histogram")
            MeanElev = rat.readColumn(ratDS, "MeanElev")

            MeanNIRRAD = rat.readColumn(ratDS, "MeanNIRRAD")
            MeanRedRAD = rat.readColumn(ratDS, "MeanRedRAD")

            radNDVI = (MeanNIRRAD - MeanRedRAD)/(MeanNIRRAD + MeanRedRAD)

            selected = Histogram * 2
            selected[...] = 0
            selected[radNDVI>0.2] = 1
            rat.writeColumn(ratDS, "Selected", selected)
            ratDS = None

            rsgislib.rastergis.spatialLocation(thresImageClumpsFinal, "Eastings", "Northings")
            rsgislib.rastergis.selectClumpsOnGrid(thresImageClumpsFinal, "Selected", "PredictAOTFor", "Eastings", "Northings", "MeanB1DOS", "min", 20, 20)

            ratDS = gdal.Open(thresImageClumpsFinal, gdal.GA_Update)
            MeanBlueDOS = rat.readColumn(ratDS, "MeanBlueDOS")
            MeanBlueDOS = MeanBlueDOS / 1000
            MeanBlueRAD = rat.readColumn(ratDS, "MeanBlueRAD")
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
                        cDist = self.run6SToOptimiseAODValue(cAOT, MeanBlueRAD[i], MeanBlueDOS[i], aeroProfile, atmosProfile, grdRefl, MeanElev[i]/1000)
                        if j == 0:
                            minAOT = cAOT
                            minDist = cDist
                        elif cDist < minDist:
                            minAOT = cAOT
                            minDist = cDist
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
            blueBand = 1
            return self.estimateSingleAOTFromDOSBandImpl(radianceImage, toaImage, inputDEMFile, tmpPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, minAOT, maxAOT, dosOutRefl, blueBand)
        except Exception as e:
            raise

    def setBandNames(self, imageFile):
        dataset = gdal.Open(imageFile, gdal.GA_Update)
        if not dataset is None:
            dataset.GetRasterBand(1).SetDescription("Blue")
            dataset.GetRasterBand(2).SetDescription("Green")
            dataset.GetRasterBand(3).SetDescription("Red")
            dataset.GetRasterBand(4).SetDescription("NIR")
            dataset = None
        else:
            print("Could not open image to set band names: ", imageFile)

    def cleanLocalFollowProcessing(self):
        if self.createdCopyKEADNImg:
            if not self.debugMode:
                gdalDriver = gdal.GetDriverByName('KEA')
                gdalDriver.Delete(self.copiedKEADNImg)
            self.fileName = self.origDNImg
        if self.tiledInputImg:
            if not self.debugMode:
                gdalDriver = gdal.GetDriverByName('KEA')
                gdalDriver.Delete(self.fileName)



