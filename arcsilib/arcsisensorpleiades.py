"""
Module that contains the ARCSIPleiadesSensor class.
"""
############################################################################
#  arcsisensorpleiades.py
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
# Purpose:  A class for read the Pleiades sensor header file and applying
#           the pre-processing operations within ARCSI to the Pleiades
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
# Import the RIOS RAT module
from rios import rat
# Import the subprocess module
import subprocess

class ARCSIPleiadesSensor (ARCSIAbstractSensor):
    """
    A class which represents the Pleiades sensor to read
    header parameters and apply data processing operations.
    """
    def __init__(self, debugMode, inputImage):
        ARCSIAbstractSensor.__init__(self, debugMode, inputImage)
        self.sensor = "PLEIADES"
        self.fileName = ""
        self.inputImgFiles = []
        self.tiledInputImg = False
        self.inputImgType = ''
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
        self.inImgNoData = 0.0
        self.inImgSatVal = 4095
        self.b0DynAdjBias = 0.0
        self.b0DynAdjSlope = 0.0
        self.b0DynAdjMinThres = 0.0
        self.b0DynAdjMaxThres = 0.0
        self.b1DynAdjBias = 0.0
        self.b1DynAdjSlope = 0.0
        self.b1DynAdjMinThres = 0.0
        self.b1DynAdjMaxThres = 0.0
        self.b2DynAdjBias = 0.0
        self.b2DynAdjSlope = 0.0
        self.b2DynAdjMinThres = 0.0
        self.b2DynAdjMaxThres = 0.0
        self.b3DynAdjBias = 0.0
        self.b3DynAdjSlope = 0.0
        self.b3DynAdjMinThres = 0.0
        self.b3DynAdjMaxThres = 0.0
        self.copiedKEADNImg = ''
        self.origDNImg = ''
        self.createdCopyKEADNImg = False
        self.pleiadesSat = ''
        self.cloudsMask = None
        self.roiMask = None

    def extractHeaderParameters(self, inputHeader, wktStr):
        """
        Understands and parses the Pleiades xml header file
        """
        try:
            self.headerFileName = os.path.split(inputHeader)[1]
            
            print("Reading header file")
            tree = ET.parse(inputHeader)
            root = tree.getroot()
            
            topLevelMetaIdent = root.find('Metadata_Identification')
            if topLevelMetaIdent == None:
                raise ARCSIException("Cannot open top level section \'Metadata_Identification\'")

            dimapVersion = topLevelMetaIdent.find('METADATA_FORMAT').attrib['version'].strip()
            if dimapVersion.split('.')[0] != '2':
                raise ARCSIException("Only DIMAP Version 2.X is supported by this reader; provided with: \'" + dimapVersion + "\'")

            metaSensorProfile = topLevelMetaIdent.find('METADATA_PROFILE').text.strip()
            if not ((metaSensorProfile == 'PHR_SENSOR') or (metaSensorProfile == 'PHR_ORTHO')):
                raise ARCSIException("Input file is not for Pleiades, \'METADATA_PROFILE\' should be \'PHR_SENSOR\' or \'PHR_ORTHO\'; provided with: \'" + metaSensorProfile + "\'")
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
                if self.inputImgType == 'PHR_SENSOR':
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
                epsgCode = int(topLevelCoordSystem.find('Projected_CRS').find('PROJECTED_CRS_NAME').text.strip())
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
                    print(child.find('BAND_ID').text)
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

            # If 8 bit product then read dynamic adjustment
            if self.prodBitRange == 8:
                b0BandAdjTag = None
                b1BandAdjTag = None
                b2BandAdjTag = None
                b3BandAdjTag = None
                bandAdjTags = topLevelRadiometricInfo.find('Dynamic_Adjustment').find('Band_Adjustment_List')
                for child in bandAdjTags:
                    if child.tag == 'Band_Adjustment':
                        bandID = child.find('BAND_ID').text.strip()
                        if bandID == 'B0':
                            b0BandAdjTag = child
                        elif bandID == 'B1':
                            b1BandAdjTag = child
                        elif bandID == 'B2':
                            b2BandAdjTag = child
                        elif bandID == 'B3':
                            b3BandAdjTag = child

                if b0BandAdjTag == None:
                    raise ARCSIException("Did not find B0 band adjustment tag")
                if b1BandAdjTag == None:
                    raise ARCSIException("Did not find B1 band adjustment tag")
                if b2BandAdjTag == None:
                    raise ARCSIException("Did not find B2 band adjustment tag")
                if b3BandAdjTag == None:
                    raise ARCSIException("Did not find B3 band adjustment tag")

                self.b0DynAdjBias = float(b0BandAdjTag.find('BIAS').text.strip())
                self.b0DynAdjSlope = float(b0BandAdjTag.find('SLOPE').text.strip())
                self.b0DynAdjMinThres = float(b0BandAdjTag.find('MIN_THRESHOLD').text.strip())
                self.b0DynAdjMaxThres = float(b0BandAdjTag.find('MAX_THRESHOLD').text.strip())
                self.b1DynAdjBias = float(b1BandAdjTag.find('BIAS').text.strip())
                self.b1DynAdjSlope = float(b1BandAdjTag.find('SLOPE').text.strip())
                self.b1DynAdjMinThres = float(b1BandAdjTag.find('MIN_THRESHOLD').text.strip())
                self.b1DynAdjMaxThres = float(b1BandAdjTag.find('MAX_THRESHOLD').text.strip())
                self.b2DynAdjBias = float(b2BandAdjTag.find('BIAS').text.strip())
                self.b2DynAdjSlope = float(b2BandAdjTag.find('SLOPE').text.strip())
                self.b2DynAdjMinThres = float(b2BandAdjTag.find('MIN_THRESHOLD').text.strip())
                self.b2DynAdjMaxThres = float(b2BandAdjTag.find('MAX_THRESHOLD').text.strip())
                self.b3DynAdjBias = float(b3BandAdjTag.find('BIAS').text.strip())
                self.b3DynAdjSlope = float(b3BandAdjTag.find('SLOPE').text.strip())
                self.b3DynAdjMinThres = float(b3BandAdjTag.find('MIN_THRESHOLD').text.strip())
                self.b3DynAdjMaxThres = float(b3BandAdjTag.find('MAX_THRESHOLD').text.strip())

            # Find the no data and saturation pixel values.
            self.inImgNoData = 0.0
            self.inImgSatVal = 4095

            rastDisTags = topLevelRasterData.find('Raster_Display')
            for child in rastDisTags:
                if child.tag == 'Special_Value':
                    if child.find('SPECIAL_VALUE_TEXT').text.strip() == 'NODATA':
                        self.inImgNoData = int(child.find('SPECIAL_VALUE_COUNT').text.strip())
                    elif child.find('SPECIAL_VALUE_TEXT').text.strip() == 'SATURATED':
                        self.inImgSatVal = int(child.find('SPECIAL_VALUE_COUNT').text.strip())

            self.pleiadesSat = topLevelDataSources.find('Source_Identification').find('Strip_Source').find('MISSION_INDEX').text.strip()

            for child in topLevelQualityAssess:
                if 'Cloud_Cotation' in child.find('MEASURE_NAME').text.strip():
                    self.cloudsMask = os.path.join(filesDIR, child.find('Quality_Mask').find('Component').find('COMPONENT_PATH').attrib['href'].strip())
                elif 'Area_Of_Interest' in child.find('MEASURE_NAME').text.strip():
                    self.roiMask = os.path.join(filesDIR, child.find('Quality_Mask').find('Component').find('COMPONENT_PATH').attrib['href'].strip())

            if not os.path.exists(self.cloudsMask):
                self.cloudsMask = None
            if not os.path.exists(self.roiMask):
                self.roiMask = None

            print("Processing Input File: ", self.fileName)

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
        Customises the generic name for the Pleiades sensor
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
        raise ARCSIException("Pleiades does not provide any image masks, do not use the MASK option.")

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
        pdsBand = collections.namedtuple('Band', ['bandName', 'bandIndex', 'bias', 'gain'])
        bandDefnSeq.append(pdsBand(bandName="Blue", bandIndex=3, bias=self.b2RadBias, gain=self.b2RadGain))
        bandDefnSeq.append(pdsBand(bandName="Green", bandIndex=2, bias=self.b1RadBias, gain=self.b1RadGain))
        bandDefnSeq.append(pdsBand(bandName="Red",  bandIndex=1, bias=self.b0RadBias, gain=self.b0RadGain))
        bandDefnSeq.append(pdsBand(bandName="NIR", bandIndex=4, bias=self.b3RadBias, gain=self.b3RadGain))
        rsgislib.imagecalibration.spot5ToRadiance(self.fileName, outputImage, outFormat, bandDefnSeq)

        return outputImage, None

    def generateImageSaturationMask(self, outputPath, outputName, outFormat):
        print("Generate Saturation Image")
        outputImage = os.path.join(outputPath, outputName)

        pdsBand = collections.namedtuple('Band', ['bandName', 'fileName', 'bandIndex', 'satVal'])
        bandDefnSeq = list()
        self.inImgSatVal = float(self.inImgSatVal)
        bandDefnSeq.append(pdsBand(bandName="Blue", fileName=self.fileName, bandIndex=3, satVal=self.inImgSatVal))
        bandDefnSeq.append(pdsBand(bandName="Green", fileName=self.fileName, bandIndex=2, satVal=self.inImgSatVal))
        bandDefnSeq.append(pdsBand(bandName="Red", fileName=self.fileName, bandIndex=1, satVal=self.inImgSatVal))
        bandDefnSeq.append(pdsBand(bandName="NIR", fileName=self.fileName, bandIndex=4, satVal=self.inImgSatVal))
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

        if self.pleiadesSat == '1A':
            # Band 1 - Blue
            s.wavelength = Py6S.Wavelength(0.43, 1, [0.001603,0.001338,0.004344,0.011072,0.017627,0.023367,0.057403,0.134555,0.223763,0.308543,0.461557,0.650821,0.755369,0.747720,0.716819,0.718538,0.756696,0.810109,0.842166,0.823437,0.775247,0.752701,0.780220,0.819927,0.851663,0.860275,0.858684,0.865826,0.882822,0.904041,0.919695,0.932596,0.950233,0.975798,0.994977,1.000000,0.995063,0.980628,0.941750,0.843622,0.671142,0.463340,0.288865,0.167078,0.090180,0.050519,0.031488,0.023814,0.021311,0.020630,0.019560,0.016794,0.011358,0.006652,0.004144,0.003030,0.002416,0.001990,0.001568,0.001136,0.001253,0.000836,0.000551,0.000420,0.000362,0.000378,0.000532,0.001109,0.001987,0.001220,0.000375,0.000147,0.000075,0.000053,0.000056,0.000057,0.000056,0.000038,0.000035,0.000021,0.000014,0.000020,0.000004,0.000011,0.000012,0.000011,0.000009,0.000012,0.000019,0.000009,0.000011,0.000017,0.000005,0.000009,0.000024,0.000039,0.000013,0.000024,0.000010,0.000011,0.000021,0.000014,0.000006,0.000003,0.000009,0.000012,0.000009,0.000009,0.000006,0.000012,0.000006,0.000014,0.000017,0.000007,0.000010,0.000027,0.000063,0.000219,0.000761,0.001119,0.000754,0.000408,0.000355,0.000406,0.000679,0.001629,0.002400,0.001032,0.000348,0.000166,0.000083,0.000097,0.000046,0.000026,0.000032,0.000041,0.000016,0.000009,0.000047,0.000079,0.000022,0.000054,0.000083,0.000105,0.000183,0.000260,0.000442,0.000710,0.000865,0.000737,0.000552,0.000395,0.000281,0.000234,0.000225,0.000192,0.000220,0.000234,0.000245,0.000245,0.000278,0.000351,0.000432,0.000533,0.000689,0.001028,0.001563,0.002609,0.004575,0.009062,0.017163,0.023793,0.021523,0.015914,0.011956,0.009482,0.007869,0.007224,0.007242,0.007568,0.008379,0.009864,0.012259,0.016267,0.022602,0.032027,0.044430,0.054669,0.056424,0.048004,0.035799,0.025834,0.018887,0.014439,0.011679,0.009911,0.008647,0.007810,0.007227,0.006996,0.006923,0.006914,0.007021,0.007253,0.007373,0.007505,0.007470,0.007067,0.006368,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000])
            s.run()
            sixsCoeffs[0,0] = float(s.outputs.values['coef_xa'])
            sixsCoeffs[0,1] = float(s.outputs.values['coef_xb'])
            sixsCoeffs[0,2] = float(s.outputs.values['coef_xc'])
            sixsCoeffs[0,3] = float(s.outputs.values['direct_solar_irradiance'])
            sixsCoeffs[0,4] = float(s.outputs.values['diffuse_solar_irradiance'])
            sixsCoeffs[0,5] = float(s.outputs.values['environmental_irradiance'])

            # Band 2 - Green
            s.wavelength = Py6S.Wavelength(0.43, 1, [0.000003,0.000003,0.000005,0.000007,0.000007,0.000003,0.000007,0.000006,0.000004,0.000003,0.000015,0.000023,0.000011,0.000013,0.000022,0.000019,0.000090,0.000115,0.000112,0.000189,0.000323,0.000485,0.000624,0.000815,0.001521,0.004130,0.009954,0.019642,0.032548,0.056793,0.123791,0.285909,0.528976,0.771625,0.883804,0.907957,0.913146,0.913728,0.922484,0.936708,0.949760,0.954499,0.958582,0.964206,0.970527,0.972265,0.967518,0.958910,0.952449,0.952466,0.956048,0.955179,0.948990,0.947145,0.954450,0.971060,0.989818,1.000000,0.995360,0.969822,0.925304,0.863324,0.794828,0.723897,0.645327,0.543852,0.417028,0.276671,0.157527,0.085607,0.049226,0.032724,0.023793,0.018197,0.014062,0.009966,0.005845,0.003038,0.001536,0.000839,0.000488,0.000312,0.000207,0.000138,0.000093,0.000070,0.000064,0.000054,0.000041,0.000070,0.000048,0.000047,0.000062,0.000067,0.000148,0.000251,0.000299,0.000230,0.000127,0.000067,0.000031,0.000032,0.000017,0.000007,0.000006,0.000018,0.000011,0.000017,0.000011,0.000003,0.000003,0.000003,0.000014,0.000013,0.000017,0.000010,0.000007,0.000024,0.000033,0.000130,0.000277,0.000189,0.000124,0.000024,0.000007,0.000007,0.000004,0.000010,0.000010,0.000003,0.000016,0.000023,0.000007,0.000010,0.000009,0.000007,0.000017,0.000023,0.000002,0.000021,0.000010,0.000012,0.000034,0.000009,0.000018,0.000017,0.000019,0.000018,0.000029,0.000029,0.000021,0.000043,0.000030,0.000053,0.000093,0.000134,0.000277,0.000705,0.003185,0.014191,0.008339,0.002244,0.000918,0.000572,0.000437,0.000403,0.000418,0.000445,0.000517,0.000602,0.000758,0.001111,0.001938,0.003691,0.006357,0.010271,0.013108,0.013260,0.011115,0.008366,0.006564,0.005685,0.005380,0.005623,0.006200,0.007239,0.008920,0.011510,0.014942,0.019269,0.023479,0.026440,0.027049,0.025545,0.023397,0.021395,0.019944,0.019253,0.019074,0.019689,0.020694,0.022011,0.023230,0.023757,0.022986,0.020660,0.017225,0.013292,0.009782,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000])
            s.run()
            sixsCoeffs[1,0] = float(s.outputs.values['coef_xa'])
            sixsCoeffs[1,1] = float(s.outputs.values['coef_xb'])
            sixsCoeffs[1,2] = float(s.outputs.values['coef_xc'])
            sixsCoeffs[1,3] = float(s.outputs.values['direct_solar_irradiance'])
            sixsCoeffs[1,4] = float(s.outputs.values['diffuse_solar_irradiance'])
            sixsCoeffs[1,5] = float(s.outputs.values['environmental_irradiance'])

            # Band 3 - Red
            s.wavelength = Py6S.Wavelength(0.43, 1, [0.000454,0.000853,0.000521,0.001195,0.005142,0.008003,0.006693,0.010859,0.024691,0.065359,0.122542,0.057009,0.021375,0.012797,0.006278,0.002354,0.000941,0.000517,0.000374,0.000382,0.000516,0.000846,0.001013,0.000643,0.000351,0.000251,0.000223,0.000252,0.000355,0.000500,0.000526,0.000386,0.000253,0.000163,0.000107,0.000088,0.000070,0.000058,0.000055,0.000048,0.000037,0.000032,0.000019,0.000018,0.000037,0.000024,0.000034,0.000015,0.000007,0.000007,0.000004,0.000002,0.000009,0.000029,0.000024,0.000029,0.000039,0.000045,0.000069,0.000081,0.000104,0.000123,0.000135,0.000154,0.000183,0.000221,0.000745,0.001244,0.002142,0.003819,0.006805,0.012333,0.022178,0.041333,0.078071,0.151934,0.277675,0.451038,0.629132,0.762549,0.832945,0.857906,0.865887,0.869263,0.875221,0.885776,0.900593,0.917488,0.934880,0.947811,0.956953,0.962330,0.964767,0.962429,0.961307,0.962025,0.969915,0.981157,0.993393,1.000000,0.980951,0.952263,0.913173,0.869401,0.825208,0.783047,0.736127,0.673489,0.587753,0.480491,0.363007,0.252303,0.162603,0.102221,0.064127,0.041916,0.028464,0.020455,0.015370,0.012117,0.009881,0.008317,0.007102,0.006095,0.005172,0.004314,0.003495,0.002771,0.003589,0.003031,0.002317,0.001784,0.001331,0.001021,0.000790,0.000639,0.000508,0.000412,0.000379,0.000359,0.000298,0.000279,0.000281,0.000262,0.000286,0.000295,0.000276,0.000316,0.000375,0.000430,0.000519,0.000575,0.000619,0.000650,0.000652,0.000604,0.000537,0.000464,0.000377,0.000322,0.000284,0.000254,0.000227,0.000223,0.000167,0.000198,0.000221,0.000209,0.000233,0.000205,0.000187,0.000251,0.000293,0.000264,0.000259,0.000449,0.000536,0.000567,0.000821,0.000986,0.001308,0.001761,0.002403,0.003268,0.004364,0.005662,0.006890,0.007822,0.008330,0.008436,0.008185,0.008138,0.008001,0.007768,0.007784,0.007998,0.008152,0.008506,0.008793,0.009312,0.009753,0.010136,0.010387,0.010406,0.010171,0.009522,0.008609,0.007337,0.005984,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000])
            s.run()
            sixsCoeffs[2,0] = float(s.outputs.values['coef_xa'])
            sixsCoeffs[2,1] = float(s.outputs.values['coef_xb'])
            sixsCoeffs[2,2] = float(s.outputs.values['coef_xc'])
            sixsCoeffs[2,3] = float(s.outputs.values['direct_solar_irradiance'])
            sixsCoeffs[2,4] = float(s.outputs.values['diffuse_solar_irradiance'])
            sixsCoeffs[2,5] = float(s.outputs.values['environmental_irradiance'])

            # Band 4 - NIR
            s.wavelength = Py6S.Wavelength(0.43, 1, [0.000196,0.000761,0.000382,0.000513,0.001687,0.002665,0.001225,0.000739,0.000741,0.000754,0.000776,0.000990,0.001549,0.002180,0.001988,0.001045,0.000492,0.000283,0.000214,0.000196,0.000213,0.000282,0.000418,0.000607,0.000641,0.000437,0.000265,0.000191,0.000162,0.000171,0.000185,0.000246,0.000383,0.000563,0.000678,0.000614,0.000484,0.000388,0.000360,0.000353,0.000397,0.000462,0.000551,0.000641,0.000731,0.000757,0.000767,0.000719,0.000636,0.000571,0.000529,0.000488,0.000476,0.000552,0.000718,0.001224,0.002266,0.002879,0.002041,0.000839,0.000379,0.000225,0.000140,0.000124,0.000100,0.000097,0.000104,0.000124,0.000172,0.000286,0.000528,0.000821,0.000660,0.000322,0.000168,0.000103,0.000068,0.000057,0.000048,0.000032,0.000036,0.000041,0.000044,0.000039,0.000058,0.000064,0.000062,0.000090,0.000101,0.000153,0.000192,0.000197,0.000186,0.000143,0.000136,0.000095,0.000090,0.000088,0.000079,0.000075,0.000076,0.000085,0.000094,0.000082,0.000082,0.000107,0.000121,0.000135,0.000169,0.000204,0.000251,0.000319,0.000397,0.000508,0.000633,0.000796,0.000969,0.001154,0.001352,0.001554,0.001756,0.001989,0.002251,0.002534,0.002903,0.005236,0.006401,0.008013,0.010147,0.013109,0.017135,0.022905,0.030978,0.042662,0.059190,0.083507,0.117888,0.166378,0.231114,0.315936,0.417216,0.528495,0.640959,0.747239,0.838552,0.908812,0.959366,0.988114,1.000000,0.999206,0.989642,0.967696,0.951436,0.937494,0.925472,0.915223,0.908783,0.902608,0.896683,0.892483,0.885491,0.877655,0.867639,0.856896,0.847441,0.836048,0.823698,0.813044,0.801627,0.792162,0.782795,0.776935,0.771465,0.763145,0.754340,0.745237,0.734921,0.725512,0.714614,0.703089,0.691948,0.681648,0.670520,0.659536,0.647596,0.629185,0.611093,0.593432,0.576059,0.559994,0.544875,0.530614,0.515483,0.497676,0.473410,0.438770,0.393422,0.336769,0.274273,0.211434,0.153554,0.106341,0.070442,0.045912,0.029535,0.019275,0.012680,0.008462,0.005769,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000])
            s.run()
            sixsCoeffs[3,0] = float(s.outputs.values['coef_xa'])
            sixsCoeffs[3,1] = float(s.outputs.values['coef_xb'])
            sixsCoeffs[3,2] = float(s.outputs.values['coef_xc'])
            sixsCoeffs[3,3] = float(s.outputs.values['direct_solar_irradiance'])
            sixsCoeffs[3,4] = float(s.outputs.values['diffuse_solar_irradiance'])
            sixsCoeffs[3,5] = float(s.outputs.values['environmental_irradiance'])
        elif self.pleiadesSat == '1B':
            # Band 1 - Blue
            s.wavelength = Py6S.Wavelength(0.43, 1, [0.001603,0.001338,0.004344,0.011072,0.017627,0.023367,0.057403,0.134555,0.223763,0.308543,0.461557,0.650821,0.755369,0.747720,0.716819,0.718538,0.756696,0.810109,0.842166,0.823437,0.775247,0.752701,0.780220,0.819927,0.851663,0.860275,0.858684,0.865826,0.882822,0.904041,0.919695,0.932596,0.950233,0.975798,0.994977,1.000000,0.995063,0.980628,0.941750,0.843622,0.671142,0.463340,0.288865,0.167078,0.090180,0.050519,0.031488,0.023814,0.021311,0.020630,0.019560,0.016794,0.011358,0.006652,0.004144,0.003030,0.002416,0.001990,0.001568,0.001136,0.001253,0.000836,0.000551,0.000420,0.000362,0.000378,0.000532,0.001109,0.001987,0.001220,0.000375,0.000147,0.000075,0.000053,0.000056,0.000057,0.000056,0.000038,0.000035,0.000021,0.000014,0.000020,0.000004,0.000011,0.000012,0.000011,0.000009,0.000012,0.000019,0.000009,0.000011,0.000017,0.000005,0.000009,0.000024,0.000039,0.000013,0.000024,0.000010,0.000011,0.000021,0.000014,0.000006,0.000003,0.000009,0.000012,0.000009,0.000009,0.000006,0.000012,0.000006,0.000014,0.000017,0.000007,0.000010,0.000027,0.000063,0.000219,0.000761,0.001119,0.000754,0.000408,0.000355,0.000406,0.000679,0.001629,0.002400,0.001032,0.000348,0.000166,0.000083,0.000097,0.000046,0.000026,0.000032,0.000041,0.000016,0.000009,0.000047,0.000079,0.000022,0.000054,0.000083,0.000105,0.000183,0.000260,0.000442,0.000710,0.000865,0.000737,0.000552,0.000395,0.000281,0.000234,0.000225,0.000192,0.000220,0.000234,0.000245,0.000245,0.000278,0.000351,0.000432,0.000533,0.000689,0.001028,0.001563,0.002609,0.004575,0.009062,0.017163,0.023793,0.021523,0.015914,0.011956,0.009482,0.007869,0.007224,0.007242,0.007568,0.008379,0.009864,0.012259,0.016267,0.022602,0.032027,0.044430,0.054669,0.056424,0.048004,0.035799,0.025834,0.018887,0.014439,0.011679,0.009911,0.008647,0.007810,0.007227,0.006996,0.006923,0.006914,0.007021,0.007253,0.007373,0.007505,0.007470,0.007067,0.006368,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000])
            s.run()
            sixsCoeffs[0,0] = float(s.outputs.values['coef_xa'])
            sixsCoeffs[0,1] = float(s.outputs.values['coef_xb'])
            sixsCoeffs[0,2] = float(s.outputs.values['coef_xc'])
            sixsCoeffs[0,3] = float(s.outputs.values['direct_solar_irradiance'])
            sixsCoeffs[0,4] = float(s.outputs.values['diffuse_solar_irradiance'])
            sixsCoeffs[0,5] = float(s.outputs.values['environmental_irradiance'])

            # Band 2 - Green
            s.wavelength = Py6S.Wavelength(0.43, 1, [0.000003,0.000003,0.000005,0.000007,0.000007,0.000003,0.000007,0.000006,0.000004,0.000003,0.000015,0.000023,0.000011,0.000013,0.000022,0.000019,0.000090,0.000115,0.000112,0.000189,0.000323,0.000485,0.000624,0.000815,0.001521,0.004130,0.009954,0.019642,0.032548,0.056793,0.123791,0.285909,0.528976,0.771625,0.883804,0.907957,0.913146,0.913728,0.922484,0.936708,0.949760,0.954499,0.958582,0.964206,0.970527,0.972265,0.967518,0.958910,0.952449,0.952466,0.956048,0.955179,0.948990,0.947145,0.954450,0.971060,0.989818,1.000000,0.995360,0.969822,0.925304,0.863324,0.794828,0.723897,0.645327,0.543852,0.417028,0.276671,0.157527,0.085607,0.049226,0.032724,0.023793,0.018197,0.014062,0.009966,0.005845,0.003038,0.001536,0.000839,0.000488,0.000312,0.000207,0.000138,0.000093,0.000070,0.000064,0.000054,0.000041,0.000070,0.000048,0.000047,0.000062,0.000067,0.000148,0.000251,0.000299,0.000230,0.000127,0.000067,0.000031,0.000032,0.000017,0.000007,0.000006,0.000018,0.000011,0.000017,0.000011,0.000003,0.000003,0.000003,0.000014,0.000013,0.000017,0.000010,0.000007,0.000024,0.000033,0.000130,0.000277,0.000189,0.000124,0.000024,0.000007,0.000007,0.000004,0.000010,0.000010,0.000003,0.000016,0.000023,0.000007,0.000010,0.000009,0.000007,0.000017,0.000023,0.000002,0.000021,0.000010,0.000012,0.000034,0.000009,0.000018,0.000017,0.000019,0.000018,0.000029,0.000029,0.000021,0.000043,0.000030,0.000053,0.000093,0.000134,0.000277,0.000705,0.003185,0.014191,0.008339,0.002244,0.000918,0.000572,0.000437,0.000403,0.000418,0.000445,0.000517,0.000602,0.000758,0.001111,0.001938,0.003691,0.006357,0.010271,0.013108,0.013260,0.011115,0.008366,0.006564,0.005685,0.005380,0.005623,0.006200,0.007239,0.008920,0.011510,0.014942,0.019269,0.023479,0.026440,0.027049,0.025545,0.023397,0.021395,0.019944,0.019253,0.019074,0.019689,0.020694,0.022011,0.023230,0.023757,0.022986,0.020660,0.017225,0.013292,0.009782,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000])
            s.run()
            sixsCoeffs[1,0] = float(s.outputs.values['coef_xa'])
            sixsCoeffs[1,1] = float(s.outputs.values['coef_xb'])
            sixsCoeffs[1,2] = float(s.outputs.values['coef_xc'])
            sixsCoeffs[1,3] = float(s.outputs.values['direct_solar_irradiance'])
            sixsCoeffs[1,4] = float(s.outputs.values['diffuse_solar_irradiance'])
            sixsCoeffs[1,5] = float(s.outputs.values['environmental_irradiance'])

            # Band 3 - Red
            s.wavelength = Py6S.Wavelength(0.43, 1, [0.000454,0.000853,0.000521,0.001195,0.005142,0.008003,0.006693,0.010859,0.024691,0.065359,0.122542,0.057009,0.021375,0.012797,0.006278,0.002354,0.000941,0.000517,0.000374,0.000382,0.000516,0.000846,0.001013,0.000643,0.000351,0.000251,0.000223,0.000252,0.000355,0.000500,0.000526,0.000386,0.000253,0.000163,0.000107,0.000088,0.000070,0.000058,0.000055,0.000048,0.000037,0.000032,0.000019,0.000018,0.000037,0.000024,0.000034,0.000015,0.000007,0.000007,0.000004,0.000002,0.000009,0.000029,0.000024,0.000029,0.000039,0.000045,0.000069,0.000081,0.000104,0.000123,0.000135,0.000154,0.000183,0.000221,0.000745,0.001244,0.002142,0.003819,0.006805,0.012333,0.022178,0.041333,0.078071,0.151934,0.277675,0.451038,0.629132,0.762549,0.832945,0.857906,0.865887,0.869263,0.875221,0.885776,0.900593,0.917488,0.934880,0.947811,0.956953,0.962330,0.964767,0.962429,0.961307,0.962025,0.969915,0.981157,0.993393,1.000000,0.980951,0.952263,0.913173,0.869401,0.825208,0.783047,0.736127,0.673489,0.587753,0.480491,0.363007,0.252303,0.162603,0.102221,0.064127,0.041916,0.028464,0.020455,0.015370,0.012117,0.009881,0.008317,0.007102,0.006095,0.005172,0.004314,0.003495,0.002771,0.003589,0.003031,0.002317,0.001784,0.001331,0.001021,0.000790,0.000639,0.000508,0.000412,0.000379,0.000359,0.000298,0.000279,0.000281,0.000262,0.000286,0.000295,0.000276,0.000316,0.000375,0.000430,0.000519,0.000575,0.000619,0.000650,0.000652,0.000604,0.000537,0.000464,0.000377,0.000322,0.000284,0.000254,0.000227,0.000223,0.000167,0.000198,0.000221,0.000209,0.000233,0.000205,0.000187,0.000251,0.000293,0.000264,0.000259,0.000449,0.000536,0.000567,0.000821,0.000986,0.001308,0.001761,0.002403,0.003268,0.004364,0.005662,0.006890,0.007822,0.008330,0.008436,0.008185,0.008138,0.008001,0.007768,0.007784,0.007998,0.008152,0.008506,0.008793,0.009312,0.009753,0.010136,0.010387,0.010406,0.010171,0.009522,0.008609,0.007337,0.005984,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000])
            s.run()
            sixsCoeffs[2,0] = float(s.outputs.values['coef_xa'])
            sixsCoeffs[2,1] = float(s.outputs.values['coef_xb'])
            sixsCoeffs[2,2] = float(s.outputs.values['coef_xc'])
            sixsCoeffs[2,3] = float(s.outputs.values['direct_solar_irradiance'])
            sixsCoeffs[2,4] = float(s.outputs.values['diffuse_solar_irradiance'])
            sixsCoeffs[2,5] = float(s.outputs.values['environmental_irradiance'])

            # Band 4 - NIR
            s.wavelength = Py6S.Wavelength(0.43, 1, [0.000196,0.000761,0.000382,0.000513,0.001687,0.002665,0.001225,0.000739,0.000741,0.000754,0.000776,0.000990,0.001549,0.002180,0.001988,0.001045,0.000492,0.000283,0.000214,0.000196,0.000213,0.000282,0.000418,0.000607,0.000641,0.000437,0.000265,0.000191,0.000162,0.000171,0.000185,0.000246,0.000383,0.000563,0.000678,0.000614,0.000484,0.000388,0.000360,0.000353,0.000397,0.000462,0.000551,0.000641,0.000731,0.000757,0.000767,0.000719,0.000636,0.000571,0.000529,0.000488,0.000476,0.000552,0.000718,0.001224,0.002266,0.002879,0.002041,0.000839,0.000379,0.000225,0.000140,0.000124,0.000100,0.000097,0.000104,0.000124,0.000172,0.000286,0.000528,0.000821,0.000660,0.000322,0.000168,0.000103,0.000068,0.000057,0.000048,0.000032,0.000036,0.000041,0.000044,0.000039,0.000058,0.000064,0.000062,0.000090,0.000101,0.000153,0.000192,0.000197,0.000186,0.000143,0.000136,0.000095,0.000090,0.000088,0.000079,0.000075,0.000076,0.000085,0.000094,0.000082,0.000082,0.000107,0.000121,0.000135,0.000169,0.000204,0.000251,0.000319,0.000397,0.000508,0.000633,0.000796,0.000969,0.001154,0.001352,0.001554,0.001756,0.001989,0.002251,0.002534,0.002903,0.005236,0.006401,0.008013,0.010147,0.013109,0.017135,0.022905,0.030978,0.042662,0.059190,0.083507,0.117888,0.166378,0.231114,0.315936,0.417216,0.528495,0.640959,0.747239,0.838552,0.908812,0.959366,0.988114,1.000000,0.999206,0.989642,0.967696,0.951436,0.937494,0.925472,0.915223,0.908783,0.902608,0.896683,0.892483,0.885491,0.877655,0.867639,0.856896,0.847441,0.836048,0.823698,0.813044,0.801627,0.792162,0.782795,0.776935,0.771465,0.763145,0.754340,0.745237,0.734921,0.725512,0.714614,0.703089,0.691948,0.681648,0.670520,0.659536,0.647596,0.629185,0.611093,0.593432,0.576059,0.559994,0.544875,0.530614,0.515483,0.497676,0.473410,0.438770,0.393422,0.336769,0.274273,0.211434,0.153554,0.106341,0.070442,0.045912,0.029535,0.019275,0.012680,0.008462,0.005769,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000])
            s.run()
            sixsCoeffs[3,0] = float(s.outputs.values['coef_xa'])
            sixsCoeffs[3,1] = float(s.outputs.values['coef_xb'])
            sixsCoeffs[3,2] = float(s.outputs.values['coef_xc'])
            sixsCoeffs[3,3] = float(s.outputs.values['direct_solar_irradiance'])
            sixsCoeffs[3,4] = float(s.outputs.values['diffuse_solar_irradiance'])
            sixsCoeffs[3,5] = float(s.outputs.values['environmental_irradiance'])
        else:
            raise ARCSIException("Do not recongise the satellite ("+self.pleiadesSat+") - don't have a spectral response function.")

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
        if self.pleiadesSat == '1A':
            s.wavelength = Py6S.Wavelength(0.43, 1, [0.001603,0.001338,0.004344,0.011072,0.017627,0.023367,0.057403,0.134555,0.223763,0.308543,0.461557,0.650821,0.755369,0.747720,0.716819,0.718538,0.756696,0.810109,0.842166,0.823437,0.775247,0.752701,0.780220,0.819927,0.851663,0.860275,0.858684,0.865826,0.882822,0.904041,0.919695,0.932596,0.950233,0.975798,0.994977,1.000000,0.995063,0.980628,0.941750,0.843622,0.671142,0.463340,0.288865,0.167078,0.090180,0.050519,0.031488,0.023814,0.021311,0.020630,0.019560,0.016794,0.011358,0.006652,0.004144,0.003030,0.002416,0.001990,0.001568,0.001136,0.001253,0.000836,0.000551,0.000420,0.000362,0.000378,0.000532,0.001109,0.001987,0.001220,0.000375,0.000147,0.000075,0.000053,0.000056,0.000057,0.000056,0.000038,0.000035,0.000021,0.000014,0.000020,0.000004,0.000011,0.000012,0.000011,0.000009,0.000012,0.000019,0.000009,0.000011,0.000017,0.000005,0.000009,0.000024,0.000039,0.000013,0.000024,0.000010,0.000011,0.000021,0.000014,0.000006,0.000003,0.000009,0.000012,0.000009,0.000009,0.000006,0.000012,0.000006,0.000014,0.000017,0.000007,0.000010,0.000027,0.000063,0.000219,0.000761,0.001119,0.000754,0.000408,0.000355,0.000406,0.000679,0.001629,0.002400,0.001032,0.000348,0.000166,0.000083,0.000097,0.000046,0.000026,0.000032,0.000041,0.000016,0.000009,0.000047,0.000079,0.000022,0.000054,0.000083,0.000105,0.000183,0.000260,0.000442,0.000710,0.000865,0.000737,0.000552,0.000395,0.000281,0.000234,0.000225,0.000192,0.000220,0.000234,0.000245,0.000245,0.000278,0.000351,0.000432,0.000533,0.000689,0.001028,0.001563,0.002609,0.004575,0.009062,0.017163,0.023793,0.021523,0.015914,0.011956,0.009482,0.007869,0.007224,0.007242,0.007568,0.008379,0.009864,0.012259,0.016267,0.022602,0.032027,0.044430,0.054669,0.056424,0.048004,0.035799,0.025834,0.018887,0.014439,0.011679,0.009911,0.008647,0.007810,0.007227,0.006996,0.006923,0.006914,0.007021,0.007253,0.007373,0.007505,0.007470,0.007067,0.006368,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000])
        elif self.pleiadesSat == '1B':
            s.wavelength = Py6S.Wavelength(0.43, 1, [0.001603,0.001338,0.004344,0.011072,0.017627,0.023367,0.057403,0.134555,0.223763,0.308543,0.461557,0.650821,0.755369,0.747720,0.716819,0.718538,0.756696,0.810109,0.842166,0.823437,0.775247,0.752701,0.780220,0.819927,0.851663,0.860275,0.858684,0.865826,0.882822,0.904041,0.919695,0.932596,0.950233,0.975798,0.994977,1.000000,0.995063,0.980628,0.941750,0.843622,0.671142,0.463340,0.288865,0.167078,0.090180,0.050519,0.031488,0.023814,0.021311,0.020630,0.019560,0.016794,0.011358,0.006652,0.004144,0.003030,0.002416,0.001990,0.001568,0.001136,0.001253,0.000836,0.000551,0.000420,0.000362,0.000378,0.000532,0.001109,0.001987,0.001220,0.000375,0.000147,0.000075,0.000053,0.000056,0.000057,0.000056,0.000038,0.000035,0.000021,0.000014,0.000020,0.000004,0.000011,0.000012,0.000011,0.000009,0.000012,0.000019,0.000009,0.000011,0.000017,0.000005,0.000009,0.000024,0.000039,0.000013,0.000024,0.000010,0.000011,0.000021,0.000014,0.000006,0.000003,0.000009,0.000012,0.000009,0.000009,0.000006,0.000012,0.000006,0.000014,0.000017,0.000007,0.000010,0.000027,0.000063,0.000219,0.000761,0.001119,0.000754,0.000408,0.000355,0.000406,0.000679,0.001629,0.002400,0.001032,0.000348,0.000166,0.000083,0.000097,0.000046,0.000026,0.000032,0.000041,0.000016,0.000009,0.000047,0.000079,0.000022,0.000054,0.000083,0.000105,0.000183,0.000260,0.000442,0.000710,0.000865,0.000737,0.000552,0.000395,0.000281,0.000234,0.000225,0.000192,0.000220,0.000234,0.000245,0.000245,0.000278,0.000351,0.000432,0.000533,0.000689,0.001028,0.001563,0.002609,0.004575,0.009062,0.017163,0.023793,0.021523,0.015914,0.011956,0.009482,0.007869,0.007224,0.007242,0.007568,0.008379,0.009864,0.012259,0.016267,0.022602,0.032027,0.044430,0.054669,0.056424,0.048004,0.035799,0.025834,0.018887,0.014439,0.011679,0.009911,0.008647,0.007810,0.007227,0.006996,0.006923,0.006914,0.007021,0.007253,0.007373,0.007505,0.007470,0.007067,0.006368,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000])
        else:
            raise ARCSIException("Do not recongise the satellite ("+self.pleiadesSat+") - don't have a spectral response function.")
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
