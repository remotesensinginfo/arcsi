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

# Import the future functionality (for Python 2)
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals
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
            if dimapVersion != '2.0':
                raise ARCSIException("Only DIMAP Version 2.0 is supported by this reader; provided with: \'" + dimapVersion + "\'")

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
            self.senorZenith = 0.0 ## TODO: Not sure what to set this too!!
            self.senorAzimuth = float(locGeoValsTag.find('Acquisition_Angles').find('AZIMUTH_ANGLE').text.strip())

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
        return (self.senorAzimuth, self.senorZenith)

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

    def mosaicImageTiles(self):
        raise ARCSIException("Image data does not need mosaicking")

    def convertImageToRadiance(self, outputPath, outputReflName, outputThermalName, outFormat):
        print("Converting to Radiance")
        outputImage = os.path.join(outputPath, outputReflName)

        bandDefnSeq = list()
        pdsBand = collections.namedtuple('Band', ['bandName', 'bandIndex', 'bias', 'gain'])
        bandDefnSeq.append(pdsBand(bandName="Blue", bandIndex=1, bias=self.b0RadBias, gain=self.b0RadGain))
        bandDefnSeq.append(pdsBand(bandName="Green", bandIndex=2, bias=self.b1RadBias, gain=self.b1RadGain))
        bandDefnSeq.append(pdsBand(bandName="Red", bandIndex=3, bias=self.b2RadBias, gain=self.b2RadGain))
        bandDefnSeq.append(pdsBand(bandName="NIR", bandIndex=4, bias=self.b3RadBias, gain=self.b3RadGain))
        rsgislib.imagecalibration.spot5ToRadiance(self.fileName, outputImage, outFormat, bandDefnSeq)

        return outputImage, None

    def generateImageSaturationMask(self, outputPath, outputName, outFormat):
        print("Generate Saturation Image")
        outputImage = os.path.join(outputPath, outputName)

        pdsBand = collections.namedtuple('Band', ['bandName', 'fileName', 'bandIndex', 'satVal'])
        bandDefnSeq = list()
        bandDefnSeq.append(pdsBand(bandName="Blue", fileName=self.fileName, bandIndex=1, satVal=self.inImgSatVal))
        bandDefnSeq.append(pdsBand(bandName="Green", fileName=self.fileName, bandIndex=2, satVal=self.inImgSatVal))
        bandDefnSeq.append(pdsBand(bandName="Red", fileName=self.fileName, bandIndex=3, satVal=self.inImgSatVal))
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
        rsgislib.imageutils.genValidMask(inimages=[self.fileName], outimage=outputImage, format=outFormat, nodata=self.inImgNoData)
        return outputImage

    def convertThermalToBrightness(self, inputRadImage, outputPath, outputName, outFormat, scaleFactor):
        raise ARCSIException("There are no thermal bands...")

    def convertImageToTOARefl(self, inputRadImage, outputPath, outputName, outFormat, scaleFactor):
        print("Converting to TOA")
        outputImage = os.path.join(outputPath, outputName)
        solarIrradianceVals = list()
        IrrVal = collections.namedtuple('SolarIrradiance', ['irradiance'])
        solarIrradianceVals.append(IrrVal(irradiance=self.b0SolarIrr))
        solarIrradianceVals.append(IrrVal(irradiance=self.b1SolarIrr))
        solarIrradianceVals.append(IrrVal(irradiance=self.b2SolarIrr))
        solarIrradianceVals.append(IrrVal(irradiance=self.b3SolarIrr))

        rsgislib.imagecalibration.radiance2TOARefl(inputRadImage, outputImage, outFormat, rsgislib.TYPE_16UINT, scaleFactor, self.acquisitionTime.year, self.acquisitionTime.month, self.acquisitionTime.day, self.solarZenith, solarIrradianceVals)
        return outputImage

    def generateCloudMask(self, inputReflImage, inputSatImage, inputThermalImage, inputValidImg, outputPath, outputName, outFormat, tmpPath, scaleFactor):
        raise ARCSIException("Cloud Masking Not Implemented for Pleiades.")

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
        s.geometry.view_z = self.senorZenith
        s.geometry.view_a = self.senorAzimuth
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
        s.wavelength = Py6S.Wavelength(0.45, 0.51, [0.6457, 0.6825, 0.7193, 0.7561, 0.7929, 0.8001, 0.8073, 0.8145, 0.8217, 0.8353, 0.8489, 0.8624, 0.8760, 0.8862, 0.8964, 0.9066, 0.9168, 0.9277, 0.9387, 0.9496, 0.9606, 0.8384, 0.7161, 0.5939, 0.4717])
        s.run()
        sixsCoeffs[0,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[0,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[0,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[0,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[0,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[0,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 2 - Green
        s.wavelength = Py6S.Wavelength(0.510, 0.580, [0.1314, 0.2879, 0.4444, 0.6009, 0.75740, 0.76330, 0.7692, 0.7751, 0.7809, 0.7959, 0.8109, 0.8258, 0.8408, 0.8498, 0.8588, 0.8678, 0.8768, 0.8722, 0.8675, 0.8628, 0.8582, 0.8914, 0.9246, 0.9579, 0.9911, 0.9614, 0.9316, 0.9018, 0.8721])
        s.run()
        sixsCoeffs[1,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[1,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[1,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[1,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[1,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[1,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 3 - Red
        s.wavelength = Py6S.Wavelength(0.350000, 1.347500, [0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000001,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000002,0.000003,0.000001,0.000001,0.000001,0.000002,0.000001,0.000001,0.000001,0.000001,0.000001,0.000002,0.000003,0.000002,0.000003,0.000002,0.000002,0.000002,0.000001,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000001,0.000001,0.000001,0.000001,0.000001,0.000007,0.000013,0.000025,0.000041,0.000063,0.000123,0.000220,0.000628,0.001380,0.004867,0.015002,0.095566,0.244827,0.596071,0.805405,0.872776,0.881036,0.927928,0.953189,0.965042,0.967019,0.965588,0.958258,0.955500,0.953447,0.955801,0.963907,0.970358,0.974797,0.987526,0.995433,0.999085,0.994134,0.983694,0.986727,0.885280,0.693925,0.291823,0.100695,0.020441,0.007152,0.001914,0.000951,0.000404,0.000232,0.000102,0.000065,0.000027,0.000018,0.000014,0.000018,0.000013,0.000014,0.000013,0.000013,0.000011,0.000008,0.000003,0.000004,0.000002,0.000001,0.000001,0.000000,0.000001,0.000001,0.000003,0.000004,0.000010,0.000009,0.000007,0.000009,0.000007,0.000005,0.000005,0.000002,0.000002,0.000001,0.000002,0.000001,0.000001,0.000001,0.000014,0.000015,0.000013,0.000016,0.000011,0.000013,0.000011,0.000008,0.000009,0.000008,0.000008,0.000004,0.000011,0.000004,0.000007,0.000007,0.000006,0.000006,0.000008,0.000007,0.000003,0.000007,0.000006,0.000007,0.000005,0.000007,0.000009,0.000011,0.000007,0.000009,0.000004,0.000004,0.000009,0.000006,0.000007,0.000007,0.000004,0.000006,0.000005,0.000005,0.000006,0.000005,0.000004,0.000004,0.000004,0.000004,0.000003,0.000004,0.000003,0.000003,0.000003,0.000003,0.000003,0.000003,0.000002,0.000002,0.000002,0.000002,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000000,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000000,0.000001,0.000001,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000])
        s.run()
        sixsCoeffs[2,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[2,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[2,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[2,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[2,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[2,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 4 - NIR
        s.wavelength = Py6S.Wavelength(0.350000, 1.347500, [0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000002,0.000004,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000001,0.000001,0.000000,0.000000,0.000001,0.000001,0.000001,0.000001,0.000002,0.000000,0.000000,0.000001,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000001,0.000000,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000002,0.000002,0.000003,0.000005,0.000006,0.000005,0.000004,0.000004,0.000004,0.000005,0.000006,0.000007,0.000010,0.000016,0.000026,0.000044,0.000065,0.000106,0.000149,0.000256,0.000375,0.000713,0.001183,0.002886,0.005210,0.012687,0.022741,0.058860,0.107139,0.260204,0.453897,0.759967,0.912456,0.995805,1.000000,0.988907,0.984506,0.972934,0.962395,0.941139,0.933343,0.919382,0.911682,0.897600,0.886046,0.871191,0.867420,0.844284,0.837202,0.838036,0.842575,0.835657,0.833975,0.842548,0.832789,0.818604,0.815395,0.806711,0.814738,0.794161,0.783077,0.767313,0.750939,0.739271,0.736832,0.734145,0.715699,0.697141,0.684957,0.658635,0.654340,0.632471,0.622861,0.609120,0.600083,0.590914,0.572770,0.507112,0.427640,0.275174,0.179263,0.073805,0.037563,0.013815,0.006654,0.002625,0.001930,0.001306,0.001077,0.000917,0.000914,0.000798,0.000755,0.000727,0.000806,0.000701,0.000666,0.000603,0.000647,0.000281,0.000273,0.000004,0.000003,0.000003,0.000002,0.000002,0.000002,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000000,0.000001,0.000001,0.000001,0.000000,0.000001,0.000000,0.000001,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000])
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

        Band6S = collections.namedtuple('Band6SCoeff', ['band', 'aX', 'bX', 'cX', 'DirIrr', 'DifIrr', 'EnvIrr'])
        imgBandCoeffs = list()

        sixsCoeffs = self.calc6SCoefficients(aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF)
        imgBandCoeffs.append(Band6S(band=1, aX=float(sixsCoeffs[0,0]), bX=float(sixsCoeffs[0,1]), cX=float(sixsCoeffs[0,2]), DirIrr=float(sixsCoeffs[0,3]), DifIrr=float(sixsCoeffs[0,4]), EnvIrr=float(sixsCoeffs[0,5])))
        imgBandCoeffs.append(Band6S(band=2, aX=float(sixsCoeffs[1,0]), bX=float(sixsCoeffs[1,1]), cX=float(sixsCoeffs[1,2]), DirIrr=float(sixsCoeffs[1,3]), DifIrr=float(sixsCoeffs[1,4]), EnvIrr=float(sixsCoeffs[1,5])))
        imgBandCoeffs.append(Band6S(band=3, aX=float(sixsCoeffs[2,0]), bX=float(sixsCoeffs[2,1]), cX=float(sixsCoeffs[2,2]), DirIrr=float(sixsCoeffs[2,3]), DifIrr=float(sixsCoeffs[2,4]), EnvIrr=float(sixsCoeffs[2,5])))
        imgBandCoeffs.append(Band6S(band=4, aX=float(sixsCoeffs[3,0]), bX=float(sixsCoeffs[3,1]), cX=float(sixsCoeffs[3,2]), DirIrr=float(sixsCoeffs[3,3]), DifIrr=float(sixsCoeffs[3,4]), EnvIrr=float(sixsCoeffs[3,5])))

        for band in imgBandCoeffs:
            print(band)
        rsgislib.imagecalibration.apply6SCoeffSingleParam(inputRadImage, outputImage, outFormat, rsgislib.TYPE_16UINT, scaleFactor, 0, True, imgBandCoeffs)
        return outputImage

    def convertImageToSurfaceReflDEMElevLUT(self, inputRadImage, inputDEMFile, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, aotVal, useBRDF, surfaceAltitudeMin, surfaceAltitudeMax, scaleFactor, elevCoeffs=None):
        print("Converting to Surface Reflectance")
        outputImage = os.path.join(outputPath, outputName)

        print("Build an LUT for elevation values.")
        elev6SCoeffsLUT = self.buildElevation6SCoeffLUT(aeroProfile, atmosProfile, grdRefl, aotVal, useBRDF, surfaceAltitudeMin, surfaceAltitudeMax)
        print("LUT has been built.")

        elevLUTFeat = collections.namedtuple('ElevLUTFeat', ['Elev', 'Coeffs'])
        Band6S = collections.namedtuple('Band6SCoeff', ['band', 'aX', 'bX', 'cX', 'DirIrr', 'DifIrr', 'EnvIrr'])

        elevCoeffs = list()
        for elevLUT in elev6SCoeffsLUT:
            imgBandCoeffs = list()
            sixsCoeffs = elevLUT.Coeffs
            elevVal = elevLUT.Elev
            imgBandCoeffs.append(Band6S(band=1, aX=float(sixsCoeffs[0,0]), bX=float(sixsCoeffs[0,1]), cX=float(sixsCoeffs[0,2]), DirIrr=float(sixsCoeffs[0,3]), DifIrr=float(sixsCoeffs[0,4]), EnvIrr=float(sixsCoeffs[0,5])))
            imgBandCoeffs.append(Band6S(band=2, aX=float(sixsCoeffs[1,0]), bX=float(sixsCoeffs[1,1]), cX=float(sixsCoeffs[1,2]), DirIrr=float(sixsCoeffs[1,3]), DifIrr=float(sixsCoeffs[1,4]), EnvIrr=float(sixsCoeffs[1,5])))
            imgBandCoeffs.append(Band6S(band=3, aX=float(sixsCoeffs[2,0]), bX=float(sixsCoeffs[2,1]), cX=float(sixsCoeffs[2,2]), DirIrr=float(sixsCoeffs[2,3]), DifIrr=float(sixsCoeffs[2,4]), EnvIrr=float(sixsCoeffs[2,5])))
            imgBandCoeffs.append(Band6S(band=4, aX=float(sixsCoeffs[3,0]), bX=float(sixsCoeffs[3,1]), cX=float(sixsCoeffs[3,2]), DirIrr=float(sixsCoeffs[3,3]), DifIrr=float(sixsCoeffs[3,4]), EnvIrr=float(sixsCoeffs[3,5])))

            elevCoeffs.append(elevLUTFeat(Elev=float(elevVal), Coeffs=imgBandCoeffs))

        rsgislib.imagecalibration.apply6SCoeffElevLUTParam(inputRadImage, inputDEMFile, outputImage, outFormat, rsgislib.TYPE_16UINT, scaleFactor, 0, True, elevCoeffs)
        return outputImage, elevCoeffs

    def convertImageToSurfaceReflAOTDEMElevLUT(self, inputRadImage, inputDEMFile, inputAOTImage, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, useBRDF, surfaceAltitudeMin, surfaceAltitudeMax, aotMin, aotMax, scaleFactor, elevAOTCoeffs=None):
        print("Converting to Surface Reflectance")
        outputImage = os.path.join(outputPath, outputName)

        if elevAOTCoeffs is None:
            print("Build an LUT for elevation and AOT values.")
            elevAOT6SCoeffsLUT = self.buildElevationAOT6SCoeffLUT(aeroProfile, atmosProfile, grdRefl, useBRDF, surfaceAltitudeMin, surfaceAltitudeMax, aotMin, aotMax)

            elevLUTFeat = collections.namedtuple('ElevLUTFeat', ['Elev', 'Coeffs'])
            aotLUTFeat = collections.namedtuple('AOTLUTFeat', ['AOT', 'Coeffs'])
            Band6S = collections.namedtuple('Band6SCoeff', ['band', 'aX', 'bX', 'cX', 'DirIrr', 'DifIrr', 'EnvIrr'])

            elevAOTCoeffs = list()
            for elevLUT in elevAOT6SCoeffsLUT:
                elevVal = elevLUT.Elev
                aotLUT = elevLUT.Coeffs
                aot6SCoeffsOut = list()
                for aotFeat in aotLUT:
                    sixsCoeffs = aotFeat.Coeffs
                    aotVal = aotFeat.AOT
                    imgBandCoeffs = list()
                    imgBandCoeffs.append(Band6S(band=1, aX=float(sixsCoeffs[0,0]), bX=float(sixsCoeffs[0,1]), cX=float(sixsCoeffs[0,2]), DirIrr=float(sixsCoeffs[0,3]), DifIrr=float(sixsCoeffs[0,4]), EnvIrr=float(sixsCoeffs[0,5])))
                    imgBandCoeffs.append(Band6S(band=2, aX=float(sixsCoeffs[1,0]), bX=float(sixsCoeffs[1,1]), cX=float(sixsCoeffs[1,2]), DirIrr=float(sixsCoeffs[1,3]), DifIrr=float(sixsCoeffs[1,4]), EnvIrr=float(sixsCoeffs[1,5])))
                    imgBandCoeffs.append(Band6S(band=3, aX=float(sixsCoeffs[2,0]), bX=float(sixsCoeffs[2,1]), cX=float(sixsCoeffs[2,2]), DirIrr=float(sixsCoeffs[2,3]), DifIrr=float(sixsCoeffs[2,4]), EnvIrr=float(sixsCoeffs[2,5])))
                    imgBandCoeffs.append(Band6S(band=4, aX=float(sixsCoeffs[3,0]), bX=float(sixsCoeffs[3,1]), cX=float(sixsCoeffs[3,2]), DirIrr=float(sixsCoeffs[3,3]), DifIrr=float(sixsCoeffs[3,4]), EnvIrr=float(sixsCoeffs[3,5])))
                    aot6SCoeffsOut.append(aotLUTFeat(AOT=float(aotVal), Coeffs=imgBandCoeffs))
                elevAOTCoeffs.append(elevLUTFeat(Elev=float(elevVal), Coeffs=aot6SCoeffsOut))

        rsgislib.imagecalibration.apply6SCoeffElevAOTLUTParam(inputRadImage, inputDEMFile, inputAOTImage, outputImage, outFormat, rsgislib.TYPE_16UINT, scaleFactor, 0, True, elevAOTCoeffs)

        return outputImage, elevAOTCoeffs

    def run6SToOptimiseAODValue(self, aotVal, radBlueVal, predBlueVal, aeroProfile, atmosProfile, grdRefl, surfaceAltitude):
        """Used as part of the optimastion for identifying values of AOD"""
        print("Testing AOD Val: ", aotVal,)
        sixsCoeffs = numpy.zeros((5, 3), dtype=numpy.float32)
        # Set up 6S model
        s = Py6S.SixS()
        s.atmos_profile = atmosProfile
        s.aero_profile = aeroProfile
        s.ground_reflectance = grdRefl
        s.geometry = Py6S.Geometry.User()
        s.geometry.solar_z = self.solarZenith
        s.geometry.solar_a = self.solarAzimuth
        s.geometry.view_z = self.senorZenith
        s.geometry.view_a = self.senorAzimuth
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
        s.wavelength = Py6S.Wavelength(0.45, 0.51, [0.6457, 0.6825, 0.7193, 0.7561, 0.7929, 0.8001, 0.8073, 0.8145, 0.8217, 0.8353, 0.8489, 0.8624, 0.8760, 0.8862, 0.8964, 0.9066, 0.9168, 0.9277, 0.9387, 0.9496, 0.9606, 0.8384, 0.7161, 0.5939, 0.4717])
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

    def cleanFollowProcessing(self):
        if self.createdCopyKEADNImg:
            if not self.debugMode:
                gdalDriver = gdal.GetDriverByName('KEA')
                gdalDriver.Delete(self.copiedKEADNImg)
            self.fileName = self.origDNImg

