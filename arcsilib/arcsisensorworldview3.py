"""
Module that contains the ARCSIWorldView3Sensor class.
"""
############################################################################
#  arcsisensorworldview3.py
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
# Purpose:  A class for read the WorldView3 sensor header file and applying
#           the pre-processing operations within ARCSI to the WorldView3
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
# Import OS path module for manipulating the file system
import os.path
# Import the RSGISLib Image Calibration Module.
import rsgislib.imagecalibration
# Import the RSGISLib Image Calculation Module
import rsgislib.imagecalc
import rsgislib.segmentation
import rsgislib.segmentation.segutils
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

class ARCSIWorldView3Sensor (ARCSIAbstractSensor):
    """
    A class which represents the WorldView3 sensor to read
    header parameters and apply data processing operations.
    """
    def __init__(self, debugMode, inputImage):
        ARCSIAbstractSensor.__init__(self, debugMode, inputImage)
        self.sensor = "WV3"
        self.fileName = ""
        self.catID = ""
        self.inProdLev = ""
        self.numBands = 0
        self.bandsProd = ""
        self.solarZenithMin = 0.0
        self.solarZenithMax = 0.0
        self.solarAzimuthMin = 0.0
        self.solarAzimuthMax = 0.0
        self.sensorZenithMin = 0.0
        self.sensorZenithMax = 0.0
        self.sensorAzimuthMin = 0.0
        self.sensorAzimuthMax = 0.0
        self.nadirViewAngle = 0.0
        # Coastal
        self.absCalFactB1 = 0.0
        self.effBandWidthB1 = 0.0
        # Blue
        self.absCalFactB2 = 0.0
        self.effBandWidthB2 = 0.0
        # Green
        self.absCalFactB3 = 0.0
        self.effBandWidthB3 = 0.0
        # Yellow
        self.absCalFactB4 = 0.0
        self.effBandWidthB4 = 0.0
        # Red
        self.absCalFactB5 = 0.0
        self.effBandWidthB5 = 0.0
        # Red-Edge
        self.absCalFactB6 = 0.0
        self.effBandWidthB6 = 0.0
        # NIR 1
        self.absCalFactB7 = 0.0
        self.effBandWidthB7 = 0.0
        # NIR 2
        self.absCalFactB8 = 0.0
        self.effBandWidthB8 = 0.0
        # SWIR 1
        self.absCalFactSWIRB1 = 0.0
        self.effBandWidthSWIRB1 = 0.0
        # SWIR 2
        self.absCalFactSWIRB2 = 0.0
        self.effBandWidthSWIRB2 = 0.0
        # SWIR 3
        self.absCalFactSWIRB3 = 0.0
        self.effBandWidthSWIRB3 = 0.0
        # SWIR 4
        self.absCalFactSWIRB4 = 0.0
        self.effBandWidthSWIRB4 = 0.0
        # SWIR 5
        self.absCalFactSWIRB5 = 0.0
        self.effBandWidthSWIRB5 = 0.0
        # SWIR 6
        self.absCalFactSWIRB6 = 0.0
        self.effBandWidthSWIRB6 = 0.0
        # SWIR 7
        self.absCalFactSWIRB7 = 0.0
        self.effBandWidthSWIRB7 = 0.0
        # SWIR 8
        self.absCalFactSWIRB8 = 0.0
        self.effBandWidthSWIRB8 = 0.0

    def extractHeaderParameters(self, inputHeader, wktStr):
        """
        Understands and parses the WorldView3 xml header file
        """
        try:
            self.headerFileName = os.path.split(inputHeader)[1]
            
            tree = ET.parse(inputHeader)
            root = tree.getroot()

            topLevelInfo = root.find('IMD')
            if topLevelInfo == None:
                raise ARCSIException("Cannot open top level section \'IMD\'")

            topLevelTileInfo = root.find('TIL')
            if topLevelTileInfo == None:
                raise ARCSIException("Cannot open top level section \'TIL\'")

            imageInfoTag = topLevelInfo.find('IMAGE')
            if imageInfoTag == None:
                raise ARCSIException("Cannot open \'IMAGE\' section within \'IMD\'")

            mapProjInfoTag = topLevelInfo.find('MAP_PROJECTED_PRODUCT')
            if mapProjInfoTag == None:
                raise ARCSIException("Cannot open \'MAP_PROJECTED_PRODUCT\' section within \'IMD\'")

            imageTileInfoTag = topLevelTileInfo.find('TILE')
            if imageTileInfoTag == None:
                raise ARCSIException("Cannot open \'TILE\' section within \'IMD\'")

            if imageInfoTag.find('SATID').text.strip() != 'WV03':
                raise ARCSIException("SATID is not \'WV03\' but \'" + imageInfoTag.find('SATID').text.strip() + "\'")

            self.inProdLev = topLevelInfo.find('PRODUCTLEVEL').text.strip()
            if self.inProdLev != 'LV2A':
                raise ARCSIException("PRODUCTLEVEL is not \'LV2A\' but \'" + self.inProdLev + "\'")
            
            self.bandsProd = topLevelInfo.find('BANDID').text.strip()
            if self.bandsProd == 'MS1':
                self.numBands = 4
            elif self.bandsProd == 'Multi':
                self.numBands = 8
            elif self.bandsProd == 'All-S':
                self.numBands = 8
            else:
                raise ARCSIException("BANDID is not \'MS1\' or \'Multi\' or \'All-S\' but \'" + self.bandsProd + "\'")
            
            self.catID = imageInfoTag.find('CATID').text.strip()
            tmpAcquasitionTime = imageInfoTag.find('FIRSTLINETIME').text.strip()
            tmpAcquasitionTime = tmpAcquasitionTime.replace('Z', '')
            self.acquisitionTime = datetime.datetime.strptime(tmpAcquasitionTime, "%Y-%m-%dT%H:%M:%S.%f")

            self.solarZenithMin = 90-float(imageInfoTag.find('MINSUNEL').text.strip())
            self.solarZenithMax = 90-float(imageInfoTag.find('MAXSUNEL').text.strip())
            self.solarZenith = 90-float(imageInfoTag.find('MEANSUNEL').text.strip())

            self.solarAzimuthMin = float(imageInfoTag.find('MINSUNAZ').text.strip())
            self.solarAzimuthMax = float(imageInfoTag.find('MAXSUNAZ').text.strip())
            self.solarAzimuth = float(imageInfoTag.find('MEANSUNAZ').text.strip())

            self.sensorZenithMin = 90-float(imageInfoTag.find('MINSATEL').text.strip())
            self.sensorZenithMax = 90-float(imageInfoTag.find('MAXSATEL').text.strip())
            self.sensorZenith = 90-float(imageInfoTag.find('MEANSATEL').text.strip())

            self.sensorAzimuthMin = float(imageInfoTag.find('MINSATAZ').text.strip())
            self.sensorAzimuthMax = float(imageInfoTag.find('MAXSATAZ').text.strip())
            self.sensorAzimuth = float(imageInfoTag.find('MEANSATAZ').text.strip())

            self.nadirViewAngle = float(imageInfoTag.find('MEANOFFNADIRVIEWANGLE').text.strip())

            if mapProjInfoTag.find('MAPPROJNAME').text.strip() != "UTM":
                raise ARCSIException("Expecting input image to be projected as UTM WGS84")

            utmZone = "WGS84" + mapProjInfoTag.find('MAPPROJNAME').text.strip() + mapProjInfoTag.find('MAPZONE').text.strip() + mapProjInfoTag.find('MAPHEMI').text.strip()
            epsgCode = self.epsgCodes[utmZone]
            inProj = osr.SpatialReference()
            inProj.ImportFromEPSG(epsgCode)
            if self.inWKT == "":
                self.inWKT = inProj.ExportToWkt()

            self.xTL = float(imageTileInfoTag.find('ULX').text.strip())
            self.yTL = float(imageTileInfoTag.find('ULY').text.strip())
            self.xTR = float(imageTileInfoTag.find('URX').text.strip())
            self.yTR = float(imageTileInfoTag.find('URY').text.strip())
            self.xBL = float(imageTileInfoTag.find('LLX').text.strip())
            self.yBL = float(imageTileInfoTag.find('LRY').text.strip())
            self.xBR = float(imageTileInfoTag.find('LLY').text.strip())
            self.yBR = float(imageTileInfoTag.find('LRY').text.strip())
            self.xCentre = self.xTL + ((self.xTR - self.xTL)/2)
            self.yCentre = self.yBR + ((self.yTL - self.yBR)/2)

            self.latTL = float(imageTileInfoTag.find('ULLAT').text.strip())
            self.lonTL = float(imageTileInfoTag.find('ULLON').text.strip())
            self.latTR = float(imageTileInfoTag.find('URLAT').text.strip())
            self.lonTR = float(imageTileInfoTag.find('URLON').text.strip())
            self.latBL = float(imageTileInfoTag.find('LLLAT').text.strip())
            self.lonBL = float(imageTileInfoTag.find('LLLON').text.strip())
            self.latBR = float(imageTileInfoTag.find('LRLAT').text.strip())
            self.lonBR = float(imageTileInfoTag.find('LRLON').text.strip())

            arcsiUtils = ARCSIUtils()
            self.lonCentre, self.latCentre = arcsiUtils.getLongLat(inProj, self.xCentre, self.yCentre)

            filesDIR = os.path.dirname(inputHeader)
            if not self.userSpInputImage is None:
                self.fileName = os.path.abspath(self.userSpInputImage)
            else:
                self.fileName = os.path.join(filesDIR, imageTileInfoTag.find('FILENAME').text.strip())

            if self.bandsProd == 'MS1':
                # Blue
                self.absCalFactB2 = float(topLevelInfo.find('BAND_B').find('ABSCALFACTOR').text.strip())
                self.effBandWidthB2 = float(topLevelInfo.find('BAND_B').find('EFFECTIVEBANDWIDTH').text.strip())
                # Green
                self.absCalFactB3 = float(topLevelInfo.find('BAND_G').find('ABSCALFACTOR').text.strip())
                self.effBandWidthB3 = float(topLevelInfo.find('BAND_G').find('EFFECTIVEBANDWIDTH').text.strip())
                # Red
                self.absCalFactB5 = float(topLevelInfo.find('BAND_R').find('ABSCALFACTOR').text.strip())
                self.effBandWidthB5 = float(topLevelInfo.find('BAND_R').find('EFFECTIVEBANDWIDTH').text.strip())
                # NIR 1
                self.absCalFactB7 = float(topLevelInfo.find('BAND_N').find('ABSCALFACTOR').text.strip())
                self.effBandWidthB7 = float(topLevelInfo.find('BAND_N').find('EFFECTIVEBANDWIDTH').text.strip())
            elif self.bandsProd == 'Multi':
                # Coastal
                self.absCalFactB1 = float(topLevelInfo.find('BAND_C').find('ABSCALFACTOR').text.strip())
                self.effBandWidthB1 = float(topLevelInfo.find('BAND_C').find('EFFECTIVEBANDWIDTH').text.strip())
                # Blue
                self.absCalFactB2 = float(topLevelInfo.find('BAND_B').find('ABSCALFACTOR').text.strip())
                self.effBandWidthB2 = float(topLevelInfo.find('BAND_B').find('EFFECTIVEBANDWIDTH').text.strip())
                # Green
                self.absCalFactB3 = float(topLevelInfo.find('BAND_G').find('ABSCALFACTOR').text.strip())
                self.effBandWidthB3 = float(topLevelInfo.find('BAND_G').find('EFFECTIVEBANDWIDTH').text.strip())
                # Yellow
                self.absCalFactB4 = float(topLevelInfo.find('BAND_Y').find('ABSCALFACTOR').text.strip())
                self.effBandWidthB4 = float(topLevelInfo.find('BAND_Y').find('EFFECTIVEBANDWIDTH').text.strip())
                # Red
                self.absCalFactB5 = float(topLevelInfo.find('BAND_R').find('ABSCALFACTOR').text.strip())
                self.effBandWidthB5 = float(topLevelInfo.find('BAND_R').find('EFFECTIVEBANDWIDTH').text.strip())
                # Red-Edge
                self.absCalFactB6 = float(topLevelInfo.find('BAND_RE').find('ABSCALFACTOR').text.strip())
                self.effBandWidthB6 = float(topLevelInfo.find('BAND_RE').find('EFFECTIVEBANDWIDTH').text.strip())
                # NIR 1
                self.absCalFactB7 = float(topLevelInfo.find('BAND_N').find('ABSCALFACTOR').text.strip())
                self.effBandWidthB7 = float(topLevelInfo.find('BAND_N').find('EFFECTIVEBANDWIDTH').text.strip())
                # NIR 2
                self.absCalFactB8 = float(topLevelInfo.find('BAND_N2').find('ABSCALFACTOR').text.strip())
                self.effBandWidthB8 = float(topLevelInfo.find('BAND_N2').find('EFFECTIVEBANDWIDTH').text.strip())
            elif self.bandsProd == 'All-S':
                # SWIR 1
                self.absCalFactSWIRB1 = float(topLevelInfo.find('BAND_S1').find('ABSCALFACTOR').text.strip())
                self.effBandWidthSWIRB1 = float(topLevelInfo.find('BAND_S1').find('EFFECTIVEBANDWIDTH').text.strip())
                # SWIR 2
                self.absCalFactSWIRB2 = float(topLevelInfo.find('BAND_S2').find('ABSCALFACTOR').text.strip())
                self.effBandWidthSWIRB2 = float(topLevelInfo.find('BAND_S2').find('EFFECTIVEBANDWIDTH').text.strip())
                # SWIR 3
                self.absCalFactSWIRB3 = float(topLevelInfo.find('BAND_S3').find('ABSCALFACTOR').text.strip())
                self.effBandWidthSWIRB3 = float(topLevelInfo.find('BAND_S3').find('EFFECTIVEBANDWIDTH').text.strip())
                # SWIR 4
                self.absCalFactSWIRB4 = float(topLevelInfo.find('BAND_S4').find('ABSCALFACTOR').text.strip())
                self.effBandWidthSWIRB4 = float(topLevelInfo.find('BAND_S4').find('EFFECTIVEBANDWIDTH').text.strip())
                # SWIR 5
                self.absCalFactSWIRB5 = float(topLevelInfo.find('BAND_S5').find('ABSCALFACTOR').text.strip())
                self.effBandWidthSWIRB5 = float(topLevelInfo.find('BAND_S5').find('EFFECTIVEBANDWIDTH').text.strip())
                # SWIR 6
                self.absCalFactSWIRB6 = float(topLevelInfo.find('BAND_S6').find('ABSCALFACTOR').text.strip())
                self.effBandWidthSWIRB6 = float(topLevelInfo.find('BAND_S6').find('EFFECTIVEBANDWIDTH').text.strip())
                # SWIR 7
                self.absCalFactSWIRB7 = float(topLevelInfo.find('BAND_S7').find('ABSCALFACTOR').text.strip())
                self.effBandWidthSWIRB7 = float(topLevelInfo.find('BAND_S7').find('EFFECTIVEBANDWIDTH').text.strip())
                # SWIR 8
                self.absCalFactSWIRB8 = float(topLevelInfo.find('BAND_S8').find('ABSCALFACTOR').text.strip())
                self.effBandWidthSWIRB8 = float(topLevelInfo.find('BAND_S8').find('EFFECTIVEBANDWIDTH').text.strip())
            else: 
                raise ARCSIException("Do not reconise the product type.")

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
        Customises the generic name for the WorldView3 sensor
        """
        outname = self.defaultGenBaseOutFileName()
        return outname

    def expectedImageDataPresent(self):
        imageDataPresent = True

        if not os.path.exists(self.fileName):
            imageDataPresent = False

        return imageDataPresent

    def applyImageDataMask(self, inputHeader, outputPath, outputMaskName, outputImgName, outFormat, outWKTFile):
        raise ARCSIException("WorldView3 does not provide any image masks, do not use the MASK option.")

    def mosaicImageTiles(self, outputPath):
        raise ARCSIException("Image data does not need mosaicking")

    def resampleImgRes(self, outputPath, resampleToLowResImg, resampleMethod='cubic', multicore=False):
        raise ARCSIException("Image data does not need resampling")

    def sharpenLowResRadImgBands(self, inputImg, outputImage, outFormat):
        raise ARCSIException("Image sharpening is not available for this sensor.")

    def convertImageToRadiance(self, outputPath, outputReflName, outputThermalName, outFormat):
        print("Converting to Radiance")
        outputImage = os.path.join(outputPath, outputReflName)

        bandDefnSeq = list()
        wv2Band = collections.namedtuple('WV2Band', ['bandName', 'bandIndex', 'absCalFact', 'effBandWidth'])
        if self.bandsProd == 'MS1':
            bandDefnSeq.append(wv2Band(bandName="Blue", bandIndex=1, absCalFact=self.absCalFactB2, effBandWidth=self.effBandWidthB2))
            bandDefnSeq.append(wv2Band(bandName="Green", bandIndex=2, absCalFact=self.absCalFactB3, effBandWidth=self.effBandWidthB3))
            bandDefnSeq.append(wv2Band(bandName="Red", bandIndex=3, absCalFact=self.absCalFactB5, effBandWidth=self.effBandWidthB5))
            bandDefnSeq.append(wv2Band(bandName="NIR", bandIndex=4, absCalFact=self.absCalFactB7, effBandWidth=self.effBandWidthB7))
        elif self.bandsProd == 'Multi':
            bandDefnSeq.append(wv2Band(bandName="Coastal", bandIndex=1, absCalFact=self.absCalFactB1, effBandWidth=self.effBandWidthB1))
            bandDefnSeq.append(wv2Band(bandName="Blue", bandIndex=2, absCalFact=self.absCalFactB2, effBandWidth=self.effBandWidthB2))
            bandDefnSeq.append(wv2Band(bandName="Green", bandIndex=3, absCalFact=self.absCalFactB3, effBandWidth=self.effBandWidthB3))
            bandDefnSeq.append(wv2Band(bandName="Yellow", bandIndex=4, absCalFact=self.absCalFactB4, effBandWidth=self.effBandWidthB4))
            bandDefnSeq.append(wv2Band(bandName="Red", bandIndex=5, absCalFact=self.absCalFactB5, effBandWidth=self.effBandWidthB5))
            bandDefnSeq.append(wv2Band(bandName="RedEdge", bandIndex=6, absCalFact=self.absCalFactB6, effBandWidth=self.effBandWidthB6))
            bandDefnSeq.append(wv2Band(bandName="NIR1", bandIndex=7, absCalFact=self.absCalFactB7, effBandWidth=self.effBandWidthB7))
            bandDefnSeq.append(wv2Band(bandName="NIR2", bandIndex=8, absCalFact=self.absCalFactB8, effBandWidth=self.effBandWidthB8))
        elif self.bandsProd == 'All-S':
            bandDefnSeq.append(wv2Band(bandName="SWIR1", bandIndex=1, absCalFact=self.absCalFactSWIRB1, effBandWidth=self.effBandWidthSWIRB1))
            bandDefnSeq.append(wv2Band(bandName="SWIR2", bandIndex=2, absCalFact=self.absCalFactSWIRB2, effBandWidth=self.effBandWidthSWIRB2))
            bandDefnSeq.append(wv2Band(bandName="SWIR3", bandIndex=3, absCalFact=self.absCalFactSWIRB3, effBandWidth=self.effBandWidthSWIRB3))
            bandDefnSeq.append(wv2Band(bandName="SWIR4", bandIndex=4, absCalFact=self.absCalFactSWIRB4, effBandWidth=self.effBandWidthSWIRB4))
            bandDefnSeq.append(wv2Band(bandName="SWIR5", bandIndex=5, absCalFact=self.absCalFactSWIRB5, effBandWidth=self.effBandWidthSWIRB5))
            bandDefnSeq.append(wv2Band(bandName="SWIR6", bandIndex=6, absCalFact=self.absCalFactSWIRB6, effBandWidth=self.effBandWidthSWIRB6))
            bandDefnSeq.append(wv2Band(bandName="SWIR7", bandIndex=7, absCalFact=self.absCalFactSWIRB7, effBandWidth=self.effBandWidthSWIRB7))
            bandDefnSeq.append(wv2Band(bandName="SWIR8", bandIndex=8, absCalFact=self.absCalFactSWIRB8, effBandWidth=self.effBandWidthSWIRB8))
        else: 
            raise ARCSIException("Do not reconise the product type.")
        rsgislib.imagecalibration.worldview2ToRadiance(self.fileName, outputImage, outFormat, bandDefnSeq)

        return outputImage, None

    def generateImageSaturationMask(self, outputPath, outputName, outFormat):
        print("Generate Saturation Image")
        outputImage = os.path.join(outputPath, outputName)

        wv2Band = collections.namedtuple('WV2Band', ['bandName', 'fileName', 'bandIndex', 'satVal'])
        bandDefnSeq = list()
        if self.bandsProd == 'MS1':
            bandDefnSeq.append(wv2Band(bandName="Blue", fileName=self.fileName, bandIndex=1, satVal=65535.0))
            bandDefnSeq.append(wv2Band(bandName="Green", fileName=self.fileName, bandIndex=2, satVal=65535.0))
            bandDefnSeq.append(wv2Band(bandName="Red", fileName=self.fileName, bandIndex=3, satVal=65535.0))
            bandDefnSeq.append(wv2Band(bandName="NIR", fileName=self.fileName, bandIndex=4, satVal=65535.0))
        elif self.bandsProd == 'Multi':
            bandDefnSeq.append(wv2Band(bandName="Coastal", fileName=self.fileName, bandIndex=1, satVal=65535.0))
            bandDefnSeq.append(wv2Band(bandName="Blue", fileName=self.fileName, bandIndex=2, satVal=65535.0))
            bandDefnSeq.append(wv2Band(bandName="Green", fileName=self.fileName, bandIndex=3, satVal=65535.0))
            bandDefnSeq.append(wv2Band(bandName="Yellow", fileName=self.fileName, bandIndex=4, satVal=65535.0))
            bandDefnSeq.append(wv2Band(bandName="Red", fileName=self.fileName, bandIndex=5, satVal=65535.0))
            bandDefnSeq.append(wv2Band(bandName="RedEdge", fileName=self.fileName, bandIndex=6, satVal=65535.0))
            bandDefnSeq.append(wv2Band(bandName="NIR1", fileName=self.fileName, bandIndex=7, satVal=65535.0))
            bandDefnSeq.append(wv2Band(bandName="NIR2", fileName=self.fileName, bandIndex=8, satVal=65535.0))
        elif self.bandsProd == 'All-S':
            bandDefnSeq.append(wv2Band(bandName="SWIR1", fileName=self.fileName, bandIndex=1, satVal=65535.0))
            bandDefnSeq.append(wv2Band(bandName="SWIR2", fileName=self.fileName, bandIndex=2, satVal=65535.0))
            bandDefnSeq.append(wv2Band(bandName="SWIR3", fileName=self.fileName, bandIndex=3, satVal=65535.0))
            bandDefnSeq.append(wv2Band(bandName="SWIR4", fileName=self.fileName, bandIndex=4, satVal=65535.0))
            bandDefnSeq.append(wv2Band(bandName="SWIR5", fileName=self.fileName, bandIndex=5, satVal=65535.0))
            bandDefnSeq.append(wv2Band(bandName="SWIR6", fileName=self.fileName, bandIndex=6, satVal=65535.0))
            bandDefnSeq.append(wv2Band(bandName="SWIR7", fileName=self.fileName, bandIndex=7, satVal=65535.0))
            bandDefnSeq.append(wv2Band(bandName="SWIR8", fileName=self.fileName, bandIndex=8, satVal=65535.0))
        else: 
            raise ARCSIException("Do not reconise the product type.")

        rsgislib.imagecalibration.saturatedPixelsMask(outputImage, outFormat, bandDefnSeq)

        return outputImage

    def generateValidImageDataMask(self, outputPath, outputMaskName, viewAngleImg, outFormat):
        print("Generate valid image mask")
        outputImage = os.path.join(outputPath, outputMaskName)
        rsgislib.imageutils.genValidMask(inimages=[self.fileName], outimage=outputImage, gdalformat=outFormat, nodata=0.0)
        return outputImage

    def convertThermalToBrightness(self, inputRadImage, outputPath, outputName, outFormat, scaleFactor):
        raise ARCSIException("There are no thermal bands...")

    def convertImageToTOARefl(self, inputRadImage, outputPath, outputName, outFormat, scaleFactor):
        print("Converting to TOA")
        outputImage = os.path.join(outputPath, outputName)
        solarIrradianceVals = list()
        IrrVal = collections.namedtuple('SolarIrradiance', ['irradiance'])
        # Spectral Irradiance from Thuillier (2003)
        if self.bandsProd == 'MS1':
            solarIrradianceVals.append(IrrVal(irradiance=2004.61))
            solarIrradianceVals.append(IrrVal(irradiance=1830.18))
            solarIrradianceVals.append(IrrVal(irradiance=1535.33))
            solarIrradianceVals.append(IrrVal(irradiance=1055.94))
        elif self.bandsProd == 'Multi':
            solarIrradianceVals.append(IrrVal(irradiance=1757.89))
            solarIrradianceVals.append(IrrVal(irradiance=2004.61))
            solarIrradianceVals.append(IrrVal(irradiance=1830.18))
            solarIrradianceVals.append(IrrVal(irradiance=1712.07))
            solarIrradianceVals.append(IrrVal(irradiance=1535.33))
            solarIrradianceVals.append(IrrVal(irradiance=1348.08))
            solarIrradianceVals.append(IrrVal(irradiance=1055.94))
            solarIrradianceVals.append(IrrVal(irradiance=858.77))
        elif self.bandsProd == 'All-S':
            solarIrradianceVals.append(IrrVal(irradiance=479.019))
            solarIrradianceVals.append(IrrVal(irradiance=263.797))
            solarIrradianceVals.append(IrrVal(irradiance=225.283))
            solarIrradianceVals.append(IrrVal(irradiance=197.552))
            solarIrradianceVals.append(IrrVal(irradiance=90.4178))
            solarIrradianceVals.append(IrrVal(irradiance=85.0642))
            solarIrradianceVals.append(IrrVal(irradiance=76.9507))
            solarIrradianceVals.append(IrrVal(irradiance=68.0988))
        else: 
            raise ARCSIException("Do not reconise the product type.")
        rsgislib.imagecalibration.radiance2TOARefl(inputRadImage, outputImage, outFormat, rsgislib.TYPE_16UINT, scaleFactor, self.acquisitionTime.year, self.acquisitionTime.month, self.acquisitionTime.day, self.solarZenith, solarIrradianceVals)
        return outputImage

    def generateCloudMask(self, inputReflImage, inputSatImage, inputThermalImage, inputViewAngleImg, inputValidImg, outputPath, outputName, outFormat, tmpPath, scaleFactor, cloud_msk_methods=None):
        raise ARCSIException("Cloud Masking Not Implemented for WorldView3.")

    def createCloudMaskDataArray(self, inImgDataArr):
        return inImgDataArr

    def defineDarkShadowImageBand(self):
        if self.bandsProd == 'MS1':
            return 4
        elif self.bandsProd == 'Multi':
            return 8
        elif self.bandsProd == 'All-S':
            return 1
        else:
            raise ARCSIException("Don't recognise product type.")


    def calc6SCoefficients(self, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF):
        sixsCoeffs = numpy.zeros((8, 6), dtype=numpy.float32)
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

        if self.bandsProd == 'MS1':
            sixsCoeffs = None
            sixsCoeffs = numpy.zeros((4, 6), dtype=numpy.float32)
            # Band 1 - Blue
            s.wavelength = Py6S.Wavelength(0.395, 0.615, [0.000000, 0.000000, 0.000001, 0.000001, 0.000002, 0.000002, 0.000003, 0.000005, 0.000007, 0.000011, 0.000016, 0.000026, 0.000046, 0.000092, 0.000187, 0.000287, 0.000757, 0.001534, 0.005096, 0.010599, 0.051236, 0.103638, 0.363869, 0.578093, 0.746604, 0.781177, 0.802868, 0.813323, 0.825085, 0.820022, 0.839875, 0.849677, 0.876742, 0.876547, 0.894562, 0.899177, 0.898469, 0.912201, 0.950024, 0.963251, 0.964915, 0.980350, 0.996347, 0.999729, 0.952372, 0.850259, 0.533870, 0.296567, 0.088810, 0.041059, 0.009926, 0.004817, 0.001574, 0.000845, 0.000389, 0.000271, 0.000162, 0.000096, 0.000044, 0.000032, 0.000022, 0.000017, 0.000012, 0.000010, 0.000008, 0.000006, 0.000004, 0.000004, 0.000003, 0.000003, 0.000003, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000003, 0.000004, 0.000010, 0.000005, 0.000000, 0.000000])
            s.run()
            sixsCoeffs[0,0] = float(s.outputs.values['coef_xa'])
            sixsCoeffs[0,1] = float(s.outputs.values['coef_xb'])
            sixsCoeffs[0,2] = float(s.outputs.values['coef_xc'])
            sixsCoeffs[0,3] = float(s.outputs.values['direct_solar_irradiance'])
            sixsCoeffs[0,4] = float(s.outputs.values['diffuse_solar_irradiance'])
            sixsCoeffs[0,5] = float(s.outputs.values['environmental_irradiance'])

            # Band 2 - Green
            s.wavelength = Py6S.Wavelength(0.46, 0.7, [0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000001, 0.000000, 0.000000, 0.000001, 0.000001, 0.000001, 0.000001, 0.000002, 0.000009, 0.000021, 0.000091, 0.000380, 0.006496, 0.046836, 0.401607, 0.664727, 0.768821, 0.769261, 0.805967, 0.790911, 0.806984, 0.813092, 0.834020, 0.842021, 0.792565, 0.802624, 0.859339, 0.874276, 0.874751, 0.877110, 0.891636, 0.913426, 0.954001, 0.971404, 0.969935, 0.973167, 0.945509, 0.937980, 0.984937, 1.000000, 0.961086, 0.835060, 0.428809, 0.209108, 0.058628, 0.024866, 0.007481, 0.003516, 0.001485, 0.000896, 0.000429, 0.000263, 0.000133, 0.000092, 0.000058, 0.000043, 0.000027, 0.000019, 0.000012, 0.000010, 0.000008, 0.000007, 0.000006, 0.000005, 0.000004, 0.000004, 0.000004, 0.000003, 0.000003, 0.000003, 0.000002, 0.000002, 0.000002, 0.000002, 0.000003, 0.000003, 0.000003, 0.000003, 0.000003, 0.000003, 0.000003, 0.000003, 0.000003, 0.000003, 0.000004, 0.000004, 0.000003, 0.000002, 0.000001, 0.000000, 0.000000])
            s.run()
            sixsCoeffs[1,0] = float(s.outputs.values['coef_xa'])
            sixsCoeffs[1,1] = float(s.outputs.values['coef_xb'])
            sixsCoeffs[1,2] = float(s.outputs.values['coef_xc'])
            sixsCoeffs[1,3] = float(s.outputs.values['direct_solar_irradiance'])
            sixsCoeffs[1,4] = float(s.outputs.values['diffuse_solar_irradiance'])
            sixsCoeffs[1,5] = float(s.outputs.values['environmental_irradiance'])

            # Band 3 - Red
            s.wavelength = Py6S.Wavelength(0.59, 0.7975, [0.000003, 0.000002, 0.000002, 0.000002, 0.000002, 0.000003, 0.000003, 0.000004, 0.000007, 0.000012, 0.000036, 0.000086, 0.000429, 0.001396, 0.020268, 0.103094, 0.486020, 0.765525, 0.937167, 0.941915, 0.934316, 0.937769, 0.946738, 0.943077, 0.953586, 0.957486, 0.959214, 0.968426, 0.977351, 0.980538, 0.991643, 0.991092, 0.993789, 0.994089, 0.987204, 0.986463, 0.995651, 1.000000, 0.970394, 0.857500, 0.469822, 0.228840, 0.060716, 0.027017, 0.009253, 0.004604, 0.001759, 0.000965, 0.000438, 0.000283, 0.000153, 0.000104, 0.000058, 0.000040, 0.000027, 0.000019, 0.000013, 0.000010, 0.000007, 0.000006, 0.000004, 0.000004, 0.000003, 0.000003, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000001, 0.000002, 0.000001, 0.000001, 0.000002, 0.000001, 0.000002, 0.000001, 0.000001, 0.000001, 0.000000])
            s.run()
            sixsCoeffs[2,0] = float(s.outputs.values['coef_xa'])
            sixsCoeffs[2,1] = float(s.outputs.values['coef_xb'])
            sixsCoeffs[2,2] = float(s.outputs.values['coef_xc'])
            sixsCoeffs[2,3] = float(s.outputs.values['direct_solar_irradiance'])
            sixsCoeffs[2,4] = float(s.outputs.values['diffuse_solar_irradiance'])
            sixsCoeffs[2,5] = float(s.outputs.values['environmental_irradiance'])

            # Band 4 - NIR
            s.wavelength = Py6S.Wavelength(0.69, 1.065, [0.000000, 0.000003, 0.000004, 0.000004, 0.000004, 0.000005, 0.000005, 0.000005, 0.000006, 0.000006, 0.000007, 0.000007, 0.000008, 0.000009, 0.000010, 0.000011, 0.000013, 0.000017, 0.000027, 0.000036, 0.000055, 0.000076, 0.000139, 0.000230, 0.000562, 0.001062, 0.002655, 0.004544, 0.010358, 0.018329, 0.058774, 0.133738, 0.394834, 0.623159, 0.869630, 0.959712, 1.000000, 0.993449, 0.980970, 0.971916, 0.952534, 0.938224, 0.924111, 0.916865, 0.911989, 0.912971, 0.913598, 0.911826, 0.906658, 0.899119, 0.889163, 0.880557, 0.870143, 0.867731, 0.865020, 0.866745, 0.864650, 0.859827, 0.852632, 0.846009, 0.834375, 0.826711, 0.809409, 0.797194, 0.776720, 0.761243, 0.739726, 0.727501, 0.709568, 0.694800, 0.674016, 0.659883, 0.641881, 0.632254, 0.620315, 0.613544, 0.605946, 0.601465, 0.601630, 0.602388, 0.561489, 0.439244, 0.201734, 0.098528, 0.030679, 0.015194, 0.006264, 0.003697, 0.002036, 0.001525, 0.001020, 0.000716, 0.000372, 0.000252, 0.000132, 0.000090, 0.000043, 0.000034, 0.000023, 0.000019, 0.000015, 0.000010, 0.000010, 0.000009, 0.000008, 0.000007, 0.000007, 0.000007, 0.000006, 0.000005, 0.000005, 0.000004, 0.000004, 0.000005, 0.000004, 0.000004, 0.000004, 0.000004, 0.000002, 0.000003, 0.000004, 0.000004, 0.000004, 0.000003, 0.000004, 0.000003, 0.000003, 0.000003, 0.000003, 0.000002, 0.000002, 0.000002, 0.000003, 0.000002, 0.000002, 0.000002, 0.000002, 0.000001, 0.000001, 0.000001, 0.000000, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000000, 0.000000])
            s.run()
            sixsCoeffs[3,0] = float(s.outputs.values['coef_xa'])
            sixsCoeffs[3,1] = float(s.outputs.values['coef_xb'])
            sixsCoeffs[3,2] = float(s.outputs.values['coef_xc'])
            sixsCoeffs[3,3] = float(s.outputs.values['direct_solar_irradiance'])
            sixsCoeffs[3,4] = float(s.outputs.values['diffuse_solar_irradiance'])
            sixsCoeffs[3,5] = float(s.outputs.values['environmental_irradiance'])
        elif self.bandsProd == 'Multi':
            # Band 1 - Coastal
            s.wavelength = Py6S.Wavelength(0.37, 0.49, [0.000000, 0.000000, 0.000000, 0.000000, 0.000002, 0.000008, 0.000030, 0.000093, 0.000475, 0.001464, 0.018576, 0.053293, 0.305115, 0.467192, 0.575919, 0.603239, 0.648781, 0.668505, 0.712588, 0.741351, 0.778487, 0.797011, 0.828532, 0.853300, 0.881852, 0.886541, 0.894989, 0.920607, 0.962121, 0.982343, 1.000000, 0.945330, 0.474345, 0.171196, 0.019747, 0.004716, 0.001079, 0.000416, 0.000096, 0.000060, 0.000016, 0.000009, 0.000004, 0.000003, 0.000002, 0.000001, 0.000000, 0.000000, 0.000000])
            s.run()
            sixsCoeffs[0,0] = float(s.outputs.values['coef_xa'])
            sixsCoeffs[0,1] = float(s.outputs.values['coef_xb'])
            sixsCoeffs[0,2] = float(s.outputs.values['coef_xc'])
            sixsCoeffs[0,3] = float(s.outputs.values['direct_solar_irradiance'])
            sixsCoeffs[0,4] = float(s.outputs.values['diffuse_solar_irradiance'])
            sixsCoeffs[0,5] = float(s.outputs.values['environmental_irradiance'])

            # Band 2 - Blue
            s.wavelength = Py6S.Wavelength(0.395, 0.615, [0.000000, 0.000000, 0.000001, 0.000001, 0.000002, 0.000002, 0.000003, 0.000005, 0.000007, 0.000011, 0.000016, 0.000026, 0.000046, 0.000092, 0.000187, 0.000287, 0.000757, 0.001534, 0.005096, 0.010599, 0.051236, 0.103638, 0.363869, 0.578093, 0.746604, 0.781177, 0.802868, 0.813323, 0.825085, 0.820022, 0.839875, 0.849677, 0.876742, 0.876547, 0.894562, 0.899177, 0.898469, 0.912201, 0.950024, 0.963251, 0.964915, 0.980350, 0.996347, 0.999729, 0.952372, 0.850259, 0.533870, 0.296567, 0.088810, 0.041059, 0.009926, 0.004817, 0.001574, 0.000845, 0.000389, 0.000271, 0.000162, 0.000096, 0.000044, 0.000032, 0.000022, 0.000017, 0.000012, 0.000010, 0.000008, 0.000006, 0.000004, 0.000004, 0.000003, 0.000003, 0.000003, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000003, 0.000004, 0.000010, 0.000005, 0.000000, 0.000000])
            s.run()
            sixsCoeffs[1,0] = float(s.outputs.values['coef_xa'])
            sixsCoeffs[1,1] = float(s.outputs.values['coef_xb'])
            sixsCoeffs[1,2] = float(s.outputs.values['coef_xc'])
            sixsCoeffs[1,3] = float(s.outputs.values['direct_solar_irradiance'])
            sixsCoeffs[1,4] = float(s.outputs.values['diffuse_solar_irradiance'])
            sixsCoeffs[1,5] = float(s.outputs.values['environmental_irradiance'])

            # Band 3 - Green
            s.wavelength = Py6S.Wavelength(0.46, 0.7, [0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000001, 0.000000, 0.000000, 0.000001, 0.000001, 0.000001, 0.000001, 0.000002, 0.000009, 0.000021, 0.000091, 0.000380, 0.006496, 0.046836, 0.401607, 0.664727, 0.768821, 0.769261, 0.805967, 0.790911, 0.806984, 0.813092, 0.834020, 0.842021, 0.792565, 0.802624, 0.859339, 0.874276, 0.874751, 0.877110, 0.891636, 0.913426, 0.954001, 0.971404, 0.969935, 0.973167, 0.945509, 0.937980, 0.984937, 1.000000, 0.961086, 0.835060, 0.428809, 0.209108, 0.058628, 0.024866, 0.007481, 0.003516, 0.001485, 0.000896, 0.000429, 0.000263, 0.000133, 0.000092, 0.000058, 0.000043, 0.000027, 0.000019, 0.000012, 0.000010, 0.000008, 0.000007, 0.000006, 0.000005, 0.000004, 0.000004, 0.000004, 0.000003, 0.000003, 0.000003, 0.000002, 0.000002, 0.000002, 0.000002, 0.000003, 0.000003, 0.000003, 0.000003, 0.000003, 0.000003, 0.000003, 0.000003, 0.000003, 0.000003, 0.000004, 0.000004, 0.000003, 0.000002, 0.000001, 0.000000, 0.000000])
            s.run()
            sixsCoeffs[2,0] = float(s.outputs.values['coef_xa'])
            sixsCoeffs[2,1] = float(s.outputs.values['coef_xb'])
            sixsCoeffs[2,2] = float(s.outputs.values['coef_xc'])
            sixsCoeffs[2,3] = float(s.outputs.values['direct_solar_irradiance'])
            sixsCoeffs[2,4] = float(s.outputs.values['diffuse_solar_irradiance'])
            sixsCoeffs[2,5] = float(s.outputs.values['environmental_irradiance'])

            # Band 4 - Yellow
            s.wavelength = Py6S.Wavelength(0.5, 0.985, [0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000001, 0.000000, 0.000000, 0.000001, 0.000001, 0.000000, 0.000000, 0.000000, 0.000000, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000002, 0.000002, 0.000006, 0.000011, 0.000038, 0.000109, 0.000676, 0.002448, 0.034177, 0.160065, 0.599934, 0.853146, 0.894733, 0.904336, 0.916896, 0.922086, 0.927731, 0.943046, 0.943101, 0.952921, 0.988084, 0.995714, 0.992549, 0.991442, 0.959479, 0.761264, 0.314076, 0.115597, 0.026138, 0.011110, 0.003647, 0.001828, 0.000721, 0.000416, 0.000202, 0.000125, 0.000058, 0.000037, 0.000021, 0.000015, 0.000009, 0.000007, 0.000004, 0.000003, 0.000003, 0.000002, 0.000002, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000000, 0.000001, 0.000001, 0.000001, 0.000000, 0.000001, 0.000000, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000002, 0.000002, 0.000002, 0.000003, 0.000003, 0.000004, 0.000005, 0.000007, 0.000007, 0.000009, 0.000008, 0.000010, 0.000012, 0.000032, 0.000076, 0.000276, 0.000466, 0.000827, 0.001524, 0.005704, 0.006749, 0.001569, 0.000589, 0.000234, 0.000159, 0.000128, 0.000098, 0.000044, 0.000030, 0.000010, 0.000010, 0.000005, 0.000006, 0.000008, 0.000008, 0.000009, 0.000010, 0.000006, 0.000005, 0.000004, 0.000004, 0.000006, 0.000011, 0.000026, 0.000027, 0.000010, 0.000004, 0.000003, 0.000004, 0.000012, 0.000017, 0.000006, 0.000001, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000001, 0.000002, 0.000001, 0.000001, 0.000000, 0.000000, 0.000000, 0.000001, 0.000012, 0.000004, 0.000001, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000001, 0.000003, 0.000001, 0.000001, 0.000001, 0.000082, 0.000211, 0.000000, 0.000000, 0.000000])
            s.run()
            sixsCoeffs[3,0] = float(s.outputs.values['coef_xa'])
            sixsCoeffs[3,1] = float(s.outputs.values['coef_xb'])
            sixsCoeffs[3,2] = float(s.outputs.values['coef_xc'])
            sixsCoeffs[3,3] = float(s.outputs.values['direct_solar_irradiance'])
            sixsCoeffs[3,4] = float(s.outputs.values['diffuse_solar_irradiance'])
            sixsCoeffs[3,5] = float(s.outputs.values['environmental_irradiance'])

            # Band 5 - Red
            s.wavelength = Py6S.Wavelength(0.59, 0.7975, [0.000003, 0.000002, 0.000002, 0.000002, 0.000002, 0.000003, 0.000003, 0.000004, 0.000007, 0.000012, 0.000036, 0.000086, 0.000429, 0.001396, 0.020268, 0.103094, 0.486020, 0.765525, 0.937167, 0.941915, 0.934316, 0.937769, 0.946738, 0.943077, 0.953586, 0.957486, 0.959214, 0.968426, 0.977351, 0.980538, 0.991643, 0.991092, 0.993789, 0.994089, 0.987204, 0.986463, 0.995651, 1.000000, 0.970394, 0.857500, 0.469822, 0.228840, 0.060716, 0.027017, 0.009253, 0.004604, 0.001759, 0.000965, 0.000438, 0.000283, 0.000153, 0.000104, 0.000058, 0.000040, 0.000027, 0.000019, 0.000013, 0.000010, 0.000007, 0.000006, 0.000004, 0.000004, 0.000003, 0.000003, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000001, 0.000002, 0.000001, 0.000001, 0.000002, 0.000001, 0.000002, 0.000001, 0.000001, 0.000001, 0.000000])
            s.run()
            sixsCoeffs[4,0] = float(s.outputs.values['coef_xa'])
            sixsCoeffs[4,1] = float(s.outputs.values['coef_xb'])
            sixsCoeffs[4,2] = float(s.outputs.values['coef_xc'])
            sixsCoeffs[4,3] = float(s.outputs.values['direct_solar_irradiance'])
            sixsCoeffs[4,4] = float(s.outputs.values['diffuse_solar_irradiance'])
            sixsCoeffs[4,5] = float(s.outputs.values['environmental_irradiance'])

            # Band 6 - Red Edge
            s.wavelength = Py6S.Wavelength(0.66, 0.83, [0.000000, 0.000001, 0.000001, 0.000002, 0.000003, 0.000005, 0.000011, 0.000019, 0.000046, 0.000088, 0.000235, 0.000450, 0.001230, 0.002539, 0.008798, 0.022401, 0.111046, 0.267335, 0.647785, 0.860870, 0.984207, 0.993436, 0.986575, 0.982586, 0.980137, 0.981445, 0.985978, 0.991431, 0.997118, 1.000000, 0.990438, 0.953738, 0.779389, 0.555312, 0.221307, 0.098916, 0.029256, 0.014336, 0.005491, 0.003042, 0.001310, 0.000802, 0.000382, 0.000245, 0.000133, 0.000092, 0.000055, 0.000039, 0.000025, 0.000019, 0.000013, 0.000010, 0.000006, 0.000005, 0.000004, 0.000003, 0.000002, 0.000002, 0.000002, 0.000002, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000000, 0.000000, 0.000000])
            s.run()
            sixsCoeffs[5,0] = float(s.outputs.values['coef_xa'])
            sixsCoeffs[5,1] = float(s.outputs.values['coef_xb'])
            sixsCoeffs[5,2] = float(s.outputs.values['coef_xc'])
            sixsCoeffs[5,3] = float(s.outputs.values['direct_solar_irradiance'])
            sixsCoeffs[5,4] = float(s.outputs.values['diffuse_solar_irradiance'])
            sixsCoeffs[5,5] = float(s.outputs.values['environmental_irradiance'])

            # Band 6 - NIR1
            s.wavelength = Py6S.Wavelength(0.69, 1.065, [0.000000, 0.000003, 0.000004, 0.000004, 0.000004, 0.000005, 0.000005, 0.000005, 0.000006, 0.000006, 0.000007, 0.000007, 0.000008, 0.000009, 0.000010, 0.000011, 0.000013, 0.000017, 0.000027, 0.000036, 0.000055, 0.000076, 0.000139, 0.000230, 0.000562, 0.001062, 0.002655, 0.004544, 0.010358, 0.018329, 0.058774, 0.133738, 0.394834, 0.623159, 0.869630, 0.959712, 1.000000, 0.993449, 0.980970, 0.971916, 0.952534, 0.938224, 0.924111, 0.916865, 0.911989, 0.912971, 0.913598, 0.911826, 0.906658, 0.899119, 0.889163, 0.880557, 0.870143, 0.867731, 0.865020, 0.866745, 0.864650, 0.859827, 0.852632, 0.846009, 0.834375, 0.826711, 0.809409, 0.797194, 0.776720, 0.761243, 0.739726, 0.727501, 0.709568, 0.694800, 0.674016, 0.659883, 0.641881, 0.632254, 0.620315, 0.613544, 0.605946, 0.601465, 0.601630, 0.602388, 0.561489, 0.439244, 0.201734, 0.098528, 0.030679, 0.015194, 0.006264, 0.003697, 0.002036, 0.001525, 0.001020, 0.000716, 0.000372, 0.000252, 0.000132, 0.000090, 0.000043, 0.000034, 0.000023, 0.000019, 0.000015, 0.000010, 0.000010, 0.000009, 0.000008, 0.000007, 0.000007, 0.000007, 0.000006, 0.000005, 0.000005, 0.000004, 0.000004, 0.000005, 0.000004, 0.000004, 0.000004, 0.000004, 0.000002, 0.000003, 0.000004, 0.000004, 0.000004, 0.000003, 0.000004, 0.000003, 0.000003, 0.000003, 0.000003, 0.000002, 0.000002, 0.000002, 0.000003, 0.000002, 0.000002, 0.000002, 0.000002, 0.000001, 0.000001, 0.000001, 0.000000, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000000, 0.000000])
            s.run()
            sixsCoeffs[6,0] = float(s.outputs.values['coef_xa'])
            sixsCoeffs[6,1] = float(s.outputs.values['coef_xb'])
            sixsCoeffs[6,2] = float(s.outputs.values['coef_xc'])
            sixsCoeffs[6,3] = float(s.outputs.values['direct_solar_irradiance'])
            sixsCoeffs[6,4] = float(s.outputs.values['diffuse_solar_irradiance'])
            sixsCoeffs[6,5] = float(s.outputs.values['environmental_irradiance'])

            # Band 8 - NIR2
            s.wavelength = Py6S.Wavelength(0.8, 1.095, [0.000000, 0.000000, 0.000000, 0.000001, 0.000000, 0.000001, 0.000000, 0.000001, 0.000001, 0.000000, 0.000000, 0.000000, 0.000001, 0.000002, 0.000003, 0.000004, 0.000017, 0.000024, 0.000086, 0.000193, 0.005004, 0.006248, 0.015480, 0.040770, 0.198451, 0.490378, 0.898905, 0.988771, 0.995955, 0.980427, 0.953102, 0.936267, 0.913444, 0.896827, 0.871787, 0.855441, 0.830337, 0.812005, 0.795488, 0.783124, 0.767600, 0.753668, 0.732714, 0.723220, 0.710427, 0.701802, 0.685938, 0.676788, 0.660788, 0.651630, 0.636879, 0.625875, 0.608110, 0.596127, 0.579630, 0.567198, 0.548800, 0.536805, 0.520057, 0.507185, 0.486073, 0.472972, 0.454661, 0.441460, 0.423699, 0.411377, 0.391904, 0.379373, 0.362515, 0.348911, 0.331009, 0.319164, 0.303764, 0.293955, 0.279563, 0.271664, 0.259319, 0.250665, 0.237601, 0.229018, 0.218230, 0.209025, 0.198132, 0.191910, 0.180931, 0.174822, 0.163719, 0.156087, 0.148322, 0.141471, 0.133976, 0.128863, 0.120766, 0.115362, 0.098007, 0.077797, 0.041083, 0.022039, 0.008047, 0.004068, 0.001705, 0.000912, 0.000493, 0.000354, 0.000170, 0.000080, 0.000055, 0.000030, 0.000019, 0.000014, 0.000007, 0.000005, 0.000002, 0.000002, 0.000001, 0.000001, 0.000001, 0.000000, 0.000000])
            s.run()
            sixsCoeffs[7,0] = float(s.outputs.values['coef_xa'])
            sixsCoeffs[7,1] = float(s.outputs.values['coef_xb'])
            sixsCoeffs[7,2] = float(s.outputs.values['coef_xc'])
            sixsCoeffs[7,3] = float(s.outputs.values['direct_solar_irradiance'])
            sixsCoeffs[7,4] = float(s.outputs.values['diffuse_solar_irradiance'])
            sixsCoeffs[7,5] = float(s.outputs.values['environmental_irradiance'])
        elif self.bandsProd == 'All-S':
            # Band 1 - SWIR1
            s.wavelength = Py6S.Wavelength(1.13, 1.3, [0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000253, 0.000293, 0.000303, 0.000390, 0.000426, 0.000535, 0.000635, 0.000843, 0.001035, 0.001458, 0.001930, 0.002951, 0.004041, 0.006816, 0.010050, 0.018794, 0.029201, 0.060660, 0.101829, 0.219839, 0.358152, 0.635030, 0.804403, 0.952379, 0.980083, 0.990362, 0.990603, 0.988023, 0.987617, 0.994926, 1.000000, 0.953605, 0.825556, 0.537288, 0.341430, 0.159654, 0.096342, 0.046975, 0.030120, 0.016281, 0.011236, 0.006716, 0.004922, 0.003226, 0.002464, 0.001736, 0.001395, 0.001029, 0.000905, 0.000702, 0.000582, 0.000504, 0.000448, 0.000375, 0.000375, 0.000320, 0.000312, 0.000283, 0.000273, 0.000294, 0.000255, 0.000000, 0.000000, 0.000000])
            s.run()
            sixsCoeffs[0,0] = float(s.outputs.values['coef_xa'])
            sixsCoeffs[0,1] = float(s.outputs.values['coef_xb'])
            sixsCoeffs[0,2] = float(s.outputs.values['coef_xc'])
            sixsCoeffs[0,3] = float(s.outputs.values['direct_solar_irradiance'])
            sixsCoeffs[0,4] = float(s.outputs.values['diffuse_solar_irradiance'])
            sixsCoeffs[0,5] = float(s.outputs.values['environmental_irradiance'])

            # Band 2 - SWIR2
            s.wavelength = Py6S.Wavelength(1.48, 1.65, [0.000000, 0.000000, 0.000000, 0.000000, 0.000081, 0.000080, 0.000088, 0.000094, 0.000080, 0.000087, 0.000092, 0.000090, 0.000099, 0.000088, 0.000090, 0.000086, 0.000075, 0.000094, 0.000128, 0.000167, 0.000218, 0.000305, 0.000623, 0.001168, 0.003351, 0.007608, 0.031037, 0.081286, 0.285810, 0.548184, 0.909531, 0.931296, 0.894925, 0.907092, 0.947587, 0.970979, 0.988443, 0.988939, 0.985879, 0.989478, 0.997065, 0.996143, 1.000000, 0.980305, 0.725741, 0.488847, 0.152377, 0.060699, 0.015334, 0.006797, 0.002333, 0.001267, 0.000554, 0.000354, 0.000232, 0.000170, 0.000137, 0.000120, 0.000105, 0.000102, 0.000130, 0.000105, 0.000090, 0.000085, 0.000114, 0.000000, 0.000000, 0.000000, 0.000000])
            s.run()
            sixsCoeffs[1,0] = float(s.outputs.values['coef_xa'])
            sixsCoeffs[1,1] = float(s.outputs.values['coef_xb'])
            sixsCoeffs[1,2] = float(s.outputs.values['coef_xc'])
            sixsCoeffs[1,3] = float(s.outputs.values['direct_solar_irradiance'])
            sixsCoeffs[1,4] = float(s.outputs.values['diffuse_solar_irradiance'])
            sixsCoeffs[1,5] = float(s.outputs.values['environmental_irradiance'])

            # Band 3 - SWIR3
            s.wavelength = Py6S.Wavelength(1.585, 1.74, [0.000000, 0.000000, 0.000036, 0.000020, 0.000000, 0.000046, 0.000024, 0.000000, 0.000000, 0.000000, 0.000000, 0.000008, 0.000128, 0.000032, 0.000155, 0.000287, 0.000542, 0.001067, 0.003346, 0.007835, 0.033883, 0.086584, 0.288351, 0.528870, 0.792966, 0.887712, 0.896967, 0.908576, 0.924562, 0.932828, 0.947428, 0.957944, 0.977268, 0.992391, 0.998861, 0.988963, 0.978623, 0.899237, 0.553352, 0.269970, 0.066549, 0.024358, 0.006187, 0.002869, 0.001056, 0.000530, 0.000294, 0.000130, 0.000027, 0.000147, 0.000019, 0.000058, 0.000034, 0.000080, 0.000040, 0.000000, 0.000000, 0.000000, 0.000009, 0.000000, 0.000088, 0.000000, 0.000000])
            s.run()
            sixsCoeffs[2,0] = float(s.outputs.values['coef_xa'])
            sixsCoeffs[2,1] = float(s.outputs.values['coef_xb'])
            sixsCoeffs[2,2] = float(s.outputs.values['coef_xc'])
            sixsCoeffs[2,3] = float(s.outputs.values['direct_solar_irradiance'])
            sixsCoeffs[2,4] = float(s.outputs.values['diffuse_solar_irradiance'])
            sixsCoeffs[2,5] = float(s.outputs.values['environmental_irradiance'])

            # Band 4 - SWIR4
            s.wavelength = Py6S.Wavelength(1.655, 1.825, [0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000067, 0.000106, 0.000161, 0.000176, 0.000262, 0.000375, 0.000698, 0.000970, 0.001851, 0.002970, 0.006362, 0.011125, 0.028352, 0.056000, 0.141859, 0.266281, 0.574651, 0.769154, 0.965049, 0.961823, 0.951606, 0.962861, 0.985168, 0.993860, 0.999488, 0.999739, 0.999678, 0.996571, 0.986755, 0.978081, 0.933916, 0.799587, 0.564041, 0.365284, 0.150568, 0.081784, 0.032797, 0.018221, 0.008213, 0.005044, 0.002552, 0.001681, 0.000911, 0.000607, 0.000392, 0.000288, 0.000149, 0.000117, 0.000110, 0.000037, 0.000064, 0.000062, 0.000077, 0.000022, 0.000025, 0.000015, 0.000000, 0.000015, 0.000059, 0.000000, 0.000000, 0.000000, 0.000000])
            s.run()
            sixsCoeffs[3,0] = float(s.outputs.values['coef_xa'])
            sixsCoeffs[3,1] = float(s.outputs.values['coef_xb'])
            sixsCoeffs[3,2] = float(s.outputs.values['coef_xc'])
            sixsCoeffs[3,3] = float(s.outputs.values['direct_solar_irradiance'])
            sixsCoeffs[3,4] = float(s.outputs.values['diffuse_solar_irradiance'])
            sixsCoeffs[3,5] = float(s.outputs.values['environmental_irradiance'])

            # Band 5 - SWIR5
            s.wavelength = Py6S.Wavelength(2.085, 2.235, [0.000000, 0.000000, 0.000000, 0.000052, 0.000000, 0.000000, 0.000000, 0.000083, 0.000019, 0.000055, 0.000022, 0.000242, 0.000118, 0.000262, 0.000531, 0.000744, 0.001409, 0.002391, 0.005357, 0.009834, 0.028200, 0.059326, 0.153355, 0.268981, 0.528209, 0.687011, 0.896047, 0.954417, 0.962773, 0.968441, 0.984622, 0.994119, 0.999917, 0.992653, 0.970129, 0.950530, 0.943482, 0.948859, 0.817985, 0.668471, 0.399201, 0.214359, 0.074121, 0.035799, 0.012435, 0.006813, 0.002946, 0.001824, 0.000910, 0.000505, 0.000264, 0.000308, 0.000073, 0.000091, 0.000162, 0.000025, 0.000021, 0.000000, 0.000000, 0.000002, 0.000000])
            s.run()
            sixsCoeffs[4,0] = float(s.outputs.values['coef_xa'])
            sixsCoeffs[4,1] = float(s.outputs.values['coef_xb'])
            sixsCoeffs[4,2] = float(s.outputs.values['coef_xc'])
            sixsCoeffs[4,3] = float(s.outputs.values['direct_solar_irradiance'])
            sixsCoeffs[4,4] = float(s.outputs.values['diffuse_solar_irradiance'])
            sixsCoeffs[4,5] = float(s.outputs.values['environmental_irradiance'])

            # Band 6 - SWIR6
            s.wavelength = Py6S.Wavelength(2.13, 2.29, [0.000000, 0.000000, 0.000007, 0.000039, 0.000059, 0.000047, 0.000073, 0.000122, 0.000142, 0.000235, 0.000379, 0.000577, 0.001098, 0.001687, 0.003320, 0.005485, 0.012517, 0.023562, 0.065047, 0.114972, 0.268342, 0.439762, 0.719972, 0.891663, 0.998379, 0.998512, 0.995415, 0.995654, 0.991513, 0.990272, 0.982527, 0.968473, 0.934127, 0.895251, 0.861456, 0.860045, 0.774778, 0.667452, 0.448839, 0.262043, 0.100543, 0.052871, 0.019373, 0.010589, 0.004691, 0.002856, 0.001507, 0.000993, 0.000540, 0.000418, 0.000260, 0.000236, 0.000075, 0.000091, 0.000039, 0.000084, 0.000037, 0.000000, 0.000001, 0.000000, 0.000064, 0.000001, 0.000000, 0.000000, 0.000000])
            s.run()
            sixsCoeffs[5,0] = float(s.outputs.values['coef_xa'])
            sixsCoeffs[5,1] = float(s.outputs.values['coef_xb'])
            sixsCoeffs[5,2] = float(s.outputs.values['coef_xc'])
            sixsCoeffs[5,3] = float(s.outputs.values['direct_solar_irradiance'])
            sixsCoeffs[5,4] = float(s.outputs.values['diffuse_solar_irradiance'])
            sixsCoeffs[5,5] = float(s.outputs.values['environmental_irradiance'])

            # Band 7 - SWIR7
            s.wavelength = Py6S.Wavelength(2.15, 2.36, [0.000000, 0.000000, 0.000000, 0.000005, 0.000033, 0.000027, 0.000008, 0.000003, 0.000004, 0.000000, 0.000051, 0.000021, 0.000001, 0.000020, 0.000005, 0.000004, 0.000000, 0.000028, 0.000039, 0.000061, 0.000028, 0.000040, 0.000086, 0.000104, 0.000147, 0.000255, 0.000557, 0.001014, 0.002501, 0.004978, 0.015550, 0.036622, 0.111236, 0.210157, 0.458375, 0.600908, 0.769271, 0.892668, 0.948770, 0.962194, 0.984210, 0.992194, 0.996299, 0.998429, 1.000000, 0.997832, 0.991417, 0.986279, 0.976906, 0.966125, 0.948347, 0.928221, 0.772989, 0.653959, 0.418250, 0.237462, 0.086496, 0.041015, 0.012498, 0.006198, 0.002371, 0.001357, 0.000649, 0.000371, 0.000204, 0.000136, 0.000104, 0.000118, 0.000037, 0.000035, 0.000033, 0.000041, 0.000048, 0.000039, 0.000019, 0.000000, 0.000000, 0.000044, 0.000000, 0.000000, 0.000012, 0.000018, 0.000052, 0.000000, 0.000000])
            s.run()
            sixsCoeffs[6,0] = float(s.outputs.values['coef_xa'])
            sixsCoeffs[6,1] = float(s.outputs.values['coef_xb'])
            sixsCoeffs[6,2] = float(s.outputs.values['coef_xc'])
            sixsCoeffs[6,3] = float(s.outputs.values['direct_solar_irradiance'])
            sixsCoeffs[6,4] = float(s.outputs.values['diffuse_solar_irradiance'])
            sixsCoeffs[6,5] = float(s.outputs.values['environmental_irradiance'])

            # Band 8 - SWIR8
            s.wavelength = Py6S.Wavelength(2.23, 2.44, [0.000000, 0.000000, 0.000219, 0.000221, 0.000198, 0.000232, 0.000225, 0.000174, 0.000216, 0.000218, 0.000250, 0.000262, 0.000290, 0.000316, 0.000369, 0.000522, 0.000862, 0.001256, 0.002483, 0.004211, 0.009952, 0.018347, 0.047885, 0.080950, 0.166367, 0.258357, 0.445946, 0.587139, 0.766937, 0.856844, 0.934044, 0.958367, 0.978409, 0.986971, 0.994644, 0.997807, 0.999978, 0.999632, 0.997006, 0.994443, 0.990477, 0.989680, 0.989941, 0.989197, 0.983878, 0.976414, 0.956834, 0.936535, 0.905523, 0.892770, 0.899677, 0.914064, 0.820145, 0.699476, 0.513352, 0.329919, 0.134818, 0.072645, 0.027000, 0.014156, 0.005814, 0.003422, 0.001681, 0.001093, 0.000634, 0.000495, 0.000342, 0.000314, 0.000245, 0.000218, 0.000216, 0.000228, 0.000184, 0.000168, 0.000203, 0.000170, 0.000153, 0.000200, 0.000189, 0.000164, 0.000175, 0.000173, 0.000208, 0.000000, 0.000000])
            s.run()
            sixsCoeffs[7,0] = float(s.outputs.values['coef_xa'])
            sixsCoeffs[7,1] = float(s.outputs.values['coef_xb'])
            sixsCoeffs[7,2] = float(s.outputs.values['coef_xc'])
            sixsCoeffs[7,3] = float(s.outputs.values['direct_solar_irradiance'])
            sixsCoeffs[7,4] = float(s.outputs.values['diffuse_solar_irradiance'])
            sixsCoeffs[7,5] = float(s.outputs.values['environmental_irradiance'])
        else: 
            raise ARCSIException("Do not reconise the product type.")
        return sixsCoeffs

    def convertImageToSurfaceReflSglParam(self, inputRadImage, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF, scaleFactor):
        print("Converting to Surface Reflectance")
        outputImage = os.path.join(outputPath, outputName)

        imgBandCoeffs = list()
        sixsCoeffs = self.calc6SCoefficients(aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF)

        if self.bandsProd == 'MS1':
            imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=1, aX=float(sixsCoeffs[0,0]), bX=float(sixsCoeffs[0,1]), cX=float(sixsCoeffs[0,2]), DirIrr=float(sixsCoeffs[0,3]), DifIrr=float(sixsCoeffs[0,4]), EnvIrr=float(sixsCoeffs[0,5])))
            imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=2, aX=float(sixsCoeffs[1,0]), bX=float(sixsCoeffs[1,1]), cX=float(sixsCoeffs[1,2]), DirIrr=float(sixsCoeffs[1,3]), DifIrr=float(sixsCoeffs[1,4]), EnvIrr=float(sixsCoeffs[1,5])))
            imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=3, aX=float(sixsCoeffs[2,0]), bX=float(sixsCoeffs[2,1]), cX=float(sixsCoeffs[2,2]), DirIrr=float(sixsCoeffs[2,3]), DifIrr=float(sixsCoeffs[2,4]), EnvIrr=float(sixsCoeffs[2,5])))
            imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=4, aX=float(sixsCoeffs[3,0]), bX=float(sixsCoeffs[3,1]), cX=float(sixsCoeffs[3,2]), DirIrr=float(sixsCoeffs[3,3]), DifIrr=float(sixsCoeffs[3,4]), EnvIrr=float(sixsCoeffs[3,5])))
        elif (self.bandsProd == 'Multi') or (self.bandsProd == 'All-S'):
            imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=1, aX=float(sixsCoeffs[0,0]), bX=float(sixsCoeffs[0,1]), cX=float(sixsCoeffs[0,2]), DirIrr=float(sixsCoeffs[0,3]), DifIrr=float(sixsCoeffs[0,4]), EnvIrr=float(sixsCoeffs[0,5])))
            imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=2, aX=float(sixsCoeffs[1,0]), bX=float(sixsCoeffs[1,1]), cX=float(sixsCoeffs[1,2]), DirIrr=float(sixsCoeffs[1,3]), DifIrr=float(sixsCoeffs[1,4]), EnvIrr=float(sixsCoeffs[1,5])))
            imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=3, aX=float(sixsCoeffs[2,0]), bX=float(sixsCoeffs[2,1]), cX=float(sixsCoeffs[2,2]), DirIrr=float(sixsCoeffs[2,3]), DifIrr=float(sixsCoeffs[2,4]), EnvIrr=float(sixsCoeffs[2,5])))
            imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=4, aX=float(sixsCoeffs[3,0]), bX=float(sixsCoeffs[3,1]), cX=float(sixsCoeffs[3,2]), DirIrr=float(sixsCoeffs[3,3]), DifIrr=float(sixsCoeffs[3,4]), EnvIrr=float(sixsCoeffs[3,5])))
            imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=5, aX=float(sixsCoeffs[4,0]), bX=float(sixsCoeffs[4,1]), cX=float(sixsCoeffs[4,2]), DirIrr=float(sixsCoeffs[4,3]), DifIrr=float(sixsCoeffs[4,4]), EnvIrr=float(sixsCoeffs[4,5])))
            imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=6, aX=float(sixsCoeffs[5,0]), bX=float(sixsCoeffs[5,1]), cX=float(sixsCoeffs[5,2]), DirIrr=float(sixsCoeffs[5,3]), DifIrr=float(sixsCoeffs[5,4]), EnvIrr=float(sixsCoeffs[5,5])))
            imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=7, aX=float(sixsCoeffs[6,0]), bX=float(sixsCoeffs[6,1]), cX=float(sixsCoeffs[6,2]), DirIrr=float(sixsCoeffs[6,3]), DifIrr=float(sixsCoeffs[6,4]), EnvIrr=float(sixsCoeffs[6,5])))
            imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=8, aX=float(sixsCoeffs[7,0]), bX=float(sixsCoeffs[7,1]), cX=float(sixsCoeffs[7,2]), DirIrr=float(sixsCoeffs[7,3]), DifIrr=float(sixsCoeffs[7,4]), EnvIrr=float(sixsCoeffs[7,5])))
        else: 
            raise ARCSIException("Do not reconise the product type.")

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
            if self.bandsProd == 'MS1':
                imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=1, aX=float(sixsCoeffs[0,0]), bX=float(sixsCoeffs[0,1]), cX=float(sixsCoeffs[0,2]), DirIrr=float(sixsCoeffs[0,3]), DifIrr=float(sixsCoeffs[0,4]), EnvIrr=float(sixsCoeffs[0,5])))
                imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=2, aX=float(sixsCoeffs[1,0]), bX=float(sixsCoeffs[1,1]), cX=float(sixsCoeffs[1,2]), DirIrr=float(sixsCoeffs[1,3]), DifIrr=float(sixsCoeffs[1,4]), EnvIrr=float(sixsCoeffs[1,5])))
                imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=3, aX=float(sixsCoeffs[2,0]), bX=float(sixsCoeffs[2,1]), cX=float(sixsCoeffs[2,2]), DirIrr=float(sixsCoeffs[2,3]), DifIrr=float(sixsCoeffs[2,4]), EnvIrr=float(sixsCoeffs[2,5])))
                imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=4, aX=float(sixsCoeffs[3,0]), bX=float(sixsCoeffs[3,1]), cX=float(sixsCoeffs[3,2]), DirIrr=float(sixsCoeffs[3,3]), DifIrr=float(sixsCoeffs[3,4]), EnvIrr=float(sixsCoeffs[3,5])))
            elif (self.bandsProd == 'Multi') or (self.bandsProd == 'All-S'):
                imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=1, aX=float(sixsCoeffs[0,0]), bX=float(sixsCoeffs[0,1]), cX=float(sixsCoeffs[0,2]), DirIrr=float(sixsCoeffs[0,3]), DifIrr=float(sixsCoeffs[0,4]), EnvIrr=float(sixsCoeffs[0,5])))
                imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=2, aX=float(sixsCoeffs[1,0]), bX=float(sixsCoeffs[1,1]), cX=float(sixsCoeffs[1,2]), DirIrr=float(sixsCoeffs[1,3]), DifIrr=float(sixsCoeffs[1,4]), EnvIrr=float(sixsCoeffs[1,5])))
                imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=3, aX=float(sixsCoeffs[2,0]), bX=float(sixsCoeffs[2,1]), cX=float(sixsCoeffs[2,2]), DirIrr=float(sixsCoeffs[2,3]), DifIrr=float(sixsCoeffs[2,4]), EnvIrr=float(sixsCoeffs[2,5])))
                imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=4, aX=float(sixsCoeffs[3,0]), bX=float(sixsCoeffs[3,1]), cX=float(sixsCoeffs[3,2]), DirIrr=float(sixsCoeffs[3,3]), DifIrr=float(sixsCoeffs[3,4]), EnvIrr=float(sixsCoeffs[3,5])))
                imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=5, aX=float(sixsCoeffs[4,0]), bX=float(sixsCoeffs[4,1]), cX=float(sixsCoeffs[4,2]), DirIrr=float(sixsCoeffs[4,3]), DifIrr=float(sixsCoeffs[4,4]), EnvIrr=float(sixsCoeffs[4,5])))
                imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=6, aX=float(sixsCoeffs[5,0]), bX=float(sixsCoeffs[5,1]), cX=float(sixsCoeffs[5,2]), DirIrr=float(sixsCoeffs[5,3]), DifIrr=float(sixsCoeffs[5,4]), EnvIrr=float(sixsCoeffs[5,5])))
                imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=7, aX=float(sixsCoeffs[6,0]), bX=float(sixsCoeffs[6,1]), cX=float(sixsCoeffs[6,2]), DirIrr=float(sixsCoeffs[6,3]), DifIrr=float(sixsCoeffs[6,4]), EnvIrr=float(sixsCoeffs[6,5])))
                imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=8, aX=float(sixsCoeffs[7,0]), bX=float(sixsCoeffs[7,1]), cX=float(sixsCoeffs[7,2]), DirIrr=float(sixsCoeffs[7,3]), DifIrr=float(sixsCoeffs[7,4]), EnvIrr=float(sixsCoeffs[7,5])))
            else: 
                raise ARCSIException("Do not reconise the product type.")

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
                    if self.bandsProd == 'MS1':
                        imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=1, aX=float(sixsCoeffs[0,0]), bX=float(sixsCoeffs[0,1]), cX=float(sixsCoeffs[0,2]), DirIrr=float(sixsCoeffs[0,3]), DifIrr=float(sixsCoeffs[0,4]), EnvIrr=float(sixsCoeffs[0,5])))
                        imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=2, aX=float(sixsCoeffs[1,0]), bX=float(sixsCoeffs[1,1]), cX=float(sixsCoeffs[1,2]), DirIrr=float(sixsCoeffs[1,3]), DifIrr=float(sixsCoeffs[1,4]), EnvIrr=float(sixsCoeffs[1,5])))
                        imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=3, aX=float(sixsCoeffs[2,0]), bX=float(sixsCoeffs[2,1]), cX=float(sixsCoeffs[2,2]), DirIrr=float(sixsCoeffs[2,3]), DifIrr=float(sixsCoeffs[2,4]), EnvIrr=float(sixsCoeffs[2,5])))
                        imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=4, aX=float(sixsCoeffs[3,0]), bX=float(sixsCoeffs[3,1]), cX=float(sixsCoeffs[3,2]), DirIrr=float(sixsCoeffs[3,3]), DifIrr=float(sixsCoeffs[3,4]), EnvIrr=float(sixsCoeffs[3,5])))
                    elif (self.bandsProd == 'Multi') or (self.bandsProd == 'All-S'):
                        imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=1, aX=float(sixsCoeffs[0,0]), bX=float(sixsCoeffs[0,1]), cX=float(sixsCoeffs[0,2]), DirIrr=float(sixsCoeffs[0,3]), DifIrr=float(sixsCoeffs[0,4]), EnvIrr=float(sixsCoeffs[0,5])))
                        imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=2, aX=float(sixsCoeffs[1,0]), bX=float(sixsCoeffs[1,1]), cX=float(sixsCoeffs[1,2]), DirIrr=float(sixsCoeffs[1,3]), DifIrr=float(sixsCoeffs[1,4]), EnvIrr=float(sixsCoeffs[1,5])))
                        imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=3, aX=float(sixsCoeffs[2,0]), bX=float(sixsCoeffs[2,1]), cX=float(sixsCoeffs[2,2]), DirIrr=float(sixsCoeffs[2,3]), DifIrr=float(sixsCoeffs[2,4]), EnvIrr=float(sixsCoeffs[2,5])))
                        imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=4, aX=float(sixsCoeffs[3,0]), bX=float(sixsCoeffs[3,1]), cX=float(sixsCoeffs[3,2]), DirIrr=float(sixsCoeffs[3,3]), DifIrr=float(sixsCoeffs[3,4]), EnvIrr=float(sixsCoeffs[3,5])))
                        imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=5, aX=float(sixsCoeffs[4,0]), bX=float(sixsCoeffs[4,1]), cX=float(sixsCoeffs[4,2]), DirIrr=float(sixsCoeffs[4,3]), DifIrr=float(sixsCoeffs[4,4]), EnvIrr=float(sixsCoeffs[4,5])))
                        imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=6, aX=float(sixsCoeffs[5,0]), bX=float(sixsCoeffs[5,1]), cX=float(sixsCoeffs[5,2]), DirIrr=float(sixsCoeffs[5,3]), DifIrr=float(sixsCoeffs[5,4]), EnvIrr=float(sixsCoeffs[5,5])))
                        imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=7, aX=float(sixsCoeffs[6,0]), bX=float(sixsCoeffs[6,1]), cX=float(sixsCoeffs[6,2]), DirIrr=float(sixsCoeffs[6,3]), DifIrr=float(sixsCoeffs[6,4]), EnvIrr=float(sixsCoeffs[6,5])))
                        imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=8, aX=float(sixsCoeffs[7,0]), bX=float(sixsCoeffs[7,1]), cX=float(sixsCoeffs[7,2]), DirIrr=float(sixsCoeffs[7,3]), DifIrr=float(sixsCoeffs[7,4]), EnvIrr=float(sixsCoeffs[7,5])))
                    else: 
                        raise ARCSIException("Do not reconise the product type.")
                    aot6SCoeffsOut.append(rsgislib.imagecalibration.AOTLUTFeat(AOT=float(aotVal), Coeffs=imgBandCoeffs))
                elevAOTCoeffs.append(rsgislib.imagecalibration.ElevLUTFeat(Elev=float(elevVal), Coeffs=aot6SCoeffsOut))

        rsgislib.imagecalibration.apply6SCoeffElevAOTLUTParam(inputRadImage, inputDEMFile, inputAOTImage, outputImage, outFormat, rsgislib.TYPE_16UINT, scaleFactor, 0, True, elevAOTCoeffs)

        return outputImage, elevAOTCoeffs

    def run6SToOptimiseAODValue(self, aotVal, radBlueVal, predBlueVal, aeroProfile, atmosProfile, grdRefl, surfaceAltitude):
        """Used as part of the optimastion for identifying values of AOD"""
        if self.bandsProd == 'All-S': 
            raise ARCSIException("Cannot estimate AOD value without a visual image band - AOT is not so important for SWIR bands so suggest you just enter a constant (e.g., 0.05).")

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

        # Band 2 - Blue
        s.wavelength = Py6S.Wavelength(0.395, 0.615, [0.000000, 0.000000, 0.000001, 0.000001, 0.000002, 0.000002, 0.000003, 0.000005, 0.000007, 0.000011, 0.000016, 0.000026, 0.000046, 0.000092, 0.000187, 0.000287, 0.000757, 0.001534, 0.005096, 0.010599, 0.051236, 0.103638, 0.363869, 0.578093, 0.746604, 0.781177, 0.802868, 0.813323, 0.825085, 0.820022, 0.839875, 0.849677, 0.876742, 0.876547, 0.894562, 0.899177, 0.898469, 0.912201, 0.950024, 0.963251, 0.964915, 0.980350, 0.996347, 0.999729, 0.952372, 0.850259, 0.533870, 0.296567, 0.088810, 0.041059, 0.009926, 0.004817, 0.001574, 0.000845, 0.000389, 0.000271, 0.000162, 0.000096, 0.000044, 0.000032, 0.000022, 0.000017, 0.000012, 0.000010, 0.000008, 0.000006, 0.000004, 0.000004, 0.000003, 0.000003, 0.000003, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000003, 0.000004, 0.000010, 0.000005, 0.000000, 0.000000])
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
            if self.bandsProd == 'All-S': 
                raise ARCSIException("Cannot estimate AOD value without a visual image band - AOT is not so important for SWIR bands so suggest you just enter a constant (e.g., 0.05).")

            print("Estimating AOD Using DOS")
            arcsiUtils = ARCSIUtils()
            outputAOTImage = os.path.join(outputPath, outputName)
            tmpBaseName = os.path.splitext(outputName)[0]
            imgExtension = arcsiUtils.getFileExtension(outFormat)

            blueBand = 2
            redBand = 5
            nirNamd = 7
            segBands = [8,6,1]
            if self.bandsProd == 'MS1':
                blueBand = 1
                redBand = 3
                nirNamd = 4
                segBands = [1,2,3,4]
            elif self.bandsProd == 'Multi':
                blueBand = 2
                redBand = 5
                nirNamd = 7
                segBands = [8,6,1]
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
            rsgislib.rastergis.selectClumpsOnGrid(thresImageClumpsFinal, "Selected", "PredictAOTFor", "Eastings", "Northings", "MeanBlueDOS", "min", 20, 20)

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

            aotVals = numpy.zeros_like(MeanBlueRAD, dtype=numpy.float)

            for i in range(len(MeanBlueRAD)):
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
            blueBand = 2
            if self.bandsProd == 'MS1':
                blueBand = 1
            elif self.bandsProd == 'Multi':
                blueBand = 2
            elif self.bandsProd == 'All-S': 
                raise ARCSIException("Cannot estimate AOD value without a visual image band - AOT is not so important for SWIR bands so suggest you just enter a constant (e.g., 0.05).")
            else: 
                raise ARCSIException("Do not reconise the product type.")
            return self.estimateSingleAOTFromDOSBandImpl(radianceImage, toaImage, inputDEMFile, tmpPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, minAOT, maxAOT, dosOutRefl, blueBand)
        except Exception as e:
            raise

    def setBandNames(self, imageFile):
        dataset = gdal.Open(imageFile, gdal.GA_Update)
        if not dataset is None:
            if self.bandsProd == 'MS1':
                dataset.GetRasterBand(1).SetDescription("Blue")
                dataset.GetRasterBand(2).SetDescription("Green")
                dataset.GetRasterBand(3).SetDescription("Red")
                dataset.GetRasterBand(4).SetDescription("NIR")
            elif self.bandsProd == 'Multi':
                dataset.GetRasterBand(1).SetDescription("Coastal")
                dataset.GetRasterBand(2).SetDescription("Blue")
                dataset.GetRasterBand(3).SetDescription("Green")
                dataset.GetRasterBand(4).SetDescription("Yellow")
                dataset.GetRasterBand(5).SetDescription("Red")
                dataset.GetRasterBand(6).SetDescription("RedEdge")
                dataset.GetRasterBand(7).SetDescription("NIR1")
                dataset.GetRasterBand(8).SetDescription("NIR2")
            elif self.bandsProd == 'All-S':
                dataset.GetRasterBand(1).SetDescription("SWIR1")
                dataset.GetRasterBand(2).SetDescription("SWIR2")
                dataset.GetRasterBand(3).SetDescription("SWIR3")
                dataset.GetRasterBand(4).SetDescription("SWIR4")
                dataset.GetRasterBand(5).SetDescription("SWIR5")
                dataset.GetRasterBand(6).SetDescription("SWIR6")
                dataset.GetRasterBand(7).SetDescription("SWIR7")
                dataset.GetRasterBand(8).SetDescription("SWIR8")
            dataset = None
        else:
            print("Could not open image to set band names: ", imageFile)

    def cleanLocalFollowProcessing(self):
        print("")


