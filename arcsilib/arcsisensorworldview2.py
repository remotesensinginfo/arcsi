"""
Module that contains the ARCSIWorldView2Sensor class.
"""
############################################################################
#  arcsisensorworldview2.py
#
#  Copyright 2014 ARCSI.
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
# Purpose:  A class for read the WorldView2 sensor header file and applying
#           the pre-processing operations within ARCSI to the WorldView2
#           datasets.
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 20/06/2014
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

class ARCSIWorldView2Sensor (ARCSIAbstractSensor):
    """
    A class which represents the WorldView2 sensor to read
    header parameters and apply data processing operations.
    """
    def __init__(self, debugMode, inputImage):
        ARCSIAbstractSensor.__init__(self, debugMode, inputImage)
        self.sensor = "WV2"
        self.fileName = ""
        self.catID = ""
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

    def extractHeaderParameters(self, inputHeader, wktStr):
        """
        Understands and parses the WorldView2 xml header file
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

            if imageInfoTag.find('SATID').text.strip() != 'WV02':
                raise ARCSIException("SATID is not \'WV02\' but \'" + imageInfoTag.find('SATID').text.strip() + "\'")

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
        Customises the generic name for the WorldView2 sensor
        """
        outname = self.defaultGenBaseOutFileName()
        return outname

    def expectedImageDataPresent(self):
        imageDataPresent = True

        if not os.path.exists(self.fileName):
            imageDataPresent = False

        return imageDataPresent

    def applyImageDataMask(self, inputHeader, outputPath, outputMaskName, outputImgName, outFormat, outWKTFile):
        raise ARCSIException("WorldView2 does not provide any image masks, do not use the MASK option.")

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
        bandDefnSeq.append(wv2Band(bandName="Coastal", bandIndex=1, absCalFact=self.absCalFactB1, effBandWidth=self.effBandWidthB1))
        bandDefnSeq.append(wv2Band(bandName="Blue", bandIndex=2, absCalFact=self.absCalFactB2, effBandWidth=self.effBandWidthB2))
        bandDefnSeq.append(wv2Band(bandName="Green", bandIndex=3, absCalFact=self.absCalFactB3, effBandWidth=self.effBandWidthB3))
        bandDefnSeq.append(wv2Band(bandName="Yellow", bandIndex=4, absCalFact=self.absCalFactB4, effBandWidth=self.effBandWidthB4))
        bandDefnSeq.append(wv2Band(bandName="Red", bandIndex=5, absCalFact=self.absCalFactB5, effBandWidth=self.effBandWidthB5))
        bandDefnSeq.append(wv2Band(bandName="RedEdge", bandIndex=6, absCalFact=self.absCalFactB6, effBandWidth=self.effBandWidthB6))
        bandDefnSeq.append(wv2Band(bandName="NIR1", bandIndex=7, absCalFact=self.absCalFactB7, effBandWidth=self.effBandWidthB7))
        bandDefnSeq.append(wv2Band(bandName="NIR2", bandIndex=8, absCalFact=self.absCalFactB8, effBandWidth=self.effBandWidthB8))

        rsgislib.imagecalibration.worldview2ToRadiance(self.fileName, outputImage, outFormat, bandDefnSeq)

        return outputImage, None

    def generateImageSaturationMask(self, outputPath, outputName, outFormat):
        print("Generate Saturation Image")
        outputImage = os.path.join(outputPath, outputName)

        wv2Band = collections.namedtuple('WV2Band', ['bandName', 'fileName', 'bandIndex', 'satVal'])
        bandDefnSeq = list()
        bandDefnSeq.append(wv2Band(bandName="Coastal", fileName=self.fileName, bandIndex=1, satVal=65535.0))
        bandDefnSeq.append(wv2Band(bandName="Blue", fileName=self.fileName, bandIndex=2, satVal=65535.0))
        bandDefnSeq.append(wv2Band(bandName="Green", fileName=self.fileName, bandIndex=3, satVal=65535.0))
        bandDefnSeq.append(wv2Band(bandName="Yellow", fileName=self.fileName, bandIndex=4, satVal=65535.0))
        bandDefnSeq.append(wv2Band(bandName="Red", fileName=self.fileName, bandIndex=5, satVal=65535.0))
        bandDefnSeq.append(wv2Band(bandName="RedEdge", fileName=self.fileName, bandIndex=6, satVal=65535.0))
        bandDefnSeq.append(wv2Band(bandName="NIR1", fileName=self.fileName, bandIndex=7, satVal=65535.0))
        bandDefnSeq.append(wv2Band(bandName="NIR2", fileName=self.fileName, bandIndex=8, satVal=65535.0))

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
        solarIrradianceVals.append(IrrVal(irradiance=1758.2229))
        solarIrradianceVals.append(IrrVal(irradiance=1974.2416))
        solarIrradianceVals.append(IrrVal(irradiance=1856.4104))
        solarIrradianceVals.append(IrrVal(irradiance=1738.4791))
        solarIrradianceVals.append(IrrVal(irradiance=1559.4555))
        solarIrradianceVals.append(IrrVal(irradiance=1342.0695))
        solarIrradianceVals.append(IrrVal(irradiance=1069.7302))
        solarIrradianceVals.append(IrrVal(irradiance=861.2866))
        rsgislib.imagecalibration.radiance2TOARefl(inputRadImage, outputImage, outFormat, rsgislib.TYPE_16UINT, scaleFactor, self.acquisitionTime.year, self.acquisitionTime.month, self.acquisitionTime.day, self.solarZenith, solarIrradianceVals)
        return outputImage

    def generateCloudMask(self, inputReflImage, inputSatImage, inputThermalImage, inputViewAngleImg, inputValidImg, outputPath, outputName, outFormat, tmpPath, scaleFactor, cloud_msk_methods=None):
        raise ARCSIException("Cloud Masking Not Implemented for WorldView2.")

    def createCloudMaskDataArray(self, inImgDataArr):
        return inImgDataArr

    def defineDarkShadowImageBand(self):
        return 8

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

        # Band 1 - Coastal
        s.wavelength = Py6S.Wavelength(0.4, 0.45, [0.4684, 0.5206, 0.5728, 0.6250, 0.6773, 0.7082, 0.7392, 0.7701, 0.8011, 0.8203, 0.8396, 0.8589, 0.8782, 0.8925, 0.9067, 0.9209, 0.9351, 0.9500, 0.9649, 0.9798, 0.9947])
        s.run()
        sixsCoeffs[0,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[0,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[0,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[0,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[0,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[0,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 2 - Blue
        s.wavelength = Py6S.Wavelength(0.45, 0.51, [0.6457, 0.6825, 0.7193, 0.7561, 0.7929, 0.8001, 0.8073, 0.8145, 0.8217, 0.8353, 0.8489, 0.8624, 0.8760, 0.8862, 0.8964, 0.9066, 0.9168, 0.9277, 0.9387, 0.9496, 0.9606, 0.8384, 0.7161, 0.5939, 0.4717])
        s.run()
        sixsCoeffs[1,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[1,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[1,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[1,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[1,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[1,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 3 - Green
        s.wavelength = Py6S.Wavelength(0.510, 0.580, [0.1314, 0.2879, 0.4444, 0.6009, 0.75740, 0.76330, 0.7692, 0.7751, 0.7809, 0.7959, 0.8109, 0.8258, 0.8408, 0.8498, 0.8588, 0.8678, 0.8768, 0.8722, 0.8675, 0.8628, 0.8582, 0.8914, 0.9246, 0.9579, 0.9911, 0.9614, 0.9316, 0.9018, 0.8721])
        s.run()
        sixsCoeffs[2,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[2,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[2,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[2,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[2,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[2,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 4 - Yellow
        s.wavelength = Py6S.Wavelength(0.350000, 1.347500, [0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000001,0.000011,0.000027,0.000001,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000001,0.000005,0.000003,0.000001,0.000001,0.000000,0.000000,0.000000,0.000001,0.000001,0.000001,0.000000,0.000000,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000000,0.000001,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000001,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000003,0.000007,0.000007,0.000004,0.000005,0.000004,0.000008,0.000027,0.000049,0.000127,0.000266,0.000895,0.001886,0.005879,0.015728,0.134652,0.329255,0.663550,0.835025,0.899806,0.917312,0.932634,0.938310,0.947026,0.951042,0.966412,0.979075,0.989723,0.992063,1.000000,0.997307,0.759433,0.501740,0.134322,0.037557,0.007385,0.002855,0.000691,0.000323,0.000158,0.000113,0.000044,0.000025,0.000013,0.000008,0.000005,0.000004,0.000002,0.000002,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000000,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000002,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000000,0.000001,0.000001,0.000001,0.000000,0.000000,0.000000,0.000000,0.000000,0.000001,0.000001,0.000023,0.000014,0.000019,0.000020,0.000011,0.000013,0.000019,0.000020,0.000019,0.000006,0.000018,0.000013,0.000009,0.000014,0.000011,0.000010,0.000012,0.000006,0.000012,0.000013,0.000017,0.000013,0.000007,0.000008,0.000009,0.000007,0.000004,0.000006,0.000005,0.000005,0.000006,0.000005,0.000002,0.000004,0.000004,0.000003,0.000003,0.000003,0.000003,0.000004,0.000004,0.000004,0.000004,0.000004,0.000003,0.000003,0.000002,0.000002,0.000002,0.000002,0.000002,0.000002,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000002,0.000002,0.000002,0.000002,0.000002,0.000001,0.000001,0.000001,0.000001,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001])
        s.run()
        sixsCoeffs[3,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[3,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[3,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[3,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[3,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[3,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 5 - Red
        s.wavelength = Py6S.Wavelength(0.350000, 1.347500, [0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000001,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000002,0.000003,0.000001,0.000001,0.000001,0.000002,0.000001,0.000001,0.000001,0.000001,0.000001,0.000002,0.000003,0.000002,0.000003,0.000002,0.000002,0.000002,0.000001,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000001,0.000001,0.000001,0.000001,0.000001,0.000007,0.000013,0.000025,0.000041,0.000063,0.000123,0.000220,0.000628,0.001380,0.004867,0.015002,0.095566,0.244827,0.596071,0.805405,0.872776,0.881036,0.927928,0.953189,0.965042,0.967019,0.965588,0.958258,0.955500,0.953447,0.955801,0.963907,0.970358,0.974797,0.987526,0.995433,0.999085,0.994134,0.983694,0.986727,0.885280,0.693925,0.291823,0.100695,0.020441,0.007152,0.001914,0.000951,0.000404,0.000232,0.000102,0.000065,0.000027,0.000018,0.000014,0.000018,0.000013,0.000014,0.000013,0.000013,0.000011,0.000008,0.000003,0.000004,0.000002,0.000001,0.000001,0.000000,0.000001,0.000001,0.000003,0.000004,0.000010,0.000009,0.000007,0.000009,0.000007,0.000005,0.000005,0.000002,0.000002,0.000001,0.000002,0.000001,0.000001,0.000001,0.000014,0.000015,0.000013,0.000016,0.000011,0.000013,0.000011,0.000008,0.000009,0.000008,0.000008,0.000004,0.000011,0.000004,0.000007,0.000007,0.000006,0.000006,0.000008,0.000007,0.000003,0.000007,0.000006,0.000007,0.000005,0.000007,0.000009,0.000011,0.000007,0.000009,0.000004,0.000004,0.000009,0.000006,0.000007,0.000007,0.000004,0.000006,0.000005,0.000005,0.000006,0.000005,0.000004,0.000004,0.000004,0.000004,0.000003,0.000004,0.000003,0.000003,0.000003,0.000003,0.000003,0.000003,0.000002,0.000002,0.000002,0.000002,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000000,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000000,0.000001,0.000001,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000])
        s.run()
        sixsCoeffs[4,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[4,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[4,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[4,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[4,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[4,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 6 - Red Edge
        s.wavelength = Py6S.Wavelength(0.350000, 1.347500, [0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000001,0.000002,0.000001,0.000001,0.000002,0.000001,0.000002,0.000001,0.000001,0.000001,0.000001,0.000002,0.000001,0.000001,0.000001,0.000001,0.000010,0.000036,0.000003,0.000000,0.000001,0.000018,0.000003,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000001,0.000001,0.000002,0.000001,0.000008,0.000014,0.000005,0.000004,0.000009,0.000015,0.000008,0.000001,0.000003,0.000004,0.000002,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000002,0.000002,0.000001,0.000001,0.000002,0.000003,0.000003,0.000004,0.000007,0.000006,0.000005,0.000007,0.000011,0.000014,0.000015,0.000024,0.000033,0.000037,0.000039,0.000054,0.000075,0.000088,0.000079,0.000090,0.000154,0.000185,0.000079,0.000009,0.000005,0.000003,0.000002,0.000001,0.000001,0.000001,0.000003,0.000005,0.000006,0.000008,0.000010,0.000012,0.000011,0.000012,0.000011,0.000011,0.000008,0.000008,0.000008,0.000010,0.000012,0.000033,0.000066,0.000207,0.000521,0.003128,0.014134,0.099610,0.289444,0.622908,0.811376,0.952792,0.975357,0.987634,0.988561,0.993901,1.000000,0.999280,0.992650,0.986813,0.980865,0.977116,0.953762,0.846399,0.667534,0.359763,0.178559,0.030544,0.008439,0.002158,0.000904,0.000270,0.000134,0.000055,0.000028,0.000010,0.000007,0.000005,0.000004,0.000004,0.000004,0.000005,0.000007,0.000008,0.000010,0.000006,0.000005,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000001,0.000001,0.000002,0.000003,0.000009,0.000011,0.000011,0.000008,0.000005,0.000004,0.000010,0.000009,0.000008,0.000007,0.000009,0.000011,0.000015,0.000018,0.000017,0.000017,0.000018,0.000018,0.000018,0.000020,0.000016,0.000014,0.000013,0.000011,0.000008,0.000007,0.000005,0.000004,0.000003,0.000003,0.000002,0.000002,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000000,0.000000,0.000001,0.000000,0.000000,0.000001,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000])
        s.run()
        sixsCoeffs[5,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[5,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[5,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[5,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[5,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[5,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 6 - NIR1
        s.wavelength = Py6S.Wavelength(0.350000, 1.347500, [0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000002,0.000004,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000001,0.000001,0.000000,0.000000,0.000001,0.000001,0.000001,0.000001,0.000002,0.000000,0.000000,0.000001,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000001,0.000000,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000002,0.000002,0.000003,0.000005,0.000006,0.000005,0.000004,0.000004,0.000004,0.000005,0.000006,0.000007,0.000010,0.000016,0.000026,0.000044,0.000065,0.000106,0.000149,0.000256,0.000375,0.000713,0.001183,0.002886,0.005210,0.012687,0.022741,0.058860,0.107139,0.260204,0.453897,0.759967,0.912456,0.995805,1.000000,0.988907,0.984506,0.972934,0.962395,0.941139,0.933343,0.919382,0.911682,0.897600,0.886046,0.871191,0.867420,0.844284,0.837202,0.838036,0.842575,0.835657,0.833975,0.842548,0.832789,0.818604,0.815395,0.806711,0.814738,0.794161,0.783077,0.767313,0.750939,0.739271,0.736832,0.734145,0.715699,0.697141,0.684957,0.658635,0.654340,0.632471,0.622861,0.609120,0.600083,0.590914,0.572770,0.507112,0.427640,0.275174,0.179263,0.073805,0.037563,0.013815,0.006654,0.002625,0.001930,0.001306,0.001077,0.000917,0.000914,0.000798,0.000755,0.000727,0.000806,0.000701,0.000666,0.000603,0.000647,0.000281,0.000273,0.000004,0.000003,0.000003,0.000002,0.000002,0.000002,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000000,0.000001,0.000001,0.000001,0.000000,0.000001,0.000000,0.000001,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000])
        s.run()
        sixsCoeffs[6,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[6,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[6,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[6,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[6,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[6,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 8 - NIR2
        s.wavelength = Py6S.Wavelength(0.350000, 1.347500, [0.000000,0.000002,0.000020,0.000064,0.000075,0.000007,0.000000,0.000001,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000001,0.000001,0.000000,0.000000,0.000002,0.000001,0.000001,0.000001,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000002,0.000004,0.000000,0.000000,0.000000,0.000001,0.000001,0.000000,0.000000,0.000000,0.000000,0.000001,0.000001,0.000000,0.000000,0.000000,0.000000,0.000001,0.000001,0.000000,0.000000,0.000000,0.000002,0.000009,0.000001,0.000000,0.000003,0.000005,0.000000,0.000000,0.000002,0.000001,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000001,0.000000,0.000000,0.000000,0.000000,0.000001,0.000001,0.000000,0.000001,0.000000,0.000001,0.000001,0.000001,0.000001,0.000000,0.000000,0.000001,0.000001,0.000001,0.000001,0.000001,0.000000,0.000000,0.000001,0.000000,0.000001,0.000001,0.000001,0.000001,0.000002,0.000002,0.000001,0.000001,0.000002,0.000000,0.000000,0.000000,0.000000,0.000000,0.000001,0.000001,0.000000,0.000000,0.000000,0.000001,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000001,0.000001,0.000001,0.000001,0.000000,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000000,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000000,0.000453,0.000487,0.000780,0.001073,0.000976,0.001071,0.001069,0.001067,0.001096,0.001097,0.001094,0.001103,0.001042,0.001105,0.001127,0.001151,0.001151,0.001162,0.001355,0.001304,0.001492,0.002052,0.004328,0.008413,0.030596,0.086262,0.313746,0.528020,0.789269,0.926443,0.987144,0.992903,0.993363,0.999381,0.996691,0.985229,0.970856,0.956271,0.933020,0.924312,0.910922,0.901617,0.892112,0.881108,0.864645,0.853138,0.841298,0.831520,0.808180,0.796630,0.774377,0.762548,0.746906,0.730666,0.708967,0.694925,0.661720,0.644251,0.620929,0.600818,0.574685,0.556569,0.530905,0.515101,0.493667,0.478132,0.454311,0.441024,0.421449,0.409005,0.389741,0.379119,0.363517,0.352547,0.336697,0.328829,0.317080,0.308137,0.294920,0.284427,0.271439,0.262990,0.250912,0.245015,0.235734,0.229439,0.217559,0.210720,0.198576,0.190246,0.178705,0.170271,0.159032,0.150452,0.137770,0.127452,0.109612,0.095607,0.071507,0.054200,0.032806,0.022246,0.009796,0.005329,0.002185,0.001196,0.000524,0.000324,0.000200,0.000157,0.000114,0.000105,0.000080,0.000073,0.000080,0.000063,0.000070,0.000068,0.000054,0.000029,0.000030,0.000028,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001])
        s.run()
        sixsCoeffs[7,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[7,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[7,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[7,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[7,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[7,5] = float(s.outputs.values['environmental_irradiance'])

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
        imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=5, aX=float(sixsCoeffs[4,0]), bX=float(sixsCoeffs[4,1]), cX=float(sixsCoeffs[4,2]), DirIrr=float(sixsCoeffs[4,3]), DifIrr=float(sixsCoeffs[4,4]), EnvIrr=float(sixsCoeffs[4,5])))
        imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=6, aX=float(sixsCoeffs[5,0]), bX=float(sixsCoeffs[5,1]), cX=float(sixsCoeffs[5,2]), DirIrr=float(sixsCoeffs[5,3]), DifIrr=float(sixsCoeffs[5,4]), EnvIrr=float(sixsCoeffs[5,5])))
        imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=7, aX=float(sixsCoeffs[6,0]), bX=float(sixsCoeffs[6,1]), cX=float(sixsCoeffs[6,2]), DirIrr=float(sixsCoeffs[6,3]), DifIrr=float(sixsCoeffs[6,4]), EnvIrr=float(sixsCoeffs[6,5])))
        imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=8, aX=float(sixsCoeffs[7,0]), bX=float(sixsCoeffs[7,1]), cX=float(sixsCoeffs[7,2]), DirIrr=float(sixsCoeffs[7,3]), DifIrr=float(sixsCoeffs[7,4]), EnvIrr=float(sixsCoeffs[7,5])))

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
            imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=5, aX=float(sixsCoeffs[4,0]), bX=float(sixsCoeffs[4,1]), cX=float(sixsCoeffs[4,2]), DirIrr=float(sixsCoeffs[4,3]), DifIrr=float(sixsCoeffs[4,4]), EnvIrr=float(sixsCoeffs[4,5])))
            imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=6, aX=float(sixsCoeffs[5,0]), bX=float(sixsCoeffs[5,1]), cX=float(sixsCoeffs[5,2]), DirIrr=float(sixsCoeffs[5,3]), DifIrr=float(sixsCoeffs[5,4]), EnvIrr=float(sixsCoeffs[5,5])))
            imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=7, aX=float(sixsCoeffs[6,0]), bX=float(sixsCoeffs[6,1]), cX=float(sixsCoeffs[6,2]), DirIrr=float(sixsCoeffs[6,3]), DifIrr=float(sixsCoeffs[6,4]), EnvIrr=float(sixsCoeffs[6,5])))
            imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=8, aX=float(sixsCoeffs[7,0]), bX=float(sixsCoeffs[7,1]), cX=float(sixsCoeffs[7,2]), DirIrr=float(sixsCoeffs[7,3]), DifIrr=float(sixsCoeffs[7,4]), EnvIrr=float(sixsCoeffs[7,5])))

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
                    imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=5, aX=float(sixsCoeffs[4,0]), bX=float(sixsCoeffs[4,1]), cX=float(sixsCoeffs[4,2]), DirIrr=float(sixsCoeffs[4,3]), DifIrr=float(sixsCoeffs[4,4]), EnvIrr=float(sixsCoeffs[4,5])))
                    imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=6, aX=float(sixsCoeffs[5,0]), bX=float(sixsCoeffs[5,1]), cX=float(sixsCoeffs[5,2]), DirIrr=float(sixsCoeffs[5,3]), DifIrr=float(sixsCoeffs[5,4]), EnvIrr=float(sixsCoeffs[5,5])))
                    imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=7, aX=float(sixsCoeffs[6,0]), bX=float(sixsCoeffs[6,1]), cX=float(sixsCoeffs[6,2]), DirIrr=float(sixsCoeffs[6,3]), DifIrr=float(sixsCoeffs[6,4]), EnvIrr=float(sixsCoeffs[6,5])))
                    imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=8, aX=float(sixsCoeffs[7,0]), bX=float(sixsCoeffs[7,1]), cX=float(sixsCoeffs[7,2]), DirIrr=float(sixsCoeffs[7,3]), DifIrr=float(sixsCoeffs[7,4]), EnvIrr=float(sixsCoeffs[7,5])))
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
        #s.ground_reflectance = Py6S.GroundReflectance.HomogeneousHapke(0.101, -0.263, 0.589, 0.046)
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

            dosBlueImage = ""
            minObjSize = 3
            darkPxlPercentile = 0.01
            blockSize = 1000
            if simpleDOS:
                outputDOSBlueName = tmpBaseName + "DOSBlue" + imgExtension
                dosBlueImage, bandOff = self.convertImageBandToReflectanceSimpleDarkSubtract(inputTOAImage, outputPath, outputDOSBlueName, outFormat, dosOutRefl, 2)
            elif globalDOS:
                dosBlueImage = self.performDOSOnSingleBand(inputTOAImage, 2, outputPath, tmpBaseName, "Blue", "KEA", tmpPath, minObjSize, darkPxlPercentile, dosOutRefl)
            else:
                dosBlueImage = self.performLocalDOSOnSingleBand(inputTOAImage, 2, outputPath, tmpBaseName, "Blue", "KEA", tmpPath, minObjSize, darkPxlPercentile, blockSize, dosOutRefl)

            thresImageClumpsFinal = os.path.join(tmpPath, tmpBaseName + "_clumps" + imgExtension)
            rsgislib.segmentation.segutils.runShepherdSegmentation(inputTOAImage, thresImageClumpsFinal, tmpath=tmpPath, gdalformat="KEA", numClusters=40, minPxls=10, bands=[8,6,1], processInMem=True)

            stats2CalcTOA = list()
            stats2CalcTOA.append(rsgislib.rastergis.BandAttStats(band=1, meanField="MeanElev"))
            rsgislib.rastergis.populateRATWithStats(inputDEMFile, thresImageClumpsFinal, stats2CalcTOA)

            stats2CalcTOA = list()
            stats2CalcTOA.append(rsgislib.rastergis.BandAttStats(band=2, meanField="MeanBlueDOS"))
            rsgislib.rastergis.populateRATWithStats(dosBlueImage, thresImageClumpsFinal, stats2CalcTOA)

            stats2CalcRad = list()
            stats2CalcRad.append(rsgislib.rastergis.BandAttStats(band=2, meanField="MeanBlueRAD"))
            stats2CalcRad.append(rsgislib.rastergis.BandAttStats(band=7, meanField="MeanNIRRAD"))
            stats2CalcRad.append(rsgislib.rastergis.BandAttStats(band=5, meanField="MeanRedRAD"))
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
            return self.estimateSingleAOTFromDOSBandImpl(radianceImage, toaImage, inputDEMFile, tmpPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, minAOT, maxAOT, dosOutRefl, 2)
        except Exception as e:
            raise

    def setBandNames(self, imageFile):
        dataset = gdal.Open(imageFile, gdal.GA_Update)
        if not dataset is None:
            dataset.GetRasterBand(1).SetDescription("Coastal")
            dataset.GetRasterBand(2).SetDescription("Blue")
            dataset.GetRasterBand(3).SetDescription("Green")
            dataset.GetRasterBand(4).SetDescription("Yellow")
            dataset.GetRasterBand(5).SetDescription("Red")
            dataset.GetRasterBand(6).SetDescription("RedEdge")
            dataset.GetRasterBand(7).SetDescription("NIR1")
            dataset.GetRasterBand(8).SetDescription("NIR2")
            dataset = None
        else:
            print("Could not open image to set band names: ", imageFile)

    def cleanLocalFollowProcessing(self):
        print("")
