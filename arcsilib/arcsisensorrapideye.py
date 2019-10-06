"""
Module that contains the ARCSIRapidEyeSensor class.
"""
############################################################################
#  arcsisensorrapideye.py
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
# Purpose:  A class for read the RapidEye sensor header file and applying
#           the pre-processing operations within ARCSI to the RapidEye
#           datasets.
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 23/08/2013
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
# Import the RSGISLib Module.
import rsgislib
# Import the RSGISLib Image Calibration Module.
import rsgislib.imagecalibration
# Import the RSGISLib Image Calculation Module
import rsgislib.imagecalc
# Import the RSGISLib Image Utilities Module
import rsgislib.imageutils
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
# Import JSON module
import json
# Import glob module
import glob

class ARCSIRapidEyeSensor (ARCSIAbstractSensor):
    """
    A class which represents the RapidEye sensor to read
    header parameters and apply data processing operations.
    """
    def __init__(self, debugMode, inputImage):
        ARCSIAbstractSensor.__init__(self, debugMode, inputImage)
        self.sensor = "RapidEye"
        self.platOrbitType = ""
        self.platSerialId = ""
        self.platShortHand = ""
        self.instShortHand = ""
        self.senrType = ""
        self.senrRes = 0.0
        self.senrScanType = ""
        self.pixelFormat = ""
        self.tileID = ""
        self.orderID = ""
        self.catalogID = ""
        self.acquIncidAngle = 0.0
        self.acquAzimuthAngle = 0.0
        self.acquCraftViewAngle = 0.0
        self.numOfBands = 0
        self.radioCorrApplied = False
        self.atmosCorrApplied = False
        self.elevCorrApplied = False
        self.geoCorrLevel = ""
        self.radioCorrVersion = ""
        self.fileName = ""
        self.origFileName = ""

    def extractHeaderParameters(self, inputHeader, wktStr):
        """
        Understands and parses the RapidEye metadata.xml header file
        """
        try:
            print("Reading header file")

            self.headerFileName = os.path.split(inputHeader)[1]
            
            arcsiUtils = ARCSIUtils()

            hdrExt = os.path.splitext(inputHeader)
            if not len(hdrExt) is 2:
                raise ARCSIException("Cannot work out what the file extention is - support either xml or json.")
            hdrExt = hdrExt[1]

            if (hdrExt.lower() == '.xml') or (hdrExt.lower() == 'xml'):
                tree = ET.parse(inputHeader)
                root = tree.getroot()

                hdrVersion = root.attrib['version'].strip() # 1.2.1 when this was implemented but the version was not changed when moved to plantlabs schema URL.
                schemaURL = root.attrib['{http://www.w3.org/2001/XMLSchema-instance}schemaLocation'].strip().split()[0]
                rapideyeUrl = '{'+schemaURL+'}'
                metaDataProperty = root.find('{http://www.opengis.net/gml}metaDataProperty')
                eoMetaData = metaDataProperty.find(rapideyeUrl+'EarthObservationMetaData')
                if eoMetaData is None:
                    rapideyeUrl = '{http://schemas.rapideye.de/products/productMetadataSensor}'
                    eoMetaData = metaDataProperty.find(rapideyeUrl+'EarthObservationMetaData')
                productType = eoMetaData.find('{http://earth.esa.int/eop}productType').text.strip()

                if (productType == "L1B") and (self.userSpInputImage is None):
                    raise ARCSIException("L1B data is supported by ARCSI only when a user defined image is provided.")
                elif (productType != "L3A") & (productType != "L1B"):
                    raise ARCSIException("Only L3A and L1B data are supported by ARCSI.")

                eoPlatform = root.find('{http://www.opengis.net/gml}using').find('{http://earth.esa.int/eop}EarthObservationEquipment').find('{http://earth.esa.int/eop}platform').find('{http://earth.esa.int/eop}Platform')
                self.platShortHand = eoPlatform.find('{http://earth.esa.int/eop}shortName').text.strip()
                self.platSerialId = eoPlatform.find('{http://earth.esa.int/eop}serialIdentifier').text.strip()
                self.platOrbitType = eoPlatform.find('{http://earth.esa.int/eop}orbitType').text.strip()

                eoInstrument = root.find('{http://www.opengis.net/gml}using').find('{http://earth.esa.int/eop}EarthObservationEquipment').find('{http://earth.esa.int/eop}instrument').find('{http://earth.esa.int/eop}Instrument')
                self.instShortHand = eoInstrument.find('{http://earth.esa.int/eop}shortName').text.strip()

                eoSensor = root.find('{http://www.opengis.net/gml}using').find('{http://earth.esa.int/eop}EarthObservationEquipment').find('{http://earth.esa.int/eop}sensor').find(rapideyeUrl+'Sensor')
                self.senrType = eoSensor.find('{http://earth.esa.int/eop}sensorType').text.strip()
                self.senrRes = float(eoSensor.find('{http://earth.esa.int/eop}resolution').text.strip())
                self.senrScanType = eoSensor.find(rapideyeUrl+'scanType').text.strip()

                eoAcquParams = root.find('{http://www.opengis.net/gml}using').find('{http://earth.esa.int/eop}EarthObservationEquipment').find('{http://earth.esa.int/eop}acquisitionParameters').find(rapideyeUrl+'Acquisition')

                self.acquIncidAngle = float(eoAcquParams.find('{http://earth.esa.int/eop}incidenceAngle').text.strip())
                self.acquAzimuthAngle = float(eoAcquParams.find(rapideyeUrl+'azimuthAngle').text.strip())
                self.acquCraftViewAngle = float(eoAcquParams.find(rapideyeUrl+'spaceCraftViewAngle').text.strip())

                self.solarZenith = 90-float(eoAcquParams.find('{http://earth.esa.int/opt}illuminationElevationAngle').text.strip())
                self.solarAzimuth = float(eoAcquParams.find('{http://earth.esa.int/opt}illuminationAzimuthAngle').text.strip())
                self.sensorZenith = self.acquCraftViewAngle
                self.sensorAzimuth = self.acquAzimuthAngle
                timeStr = eoAcquParams.find(rapideyeUrl+'acquisitionDateTime').text.strip()
                timeStr = timeStr.replace('Z', '')
                try:
                    self.acquisitionTime = datetime.datetime.strptime(timeStr, "%Y-%m-%dT%H:%M:%S.%f")
                except Exception as e:
                    try:
                        self.acquisitionTime = datetime.datetime.strptime(timeStr, "%Y-%m-%dT%H:%M:%S")
                    except Exception as e:
                        try:
                            timeStrTmp = timeStr.split('+')[0]
                            self.acquisitionTime = datetime.datetime.strptime(timeStrTmp, "%Y-%m-%dT%H:%M:%S.%f")
                        except Exception as e:
                            raise e

                metadata = root.find('{http://www.opengis.net/gml}metaDataProperty').find(rapideyeUrl+'EarthObservationMetaData')
                if not  metadata.find(rapideyeUrl+'tileId') is None:
                    self.tileID = metadata.find(rapideyeUrl+'tileId').text.strip()
                else:
                    self.tileID = ""
                self.orderID = metadata.find(rapideyeUrl+'orderId').text.strip()
                self.pixelFormat = metadata.find(rapideyeUrl+'pixelFormat').text.strip()

                centrePt = root.find('{http://www.opengis.net/gml}target').find(rapideyeUrl+'Footprint').find('{http://www.opengis.net/gml}centerOf').find('{http://www.opengis.net/gml}Point').find('{http://www.opengis.net/gml}pos').text.strip()
                centrePtSplit = centrePt.split(' ')
                self.latCentre = float(centrePtSplit[0])
                self.lonCentre = float(centrePtSplit[1])

                imgBounds = root.find('{http://www.opengis.net/gml}target').find(rapideyeUrl+'Footprint').find(rapideyeUrl+'geographicLocation')
                tlPoint = imgBounds.find(rapideyeUrl+'topLeft')
                self.latTL = float(tlPoint.find(rapideyeUrl+'latitude').text)
                self.lonTL = float(tlPoint.find(rapideyeUrl+'longitude').text)
                trPoint = imgBounds.find(rapideyeUrl+'topRight')
                self.latTR = float(trPoint.find(rapideyeUrl+'latitude').text)
                self.lonTR = float(trPoint.find(rapideyeUrl+'longitude').text)
                brPoint = imgBounds.find(rapideyeUrl+'bottomRight')
                self.latBR = float(brPoint.find(rapideyeUrl+'latitude').text)
                self.lonBR = float(brPoint.find(rapideyeUrl+'longitude').text)
                blPoint = imgBounds.find(rapideyeUrl+'bottomLeft')
                self.latBL = float(blPoint.find(rapideyeUrl+'latitude').text)
                self.lonBL = float(blPoint.find(rapideyeUrl+'longitude').text)

                productInfo = root.find('{http://www.opengis.net/gml}resultOf').find(rapideyeUrl+'EarthObservationResult').find('{http://earth.esa.int/eop}product').find(rapideyeUrl+'ProductInformation')

                spatialRef = productInfo.find(rapideyeUrl+'spatialReferenceSystem')

                epsgCode = int(spatialRef.find(rapideyeUrl+'epsgCode').text)
                inProj = osr.SpatialReference()
                inProj.ImportFromEPSG(epsgCode)
                if self.inWKT is "":
                    self.inWKT = inProj.ExportToWkt()

                self.numOfBands = int(productInfo.find(rapideyeUrl+'numBands').text.strip())
                if self.numOfBands != 5:
                    raise ARCSIException("The number of image band is not equal to 5 according to XML header.")

                radioCorrAppliedStr = productInfo.find(rapideyeUrl+'radiometricCorrectionApplied').text.strip()
                if radioCorrAppliedStr == "true":
                    self.radioCorrApplied = True
                else:
                    self.radioCorrApplied = False

                if self.radioCorrApplied:
                    try:
                        self.radioCorrVersion = productInfo.find(rapideyeUrl+'radiometricCalibrationVersion').text.strip()
                    except Exception:
                        self.radioCorrVersion = 'Unknown'
                else:
                    self.radioCorrVersion = 'Not Applied'

                atmosCorrAppliedStr = productInfo.find(rapideyeUrl+'atmosphericCorrectionApplied').text.strip()
                if atmosCorrAppliedStr == "true":
                    self.atmosCorrApplied = True
                else:
                    self.atmosCorrApplied = False

                if self.atmosCorrApplied:
                    raise ARCSIException("An atmosheric correction has already been applied according to the metadata.")

                elevCorrAppliedStr = productInfo.find(rapideyeUrl+'elevationCorrectionApplied').text.strip()
                if elevCorrAppliedStr == "true":
                    self.elevCorrApplied = True
                else:
                    self.elevCorrApplied = False

                self.geoCorrLevel = productInfo.find(rapideyeUrl+'geoCorrectionLevel').text.strip()

                filesDIR = os.path.dirname(inputHeader)
                if not self.userSpInputImage is None:
                    self.fileName = os.path.abspath(self.userSpInputImage)
                else:
                    self.fileName = os.path.join(filesDIR, productInfo.find('{http://earth.esa.int/eop}fileName').text.strip())
                print('self.fileName = ', self.fileName)

                rsgisUtils = rsgislib.RSGISPyUtils()
                minX, maxX, minY, maxY = rsgisUtils.getImageBBOX(self.fileName)

                self.xTL = minX
                self.yTL = maxY
                self.xTR = maxX
                self.yTR = maxY
                self.xBL = minX
                self.yBL = minY
                self.xBR = maxX
                self.yBR = minY
                self.xCentre = minX + ((maxX - minX) / 2)
                self.yCentre = minY + ((maxY - minY) / 2)

            elif (hdrExt.lower() == '.json') or (hdrExt.lower() == 'json'):
                with open(inputHeader, 'r') as f:
                    jsonStrData = f.read()
                reHdrInfo = json.loads(jsonStrData)

                if 'properties' in reHdrInfo:
                    if 'provider' in reHdrInfo['properties']:
                        if reHdrInfo['properties']['provider'].lower() != 'rapideye':
                            raise ARCSIException("JSON Header is not expect format for RapidEye.")
                    else:
                        raise ARCSIException("JSON Header is not expect format for RapidEye.")
                else:
                    raise ARCSIException("JSON Header is not expect format for RapidEye.")


                if 'sat' in reHdrInfo['properties']:
                    self.acquIncidAngle = arcsiUtils.str2Float(reHdrInfo['properties']['sat']['off_nadir'])
                    self.acquAzimuthAngle = arcsiUtils.str2Float(reHdrInfo['properties']['sat']['azimuth_angle'])
                    self.acquCraftViewAngle = arcsiUtils.str2Float(reHdrInfo['properties']['sat']['view_angle'])

                    self.solarZenith = 90-arcsiUtils.str2Float(reHdrInfo['properties']['sun']['altitude'])
                    self.solarAzimuth = arcsiUtils.str2Float(reHdrInfo['properties']['sun']['azimuth'])
                    self.sensorZenith = self.acquCraftViewAngle
                    self.sensorAzimuth = self.acquAzimuthAngle
                else:
                    raise ARCSIException("JSON Header is not expect format for RapidEye.")

                if 'acquired' in reHdrInfo['properties']:
                    timeStr = reHdrInfo['properties']['acquired']
                    timeStr = timeStr.replace('Z', '')
                    try:
                        self.acquisitionTime = datetime.datetime.strptime(timeStr, "%Y-%m-%dT%H:%M:%S.%f")
                    except Exception as e:
                        try:
                            self.acquisitionTime = datetime.datetime.strptime(timeStr, "%Y-%m-%dT%H:%M:%S")
                        except Exception as e:
                            raise e
                else:
                    raise ARCSIException("JSON Header is not expect format for RapidEye.")

                if 'rapideye' in reHdrInfo['properties']:
                    self.tileID = reHdrInfo['properties']['rapideye']['tile_id']
                    self.catalogID = reHdrInfo['properties']['rapideye']['catalog_id']
                else:
                    raise ARCSIException("JSON Header is not expect format for RapidEye.")

                if 'geometry' in reHdrInfo:
                    if 'coordinates' in reHdrInfo['geometry']:
                        self.latTL = arcsiUtils.str2Float(reHdrInfo['geometry']['coordinates'][0][0][1])
                        self.lonTL = arcsiUtils.str2Float(reHdrInfo['geometry']['coordinates'][0][0][0])
                        self.latTR = arcsiUtils.str2Float(reHdrInfo['geometry']['coordinates'][0][1][1])
                        self.lonTR = arcsiUtils.str2Float(reHdrInfo['geometry']['coordinates'][0][1][0])
                        self.latBR = arcsiUtils.str2Float(reHdrInfo['geometry']['coordinates'][0][2][1])
                        self.lonBR = arcsiUtils.str2Float(reHdrInfo['geometry']['coordinates'][0][2][0])
                        self.latBL = arcsiUtils.str2Float(reHdrInfo['geometry']['coordinates'][0][3][1])
                        self.lonBL = arcsiUtils.str2Float(reHdrInfo['geometry']['coordinates'][0][3][0])

                        self.latCentre = self.latTL + (self.latBR - self.latTL)/2
                        self.lonCentre = self.lonTL + (self.lonTR - self.lonTL)/2
                    else:
                        raise ARCSIException("JSON Header is not expect format for RapidEye.")
                else:
                    raise ARCSIException("JSON Header is not expect format for RapidEye.")

                if not self.userSpInputImage is None:
                    self.fileName = os.path.abspath(self.userSpInputImage)
                else:
                    baseHdrName = os.path.splitext(inputHeader)[0]
                    fileBaseName = baseHdrName.replace('_metadata', '')
                    imgFiles = glob.glob(fileBaseName+'*analytic.tif')
                    if len(imgFiles) == 0:
                        raise ARCSIException("Could not find input image file.")
                    if len(imgFiles) > 1:
                        raise ARCSIException("Found multiple potential input image files - don't know which one is correct specify input image using arcsi.py.")
                    self.fileName = imgFiles[0]
                print('self.fileName = ', self.fileName)

                self.radioCorrApplied = True # JSON doesn't specify this!! Assume true :s

                rsgisUtils = rsgislib.RSGISPyUtils()
                minX, maxX, minY, maxY = rsgisUtils.getImageBBOX(self.fileName)

                self.xTL = minX
                self.yTL = maxY
                self.xTR = maxX
                self.yTR = maxY
                self.xBL = minX
                self.yBL = minY
                self.xBR = maxX
                self.yBR = minY
                self.xCentre = minX + ((maxX - minX) / 2)
                self.yCentre = minY + ((maxY - minY) / 2)
            else:
                raise ARCSIException("Header file extention is not recognised - support either xml or json.")

        except Exception as e:
            raise e

    def checkInputImageValid(self):
        if not self.expectedImageDataPresent():
            raise ARCSIException("Error image image was not present.")
        rasterDS = gdal.Open(self.fileName, gdal.GA_ReadOnly)
        if rasterDS == None:
            raise ARCSIException('Could not open raster image: ' + self.fileName)
        nBands = rasterDS.RasterCount
        if nBands == 6:
            # Subset bands...
            rsgisUtils = rsgislib.RSGISPyUtils()
            self.origFileName = self.fileName
            self.fileName = self.generateOutputBaseName()+'BandSubDNImg.kea'
            rsgislib.imageutils.selectImageBands(self.origFileName, self.fileName, 'KEA', rsgisUtils.getRSGISLibDataTypeFromImg(self.origFileName), [1,2,3,4,5])
        elif nBands != 5:
            raise ARCSIException('Input Image \'' + self.fileName + '\' does not have the expect 5 image bands.')

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
        Customises the generic name for the RapidEye sensor
        """
        reTileID = ""
        if self.tileID != "":
            reTileID = "_tid" + str(self.tileID)
        reOrderID = ""
        if self.orderID != "":
            reOrderID = "_oid" + str(self.orderID)
        outname = self.defaultGenBaseOutFileName()
        outname = outname + reTileID + reOrderID
        return outname

    def expectedImageDataPresent(self):
        imageDataPresent = True

        if not os.path.exists(self.fileName):
            imageDataPresent = False

        return imageDataPresent

    def applyImageDataMask(self, inputHeader, outputPath, outputMaskName, outputImgName, outFormat, outWKTFile):
        raise ARCSIException("RapidEye does not provide any image masks, do not use the MASK option.")

    def mosaicImageTiles(self, outputPath):
        raise ARCSIException("Image data does not need mosaicking")

    def resampleImgRes(self, outputPath, resampleToLowResImg, resampleMethod='cubic', multicore=False):
        raise ARCSIException("Image data does not need resampling")

    def sharpenLowResRadImgBands(self, inputImg, outputImage, outFormat):
        raise ARCSIException("Image sharpening is not available for this sensor.")

    def convertImageToRadiance(self, outputPath, outputReflName, outputThermalName, outFormat):
        print("Converting to Radiance")
        outputImage = os.path.join(outputPath, outputReflName)
        if self.radioCorrApplied:
            # Rescale the data to be between 0 and 1.
            rsgislib.imagecalc.imageMath(self.fileName, outputImage, "b1/100", outFormat, rsgislib.TYPE_32FLOAT)
        else:
            raise ARCSIException("Radiometric correction has not been applied - this is not implemented within ARCSI yet. Check your data version.")

        return outputImage, None

    def generateImageSaturationMask(self, outputPath, outputName, outFormat):
        print("Generate Saturation Image")
        outputImage = os.path.join(outputPath, outputName)

        reBand = collections.namedtuple('REBand', ['bandName', 'fileName', 'bandIndex', 'satVal'])
        bandDefnSeq = list()
        bandDefnSeq.append(reBand(bandName="Blue", fileName=self.fileName, bandIndex=1, satVal=65535.0))
        bandDefnSeq.append(reBand(bandName="Green", fileName=self.fileName, bandIndex=2, satVal=65535.0))
        bandDefnSeq.append(reBand(bandName="Red", fileName=self.fileName, bandIndex=3, satVal=65535.0))
        bandDefnSeq.append(reBand(bandName="RedEdge", fileName=self.fileName, bandIndex=4, satVal=65535.0))
        bandDefnSeq.append(reBand(bandName="NIR", fileName=self.fileName, bandIndex=5, satVal=65535.0))

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
        solarIrradianceVals.append(IrrVal(irradiance=1997.8))
        solarIrradianceVals.append(IrrVal(irradiance=1863.5))
        solarIrradianceVals.append(IrrVal(irradiance=1560.4))
        solarIrradianceVals.append(IrrVal(irradiance=1395.0))
        solarIrradianceVals.append(IrrVal(irradiance=1124.4))
        rsgislib.imagecalibration.radiance2TOARefl(inputRadImage, outputImage, outFormat, rsgislib.TYPE_16UINT, scaleFactor, self.acquisitionTime.year, self.acquisitionTime.month, self.acquisitionTime.day, self.solarZenith, solarIrradianceVals)
        return outputImage

    def generateCloudMask(self, inputReflImage, inputSatImage, inputThermalImage, inputViewAngleImg, inputValidImg, outputPath, outputName, outFormat, tmpPath, scaleFactor, cloud_msk_methods=None):
        raise ARCSIException("Cloud Masking Not Implemented for Rapideye.")

    def createCloudMaskDataArray(self, inImgDataArr):
        # Calc Whiteness
        meanArr = numpy.mean(inImgDataArr, axis=1)
        whitenessArr = numpy.absolute((inImgDataArr[...,0] - meanArr)/meanArr) + numpy.absolute((inImgDataArr[...,1] - meanArr)/meanArr) + numpy.absolute((inImgDataArr[...,2] - meanArr)/meanArr) + numpy.absolute((inImgDataArr[...,3] - meanArr)/meanArr) + numpy.absolute((inImgDataArr[...,4] - meanArr)/meanArr)
        # Calc NDVI
        ndvi = (inImgDataArr[...,4] - inImgDataArr[...,2]) / (inImgDataArr[...,4] + inImgDataArr[...,2])
        
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
        return 5

    def calc6SCoefficients(self, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF):
        sixsCoeffs = numpy.zeros((5, 6), dtype=numpy.float32)
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
        if useBRDF:
            s.atmos_corr = Py6S.AtmosCorr.AtmosCorrBRDFFromRadiance(200)
        else:
            s.atmos_corr = Py6S.AtmosCorr.AtmosCorrLambertianFromRadiance(200)
        s.aot550 = aotVal

        # Band 1
        s.wavelength = Py6S.Wavelength(0.435, 0.515, [0.001, 0.004, 0.321, 0.719, 0.74, 0.756, 0.77, 0.78, 0.784, 0.792, 0.796, 0.799, 0.806, 0.804, 0.807, 0.816, 0.82, 0.825, 0.84, 0.845, 0.862, 0.875, 0.886, 0.905, 0.928, 0.936, 0.969, 0.967, 1, 0.976, 0.437, 0.029, 0.001])
        s.run()
        sixsCoeffs[0,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[0,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[0,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[0,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[0,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[0,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 2
        s.wavelength = Py6S.Wavelength(0.510, 0.5975, [0.001, 0.002, 0.013, 0.054, 0.539, 0.868, 0.868, 0.877, 0.871, 0.874, 0.882, 0.882, 0.881, 0.886, 0.897, 0.899, 0.901, 0.91, 0.924, 0.928, 0.936, 0.946, 0.953, 0.96, 0.974, 0.976, 0.976, 0.989, 0.988, 0.984, 0.994, 0.97, 0.417, 0.039, 0.002, 0.001])
        s.run()
        sixsCoeffs[1,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[1,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[1,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[1,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[1,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[1,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 3
        s.wavelength = Py6S.Wavelength(0.620, 0.6925, [0.001, 0.002, 0.009, 0.038, 0.437, 0.856, 0.854, 0.876, 0.881, 0.885, 0.902, 0.909, 0.915, 0.923, 0.939, 0.947, 0.958, 0.963, 0.97, 0.976, 0.989, 0.991, 0.985, 0.994, 0.989, 0.989, 0.463, 0.062, 0.005, 0.001])
        s.run()
        sixsCoeffs[2,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[2,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[2,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[2,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[2,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[2,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 4
        s.wavelength = Py6S.Wavelength(0.6775, 0.7425, [0.001, 0.002, 0.004, 0.021, 0.074, 0.491, 0.914, 0.998, 0.999, 0.998, 0.993, 0.987, 0.986, 0.982, 0.976, 0.966, 0.964, 0.961, 0.949, 0.939, 0.936, 0.425, 0.123, 0.02, 0.007, 0.002, 0.001])
        s.run()
        sixsCoeffs[3,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[3,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[3,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[3,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[3,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[3,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 5
        s.wavelength = Py6S.Wavelength(0.740, 0.870, [0.001, 0.001, 0.003, 0.005, 0.012, 0.023, 0.068, 0.153, 0.497, 0.828, 1, 0.982, 0.967, 0.974, 0.983, 0.981, 0.97, 0.963, 0.958, 0.957, 0.958, 0.959, 0.956, 0.954, 0.948, 0.944, 0.937, 0.933, 0.928, 0.927, 0.926, 0.926, 0.923, 0.918, 0.906, 0.898, 0.889, 0.885, 0.882, 0.876, 0.857, 0.842, 0.84, 0.832, 0.582, 0.295, 0.08, 0.034, 0.011, 0.006, 0.002, 0.001, 0.001])
        s.run()
        sixsCoeffs[4,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[4,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[4,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[4,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[4,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[4,5] = float(s.outputs.values['environmental_irradiance'])

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
                imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=5, aX=float(sixsCoeffs[4,0]), bX=float(sixsCoeffs[4,1]), cX=float(sixsCoeffs[4,2]), DirIrr=float(sixsCoeffs[4,3]), DifIrr=float(sixsCoeffs[4,4]), EnvIrr=float(sixsCoeffs[4,5])))

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

        # Band 1 (Blue!)
        s.wavelength = Py6S.Wavelength(0.435, 0.515, [0.001, 0.004, 0.321, 0.719, 0.74, 0.756, 0.77, 0.78, 0.784, 0.792, 0.796, 0.799, 0.806, 0.804, 0.807, 0.816, 0.82, 0.825, 0.84, 0.845, 0.862, 0.875, 0.886, 0.905, 0.928, 0.936, 0.969, 0.967, 1, 0.976, 0.437, 0.029, 0.001])
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
        raise ARCSIException("findDDVTargets is not implemented.")

    def estimateImageToAODUsingDDV(self, inputRADImage, inputTOAImage, inputDEMFile, shadowMask, outputPath, outputName, outFormat, tmpPath, aeroProfile, atmosProfile, grdRefl, aotValMin, aotValMax):
        raise ARCSIException("estimateImageToAODUsingDDV is not implemented.")

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
            rsgislib.segmentation.segutils.runShepherdSegmentation(inputTOAImage, thresImageClumpsFinal, tmpath=tmpPath, gdalformat="KEA", numClusters=40, minPxls=10, bands=[5,4,1], processInMem=True)

            stats2CalcTOA = list()
            stats2CalcTOA.append(rsgislib.rastergis.BandAttStats(band=1, meanField="MeanElev"))
            rsgislib.rastergis.populateRATWithStats(inputDEMFile, thresImageClumpsFinal, stats2CalcTOA)

            stats2CalcTOA = list()
            stats2CalcTOA.append(rsgislib.rastergis.BandAttStats(band=1, meanField="MeanB1DOS"))
            rsgislib.rastergis.populateRATWithStats(dosBlueImage, thresImageClumpsFinal, stats2CalcTOA)

            stats2CalcRad = list()
            stats2CalcRad.append(rsgislib.rastergis.BandAttStats(band=1, meanField="MeanB1RAD"))
            stats2CalcRad.append(rsgislib.rastergis.BandAttStats(band=5, meanField="MeanB5RAD"))
            stats2CalcRad.append(rsgislib.rastergis.BandAttStats(band=3, meanField="MeanB3RAD"))
            rsgislib.rastergis.populateRATWithStats(inputRADImage, thresImageClumpsFinal, stats2CalcRad)

            ratDS = gdal.Open(thresImageClumpsFinal, gdal.GA_Update)
            Histogram = rat.readColumn(ratDS, "Histogram")
            MeanElev = rat.readColumn(ratDS, "MeanElev")

            MeanB5RAD = rat.readColumn(ratDS, "MeanB5RAD")
            MeanB3RAD = rat.readColumn(ratDS, "MeanB3RAD")

            radNDVI = (MeanB5RAD - MeanB3RAD)/(MeanB5RAD + MeanB3RAD)

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
            raise e

    def setBandNames(self, imageFile):
        dataset = gdal.Open(imageFile, gdal.GA_Update)
        if not dataset is None:
            dataset.GetRasterBand(1).SetDescription("Blue")
            dataset.GetRasterBand(2).SetDescription("Green")
            dataset.GetRasterBand(3).SetDescription("Red")
            dataset.GetRasterBand(4).SetDescription("RedEdge")
            dataset.GetRasterBand(5).SetDescription("NIR")
            dataset = None
        else:
            print("Could not open image to set band names: ", imageFile)

    def cleanLocalFollowProcessing(self):
        if not self.origFileName is '':
            rsgisUtils = rsgislib.RSGISPyUtils()
            rsgisUtils.deleteFileWithBasename(self.fileName)
            self.fileName = self.origFileName
            self.origFileName = ''


