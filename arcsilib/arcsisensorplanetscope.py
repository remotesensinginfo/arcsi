"""
Module that contains the ARCSIRapidEyeSensor class.
"""
############################################################################
#  arcsisensorrapideye.py
#
#  Copyright 2018 ARCSI.
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
# Date: 13/10/2018
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
# Import python XML Parser
import xml.etree.ElementTree as ET
# Import the numpy module
import numpy
# Import the GDAL python module
import osgeo.gdal as gdal


class _RadPlanetScopeBandInfo():

    def __init__(self, band, radCoef, reflCoef):
        self.band = band
        self.radCoef = radCoef
        self.reflCoef = reflCoef


class ARCSIPlanetScopeSensor (ARCSIAbstractSensor):
    """
    A class which represents the RapidEye sensor to read
    header parameters and apply data processing operations.
    """
    def __init__(self, debugMode, inputImage):
        ARCSIAbstractSensor.__init__(self, debugMode, inputImage)
        self.sensor = "PlanetScope"
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
        self.radioCoefficentsDict = dict()

    def extractHeaderParameters(self, inputHeader, wktStr):
        """
        Understands and parses the RapidEye metadata.xml header file
        """
        try:
            print("Reading header file")
            self.headerFileName = os.path.split(inputHeader)[1]

            hdrExt = os.path.splitext(inputHeader)
            if not len(hdrExt) is 2:
                raise ARCSIException("Cannot work out what the file extension is - supports xml.")
            hdrExt = hdrExt[1]

            if (hdrExt.lower() == '.xml') or (hdrExt.lower() == 'xml'):
                tree = ET.parse(inputHeader)
                root = tree.getroot()
                hdrVersion = root.attrib['version'].strip() # 1.2.1 when this was implemented.
                schemaURL = root.attrib['{http://www.w3.org/2001/XMLSchema-instance}schemaLocation'].strip().split()[0]
                planetScopeUrl = '{'+schemaURL+'}'
                metaDataProperty = root.find('{http://www.opengis.net/gml}metaDataProperty')
                eoMetaData = metaDataProperty.find(planetScopeUrl+'EarthObservationMetaData')
                if eoMetaData is None:
                    planetScopeUrl = '{http://schemas.rapideye.de/products/productMetadataSensor}'
                    eoMetaData = metaDataProperty.find(planetScopeUrl+'EarthObservationMetaData')
                productType = eoMetaData.find('{http://earth.esa.int/eop}productType').text.strip()

                if (productType != "L3A"):
                    raise ARCSIException("Only L3A data is supported by ARCSI.")

                eoPlatform = root.find('{http://www.opengis.net/gml}using').find('{http://earth.esa.int/eop}EarthObservationEquipment').find('{http://earth.esa.int/eop}platform').find('{http://earth.esa.int/eop}Platform')
                self.platShortHand = eoPlatform.find('{http://earth.esa.int/eop}shortName').text.strip()
                self.platSerialId = eoPlatform.find('{http://earth.esa.int/eop}serialIdentifier').text.strip()
                self.platOrbitType = eoPlatform.find('{http://earth.esa.int/eop}orbitType').text.strip()

                eoInstrument = root.find('{http://www.opengis.net/gml}using').find('{http://earth.esa.int/eop}EarthObservationEquipment').find('{http://earth.esa.int/eop}instrument').find('{http://earth.esa.int/eop}Instrument')
                self.instShortHand = eoInstrument.find('{http://earth.esa.int/eop}shortName').text.strip()

                eoSensor = root.find('{http://www.opengis.net/gml}using').find('{http://earth.esa.int/eop}EarthObservationEquipment').find('{http://earth.esa.int/eop}sensor').find(planetScopeUrl+'Sensor')
                self.senrType = eoSensor.find('{http://earth.esa.int/eop}sensorType').text.strip()
                self.senrRes = float(eoSensor.find('{http://earth.esa.int/eop}resolution').text.strip())
                self.senrScanType = eoSensor.find(planetScopeUrl+'scanType').text.strip()

                eoAcquParams = root.find('{http://www.opengis.net/gml}using').find('{http://earth.esa.int/eop}EarthObservationEquipment').find('{http://earth.esa.int/eop}acquisitionParameters').find(planetScopeUrl+'Acquisition')

                self.acquIncidAngle = float(eoAcquParams.find('{http://earth.esa.int/eop}incidenceAngle').text.strip())
                self.acquAzimuthAngle = float(eoAcquParams.find(planetScopeUrl+'azimuthAngle').text.strip())
                self.acquCraftViewAngle = float(eoAcquParams.find(planetScopeUrl+'spaceCraftViewAngle').text.strip())

                self.solarZenith = 90-float(eoAcquParams.find('{http://earth.esa.int/opt}illuminationElevationAngle').text.strip())
                self.solarAzimuth = float(eoAcquParams.find('{http://earth.esa.int/opt}illuminationAzimuthAngle').text.strip())
                self.sensorZenith = self.acquCraftViewAngle
                self.sensorAzimuth = self.acquAzimuthAngle
                timeStr = eoAcquParams.find(planetScopeUrl+'acquisitionDateTime').text.strip()
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

                metadata = root.find('{http://www.opengis.net/gml}metaDataProperty').find(planetScopeUrl+'EarthObservationMetaData')
                if not  metadata.find(planetScopeUrl+'tileId') is None:
                    self.tileID = metadata.find(planetScopeUrl+'tileId').text.strip()
                else:
                    self.tileID = ""
                self.pixelFormat = metadata.find(planetScopeUrl+'pixelFormat').text.strip()

                centrePt = root.find('{http://www.opengis.net/gml}target').find(planetScopeUrl+'Footprint').find('{http://www.opengis.net/gml}centerOf').find('{http://www.opengis.net/gml}Point').find('{http://www.opengis.net/gml}pos').text.strip()
                centrePtSplit = centrePt.split(' ')
                self.latCentre = float(centrePtSplit[0])
                self.lonCentre = float(centrePtSplit[1])

                imgBounds = root.find('{http://www.opengis.net/gml}target').find(planetScopeUrl+'Footprint').find(planetScopeUrl+'geographicLocation')
                tlPoint = imgBounds.find(planetScopeUrl+'topLeft')
                self.latTL = float(tlPoint.find(planetScopeUrl+'latitude').text)
                self.lonTL = float(tlPoint.find(planetScopeUrl+'longitude').text)
                trPoint = imgBounds.find(planetScopeUrl+'topRight')
                self.latTR = float(trPoint.find(planetScopeUrl+'latitude').text)
                self.lonTR = float(trPoint.find(planetScopeUrl+'longitude').text)
                brPoint = imgBounds.find(planetScopeUrl+'bottomRight')
                self.latBR = float(brPoint.find(planetScopeUrl+'latitude').text)
                self.lonBR = float(brPoint.find(planetScopeUrl+'longitude').text)
                blPoint = imgBounds.find(planetScopeUrl+'bottomLeft')
                self.latBL = float(blPoint.find(planetScopeUrl+'latitude').text)
                self.lonBL = float(blPoint.find(planetScopeUrl+'longitude').text)

                productInfo = root.find('{http://www.opengis.net/gml}resultOf').find(planetScopeUrl+'EarthObservationResult').find('{http://earth.esa.int/eop}product').find(planetScopeUrl+'ProductInformation')

                spatialRef = productInfo.find(planetScopeUrl+'spatialReferenceSystem')

                epsgCode = int(spatialRef.find(planetScopeUrl+'epsgCode').text)
                inProj = osr.SpatialReference()
                inProj.ImportFromEPSG(epsgCode)
                if self.inWKT is "":
                    self.inWKT = inProj.ExportToWkt()

                self.numOfBands = int(productInfo.find(planetScopeUrl+'numBands').text.strip())
                if not ((self.numOfBands == 3) or (self.numOfBands == 4)):
                    raise ARCSIException("The number of image band is not equal to 3 (RGB) or 4 (RGBNIR) according to XML header.")

                radioCorrAppliedStr = productInfo.find(planetScopeUrl+'radiometricCorrectionApplied').text.strip()
                if radioCorrAppliedStr == "true":
                    self.radioCorrApplied = True
                else:
                    self.radioCorrApplied = False

                if self.radioCorrApplied:
                    try:
                        self.radioCorrVersion = productInfo.find(planetScopeUrl+'radiometricCalibrationVersion').text.strip()
                    except Exception:
                        self.radioCorrVersion = 'Unknown'
                else:
                    self.radioCorrVersion = 'Not Applied'

                atmosCorrAppliedStr = productInfo.find(planetScopeUrl+'atmosphericCorrectionApplied').text.strip()
                if atmosCorrAppliedStr == "true":
                    self.atmosCorrApplied = True
                else:
                    self.atmosCorrApplied = False

                if self.atmosCorrApplied:
                    raise ARCSIException("An atmosheric correction has already been applied according to the metadata.")

                elevCorrAppliedStr = productInfo.find(planetScopeUrl+'elevationCorrectionApplied').text.strip()
                if elevCorrAppliedStr == "true":
                    self.elevCorrApplied = True
                else:
                    self.elevCorrApplied = False

                self.geoCorrLevel = productInfo.find(planetScopeUrl+'geoCorrectionLevel').text.strip()

                filesDIR = os.path.dirname(inputHeader)
                if not self.userSpInputImage is None:
                    self.fileName = os.path.abspath(self.userSpInputImage)
                else:
                    self.fileName = os.path.join(filesDIR, productInfo.find('{http://earth.esa.int/eop}fileName').text.strip())

                bandSpecificMetadata = root.find('{http://www.opengis.net/gml}resultOf').find(planetScopeUrl + 'EarthObservationResult').findall(planetScopeUrl + 'bandSpecificMetadata')
                for bandMetadata in bandSpecificMetadata:
                    bandN = int(bandMetadata.find(planetScopeUrl + 'bandNumber').text.strip())
                    radianceCoeff = float(bandMetadata.find(planetScopeUrl + 'radiometricScaleFactor').text.strip())
                    toaReflCoeff = float(bandMetadata.find(planetScopeUrl + 'reflectanceCoefficient').text.strip())
                    self.radioCoefficentsDict[bandN] = _RadPlanetScopeBandInfo(bandN, radianceCoeff, toaReflCoeff)

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
                self.xCentre = minX + ((maxX-minX)/2)
                self.yCentre = minY + ((maxY-minY)/2)
            else:
                raise ARCSIException("Header file extention is not recognised - supports xml.")



        except Exception as e:
            raise e

    def checkInputImageValid(self):
        if not self.expectedImageDataPresent():
            raise ARCSIException("Error image image was not present.")
        rasterDS = gdal.Open(self.fileName, gdal.GA_ReadOnly)
        if rasterDS == None:
            raise ARCSIException('Could not open raster image: ' + self.fileName)
        nBands = rasterDS.RasterCount
        if not ((nBands == 3) or (nBands == 4)):
            raise ARCSIException('Input Image \'' + self.fileName + '\' does not have the expected 3 or 4 image bands.')

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
            rsgisUtils = rsgislib.RSGISPyUtils()
            nBands = rsgisUtils.getImageBandCount(self.fileName)
            bandDefnSeq = list()
            imgBand = collections.namedtuple('ImgBand', ['bandName', 'bandIndex', 'bias', 'gain'])
            bandDefnSeq.append(imgBand(bandName="Blue", bandIndex=1, bias=0.0, gain=self.radioCoefficentsDict[1].radCoef))
            bandDefnSeq.append(imgBand(bandName="Green", bandIndex=2, bias=0.0, gain=self.radioCoefficentsDict[2].radCoef))
            bandDefnSeq.append(imgBand(bandName="Red", bandIndex=3, bias=0.0, gain=self.radioCoefficentsDict[3].radCoef))
            if nBands >= 4:
                bandDefnSeq.append(imgBand(bandName="NIR", bandIndex=4, bias=0.0, gain=self.radioCoefficentsDict[4].radCoef))
            rsgislib.imagecalibration.spot5ToRadiance(self.fileName, outputImage, outFormat, bandDefnSeq)


            # Rescale the data to be between 0 and 1.
            rsgislib.imagecalc.imageMath(self.fileName, outputImage, "b1/100", outFormat, rsgislib.TYPE_32FLOAT)
        else:
            raise ARCSIException("Radiometric correction has not been applied - this is not implemented within ARCSI yet. Check your data version.")

        return outputImage, None

    def generateImageSaturationMask(self, outputPath, outputName, outFormat):
        print("Generate Saturation Image")
        outputImage = os.path.join(outputPath, outputName)
        rsgisUtils = rsgislib.RSGISPyUtils()
        nBands = rsgisUtils.getImageBandCount(self.fileName)

        reBand = collections.namedtuple('REBand', ['bandName', 'fileName', 'bandIndex', 'satVal'])
        bandDefnSeq = list()
        if nBands >= 3:
            bandDefnSeq.append(reBand(bandName="Blue", fileName=self.fileName, bandIndex=1, satVal=65535.0))
            bandDefnSeq.append(reBand(bandName="Green", fileName=self.fileName, bandIndex=2, satVal=65535.0))
            bandDefnSeq.append(reBand(bandName="Red", fileName=self.fileName, bandIndex=3, satVal=65535.0))
            if nBands >= 4:
                bandDefnSeq.append(reBand(bandName="NIR", fileName=self.fileName, bandIndex=4, satVal=65535.0))
        else:
            raise Exception("The image must have at least 3 bands.")

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
        rsgisUtils = rsgislib.RSGISPyUtils()
        nBands = rsgisUtils.getImageBandCount(self.fileName)

        arcsiUtils = ARCSIUtils()
        dn = 1000.0
        b1_rad = dn * self.radioCoefficentsDict[1].radCoef
        b1_toa = dn * self.radioCoefficentsDict[1].reflCoef
        b1_solar_irr = arcsiUtils.getESUNValue(b1_rad, b1_toa, self.acquisitionTime.day, self.acquisitionTime.month,
                                               self.acquisitionTime.year, self.solarZenith)
        print("b1_solar_irr = {}".format(b1_solar_irr))

        b2_rad = dn * self.radioCoefficentsDict[2].radCoef
        b2_toa = dn * self.radioCoefficentsDict[2].reflCoef
        b2_solar_irr = arcsiUtils.getESUNValue(b2_rad, b2_toa, self.acquisitionTime.day, self.acquisitionTime.month,
                                               self.acquisitionTime.year, self.solarZenith)
        print("b2_solar_irr = {}".format(b2_solar_irr))

        b3_rad = dn * self.radioCoefficentsDict[3].radCoef
        b3_toa = dn * self.radioCoefficentsDict[3].reflCoef
        b3_solar_irr = arcsiUtils.getESUNValue(b3_rad, b3_toa, self.acquisitionTime.day, self.acquisitionTime.month,
                                               self.acquisitionTime.year, self.solarZenith)
        print("b3_solar_irr = {}".format(b3_solar_irr))

        outputImage = os.path.join(outputPath, outputName)
        solarIrradianceVals = list()
        IrrVal = collections.namedtuple('SolarIrradiance', ['irradiance'])
        solarIrradianceVals.append(IrrVal(irradiance=b1_solar_irr))
        solarIrradianceVals.append(IrrVal(irradiance=b2_solar_irr))
        solarIrradianceVals.append(IrrVal(irradiance=b3_solar_irr))
        if nBands >= 4:
            b4_rad = dn * self.radioCoefficentsDict[4].radCoef
            b4_toa = dn * self.radioCoefficentsDict[4].reflCoef
            b4_solar_irr = arcsiUtils.getESUNValue(b4_rad, b4_toa, self.acquisitionTime.day, self.acquisitionTime.month,
                                                   self.acquisitionTime.year, self.solarZenith)
            print("b4_solar_irr = {}".format(b4_solar_irr))
            solarIrradianceVals.append(IrrVal(irradiance=b4_solar_irr))

        rsgislib.imagecalibration.radiance2TOARefl(inputRadImage, outputImage, outFormat, rsgislib.TYPE_16UINT, scaleFactor, self.acquisitionTime.year, self.acquisitionTime.month, self.acquisitionTime.day, self.solarZenith, solarIrradianceVals)
        return outputImage

    def generateCloudMask(self, inputReflImage, inputSatImage, inputThermalImage, inputViewAngleImg, inputValidImg, outputPath, outputName, outFormat, tmpPath, scaleFactor):
        raise ARCSIException("Cloud Masking Not Implemented for PlanetScope.")

    def createCloudMaskDataArray(self, inImgDataArr):
        # Calc Whiteness
        meanArr = numpy.mean(inImgDataArr, axis=1)
        whitenessArr = numpy.absolute((inImgDataArr[...,0] - meanArr)/meanArr) + numpy.absolute((inImgDataArr[...,1] - meanArr)/meanArr) + numpy.absolute((inImgDataArr[...,2] - meanArr)/meanArr) + numpy.absolute((inImgDataArr[...,3] - meanArr)/meanArr)
        # Calc NDVI
        #ndvi = (inImgDataArr[...,4] - inImgDataArr[...,2]) / (inImgDataArr[...,4] + inImgDataArr[...,2])
        
        # Create and populate the output array.
        inShape = inImgDataArr.shape
        outShape = [inShape[0], inShape[1]+2]
        outArr = numpy.zeros(outShape, dtype=float)
        
        for i in range(inShape[1]):
            outArr[...,i] = inImgDataArr[...,i]
        
        idx = inShape[1]
        outArr[...,idx] = meanArr
        outArr[...,idx+1] = whitenessArr
        #outArr[...,idx+2] = ndvi
        
        return outArr

    def defineDarkShadowImageBand(self):
        raise ARCSIException("For the PlanetScope Sensor there is not implementation of defineDarkShadowImageBand")

    def calc6SCoefficients(self, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF):
        raise ARCSIException("For the PlanetScope Sensor there is not implementation of calc6SCoefficients")

        return sixsCoeffs

    def convertImageToSurfaceReflSglParam(self, inputRadImage, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF, scaleFactor):
        raise ARCSIException("For the PlanetScope Sensor there is not implementation of convertImageToSurfaceReflSglParam")

    def convertImageToSurfaceReflDEMElevLUT(self, inputRadImage, inputDEMFile, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, aotVal, useBRDF, surfaceAltitudeMin, surfaceAltitudeMax, scaleFactor, elevCoeffs=None):
        raise ARCSIException("For the PlanetScope Sensor there is not implementation of convertImageToSurfaceReflDEMElevLUT")

    def convertImageToSurfaceReflAOTDEMElevLUT(self, inputRadImage, inputDEMFile, inputAOTImage, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, useBRDF, surfaceAltitudeMin, surfaceAltitudeMax, aotMin, aotMax, scaleFactor, elevAOTCoeffs=None):
        raise ARCSIException("For the PlanetScope Sensor there is not implementation of convertImageToSurfaceReflAOTDEMElevLUT")

    def run6SToOptimiseAODValue(self, aotVal, radBlueVal, predBlueVal, aeroProfile, atmosProfile, grdRefl, surfaceAltitude):
        """Used as part of the optimastion for identifying values of AOD"""
        raise ARCSIException("For the PlanetScope Sensor there is not implementation of run6SToOptimiseAODValue")

    def findDDVTargets(self, inputTOAImage, outputPath, outputName, outFormat, tmpPath):
        raise ARCSIException("For the PlanetScope Sensor there is not implementation of findDDVTargets")

    def estimateImageToAODUsingDDV(self, inputRADImage, inputTOAImage, inputDEMFile, shadowMask, outputPath, outputName, outFormat, tmpPath, aeroProfile, atmosProfile, grdRefl, aotValMin, aotValMax):
        raise ARCSIException("For the PlanetScope Sensor there is not implementation of estimateImageToAODUsingDDV")

    def estimateImageToAODUsingDOS(self, inputRADImage, inputTOAImage, inputDEMFile, shadowMask, outputPath, outputName, outFormat, tmpPath, aeroProfile, atmosProfile, grdRefl, aotValMin, aotValMax, globalDOS, simpleDOS, dosOutRefl):
        raise ARCSIException("For the PlanetScope Sensor there is not implementation of estimateImageToAODUsingDOS")

    def estimateSingleAOTFromDOS(self, radianceImage, toaImage, inputDEMFile, tmpPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, minAOT, maxAOT, dosOutRefl):
        raise ARCSIException("For the PlanetScope Sensor there is not implementation of estimateSingleAOTFromDOS")

    def setBandNames(self, imageFile):
        rsgisUtils = rsgislib.RSGISPyUtils()
        nBands = rsgisUtils.getImageBandCount(imageFile)
        dataset = gdal.Open(imageFile, gdal.GA_Update)
        if not dataset is None:
            dataset.GetRasterBand(1).SetDescription("Blue")
            dataset.GetRasterBand(2).SetDescription("Green")
            dataset.GetRasterBand(3).SetDescription("Red")
            if nBands == 4:
                dataset.GetRasterBand(4).SetDescription("NIR")
            dataset = None
        else:
            print("Could not open image to set band names: ", imageFile)

    def cleanLocalFollowProcessing(self):
        if not self.origFileName is '':
            rsgisUtils = rsgislib.RSGISPyUtils()
            rsgisUtils.deleteFileWithBasename(self.fileName)
            self.fileName = self.origFileName
            self.origFileName = ''


