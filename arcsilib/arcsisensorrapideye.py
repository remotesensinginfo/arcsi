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

class ARCSIRapidEyeSensor (ARCSIAbstractSensor):
    """
    A class which represents the RapidEye sensor to read
    header parameters and apply data processing operations.
    """
    def __init__(self):
        ARCSIAbstractSensor.__init__(self)
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
    
    def extractHeaderParameters(self, inputHeader, wktStr):
        """
        Understands and parses the RapidEye metadata.xml header file
        """
        try:
            print("Reading header file")
            tree = ET.parse(inputHeader)
            root = tree.getroot()
                
            eoPlatform = root.find('{http://www.opengis.net/gml}using').find('{http://earth.esa.int/eop}EarthObservationEquipment').find('{http://earth.esa.int/eop}platform').find('{http://earth.esa.int/eop}Platform')               
            self.platShortHand = eoPlatform.find('{http://earth.esa.int/eop}shortName').text.strip()
            print("self.platShortHand = ", self.platShortHand)
            self.platSerialId = eoPlatform.find('{http://earth.esa.int/eop}serialIdentifier').text.strip()
            print("self.platSerialId = ", self.platSerialId)
            self.platOrbitType = eoPlatform.find('{http://earth.esa.int/eop}orbitType').text.strip()
            print("self.platOrbitType = ", self.platOrbitType)
            
            #if (self.platSerialId != "RE-1") or (self.platSerialId != "RE-2") or (self.platSerialId != "RE-3") or (self.platSerialId != "RE-4"):
            #    raise ARCSIException("Do no recognise the spacecraft needs to be RE-1, RE-2, RE-3 or RE-4.")
            
            eoInstrument = root.find('{http://www.opengis.net/gml}using').find('{http://earth.esa.int/eop}EarthObservationEquipment').find('{http://earth.esa.int/eop}instrument').find('{http://earth.esa.int/eop}Instrument')
            self.instShortHand = eoInstrument.find('{http://earth.esa.int/eop}shortName').text.strip()
            print("self.instShortHand = ", self.instShortHand)
                
            eoSensor = root.find('{http://www.opengis.net/gml}using').find('{http://earth.esa.int/eop}EarthObservationEquipment').find('{http://earth.esa.int/eop}sensor').find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}Sensor')
            self.senrType = eoSensor.find('{http://earth.esa.int/eop}sensorType').text.strip()
            print("self.senrType = ", self.senrType)
            self.senrRes = float(eoSensor.find('{http://earth.esa.int/eop}resolution').text.strip())
            print("self.senrRes = ", self.senrRes)
            self.senrScanType = eoSensor.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}scanType').text.strip()
            print("self.senrScanType = ", self.senrScanType)
                
            eoAcquParams = root.find('{http://www.opengis.net/gml}using').find('{http://earth.esa.int/eop}EarthObservationEquipment').find('{http://earth.esa.int/eop}acquisitionParameters').find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}Acquisition')
            
            self.acquIncidAngle = float(eoAcquParams.find('{http://earth.esa.int/eop}incidenceAngle').text.strip())
            print("self.acquIncidAngle: ", self.acquIncidAngle)
            self.acquAzimuthAngle = float(eoAcquParams.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}azimuthAngle').text.strip())
            print("self.acquAzimuthAngle: ", self.acquAzimuthAngle)
            self.acquCraftViewAngle = float(eoAcquParams.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}spaceCraftViewAngle').text.strip())
            print("self.acquCraftViewAngle: ", self.acquCraftViewAngle)
            
            self.solarZenith = 90-float(eoAcquParams.find('{http://earth.esa.int/opt}illuminationElevationAngle').text.strip())
            print("self.solarZenith: ", self.solarZenith)
            self.solarAzimuth = float(eoAcquParams.find('{http://earth.esa.int/opt}illuminationAzimuthAngle').text.strip())
            print("self.solarAzimuth: ", self.solarAzimuth)
            self.senorZenith = self.acquCraftViewAngle
            print("self.senorZenith: ", self.senorZenith)
            self.senorAzimuth = self.acquAzimuthAngle
            print("self.senorAzimuth: ", self.senorAzimuth)
            timeStr = eoAcquParams.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}acquisitionDateTime').text.strip()
            timeStr = timeStr.replace('Z', '')
            try:
                self.acquisitionTime = datetime.datetime.strptime(timeStr, "%Y-%m-%dT%H:%M:%S.%f")
            except Exception as e:
                try:
                    self.acquisitionTime = datetime.datetime.strptime(timeStr, "%Y-%m-%dT%H:%M:%S")
                except Exception as e:
                    raise e
            print("self.acquisitionTime: ", self.acquisitionTime)
            
            metadata = root.find('{http://www.opengis.net/gml}metaDataProperty').find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}EarthObservationMetaData')
            self.tileID = metadata.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}tileId').text.strip()
            self.pixelFormat = metadata.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}pixelFormat').text.strip()
            print("self.tileID = ", self.tileID)
            print("self.pixelFormat = ", self.pixelFormat)
            
            
            centrePt = root.find('{http://www.opengis.net/gml}target').find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}Footprint').find('{http://www.opengis.net/gml}centerOf').find('{http://www.opengis.net/gml}Point').find('{http://www.opengis.net/gml}pos').text.strip()
            centrePtSplit = centrePt.split(' ')
            self.latCentre = float(centrePtSplit[0])
            self.lonCentre = float(centrePtSplit[1])
            print("self.latCentre = ", self.latCentre)
            print("self.lonCentre = ", self.lonCentre)
            
            imgBounds = root.find('{http://www.opengis.net/gml}target').find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}Footprint').find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}geographicLocation')
            tlPoint = imgBounds.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}topLeft')
            self.latTL = float(tlPoint.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}latitude').text)
            self.lonTL = float(tlPoint.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}longitude').text)
            trPoint = imgBounds.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}topRight')
            self.latTR = float(trPoint.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}latitude').text)
            self.lonTR = float(trPoint.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}longitude').text)
            brPoint = imgBounds.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}bottomRight')
            self.latBR = float(brPoint.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}latitude').text)
            self.lonBR = float(brPoint.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}longitude').text)
            blPoint = imgBounds.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}bottomLeft')
            self.latBL = float(blPoint.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}latitude').text)
            self.lonBL = float(blPoint.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}longitude').text)
                        
            print("self.latTL = ", self.latTL)
            print("self.lonTL = ", self.lonTL)
            print("self.latTR = ", self.latTR)
            print("self.lonTR = ", self.lonTR)           
            print("self.latBR = ", self.latBR)
            print("self.lonBR = ", self.lonBR)
            print("self.latBL = ", self.latBL)
            print("self.lonBL = ", self.lonBL)
            
            
            productInfo = root.find('{http://www.opengis.net/gml}resultOf').find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}EarthObservationResult').find('{http://earth.esa.int/eop}product').find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}ProductInformation')

            spatialRef = productInfo.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}spatialReferenceSystem')
            
            epsgCode = int(spatialRef.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}epsgCode').text)
            inProj = osr.SpatialReference()
            inProj.ImportFromEPSG(epsgCode)
            if self.inWKT == "":
                self.inWKT = inProj.ExportToWkt()
            print("WKT: ", self.inWKT)
            
            self.numOfBands = int(productInfo.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}numBands').text.strip())
            print('self.numOfBands = ', self.numOfBands)
            if self.numOfBands != 5:
                raise ARCSIException("The number of image band is not equal to 5 according to XML header.")
            
            radioCorrAppliedStr = productInfo.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}radiometricCorrectionApplied').text.strip()
            if radioCorrAppliedStr == "true":
                self.radioCorrApplied = True
            else:
                self.radioCorrApplied = False
            
            if self.radioCorrApplied:
                try:
                    self.radioCorrVersion = productInfo.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}radiometricCalibrationVersion').text.strip()
                except Exception:
                    self.radioCorrVersion = 'Unknown'
            else:
                self.radioCorrVersion = 'Not Applied'
            print('self.radioCorrVersion = ', self.radioCorrVersion)
            
            atmosCorrAppliedStr = productInfo.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}atmosphericCorrectionApplied').text.strip()
            if atmosCorrAppliedStr == "true":
                self.atmosCorrApplied = True
            else:
                self.atmosCorrApplied = False
            print('self.atmosCorrApplied = ', self.atmosCorrApplied)
            
            if self.atmosCorrApplied:
                raise ARCSIException("An atmosheric correction has already been applied according to the metadata.")
            
            elevCorrAppliedStr = productInfo.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}elevationCorrectionApplied').text.strip()
            if elevCorrAppliedStr == "true":
                self.elevCorrApplied = True
            else:
                self.elevCorrApplied = False
            print('self.elevCorrApplied = ', self.elevCorrApplied)
            
            self.geoCorrLevel = productInfo.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}geoCorrectionLevel').text.strip()
            print('self.geoCorrLevel = ', self.geoCorrLevel)
            
            filesDIR = os.path.dirname(inputHeader)
            self.fileName = os.path.join(filesDIR, productInfo.find('{http://earth.esa.int/eop}fileName').text.strip())
            print('self.fileName = ', self.fileName)
            
            # Haven't been defined yet!!
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
                        
        except Exception as e:
            raise e
        
    def generateOutputBaseName(self):
        """
        Provides an implementation for the landsat sensor
        """
        reTileID = "tid" + str(self.tileID)
        outname = self.defaultGenBaseOutFileName()
        outname = outname + str("_") + reTileID
        return outname
    
    def applyImageDataMask(self, inputHeader, outputPath, outputMaskName, outputImgName, outFormat):
    	raise ARCSIException("RapidEye does not provide any image masks, do not use the MASK option.")
        
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
        
        lsBand = collections.namedtuple('LSBand', ['bandName', 'fileName', 'bandIndex', 'satVal'])
        bandDefnSeq = list()
        bandDefnSeq.append(lsBand(bandName="Blue", fileName=self.fileName, bandIndex=1, satVal=65535))
        bandDefnSeq.append(lsBand(bandName="Green", fileName=self.fileName, bandIndex=2, satVal=65535))
        bandDefnSeq.append(lsBand(bandName="Red", fileName=self.fileName, bandIndex=3, satVal=65535))
        bandDefnSeq.append(lsBand(bandName="RedEdge", fileName=self.fileName, bandIndex=4, satVal=65535))
        bandDefnSeq.append(lsBand(bandName="NIR", fileName=self.fileName, bandIndex=5, satVal=65535))
        
        rsgislib.imagecalibration.saturatedPixelsMask(outputImage, outFormat, bandDefnSeq)
        
        return outputImage
    
    def convertThermalToBrightness(self, inputRadImage, outputPath, outputName, outFormat):
        raise ARCSIException("There are no thermal bands...")
    
    def convertImageToTOARefl(self, inputRadImage, outputPath, outputName, outFormat):
        print("Converting to TOA")
        outputImage = os.path.join(outputPath, outputName)
        solarIrradianceVals = list()
        IrrVal = collections.namedtuple('SolarIrradiance', ['irradiance'])
        solarIrradianceVals.append(IrrVal(irradiance=1997.8))
        solarIrradianceVals.append(IrrVal(irradiance=1863.5))
        solarIrradianceVals.append(IrrVal(irradiance=1560.4))
        solarIrradianceVals.append(IrrVal(irradiance=1395.0))
        solarIrradianceVals.append(IrrVal(irradiance=1124.4))
        rsgislib.imagecalibration.radiance2TOARefl(inputRadImage, outputImage, outFormat, rsgislib.TYPE_16UINT, 1000, self.acquisitionTime.year, self.acquisitionTime.month, self.acquisitionTime.day, self.solarZenith, solarIrradianceVals)
        return outputImage
    
    def generateCloudMask(self, inputReflImage, inputSatImage, inputThermalImage, outputPath, outputName, outFormat, tmpPath):
        print("Generating Cloud Mask")
        try:
            arcsiUtils = ARCSIUtils()
            tmpBaseName = os.path.splitext(outputName)[0]
            imgExtension = arcsiUtils.getFileExtension(outFormat)        
            stretchedImg = os.path.join(tmpPath, tmpBaseName + "_stchd" + imgExtension)
            outputCloudsImage = os.path.join(outputPath, outputName)
            outputCloudsImage = outputCloudsImage.replace(imgExtension, ".tif")
        
            rsgislib.imageutils.stretchImage(inputReflImage, stretchedImg, False, "", True, False, outFormat, rsgislib.TYPE_8UINT, rsgislib.imageutils.STRETCH_LINEARSTDDEV, 2)
            
            cloudsCmd = "recloud"
            cloudsOpts = "-i " + inputReflImage + " -s " + stretchedImg + " -o " + outputCloudsImage
            print(cloudsCmd + " " + cloudsOpts)
            
            os.system(cloudsCmd + " " + cloudsOpts)
            #subprocess.call([cloudsCmd, cloudsOpts]) # TODO: WHY DID THIS NOT WORK!?!?!?
            
            arcsiUtils.setImgThematic(outputCloudsImage)
            
            gdalDriver = gdal.GetDriverByName(outFormat)
            gdalDriver.Delete(stretchedImg)
            
            return outputCloudsImage
        except Exception as e:
            raise e    
        
        
    def calc6SCoefficients(self, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF):
        sixsCoeffs = numpy.zeros((5, 3), dtype=numpy.float32)    
        # Set up 6S model
        s = Py6S.SixS()
        s.atmos_profile = atmosProfile
        s.aero_profile = aeroProfile
        #s.ground_reflectance = Py6S.GroundReflectance.HomogeneousHapke(0.101, -0.263, 0.589, 0.046)
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
        
        # Band 1
        s.wavelength = Py6S.Wavelength(0.435, 0.515, [0.001, 0.004, 0.321, 0.719, 0.74, 0.756, 0.77, 0.78, 0.784, 0.792, 0.796, 0.799, 0.806, 0.804, 0.807, 0.816, 0.82, 0.825, 0.84, 0.845, 0.862, 0.875, 0.886, 0.905, 0.928, 0.936, 0.969, 0.967, 1, 0.976, 0.437, 0.029, 0.001])
        s.run()
        sixsCoeffs[0,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[0,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[0,2] = float(s.outputs.values['coef_xc'])
        
        # Band 2
        s.wavelength = Py6S.Wavelength(0.510, 0.5975, [0.001, 0.002, 0.013, 0.054, 0.539, 0.868, 0.868, 0.877, 0.871, 0.874, 0.882, 0.882, 0.881, 0.886, 0.897, 0.899, 0.901, 0.91, 0.924, 0.928, 0.936, 0.946, 0.953, 0.96, 0.974, 0.976, 0.976, 0.989, 0.988, 0.984, 0.994, 0.97, 0.417, 0.039, 0.002, 0.001])
        s.run()
        sixsCoeffs[1,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[1,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[1,2] = float(s.outputs.values['coef_xc'])
        
        # Band 3
        s.wavelength = Py6S.Wavelength(0.620, 0.6925, [0.001, 0.002, 0.009, 0.038, 0.437, 0.856, 0.854, 0.876, 0.881, 0.885, 0.902, 0.909, 0.915, 0.923, 0.939, 0.947, 0.958, 0.963, 0.97, 0.976, 0.989, 0.991, 0.985, 0.994, 0.989, 0.989, 0.463, 0.062, 0.005, 0.001])
        s.run()
        sixsCoeffs[2,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[2,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[2,2] = float(s.outputs.values['coef_xc'])
        
        # Band 4
        s.wavelength = Py6S.Wavelength(0.6775, 0.7425, [0.001, 0.002, 0.004, 0.021, 0.074, 0.491, 0.914, 0.998, 0.999, 0.998, 0.993, 0.987, 0.986, 0.982, 0.976, 0.966, 0.964, 0.961, 0.949, 0.939, 0.936, 0.425, 0.123, 0.02, 0.007, 0.002, 0.001])
        s.run()
        sixsCoeffs[3,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[3,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[3,2] = float(s.outputs.values['coef_xc'])
        
        # Band 5
        s.wavelength = Py6S.Wavelength(0.740, 0.870, [0.001, 0.001, 0.003, 0.005, 0.012, 0.023, 0.068, 0.153, 0.497, 0.828, 1, 0.982, 0.967, 0.974, 0.983, 0.981, 0.97, 0.963, 0.958, 0.957, 0.958, 0.959, 0.956, 0.954, 0.948, 0.944, 0.937, 0.933, 0.928, 0.927, 0.926, 0.926, 0.923, 0.918, 0.906, 0.898, 0.889, 0.885, 0.882, 0.876, 0.857, 0.842, 0.84, 0.832, 0.582, 0.295, 0.08, 0.034, 0.011, 0.006, 0.002, 0.001, 0.001])
        s.run()
        sixsCoeffs[4,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[4,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[4,2] = float(s.outputs.values['coef_xc'])
        
        return sixsCoeffs
    
    def convertImageToSurfaceReflSglParam(self, inputRadImage, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF):
        print("Converting to Surface Reflectance")
        outputImage = os.path.join(outputPath, outputName)
        
        Band6S = collections.namedtuple('Band6SCoeff', ['band', 'aX', 'bX', 'cX'])
        imgBandCoeffs = list()
        
        sixsCoeffs = self.calc6SCoefficients(aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF)
        
        imgBandCoeffs.append(Band6S(band=1, aX=float(sixsCoeffs[0,0]), bX=float(sixsCoeffs[0,1]), cX=float(sixsCoeffs[0,2])))
        imgBandCoeffs.append(Band6S(band=2, aX=float(sixsCoeffs[1,0]), bX=float(sixsCoeffs[1,1]), cX=float(sixsCoeffs[1,2])))
        imgBandCoeffs.append(Band6S(band=3, aX=float(sixsCoeffs[2,0]), bX=float(sixsCoeffs[2,1]), cX=float(sixsCoeffs[2,2])))
        imgBandCoeffs.append(Band6S(band=4, aX=float(sixsCoeffs[3,0]), bX=float(sixsCoeffs[3,1]), cX=float(sixsCoeffs[3,2])))
        imgBandCoeffs.append(Band6S(band=5, aX=float(sixsCoeffs[4,0]), bX=float(sixsCoeffs[4,1]), cX=float(sixsCoeffs[4,2])))
        
        for band in imgBandCoeffs:
            print(band)
        rsgislib.imagecalibration.apply6SCoeffSingleParam(inputRadImage, outputImage, outFormat, rsgislib.TYPE_16UINT, 1000, 0, True, imgBandCoeffs)
        return outputImage
        
    def convertImageToSurfaceReflDEMElevLUT(self, inputRadImage, inputDEMFile, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, aotVal, useBRDF, surfaceAltitudeMin, surfaceAltitudeMax):
        print("Converting to Surface Reflectance")
        outputImage = os.path.join(outputPath, outputName)        
        
        print("Build an LUT for elevation values.")    
        elev6SCoeffsLUT = self.buildElevation6SCoeffLUT(aeroProfile, atmosProfile, grdRefl, aotVal, useBRDF, surfaceAltitudeMin, surfaceAltitudeMax)
        print("LUT has been built.")
        
        elevLUTFeat = collections.namedtuple('ElevLUTFeat', ['Elev', 'Coeffs'])
        Band6S = collections.namedtuple('Band6SCoeff', ['band', 'aX', 'bX', 'cX'])
        
        elevCoeffs = list()
        for elevLUT in elev6SCoeffsLUT:
            imgBandCoeffs = list()
            sixsCoeffs = elevLUT.Coeffs
            elevVal = elevLUT.Elev
            imgBandCoeffs.append(Band6S(band=1, aX=float(sixsCoeffs[0,0]), bX=float(sixsCoeffs[0,1]), cX=float(sixsCoeffs[0,2])))
            imgBandCoeffs.append(Band6S(band=2, aX=float(sixsCoeffs[1,0]), bX=float(sixsCoeffs[1,1]), cX=float(sixsCoeffs[1,2])))
            imgBandCoeffs.append(Band6S(band=3, aX=float(sixsCoeffs[2,0]), bX=float(sixsCoeffs[2,1]), cX=float(sixsCoeffs[2,2])))
            imgBandCoeffs.append(Band6S(band=4, aX=float(sixsCoeffs[3,0]), bX=float(sixsCoeffs[3,1]), cX=float(sixsCoeffs[3,2])))
            imgBandCoeffs.append(Band6S(band=5, aX=float(sixsCoeffs[4,0]), bX=float(sixsCoeffs[4,1]), cX=float(sixsCoeffs[4,2])))
            
            elevCoeffs.append(elevLUTFeat(Elev=float(elevVal), Coeffs=imgBandCoeffs))
            
        rsgislib.imagecalibration.apply6SCoeffElevLUTParam(inputRadImage, inputDEMFile, outputImage, outFormat, rsgislib.TYPE_16UINT, 1000, 0, True, elevCoeffs)
        return outputImage
        
    
    def convertImageToSurfaceReflAOTDEMElevLUT(self, inputRadImage, inputDEMFile, inputAOTImage, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, useBRDF, surfaceAltitudeMin, surfaceAltitudeMax, aotMin, aotMax):
        print("Converting to Surface Reflectance")
        outputImage = os.path.join(outputPath, outputName) 
    
        print("Build an LUT for elevation and AOT values.")
        elevAOT6SCoeffsLUT = self.buildElevationAOT6SCoeffLUT(aeroProfile, atmosProfile, grdRefl, useBRDF, surfaceAltitudeMin, surfaceAltitudeMax, aotMin, aotMax)
                
        elevLUTFeat = collections.namedtuple('ElevLUTFeat', ['Elev', 'Coeffs'])
        aotLUTFeat = collections.namedtuple('AOTLUTFeat', ['AOT', 'Coeffs'])
        Band6S = collections.namedtuple('Band6SCoeff', ['band', 'aX', 'bX', 'cX'])
        
        elevAOTCoeffs = list()
        for elevLUT in elevAOT6SCoeffsLUT:
            elevVal = elevLUT.Elev
            aotLUT = elevLUT.Coeffs
            aot6SCoeffsOut = list()
            for aotFeat in aotLUT: 
                sixsCoeffs = aotFeat.Coeffs
                aotVal = aotFeat.AOT
                imgBandCoeffs = list()
                imgBandCoeffs.append(Band6S(band=1, aX=float(sixsCoeffs[0,0]), bX=float(sixsCoeffs[0,1]), cX=float(sixsCoeffs[0,2])))
                imgBandCoeffs.append(Band6S(band=2, aX=float(sixsCoeffs[1,0]), bX=float(sixsCoeffs[1,1]), cX=float(sixsCoeffs[1,2])))
                imgBandCoeffs.append(Band6S(band=3, aX=float(sixsCoeffs[2,0]), bX=float(sixsCoeffs[2,1]), cX=float(sixsCoeffs[2,2])))
                imgBandCoeffs.append(Band6S(band=4, aX=float(sixsCoeffs[3,0]), bX=float(sixsCoeffs[3,1]), cX=float(sixsCoeffs[3,2])))
                imgBandCoeffs.append(Band6S(band=5, aX=float(sixsCoeffs[4,0]), bX=float(sixsCoeffs[4,1]), cX=float(sixsCoeffs[4,2])))
                aot6SCoeffsOut.append(aotLUTFeat(AOT=float(aotVal), Coeffs=imgBandCoeffs))
            elevAOTCoeffs.append(elevLUTFeat(Elev=float(elevVal), Coeffs=aot6SCoeffsOut))
                        
        rsgislib.imagecalibration.apply6SCoeffElevAOTLUTParam(inputRadImage, inputDEMFile, inputAOTImage, outputImage, outFormat, rsgislib.TYPE_16UINT, 1000, 0, True, elevAOTCoeffs)
            
        return outputImage
    
    def convertImageToReflectanceDarkSubstract(self, inputTOAImage, outputPath, outputName, outFormat, tmpPath):
        try:
            print("Opening: ", inputTOAImage)
            toaDataset = gdal.Open(inputTOAImage, gdal.GA_ReadOnly)
            if toaDataset == None:
                raise Exception('Could not open the image dataset \'' + inputTOAImage + '\'')
            
            numBands = toaDataset.RasterCount
            toaDataset = None 
            
            print("Number of bands = ", numBands)
            
            darkPxlPercentile = 0.01
            minObjSize = 5
            
            offsetsImage = self.findPerBandDarkTargetsOffsets(inputTOAImage, numBands, outputPath, outputName, outFormat, tmpPath, minObjSize, darkPxlPercentile)
                       
            # TOA Image - Offset Image (if data and < 1 then set min value as 1)... 
            outputImage = os.path.join(outputPath, outputName)
            rsgislib.imagecalibration.applySubtractOffsets(inputTOAImage, offsetsImage, outputImage, outFormat, rsgislib.TYPE_16UINT, True, True, 0.0)
            
            return outputImage
            
        except Exception as e:
            raise e
        
    
    def findDDVTargets(self, inputTOAImage, outputPath, outputName, outFormat, tmpPath):
        print("Not implemented\n")
        sys.exit()
    
    def estimateImageToAOD(self, inputTOAImage, outputPath, outputName, outFormat, tmpPath, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotValMin, aotValMax):
        print("Not implemented\n")
        sys.exit()

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
            
            
            