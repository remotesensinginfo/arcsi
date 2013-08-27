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
from arcsisensor import ARCSIAbstractSensor
# Import the ARCSI exception class
from arcsiexception import ARCSIException
# Import the ARCSI utilities class
from arcsiutils import ARCSIUtils
# Import the datetime module
import datetime
# Import the GDAL/OGR spatial reference library
from osgeo import osr
from osgeo import ogr
# Import OS path module for manipulating the file system 
import os.path
# Import the RSGISLib Image Calibration Module.
import rsgislib.imagecalibration
# Import the collections module
import collections
# Import the py6s module for running 6S from python.
import Py6S
# Import the python maths library
import math
# Import python XML Parser
import xml.etree.ElementTree as ET

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
            for child in root:
                print child.tag
                
            eoPlatform = root.find('{http://www.opengis.net/gml}using').find('{http://earth.esa.int/eop}EarthObservationEquipment').find('{http://earth.esa.int/eop}platform').find('{http://earth.esa.int/eop}Platform')               
            self.platShortHand = eoPlatform.find('{http://earth.esa.int/eop}shortName').text.strip()
            print "self.platShortHand = ", self.platShortHand
            self.platSerialId = eoPlatform.find('{http://earth.esa.int/eop}serialIdentifier').text.strip()
            print "self.platSerialId = ", self.platSerialId
            self.platOrbitType = eoPlatform.find('{http://earth.esa.int/eop}orbitType').text.strip()
            print "self.platOrbitType = ", self.platOrbitType
            
            #if (self.platSerialId != "RE-1") or (self.platSerialId != "RE-2") or (self.platSerialId != "RE-3") or (self.platSerialId != "RE-4"):
            #    raise ARCSIException("Do no recognise the spacecraft needs to be RE-1, RE-2, RE-3 or RE-4.")
            
            eoInstrument = root.find('{http://www.opengis.net/gml}using').find('{http://earth.esa.int/eop}EarthObservationEquipment').find('{http://earth.esa.int/eop}instrument').find('{http://earth.esa.int/eop}Instrument')
            self.instShortHand = eoInstrument.find('{http://earth.esa.int/eop}shortName').text.strip()
            print "self.instShortHand = ", self.instShortHand
                
            eoSensor = root.find('{http://www.opengis.net/gml}using').find('{http://earth.esa.int/eop}EarthObservationEquipment').find('{http://earth.esa.int/eop}sensor').find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}Sensor')
            self.senrType = eoSensor.find('{http://earth.esa.int/eop}sensorType').text.strip()
            print "self.senrType = ", self.senrType
            self.senrRes = float(eoSensor.find('{http://earth.esa.int/eop}resolution').text.strip())
            print "self.senrRes = ", self.senrRes
            self.senrScanType = eoSensor.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}scanType').text.strip()
            print "self.senrScanType = ", self.senrScanType
                
            eoAcquParams = root.find('{http://www.opengis.net/gml}using').find('{http://earth.esa.int/eop}EarthObservationEquipment').find('{http://earth.esa.int/eop}acquisitionParameters').find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}Acquisition')
            
            self.acquIncidAngle = float(eoAcquParams.find('{http://earth.esa.int/eop}incidenceAngle').text.strip())
            print "self.acquIncidAngle: ", self.acquIncidAngle
            self.acquAzimuthAngle = float(eoAcquParams.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}azimuthAngle').text.strip())
            print "self.acquAzimuthAngle: ", self.acquAzimuthAngle
            self.acquCraftViewAngle = float(eoAcquParams.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}spaceCraftViewAngle').text.strip())
            print "self.acquCraftViewAngle: ", self.acquCraftViewAngle
            
            self.solarZenith = 90-float(eoAcquParams.find('{http://earth.esa.int/opt}illuminationElevationAngle').text.strip())
            print "self.solarZenith: ", self.solarZenith
            self.solarAzimuth = float(eoAcquParams.find('{http://earth.esa.int/opt}illuminationAzimuthAngle').text.strip())
            print "self.solarAzimuth: ", self.solarAzimuth
            self.senorZenith = self.acquCraftViewAngle
            print "self.senorZenith: ", self.senorZenith
            self.senorAzimuth = self.acquAzimuthAngle
            print "self.senorAzimuth: ", self.senorAzimuth
            timeStr = eoAcquParams.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}acquisitionDateTime').text.strip()
            timeStr = timeStr.replace('Z', '')
            self.acquisitionTime = datetime.datetime.strptime(timeStr, "%Y-%m-%dT%H:%M:%S.%f")
            print "self.acquisitionTime: ", self.acquisitionTime
            
            metadata = root.find('{http://www.opengis.net/gml}metaDataProperty').find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}EarthObservationMetaData')
            self.tileID = metadata.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}tileId').text.strip()
            self.pixelFormat = metadata.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}pixelFormat').text.strip()
            print "self.tileID = ", self.tileID
            print "self.pixelFormat = ", self.pixelFormat
            
            
            centrePt = root.find('{http://www.opengis.net/gml}target').find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}Footprint').find('{http://www.opengis.net/gml}centerOf').find('{http://www.opengis.net/gml}Point').find('{http://www.opengis.net/gml}pos').text.strip()
            centrePtSplit = centrePt.split(' ')
            self.latCentre = float(centrePtSplit[0])
            self.lonCentre = float(centrePtSplit[1])
            print "self.latCentre = ", self.latCentre
            print "self.lonCentre = ", self.lonCentre
            
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
                        
            print "self.latTL = ", self.latTL
            print "self.lonTL = ", self.lonTL
            print "self.latTR = ", self.latTR
            print "self.lonTR = ", self.lonTR            
            print "self.latBR = ", self.latBR
            print "self.lonBR = ", self.lonBR
            print "self.latBL = ", self.latBL
            print "self.lonBL = ", self.lonBL
            
            
            productInfo = root.find('{http://www.opengis.net/gml}resultOf').find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}EarthObservationResult').find('{http://earth.esa.int/eop}product').find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}ProductInformation')
            
            for child in productInfo:
                print child.tag
            
            spatialRef = productInfo.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}spatialReferenceSystem')
            
            epsgCode = int(spatialRef.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}epsgCode').text)
            inProj = osr.SpatialReference()
            inProj.ImportFromEPSG(epsgCode)
            if self.inWKT == "":
                self.inWKT = inProj.ExportToWkt()
            print "WKT: ", self.inWKT
            
            self.numOfBands = int(productInfo.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}numBands').text.strip())
            print 'self.numOfBands = ', self.numOfBands
            
            radioCorrAppliedStr = productInfo.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}radiometricCorrectionApplied').text.strip()
            if radioCorrAppliedStr == "true":
                self.radioCorrApplied = True
            self.radioCorrVersion = productInfo.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}radiometricCalibrationVersion').text.strip()
            print 'self.radioCorrApplied = ', self.radioCorrApplied
            print 'self.radioCorrVersion = ', self.radioCorrVersion
            
            atmosCorrAppliedStr = productInfo.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}atmosphericCorrectionApplied').text.strip()
            if atmosCorrAppliedStr == "true":
                self.atmosCorrApplied = True
            print 'self.atmosCorrApplied = ', self.atmosCorrApplied
            
            if self.atmosCorrApplied:
                raise ARCSIException("An atmosheric correction has already been applied according to the metadata.")
            
            elevCorrAppliedStr = productInfo.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}elevationCorrectionApplied').text.strip()
            if elevCorrAppliedStr == "true":
                self.elevCorrApplied = True
            print 'self.elevCorrApplied = ', self.elevCorrApplied
            
            self.geoCorrLevel = productInfo.find('{http://schemas.rapideye.de/products/productMetadataGeocorrected}geoCorrectionLevel').text.strip()
            print 'self.geoCorrLevel = ', self.geoCorrLevel
            
            self.fileName = productInfo.find('{http://earth.esa.int/eop}fileName').text.strip()
            print 'self.fileName = ', self.fileName
            
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
        
    def convertImageToRadiance(self, outputPath, outputName, outFormat):
        print("Converting to Radiance")
        outputImage = os.path.join(outputPath, outputName)
        if self.radioCorrApplied:
            # Rescale the data to be between 0 and 1.
        else:
            raise ARCSIException("Radiometric correction has not been applied - this is not implemented within ARCSI yet. Check your data version.")
        
        bandDefnSeq = list()
        lsBand = collections.namedtuple('LSBand', ['bandName', 'fileName', 'bandIndex', 'lMin', 'lMax', 'qCalMin', 'qCalMax'])
        bandDefnSeq.append(lsBand(bandName="Green", fileName=self.band4File, bandIndex=1, lMin=self.b4MinRad, lMax=self.b4MaxRad, qCalMin=self.b4CalMin, qCalMax=self.b4CalMax))
        bandDefnSeq.append(lsBand(bandName="Red", fileName=self.band5File, bandIndex=1, lMin=self.b5MinRad, lMax=self.b5MaxRad, qCalMin=self.b5CalMin, qCalMax=self.b5CalMax))
        bandDefnSeq.append(lsBand(bandName="NIR1", fileName=self.band6File, bandIndex=1, lMin=self.b6MinRad, lMax=self.b6MaxRad, qCalMin=self.b6CalMin, qCalMax=self.b6CalMax))
        bandDefnSeq.append(lsBand(bandName="NIR2", fileName=self.band7File, bandIndex=1, lMin=self.b7MinRad, lMax=self.b7MaxRad, qCalMin=self.b7CalMin, qCalMax=self.b7CalMax))
        rsgislib.imagecalibration.landsat2Radiance(outputImage, outFormat, bandDefnSeq)
        return outputImage
    
    def convertImageToTOARefl(self, inputRadImage, outputPath, outputName, outFormat):
        print("Converting to TOA")
        outputImage = os.path.join(outputPath, outputName)
        solarIrradianceVals = list()
        IrrVal = collections.namedtuple('SolarIrradiance', ['irradiance'])
        solarIrradianceVals.append(IrrVal(irradiance=1823.0))
        solarIrradianceVals.append(IrrVal(irradiance=1559.0))
        solarIrradianceVals.append(IrrVal(irradiance=1276.0))
        solarIrradianceVals.append(IrrVal(irradiance=880.1))
        rsgislib.imagecalibration.radiance2TOARefl(inputRadImage, outputImage, outFormat, rsgislib.TYPE_16UINT, 1000, self.acquisitionTime.year, self.acquisitionTime.month, self.acquisitionTime.day, self.solarZenith, solarIrradianceVals)
        return outputImage
        
    def convertImageToSurfaceReflSglParam(self, inputRadImage, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal):
        print("Converting to Surface Reflectance")
        outputImage = os.path.join(outputPath, outputName)
        
        Band6S = collections.namedtuple('Band6SCoeff', ['band', 'aX', 'bX', 'cX'])
        imgBandCoeffs = list()
        # Set up 6S model
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
        s.atmos_corr = Py6S.AtmosCorr.AtmosCorrLambertianFromReflectance(0.40)
        s.aot550 = aotVal
        
        # Band 1
        s.wavelength = Py6S.Wavelength(Py6S.SixSHelpers.PredefinedWavelengths.LANDSAT_MSS_B1)
        s.run()
        imgBandCoeffs.append(Band6S(band=1, aX=s.outputs.values['coef_xa'], bX=s.outputs.values['coef_xb'], cX=s.outputs.values['coef_xc']))
        
        # Band 2
        s.wavelength = Py6S.Wavelength(Py6S.SixSHelpers.PredefinedWavelengths.LANDSAT_MSS_B2)
        s.run()
        imgBandCoeffs.append(Band6S(band=2, aX=s.outputs.values['coef_xa'], bX=s.outputs.values['coef_xb'], cX=s.outputs.values['coef_xc']))
        
        # Band 3
        s.wavelength = Py6S.Wavelength(Py6S.SixSHelpers.PredefinedWavelengths.LANDSAT_MSS_B3)
        s.run()
        imgBandCoeffs.append(Band6S(band=3, aX=s.outputs.values['coef_xa'], bX=s.outputs.values['coef_xb'], cX=s.outputs.values['coef_xc']))
        
        # Band 4
        s.wavelength = Py6S.Wavelength(Py6S.SixSHelpers.PredefinedWavelengths.LANDSAT_MSS_B4)
        s.run()
        imgBandCoeffs.append(Band6S(band=4, aX=s.outputs.values['coef_xa'], bX=s.outputs.values['coef_xb'], cX=s.outputs.values['coef_xc']))
        
        for band in imgBandCoeffs:
            print(band)
        
        rsgislib.imagecalibration.apply6SCoeffSingleParam(inputRadImage, outputImage, outFormat, rsgislib.TYPE_16UINT, 1000, 0, True, imgBandCoeffs)
        
        return outputImage


    def estimateImageToAOD(self, inputTOAImage, outputPath, outputName, outFormat, tmpPath, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotValMin, aotValMax):
        print("Not implemented\n")
        sys.exit()


