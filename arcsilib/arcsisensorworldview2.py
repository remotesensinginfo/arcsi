"""
Module that contains the ARCSIWorldView2Sensor class.
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
# Date: 20/06/2014
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
# Import the RIOS RAT module
from rios import rat

class ARCSIWorldView2Sensor (ARCSIAbstractSensor):
    """
    A class which represents the WorldView2 sensor to read
    header parameters and apply data processing operations.
    """
    def __init__(self):
        ARCSIAbstractSensor.__init__(self)
        self.sensor = "WV2"
        self.fileName = ""
        self.catID = ""
        self.solarZenithMin = 0.0
        self.solarZenithMax = 0.0
        self.solarAzimuthMin = 0.0
        self.solarAzimuthMax = 0.0
        self.senorZenithMin = 0.0
        self.senorZenithMax = 0.0
        self.senorAzimuthMin = 0.0
        self.senorAzimuthMax = 0.0
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
            print("Reading header file")
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
                 
            #for child in topLevelInfo:
            #   print(child.tag, child.attrib)
                 
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
            print("Acqusition Time: ", self.acquisitionTime)
            
            self.solarZenithMin = 90-float(imageInfoTag.find('MINSUNEL').text.strip())
            self.solarZenithMax = 90-float(imageInfoTag.find('MAXSUNEL').text.strip())
            self.solarZenith = 90-float(imageInfoTag.find('MEANSUNEL').text.strip())
                        
            self.solarAzimuthMin = float(imageInfoTag.find('MINSUNAZ').text.strip())
            self.solarAzimuthMax = float(imageInfoTag.find('MAXSUNAZ').text.strip())
            self.solarAzimuth = float(imageInfoTag.find('MEANSATEL').text.strip())
            
            self.senorZenithMin = float(imageInfoTag.find('MINSATEL').text.strip())
            self.senorZenithMax = float(imageInfoTag.find('MAXSATEL').text.strip())
            self.senorZenith = float(imageInfoTag.find('MINSUNAZ').text.strip())
            
            self.senorAzimuthMin = float(imageInfoTag.find('MINSATAZ').text.strip())
            self.senorAzimuthMax = float(imageInfoTag.find('MAXSATAZ').text.strip())
            self.senorAzimuth = float(imageInfoTag.find('MEANSATAZ').text.strip())
            
            self.nadirViewAngle = float(imageInfoTag.find('MEANOFFNADIRVIEWANGLE').text.strip())


            if mapProjInfoTag.find('MAPPROJNAME').text.strip() != "UTM":
                raise ARCSIException("Expecting input image to be projected as UTM WGS84")
            
            utmZone = "WGS84" + mapProjInfoTag.find('MAPPROJNAME').text.strip() + mapProjInfoTag.find('MAPZONE').text.strip() + mapProjInfoTag.find('MAPHEMI').text.strip()            
            epsgCode = self.epsgCodes[utmZone]
            inProj = osr.SpatialReference()
            inProj.ImportFromEPSG(epsgCode)
            if self.inWKT == "":
                self.inWKT = inProj.ExportToWkt()
            #print("WKT: ", self.inWKT)
            
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
            
            wgs84latlonProj = osr.SpatialReference()
            wgs84latlonProj.ImportFromEPSG(4326)
            
            wktPt = 'POINT(%s %s)' % (self.xCentre, self.yCentre)
            #print(wktPt)
            point = ogr.CreateGeometryFromWkt(wktPt)
            point.AssignSpatialReference(inProj)
            point.TransformTo(wgs84latlonProj)            
            
            self.latCentre = point.GetY()
            self.lonCentre = point.GetX()
            #print("Lat: " + str(self.latCentre) + " Long: " + str(self.lonCentre))
            
            self.fileName = imageTileInfoTag.find('FILENAME').text.strip()
            
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
                        
            print("Processing Input File: ", self.fileName)
                        
        except Exception as e:
            raise e
        
    def generateOutputBaseName(self):
        """
        Provides an implementation for the landsat sensor
        """
        outname = self.defaultGenBaseOutFileName()
        return outname
    
    def applyImageDataMask(self, inputHeader, outputPath, outputMaskName, outputImgName, outFormat):
        raise ARCSIException("WorldView2 does not provide any image masks, do not use the MASK option.")
        
    def convertImageToRadiance(self, outputPath, outputReflName, outputThermalName, outFormat):
        raise ARCSIException("Radiance conversion not yet implemented for WV2...")
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
    
    def run6SToOptimiseAODValue(self, aotVal, radBlueVal, predBlueVal, aeroProfile, atmosProfile, grdRefl, surfaceAltitude):
        """Used as part of the optimastion for identifying values of AOD"""
        print("Testing AOD Val: ", aotVal,)
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
    
    def convertImageToReflectanceDarkSubstract(self, inputTOAImage, outputPath, outputName, outFormat, tmpPath, globalDOS, dosOutRefl):
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
            offsetsImage = ""
            
            if globalDOS:
            	offsetsImage = self.findPerBandDarkTargetsOffsets(inputTOAImage, numBands, outputPath, outputName, outFormat, tmpPath, minObjSize, darkPxlPercentile)
            else:
            	blockSize = 200
            	offsetsImage = self.findPerBandLocalDarkTargetsOffsets(inputTOAImage, numBands, outputPath, outputName, outFormat, tmpPath, blockSize, minObjSize, darkPxlPercentile)
            
                       
            # TOA Image - Offset Image (if data and < 1 then set min value as 1)... 
            outputImage = os.path.join(outputPath, outputName)
            rsgislib.imagecalibration.applySubtractOffsets(inputTOAImage, offsetsImage, outputImage, outFormat, rsgislib.TYPE_16UINT, True, True, 0.0, dosOutRefl)
            
            return outputImage
            
        except Exception as e:
            raise e
        
    def findDDVTargets(self, inputTOAImage, outputPath, outputName, outFormat, tmpPath):
        print("Not implemented\n")
        sys.exit()
    
    def estimateImageToAOD(self, inputRADImage, inputTOAImage, inputDEMFile, outputPath, outputName, outFormat, tmpPath, aeroProfile, atmosProfile, grdRefl, aotValMin, aotValMax):
        print("Not implemented\n")
        sys.exit()
        
    def estimateImageToAODUsingDOS(self, inputRADImage, inputTOAImage, inputDEMFile, outputPath, outputName, outFormat, tmpPath, aeroProfile, atmosProfile, grdRefl, aotValMin, aotValMax, globalDOS, dosOutRefl):
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
            if globalDOS:
            	dosBlueImage = self.performDOSOnSingleBand(inputTOAImage, 1, outputPath, tmpBaseName, "Blue", outFormat, tmpPath, minObjSize, darkPxlPercentile, dosOutRefl)
            else:
	            dosBlueImage = self.performLocalDOSOnSingleBand(inputTOAImage, 1, outputPath, tmpBaseName, "Blue", outFormat, tmpPath, minObjSize, darkPxlPercentile, blockSize, dosOutRefl)
                        
            thresImageClumpsFinal = os.path.join(tmpPath, tmpBaseName + "_clumps" + imgExtension)
            rsgislib.segmentation.segutils.runShepherdSegmentation(inputTOAImage, thresImageClumpsFinal, tmpath=tmpPath, gdalFormat=outFormat, numClusters=40, minPxls=10, bands=[5,4,1])
            
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
            self.interpolateImageFromPointData(inputTOAImage, Eastings, Northings, aotVals, outputAOTImage, outFormat, interpSmoothing)
                    
            gdalDriver = gdal.GetDriverByName(outFormat)
            gdalDriver.Delete(thresImageClumpsFinal)
            gdalDriver.Delete(dosBlueImage)        
        
            return outputAOTImage
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
            
            
            