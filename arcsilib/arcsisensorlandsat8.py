"""
Module that contains the ARCSILandsat8Sensor class.
"""
############################################################################
#  arcsisensorlandsat.py
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
# Purpose:  A class for read the landsat sensor header file and applying
#           the pre-processing operations within ARCSI to the landsat 8
#           datasets.
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 05/07/2013
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
# Import the RSGISLib Module.
import rsgislib
# Import the RSGISLib Image Calibration Module.
import rsgislib.imagecalibration
# Import the RSGISLib Image Calculations Module.
import rsgislib.imagecalc
# Import the RSGISLib segmentation Module
import rsgislib.segmentation
# Import the RSGISLib Raster GIS Module
import rsgislib.rastergis
# Import the collections module
import collections
# Import the py6s module for running 6S from python.
import Py6S
# Import the python maths library
import math
# Import the RIOS RAT library
from rios import rat
# Import the GDAL python library
import osgeo.gdal as gdal
# Import the scipy optimisation library - used for finding AOD values form the imagery.
from scipy.optimize import minimize
import numpy

class ARCSILandsat8Sensor (ARCSIAbstractSensor):
    """
    A class which represents the landsat 8 sensor to read
    header parameters and apply data processing operations.
    """
    def __init__(self):
        ARCSIAbstractSensor.__init__(self)
        self.sensor = "LS8"
        self.band1File = ""
        self.band2File = ""
        self.band3File = ""
        self.band4File = ""
        self.band5File = ""
        self.band6File = ""
        self.band7File = ""
        self.band8File = ""
        self.band9File = ""
        self.band10File = ""
        self.band11File = ""
        self.bandQAFile = ""
        self.row = 0
        self.path = 0
        
        self.b1RadMulti = 0
        self.b1CalMax = 0
        self.b2RadMulti = 0
        self.b2CalMax = 0
        self.b3RadMulti = 0
        self.b3CalMax = 0
        self.b4RadMulti = 0
        self.b4CalMax = 0
        self.b5RadMulti = 0
        self.b5CalMax = 0
        self.b6aRadMulti = 0
        self.b6aCalMax = 0
        self.b6bRadMulti = 0
        self.b6bCalMax = 0
        self.b7RadMulti = 0
        self.b7CalMax = 0
        self.b8RadMulti = 0
        self.b8CalMax = 0
        
        self.b1RadAdd = 0.0
        self.b1MaxRad = 0.0
        self.b2RadAdd = 0.0
        self.b2MaxRad = 0.0
        self.b3RadAdd = 0.0
        self.b3MaxRad = 0.0
        self.b4RadAdd = 0.0
        self.b4MaxRad = 0.0
        self.b5RadAdd = 0.0
        self.b5MaxRad = 0.0
        self.b6aRadAdd = 0.0
        self.b6aMaxRad = 0.0
        self.b6bRadAdd = 0.0
        self.b6bMaxRad = 0.0
        self.b7RadAdd = 0.0
        self.b7MaxRad = 0.0
        self.b8RadAdd = 0.0
        self.b8MaxRad = 0.0
    
    def extractHeaderParameters(self, inputHeader, wktStr):
        """
        Understands and parses the Landsat MTL header files
        """
        try:
            print("Reading header file")
            hFile = open(inputHeader, 'r')
            headerParams = dict()
            for line in hFile:
                line = line.strip()
                if line:
                    lineVals = line.split('=')
                    if len(lineVals) == 2:
                        if (lineVals[0].strip() != "GROUP") or (lineVals[0].strip() != "END_GROUP"):
                            headerParams[lineVals[0].strip()] = lineVals[1].strip().replace('"','')
            hFile.close()
            print("Extracting Header Values")
            # Get the sensor info.
            if (headerParams["SPACECRAFT_ID"] == "LANDSAT_8") and (headerParams["SENSOR_ID"] == "OLI_TIRS"):
                self.sensor = "LS8"
            else:
                raise ARCSIException("Do no recognise the spacecraft and sensor or combination.")
            
            # Get row/path
            self.row = int(headerParams["WRS_ROW"])
            self.path = int(headerParams["WRS_PATH"])
            
            # Get date and time of the acquisition
            acData = headerParams["DATE_ACQUIRED"].split('-')
            acTime = headerParams["SCENE_CENTER_TIME"].split(':')
            secsTime = acTime[2].split('.')
            self.acquisitionTime = datetime.datetime(int(acData[0]), int(acData[1]), int(acData[2]), int(acTime[0]), int(acTime[1]), int(secsTime[0]))
            
            self.solarZenith = 90-float(headerParams["SUN_ELEVATION"])
            self.solarAzimuth = float(headerParams["SUN_AZIMUTH"])
            
            # Get the geographic lat/long corners of the image.
            self.latTL = float(headerParams["CORNER_UL_LAT_PRODUCT"])
            self.lonTL = float(headerParams["CORNER_UL_LON_PRODUCT"])
            self.latTR = float(headerParams["CORNER_UR_LAT_PRODUCT"])
            self.lonTR = float(headerParams["CORNER_UR_LON_PRODUCT"])
            self.latBL = float(headerParams["CORNER_LL_LAT_PRODUCT"])
            self.lonBL = float(headerParams["CORNER_LL_LON_PRODUCT"])
            self.latBR = float(headerParams["CORNER_LR_LAT_PRODUCT"])
            self.lonBR = float(headerParams["CORNER_LR_LON_PRODUCT"])
            
            # Get the projected X/Y corners of the image
            self.xTL = float(headerParams["CORNER_UL_PROJECTION_X_PRODUCT"])
            self.yTL = float(headerParams["CORNER_UL_PROJECTION_Y_PRODUCT"])
            self.xTR = float(headerParams["CORNER_UR_PROJECTION_X_PRODUCT"])
            self.yTR = float(headerParams["CORNER_UR_PROJECTION_Y_PRODUCT"])
            self.xBL = float(headerParams["CORNER_LL_PROJECTION_X_PRODUCT"])
            self.yBL = float(headerParams["CORNER_LL_PROJECTION_Y_PRODUCT"])
            self.xBR = float(headerParams["CORNER_LR_PROJECTION_X_PRODUCT"])
            self.yBR = float(headerParams["CORNER_LR_PROJECTION_Y_PRODUCT"])
            
            # Get projection
            inProj = osr.SpatialReference()
            if (headerParams["MAP_PROJECTION"] == "UTM") and (headerParams["DATUM"] == "WGS84") and (headerParams["ELLIPSOID"] == "WGS84"):
                utmZone = int(headerParams["UTM_ZONE"])
                utmCode = "WGS84UTM" + str(utmZone) + str("N")
                #print("UTM: ", utmCode)
                inProj.ImportFromEPSG(self.epsgCodes[utmCode])
            else:
                raise ARCSIException("Expecting Landsat to be projected in UTM with datum=WGS84 and ellipsoid=WGS84.")
            
            # Check image is square!
            if not ((self.xTL == self.xBL) and (self.yTL == self.yTR) and (self.xTR == self.xBR) and (self.yBL == self.yBR)):
                raise ARCSIException("Image is not square in projected coordinates.")
            
            self.xCentre = self.xTL + ((self.xTR - self.xTL)/2)
            self.yCentre = self.yBR + ((self.yTL - self.yBR)/2)
            
            wgs84latlonProj = osr.SpatialReference()
            wgs84latlonProj.ImportFromEPSG(4326)
            
            wktPt = 'POINT(%s %s)' % (self.xCentre, self.yCentre)
            #print(wktPt)
            point = ogr.CreateGeometryFromWkt(wktPt)
            point.AssignSpatialReference(inProj)
            point.TransformTo(wgs84latlonProj)
            #print(point)
            
            self.latCentre = point.GetY()
            self.lonCentre = point.GetX()
            
            #print("Lat: " + str(self.latCentre) + " Long: " + str(self.lonCentre))
            
            filesDIR = os.path.dirname(inputHeader)
            
            self.band1File = os.path.join(filesDIR, headerParams["FILE_NAME_BAND_1"])
            self.band2File = os.path.join(filesDIR, headerParams["FILE_NAME_BAND_2"])
            self.band3File = os.path.join(filesDIR, headerParams["FILE_NAME_BAND_3"])
            self.band4File = os.path.join(filesDIR, headerParams["FILE_NAME_BAND_4"])
            self.band5File = os.path.join(filesDIR, headerParams["FILE_NAME_BAND_5"])
            self.band6File = os.path.join(filesDIR, headerParams["FILE_NAME_BAND_6"])
            self.band7File = os.path.join(filesDIR, headerParams["FILE_NAME_BAND_7"])
            self.band8File = os.path.join(filesDIR, headerParams["FILE_NAME_BAND_8"])
            self.band9File = os.path.join(filesDIR, headerParams["FILE_NAME_BAND_9"])
            self.band10File = os.path.join(filesDIR, headerParams["FILE_NAME_BAND_10"])
            self.band11File = os.path.join(filesDIR, headerParams["FILE_NAME_BAND_11"])
            self.bandQAFile = os.path.join(filesDIR, headerParams["FILE_NAME_BAND_QUALITY"])
            
            self.b1RadMulti = float(headerParams["RADIANCE_MULT_BAND_1"])
            self.b2RadMulti = float(headerParams["RADIANCE_MULT_BAND_2"])
            self.b3RadMulti = float(headerParams["RADIANCE_MULT_BAND_3"])
            self.b4RadMulti = float(headerParams["RADIANCE_MULT_BAND_4"])
            self.b5RadMulti = float(headerParams["RADIANCE_MULT_BAND_5"])
            self.b6RadMulti = float(headerParams["RADIANCE_MULT_BAND_6"])
            self.b7RadMulti = float(headerParams["RADIANCE_MULT_BAND_7"])
            self.b8RadMulti = float(headerParams["RADIANCE_MULT_BAND_8"])
            self.b9RadMulti = float(headerParams["RADIANCE_MULT_BAND_9"])
            self.b10RadMulti = float(headerParams["RADIANCE_MULT_BAND_10"])
            self.b11RadMulti = float(headerParams["RADIANCE_MULT_BAND_11"])
            
            self.b1RadAdd = float(headerParams["RADIANCE_ADD_BAND_1"])
            self.b2RadAdd = float(headerParams["RADIANCE_ADD_BAND_2"])
            self.b3RadAdd = float(headerParams["RADIANCE_ADD_BAND_3"])
            self.b4RadAdd = float(headerParams["RADIANCE_ADD_BAND_4"])
            self.b5RadAdd = float(headerParams["RADIANCE_ADD_BAND_5"])
            self.b6RadAdd = float(headerParams["RADIANCE_ADD_BAND_6"])
            self.b7RadAdd = float(headerParams["RADIANCE_ADD_BAND_7"])
            self.b8RadAdd = float(headerParams["RADIANCE_ADD_BAND_8"])
            self.b9RadAdd = float(headerParams["RADIANCE_ADD_BAND_9"])
            self.b10RadAdd = float(headerParams["RADIANCE_ADD_BAND_10"])
            self.b11RadAdd = float(headerParams["RADIANCE_ADD_BAND_11"])
            
        except Exception as e:
            raise e
        
    def generateOutputBaseName(self):
        """
        Provides an implementation for the landsat sensor
        """
        rowpath = "r" + str(self.row) + "p" + str(self.path)
        outname = self.defaultGenBaseOutFileName()
        outname = outname + str("_") + rowpath
        return outname
        
    def convertImageToRadiance(self, outputPath, outputName, outFormat):
        print("Converting to Radiance")
        outputImage = os.path.join(outputPath, outputName)
        bandDefnSeq = list()
        lsBand = collections.namedtuple('LSBand', ['bandName', 'fileName', 'bandIndex', 'addVal', 'multiVal'])
        bandDefnSeq.append(lsBand(bandName="Coastal", fileName=self.band1File, bandIndex=1, addVal=self.b1RadAdd, multiVal=self.b1RadMulti))
        bandDefnSeq.append(lsBand(bandName="Blue", fileName=self.band2File, bandIndex=1, addVal=self.b2RadAdd, multiVal=self.b2RadMulti))
        bandDefnSeq.append(lsBand(bandName="Green", fileName=self.band3File, bandIndex=1, addVal=self.b3RadAdd, multiVal=self.b3RadMulti))
        bandDefnSeq.append(lsBand(bandName="Red", fileName=self.band4File, bandIndex=1, addVal=self.b4RadAdd, multiVal=self.b4RadMulti))
        bandDefnSeq.append(lsBand(bandName="NIR", fileName=self.band5File, bandIndex=1, addVal=self.b5RadAdd, multiVal=self.b5RadMulti))
        bandDefnSeq.append(lsBand(bandName="SWIR1", fileName=self.band6File, bandIndex=1, addVal=self.b6RadAdd, multiVal=self.b6RadMulti))
        bandDefnSeq.append(lsBand(bandName="SWIR2", fileName=self.band7File, bandIndex=1, addVal=self.b7RadAdd, multiVal=self.b7RadMulti))
        rsgislib.imagecalibration.landsat2RadianceMultiAdd(outputImage, outFormat, bandDefnSeq)
        return outputImage
    
    def convertImageToTOARefl(self, inputRadImage, outputPath, outputName, outFormat):
        print("Converting to TOA")
        outputImage = os.path.join(outputPath, outputName)
        solarIrradianceVals = list()
        IrrVal = collections.namedtuple('SolarIrradiance', ['irradiance'])
        solarIrradianceVals.append(IrrVal(irradiance=1876.61))
        solarIrradianceVals.append(IrrVal(irradiance=1970.03))
        solarIrradianceVals.append(IrrVal(irradiance=1848.9))
        solarIrradianceVals.append(IrrVal(irradiance=1571.3))
        solarIrradianceVals.append(IrrVal(irradiance=967.66))
        solarIrradianceVals.append(IrrVal(irradiance=245.73))
        solarIrradianceVals.append(IrrVal(irradiance=82.03))
        rsgislib.imagecalibration.radiance2TOARefl(inputRadImage, outputImage, outFormat, rsgislib.TYPE_16UINT, 1000, self.acquisitionTime.year, self.acquisitionTime.month, self.acquisitionTime.day, self.solarZenith, solarIrradianceVals)
        return outputImage
        
    def convertImageToSurfaceReflSglParam(self, inputRadImage, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF):
        print("Converting to Surface Reflectance")
        outputImage = os.path.join(outputPath, outputName)
        
        Band6S = collections.namedtuple('Band6SCoeff', ['band', 'aX', 'bX', 'cX'])
        imgBandCoeffs = list()
        # Set up 6S model
        s = Py6S.SixS()
        s.atmos_profile = atmosProfile
        s.aero_profile = aeroProfile
        #s.ground_reflectance = Py6S.GroundReflectance.HomogeneousHapke(0.101, -0.263, 0.589, 0.046)
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
        if useBRDF:
            s.atmos_corr = Py6S.AtmosCorr.AtmosCorrBRDFFromRadiance(200)
        else:
            s.atmos_corr = Py6S.AtmosCorr.AtmosCorrLambertianFromRadiance(200)
        s.aot550 = aotVal
        
        # Band 1
        s.wavelength = Py6S.Wavelength(0.427, 0.4595, [0.000073, 0.001628, 0.024767, 0.254149, 0.908749, 0.977393, 0.986713, 0.993137, 0.982780, 0.905808, 0.226412, 0.036603, 0.002414, 0.000255])
        s.run()
        imgBandCoeffs.append(Band6S(band=1, aX=s.outputs.values['coef_xa'], bX=s.outputs.values['coef_xb'], cX=s.outputs.values['coef_xc']))
        
        # Band 2
        s.wavelength = Py6S.Wavelength(0.436, 0.5285, [0.000010, 0.000117, 0.000455, 0.001197, 0.006869, 0.027170, 0.271370, 0.723971, 0.903034, 0.909880, 0.889667, 0.877453, 0.879688, 0.891913, 0.848533, 0.828339, 0.868497, 0.912538, 0.931726, 0.954248, 0.956424, 0.978564, 0.989469, 0.968801, 0.988729, 0.967361, 0.966125, 0.981834, 0.963135, 0.996498, 0.844893, 0.190738, 0.005328, 0.001557, 0.000516, 0.000162, 0.000023, -0.000016])
        s.run()
        imgBandCoeffs.append(Band6S(band=2, aX=s.outputs.values['coef_xa'], bX=s.outputs.values['coef_xb'], cX=s.outputs.values['coef_xc']))
        
        # Band 3
        s.wavelength = Py6S.Wavelength(0.512, 0.6095, [-0.000046, 0.00011, 0.000648, 0.001332, 0.003446, 0.007024, 0.025513, 0.070551, 0.353885, 0.741205, 0.954627, 0.959215, 0.969873, 0.961397, 0.977001, 0.990784, 0.982642, 0.977765, 0.946245, 0.959038, 0.966447, 0.958314, 0.983397, 0.974522, 0.978208, 0.974392, 0.969181, 0.982956, 0.968886, 0.986657, 0.904478, 0.684974, 0.190467, 0.035393, 0.002574, 0.000394, -0.000194, -0.000292, -0.000348, -0.000351])
        s.run()
        imgBandCoeffs.append(Band6S(band=3, aX=s.outputs.values['coef_xa'], bX=s.outputs.values['coef_xb'], cX=s.outputs.values['coef_xc']))
        
        # Band 4
        s.wavelength = Py6S.Wavelength(0.625, 0.690, [-0.000342, 0.000895, 0.007197, 0.030432, 0.299778, 0.764443, 0.950823, 0.951831, 0.984173, 0.983434, 0.959441, 0.955548, 0.981688, 0.992388, 0.97696, 0.98108, 0.980678, 0.962154, 0.966928, 0.848855, 0.123946, 0.017702, 0.001402, 0.000117, -0.000376, -0.000458, -0.000429])
        s.run()
        imgBandCoeffs.append(Band6S(band=4, aX=s.outputs.values['coef_xa'], bX=s.outputs.values['coef_xb'], cX=s.outputs.values['coef_xc']))
        
        # Band 5
        s.wavelength = Py6S.Wavelength(0.829, 0.899, [-0.000034, 0.000050, 0.000314, 0.000719, 0.002107, 0.004744, 0.017346, 0.048191, 0.249733, 0.582623, 0.960215, 0.973133, 1.000000, 0.980733, 0.957357, 0.947044, 0.948450, 0.950632, 0.969821, 0.891066, 0.448364, 0.174619, 0.034532, 0.012440, 0.002944, 0.001192, 0.000241, 0.000044, -0.000084])
        s.run()
        imgBandCoeffs.append(Band6S(band=5, aX=s.outputs.values['coef_xa'], bX=s.outputs.values['coef_xb'], cX=s.outputs.values['coef_xc']))
        
        # Band 6
        s.wavelength = Py6S.Wavelength(1.515, 1.6975, [-0.00002, 0.00015, 0.00047, 0.00076, 0.00137, 0.00186, 0.00288, 0.00377, 0.00553, 0.00732, 0.01099, 0.01430, 0.02183, 0.02995, 0.04786, 0.06573, 0.10189, 0.13864, 0.22026, 0.29136, 0.42147, 0.52568, 0.67668, 0.75477, 0.85407, 0.89183, 0.91301, 0.92295, 0.92641, 0.92368, 0.92283, 0.92206, 0.92661, 0.94253, 0.94618, 0.94701, 0.95286, 0.94967, 0.95905, 0.96005, 0.96147, 0.96018, 0.96470, 0.96931, 0.97691, 0.98126, 0.98861, 0.99802, 0.99964, 0.99344, 0.96713, 0.93620, 0.84097, 0.75189, 0.57323, 0.45197, 0.29175, 0.21115, 0.12846, 0.09074, 0.05275, 0.03731, 0.02250, 0.01605, 0.00959, 0.00688, 0.00426, 0.00306, 0.00178, 0.00124, 0.00068, 0.00041, 0.00011, -0.00003])
        s.run()
        imgBandCoeffs.append(Band6S(band=6, aX=s.outputs.values['coef_xa'], bX=s.outputs.values['coef_xb'], cX=s.outputs.values['coef_xc']))
        
         # Band 7
        s.wavelength = Py6S.Wavelength(2.037, 2.3545, [-0.000010, 0.000083, 0.000240, 0.000368, 0.000599, 0.000814, 0.001222, 0.001546, 0.002187, 0.002696, 0.003733, 0.004627, 0.006337, 0.007996, 0.011005, 0.013610, 0.018899, 0.023121, 0.032071, 0.040206, 0.056429, 0.070409, 0.100640, 0.128292, 0.179714, 0.227234, 0.311347, 0.377044, 0.488816, 0.554715, 0.663067, 0.722284, 0.792667, 0.836001, 0.867845, 0.886411, 0.906527, 0.911091, 0.929693, 0.936544, 0.942952, 0.943194, 0.948776, 0.949643, 0.956635, 0.947423, 0.950874, 0.947014, 0.957717, 0.946412, 0.951641, 0.948644, 0.940311, 0.947923, 0.938737, 0.941859, 0.944482, 0.951661, 0.939939, 0.935493, 0.938955, 0.929162, 0.930508, 0.933908, 0.936472, 0.933523, 0.946217, 0.955661, 0.963135, 0.964365, 0.962905, 0.962473, 0.957814, 0.958041, 0.951706, 0.960212, 0.947696, 0.959060, 0.955750, 0.953245, 0.966786, 0.960173, 0.977637, 0.982760, 0.985056, 0.999600, 0.992469, 0.995894, 0.997261, 0.991127, 0.986037, 0.984536, 0.972794, 0.976540, 0.974409, 0.967502, 0.955095, 0.955588, 0.922405, 0.894940, 0.823876, 0.744025, 0.602539, 0.502693, 0.355569, 0.278260, 0.186151, 0.141435, 0.092029, 0.069276, 0.046332, 0.035634, 0.024000, 0.018688, 0.012930, 0.010155, 0.007088, 0.005643, 0.003903, 0.003025, 0.002047, 0.001554, 0.000974, 0.000680, 0.000320, 0.000119, -0.000134, -0.000263])
        s.run()
        imgBandCoeffs.append(Band6S(band=7, aX=s.outputs.values['coef_xa'], bX=s.outputs.values['coef_xb'], cX=s.outputs.values['coef_xc']))
        
        for band in imgBandCoeffs:
            print(band)
        
        rsgislib.imagecalibration.apply6SCoeffSingleParam(inputRadImage, outputImage, outFormat, rsgislib.TYPE_16UINT, 1000, 0, True, imgBandCoeffs)
        
        return outputImage

    def run6SToOptimiseAODValue(self, aotVal, minB2Val, minB3Val, minB4Val, predB2Val, predB3Val, predB4Val, aeroProfile, atmosProfile, grdRefl, surfaceAltitude):
        """Used as part of the optimastion for identifying values of AOD"""
        print "Testing AOD Val: ", aotVal
        s = Py6S.SixS()
        s.atmos_profile = Py6S.AtmosProfile.PredefinedType(Py6S.AtmosProfile.MidlatitudeSummer)
        s.aero_profile = Py6S.AeroProfile.PredefinedType(Py6S.AeroProfile.Maritime)
        s.ground_reflectance = Py6S.GroundReflectance.HomogeneousLambertian(Py6S.GroundReflectance.GreenVegetation)
        s.geometry = Py6S.Geometry.Landsat_TM()
        s.geometry.month = self.acquisitionTime.month
        s.geometry.day = self.acquisitionTime.day
        s.geometry.gmt_decimal_hour = float(self.acquisitionTime.hour) + float(self.acquisitionTime.minute)/60.0
        s.geometry.latitude = self.latCentre
        s.geometry.longitude = self.lonCentre
        s.altitudes = Py6S.Altitudes()
        s.altitudes.set_target_custom_altitude(200)
        s.altitudes.set_sensor_satellite_level()
        s.atmos_corr = Py6S.AtmosCorr.AtmosCorrLambertianFromReflectance(0.40)
        s.aot550 = aotVal
        
        
        # Band 2
        s.wavelength = Py6S.Wavelength(0.436, 0.5285, [0.000010, 0.000117, 0.000455, 0.001197, 0.006869, 0.027170, 0.271370, 0.723971, 0.903034, 0.909880, 0.889667, 0.877453, 0.879688, 0.891913, 0.848533, 0.828339, 0.868497, 0.912538, 0.931726, 0.954248, 0.956424, 0.978564, 0.989469, 0.968801, 0.988729, 0.967361, 0.966125, 0.981834, 0.963135, 0.996498, 0.844893, 0.190738, 0.005328, 0.001557, 0.000516, 0.000162, 0.000023, -0.000016])
        s.run()
        aX = float(s.outputs.values['coef_xa'])
        bX = float(s.outputs.values['coef_xb'])
        cX = float(s.outputs.values['coef_xc'])
        #print "B2 aX: ", aX
        #print "B2 bX: ", bX
        #print "B2 cX: ", cX
        tmpVal = (aX*minB2Val)-bX;
        b2Val = tmpVal/(1.0+cX*tmpVal)
        
        # Band 3
        s.wavelength = Py6S.Wavelength(0.512, 0.6095, [-0.000046, 0.00011, 0.000648, 0.001332, 0.003446, 0.007024, 0.025513, 0.070551, 0.353885, 0.741205, 0.954627, 0.959215, 0.969873, 0.961397, 0.977001, 0.990784, 0.982642, 0.977765, 0.946245, 0.959038, 0.966447, 0.958314, 0.983397, 0.974522, 0.978208, 0.974392, 0.969181, 0.982956, 0.968886, 0.986657, 0.904478, 0.684974, 0.190467, 0.035393, 0.002574, 0.000394, -0.000194, -0.000292, -0.000348, -0.000351])
        s.run()
        aX = float(s.outputs.values['coef_xa'])
        bX = float(s.outputs.values['coef_xb'])
        cX = float(s.outputs.values['coef_xc'])
        #print "B3 aX: ", aX
        #print "B3 bX: ", bX
        #print "B3 cX: ", cX
        tmpVal = (aX*minB3Val)-bX;
        b3Val = tmpVal/(1.0+cX*tmpVal)
        
        # Band 4
        s.wavelength = Py6S.Wavelength(0.625, 0.690, [-0.000342, 0.000895, 0.007197, 0.030432, 0.299778, 0.764443, 0.950823, 0.951831, 0.984173, 0.983434, 0.959441, 0.955548, 0.981688, 0.992388, 0.97696, 0.98108, 0.980678, 0.962154, 0.966928, 0.848855, 0.123946, 0.017702, 0.001402, 0.000117, -0.000376, -0.000458, -0.000429])
        s.run()
        aX = float(s.outputs.values['coef_xa'])
        bX = float(s.outputs.values['coef_xb'])
        cX = float(s.outputs.values['coef_xc'])
        #print "B4 aX: ", aX
        #print "B4 bX: ", bX
        #print "B4 cX: ", cX
        tmpVal = (aX*minB4Val)-bX;
        b4Val = tmpVal/(1.0+cX*tmpVal)
        
        outDist = math.sqrt(math.pow((b2Val - predB2Val),2) + math.pow((b3Val - predB3Val),2) + math.pow((b4Val - predB4Val),2))
        print "Dist ", outDist
        return outDist
    
    def estimateImageToAOD(self, inputTOAImage, outputPath, outputName, outFormat, tmpPath, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotValMin, aotValMax):
        try:
            arcsiUtils = ARCSIUtils()
            tmpBaseName = os.path.splitext(outputName)[0]
            thresImage = os.path.join(tmpPath, tmpBaseName+"_thresd"+arcsiUtils.getFileExtension(outFormat))
            thresImageClumps = os.path.join(tmpPath, tmpBaseName+"_thresdclumps"+arcsiUtils.getFileExtension(outFormat))
            thresImageClumpsRMSmall = os.path.join(tmpPath, tmpBaseName+"_thresdclumpsgt10"+arcsiUtils.getFileExtension(outFormat))
            thresImageClumpsFinal = os.path.join(tmpPath, tmpBaseName+"_thresdclumpsFinal"+arcsiUtils.getFileExtension(outFormat))
            
            thresMathBands = list()
            thresMathBands.append(rsgislib.imagecalc.BandDefn(bandName='b2', fileName=inputTOAImage, bandIndex=2))
            thresMathBands.append(rsgislib.imagecalc.BandDefn(bandName='b3', fileName=inputTOAImage, bandIndex=3))
            thresMathBands.append(rsgislib.imagecalc.BandDefn(bandName='b4', fileName=inputTOAImage, bandIndex=4))
            thresMathBands.append(rsgislib.imagecalc.BandDefn(bandName='b5', fileName=inputTOAImage, bandIndex=5))
            thresMathBands.append(rsgislib.imagecalc.BandDefn(bandName='b6', fileName=inputTOAImage, bandIndex=5))
            thresMathBands.append(rsgislib.imagecalc.BandDefn(bandName='b7', fileName=inputTOAImage, bandIndex=7))
            #rsgislib.imagecalc.bandMath(thresImage, "((((b2+b3+b4+b5+b6+b7)/6)<100)&&((b7>15)&&(b7<50))&&(((b5-b4)/(b5+b4))>0.2))?1:0", outFormat, rsgislib.TYPE_8UINT, thresMathBands)
            #rsgislib.segmentation.clump(thresImage, thresImageClumps, outFormat, False, 0.0)
            #rsgislib.rastergis.populateStats(thresImageClumps, True, True)
            rsgislib.segmentation.rmSmallClumps(thresImageClumps, thresImageClumpsRMSmall, 10, outFormat)
            rsgislib.segmentation.relabelClumps(thresImageClumpsRMSmall, thresImageClumpsFinal, outFormat, False)
            rsgislib.rastergis.populateStats(thresImageClumpsFinal, True, True)
            stats2Calc = list()
            stats2Calc.append(rsgislib.rastergis.BandAttStats(band=1, minField="MinB1", meanField="MeanB1"))
            stats2Calc.append(rsgislib.rastergis.BandAttStats(band=2, minField="MinB2", meanField="MeanB2"))
            stats2Calc.append(rsgislib.rastergis.BandAttStats(band=3, minField="MinB3", meanField="MeanB3"))
            stats2Calc.append(rsgislib.rastergis.BandAttStats(band=4, minField="MinB4", meanField="MeanB4"))
            stats2Calc.append(rsgislib.rastergis.BandAttStats(band=5, minField="MinB5", meanField="MeanB5"))
            stats2Calc.append(rsgislib.rastergis.BandAttStats(band=6, minField="MinB6", meanField="MeanB6"))
            stats2Calc.append(rsgislib.rastergis.BandAttStats(band=7, minField="MinB7", meanField="MeanB7"))
            rsgislib.rastergis.populateRATWithStats(inputTOAImage, thresImageClumpsFinal, stats2Calc)
            ratDS = gdal.Open(thresImageClumpsFinal, gdal.GA_Update)
            MinB2 = rat.readColumn(ratDS, "MinB2")
            MinB3 = rat.readColumn(ratDS, "MinB3")
            MinB4 = rat.readColumn(ratDS, "MinB4")
            MinB7 = rat.readColumn(ratDS, "MinB7")
            
            PredB2 = MinB7/4.3
            PredB3 = MinB7/1.5
            PredB4 = MinB7/2.0
            
            rat.writeColumn(ratDS, "PredB2", PredB2)
            rat.writeColumn(ratDS, "PredB3", PredB3)
            rat.writeColumn(ratDS, "PredB4", PredB4)
            
            numAOTValTests = int(math.ceil((aotValMax - aotValMin)/0.05))
            
            if not numAOTValTests >= 1:
                raise ARCSIException("min and max AOT range are too close together, they need to be at least 0.05 apart.")
            
            cAOT = aotValMin
            cDist = 0.0
            minAOT = 0.0
            minDist = 0.0
            
            for i in range(len(PredB2)):
                if i != 0:
                    print "Predicting AOD for Segment ", i
                    predAOTArgs = list()
                    predAOTArgs.append(MinB2[i])
                    predAOTArgs.append(MinB3[i])
                    predAOTArgs.append(MinB4[i])
                    predAOTArgs.append(PredB2[i])
                    predAOTArgs.append(PredB3[i])
                    predAOTArgs.append(PredB4[i])
                    predAOTArgs.append(aeroProfile)
                    predAOTArgs.append(atmosProfile)
                    predAOTArgs.append(grdRefl)
                    predAOTArgs.append(surfaceAltitude)
                    for j in range(numAOTValTests):
                        cAOT = aotValMin + (0.05 * j)
                        cDist = self.run6SToOptimiseAODValue(cAOT, MinB2[i], MinB3[i], MinB4[i], PredB2[i], PredB3[i], PredB4[i], aeroProfile, atmosProfile, grdRefl, surfaceAltitude)
                        if j == 0:
                            minAOT = cAOT
                            minDist = cDist
                        elif cDist < minDist:
                            minAOT = cAOT
                            minDist = cDist
                    
                    res = minimize(self.run6SToOptimiseAODValue, minAOT, method='nelder-mead', options={'maxiter': 20, 'xtol': 0.001, 'disp': True}, args=predAOTArgs)
                    print res
        except Exception as e:
            raise e

