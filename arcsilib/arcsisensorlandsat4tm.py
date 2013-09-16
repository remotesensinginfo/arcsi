"""
Module that contains the ARCSILandsat4TMSensor class.
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
#           the pre-processing operations within ARCSI to the landsat 4 TM
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
# Import the RSGISLib Image Calibration Module.
import rsgislib.imagecalibration
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

class ARCSILandsat4TMSensor (ARCSIAbstractSensor):
    """
    A class which represents the landsat 4 TM sensor to read
    header parameters and apply data processing operations.
    """
    def __init__(self):
        ARCSIAbstractSensor.__init__(self)
        self.sensor = "LS4TM"
        self.band1File = ""
        self.band2File = ""
        self.band3File = ""
        self.band4File = ""
        self.band5File = ""
        self.band6File = ""
        self.band7File = ""
        self.row = 0
        self.path = 0
        
        self.b1CalMin = 0
        self.b1CalMax = 0
        self.b2CalMin = 0
        self.b2CalMax = 0
        self.b3CalMin = 0
        self.b3CalMax = 0
        self.b4CalMin = 0
        self.b4CalMax = 0
        self.b5CalMin = 0
        self.b5CalMax = 0
        self.b6CalMin = 0
        self.b6CalMax = 0
        self.b7CalMin = 0
        self.b7CalMax = 0
        
        self.b1MinRad = 0.0
        self.b1MaxRad = 0.0
        self.b2MinRad = 0.0
        self.b2MaxRad = 0.0
        self.b3MinRad = 0.0
        self.b3MaxRad = 0.0
        self.b4MinRad = 0.0
        self.b4MaxRad = 0.0
        self.b5MinRad = 0.0
        self.b5MaxRad = 0.0
        self.b6MinRad = 0.0
        self.b6MaxRad = 0.0
        self.b7MinRad = 0.0
        self.b7MaxRad = 0.0
    
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
            if (headerParams["SPACECRAFT_ID"] == "LANDSAT_4") and (headerParams["SENSOR_ID"] == "TM"):
                self.sensor = "LS4TM"
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
            
            self.b1CalMin = float(headerParams["QUANTIZE_CAL_MIN_BAND_1"])
            self.b1CalMax = float(headerParams["QUANTIZE_CAL_MAX_BAND_1"])
            self.b2CalMin = float(headerParams["QUANTIZE_CAL_MIN_BAND_2"])
            self.b2CalMax = float(headerParams["QUANTIZE_CAL_MAX_BAND_2"])
            self.b3CalMin = float(headerParams["QUANTIZE_CAL_MIN_BAND_3"])
            self.b3CalMax = float(headerParams["QUANTIZE_CAL_MAX_BAND_3"])
            self.b4CalMin = float(headerParams["QUANTIZE_CAL_MIN_BAND_4"])
            self.b4CalMax = float(headerParams["QUANTIZE_CAL_MAX_BAND_4"])
            self.b5CalMin = float(headerParams["QUANTIZE_CAL_MIN_BAND_5"])
            self.b5CalMax = float(headerParams["QUANTIZE_CAL_MAX_BAND_5"])
            self.b6CalMin = float(headerParams["QUANTIZE_CAL_MIN_BAND_6"])
            self.b6CalMax = float(headerParams["QUANTIZE_CAL_MAX_BAND_6"])
            self.b7CalMin = float(headerParams["QUANTIZE_CAL_MIN_BAND_7"])
            self.b7CalMax = float(headerParams["QUANTIZE_CAL_MAX_BAND_7"])
            
            self.b1MinRad = float(headerParams["RADIANCE_MINIMUM_BAND_1"])
            self.b1MaxRad = float(headerParams["RADIANCE_MAXIMUM_BAND_1"])
            self.b2MinRad = float(headerParams["RADIANCE_MINIMUM_BAND_2"])
            self.b2MaxRad = float(headerParams["RADIANCE_MAXIMUM_BAND_2"])
            self.b3MinRad = float(headerParams["RADIANCE_MINIMUM_BAND_3"])
            self.b3MaxRad = float(headerParams["RADIANCE_MAXIMUM_BAND_3"])
            self.b4MinRad = float(headerParams["RADIANCE_MINIMUM_BAND_4"])
            self.b4MaxRad = float(headerParams["RADIANCE_MAXIMUM_BAND_4"])
            self.b5MinRad = float(headerParams["RADIANCE_MINIMUM_BAND_5"])
            self.b5MaxRad = float(headerParams["RADIANCE_MAXIMUM_BAND_5"])
            self.b6MinRad = float(headerParams["RADIANCE_MINIMUM_BAND_6"])
            self.b6MaxRad = float(headerParams["RADIANCE_MAXIMUM_BAND_6"])
            self.b7MinRad = float(headerParams["RADIANCE_MINIMUM_BAND_7"])
            self.b7MaxRad = float(headerParams["RADIANCE_MAXIMUM_BAND_7"])
            
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
        lsBand = collections.namedtuple('LSBand', ['bandName', 'fileName', 'bandIndex', 'lMin', 'lMax', 'qCalMin', 'qCalMax'])
        bandDefnSeq.append(lsBand(bandName="Blue", fileName=self.band1File, bandIndex=1, lMin=self.b1MinRad, lMax=self.b1MaxRad, qCalMin=self.b1CalMin, qCalMax=self.b1CalMax))
        bandDefnSeq.append(lsBand(bandName="Green", fileName=self.band2File, bandIndex=1, lMin=self.b2MinRad, lMax=self.b2MaxRad, qCalMin=self.b2CalMin, qCalMax=self.b2CalMax))
        bandDefnSeq.append(lsBand(bandName="Red", fileName=self.band3File, bandIndex=1, lMin=self.b3MinRad, lMax=self.b3MaxRad, qCalMin=self.b3CalMin, qCalMax=self.b3CalMax))
        bandDefnSeq.append(lsBand(bandName="NIR", fileName=self.band4File, bandIndex=1, lMin=self.b4MinRad, lMax=self.b4MaxRad, qCalMin=self.b4CalMin, qCalMax=self.b4CalMax))
        bandDefnSeq.append(lsBand(bandName="SWIR1", fileName=self.band5File, bandIndex=1, lMin=self.b5MinRad, lMax=self.b5MaxRad, qCalMin=self.b5CalMin, qCalMax=self.b5CalMax))
        bandDefnSeq.append(lsBand(bandName="SWIR2", fileName=self.band7File, bandIndex=1, lMin=self.b7MinRad, lMax=self.b7MaxRad, qCalMin=self.b7CalMin, qCalMax=self.b7CalMax))
        rsgislib.imagecalibration.landsat2Radiance(outputImage, outFormat, bandDefnSeq)
        return outputImage
    
    def convertImageToTOARefl(self, inputRadImage, outputPath, outputName, outFormat):
        print("Converting to TOA")
        outputImage = os.path.join(outputPath, outputName)
        solarIrradianceVals = list()
        IrrVal = collections.namedtuple('SolarIrradiance', ['irradiance'])
        solarIrradianceVals.append(IrrVal(irradiance=1957.0))
        solarIrradianceVals.append(IrrVal(irradiance=1825.0))
        solarIrradianceVals.append(IrrVal(irradiance=1557.0))
        solarIrradianceVals.append(IrrVal(irradiance=1033.0))
        solarIrradianceVals.append(IrrVal(irradiance=214.9))
        solarIrradianceVals.append(IrrVal(irradiance=80.72))
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
        s.wavelength = Py6S.Wavelength(Py6S.SixSHelpers.PredefinedWavelengths.LANDSAT_TM_B1)
        s.run()
        imgBandCoeffs.append(Band6S(band=1, aX=s.outputs.values['coef_xa'], bX=s.outputs.values['coef_xb'], cX=s.outputs.values['coef_xc']))
        
        # Band 2
        s.wavelength = Py6S.Wavelength(Py6S.SixSHelpers.PredefinedWavelengths.LANDSAT_TM_B2)
        s.run()
        imgBandCoeffs.append(Band6S(band=2, aX=s.outputs.values['coef_xa'], bX=s.outputs.values['coef_xb'], cX=s.outputs.values['coef_xc']))
        
        # Band 3
        s.wavelength = Py6S.Wavelength(Py6S.SixSHelpers.PredefinedWavelengths.LANDSAT_TM_B3)
        s.run()
        imgBandCoeffs.append(Band6S(band=3, aX=s.outputs.values['coef_xa'], bX=s.outputs.values['coef_xb'], cX=s.outputs.values['coef_xc']))
        
        # Band 4
        s.wavelength = Py6S.Wavelength(Py6S.SixSHelpers.PredefinedWavelengths.LANDSAT_TM_B4)
        s.run()
        imgBandCoeffs.append(Band6S(band=4, aX=s.outputs.values['coef_xa'], bX=s.outputs.values['coef_xb'], cX=s.outputs.values['coef_xc']))
        
        # Band 5
        s.wavelength = Py6S.Wavelength(Py6S.SixSHelpers.PredefinedWavelengths.LANDSAT_TM_B5)
        s.run()
        imgBandCoeffs.append(Band6S(band=5, aX=s.outputs.values['coef_xa'], bX=s.outputs.values['coef_xb'], cX=s.outputs.values['coef_xc']))
        
        # Band 7
        s.wavelength = Py6S.Wavelength(Py6S.SixSHelpers.PredefinedWavelengths.LANDSAT_TM_B7)
        s.run()
        imgBandCoeffs.append(Band6S(band=6, aX=s.outputs.values['coef_xa'], bX=s.outputs.values['coef_xb'], cX=s.outputs.values['coef_xc']))
        
        for band in imgBandCoeffs:
            print(band)
        
        rsgislib.imagecalibration.apply6SCoeffSingleParam(inputRadImage, outputImage, outFormat, rsgislib.TYPE_16UINT, 1000, 0, True, imgBandCoeffs)
        
        return outputImage



    def run6SToOptimiseAODValue(self, aotVal, radBlueVal, predBlueVal, aeroProfile, atmosProfile, grdRefl, surfaceAltitude):
        """Used as part of the optimastion for identifying values of AOD"""
        print "Testing AOD Val: ", aotVal
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
        s.atmos_corr = Py6S.AtmosCorr.AtmosCorrLambertianFromRadiance(200)
        s.aot550 = aotVal
                
        # Band 2
        s.wavelength = Py6S.Wavelength(Py6S.SixSHelpers.PredefinedWavelengths.LANDSAT_TM_B2)
        s.run()
        aX = float(s.outputs.values['coef_xa'])
        bX = float(s.outputs.values['coef_xb'])
        cX = float(s.outputs.values['coef_xc'])
        print "\taX: ", aX
        print "\tbX: ", bX
        print "\tcX: ", cX
        tmpVal = (aX*radBlueVal)-bX;
        reflBlueVal = tmpVal/(1.0+cX*tmpVal)
        
        
        outDist = math.sqrt(math.pow((reflBlueVal - predBlueVal),2))
        print "\tDist ", outDist
        return outDist
    
    def estimateImageToAOD(self, inputRADImage, inputTOAImage, outputPath, outputName, outFormat, tmpPath, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotValMin, aotValMax):
        try:
            arcsiUtils = ARCSIUtils()
            tmpBaseName = os.path.splitext(outputName)[0]
            thresImage = os.path.join(tmpPath, tmpBaseName+"_thresd"+arcsiUtils.getFileExtension(outFormat))
            thresImageClumps = os.path.join(tmpPath, tmpBaseName+"_thresdclumps"+arcsiUtils.getFileExtension(outFormat))
            thresImageClumpsRMSmall = os.path.join(tmpPath, tmpBaseName+"_thresdclumpsgt10"+arcsiUtils.getFileExtension(outFormat))
            thresImageClumpsFinal = os.path.join(tmpPath, tmpBaseName+"_thresdclumpsFinal"+arcsiUtils.getFileExtension(outFormat))
            
            outputAOTImage = os.path.join(outputPath, outputName)
            
            thresMathBands = list()
            thresMathBands.append(rsgislib.imagecalc.BandDefn(bandName='b1', fileName=inputTOAImage, bandIndex=1))
            thresMathBands.append(rsgislib.imagecalc.BandDefn(bandName='b2', fileName=inputTOAImage, bandIndex=2))
            thresMathBands.append(rsgislib.imagecalc.BandDefn(bandName='b3', fileName=inputTOAImage, bandIndex=3))
            thresMathBands.append(rsgislib.imagecalc.BandDefn(bandName='b4', fileName=inputTOAImage, bandIndex=4))
            thresMathBands.append(rsgislib.imagecalc.BandDefn(bandName='b5', fileName=inputTOAImage, bandIndex=5))
            thresMathBands.append(rsgislib.imagecalc.BandDefn(bandName='b7', fileName=inputTOAImage, bandIndex=6))
            rsgislib.imagecalc.bandMath(thresImage, "((((b1+b2+b3+b4+b5+b7)/6)<100)&&((b7>15)&&(b7<50))&&(((b4-b3)/(b4+b3))>0.2))?1:0", outFormat, rsgislib.TYPE_8UINT, thresMathBands)
            rsgislib.segmentation.clump(thresImage, thresImageClumps, outFormat, False, 0.0)
            rsgislib.rastergis.populateStats(thresImageClumps, True, True)
            rsgislib.segmentation.rmSmallClumps(thresImageClumps, thresImageClumpsRMSmall, 100, outFormat)
            rsgislib.segmentation.relabelClumps(thresImageClumpsRMSmall, thresImageClumpsFinal, outFormat, False)
            rsgislib.rastergis.populateStats(thresImageClumpsFinal, True, True)
            stats2CalcTOA = list()
            stats2CalcTOA.append(rsgislib.rastergis.BandAttStats(band=1, minField="MinB1TOA", meanField="MeanB1TOA"))
            stats2CalcTOA.append(rsgislib.rastergis.BandAttStats(band=6, minField="MinB7TOA", meanField="MeanB7TOA"))
            rsgislib.rastergis.populateRATWithStats(inputTOAImage, thresImageClumpsFinal, stats2CalcTOA)
            stats2CalcRad = list()
            stats2CalcRad.append(rsgislib.rastergis.BandAttStats(band=1, minField="MinB1RAD", meanField="MeanB2RAD"))
            rsgislib.rastergis.populateRATWithStats(inputRADImage, thresImageClumpsFinal, stats2CalcRad)

            ratDS = gdal.Open(thresImageClumpsFinal, gdal.GA_Update)
            Histogram = rat.readColumn(ratDS, "Histogram")
            selected = Histogram * 2
            selected[...] = 1
            selected[0] = 0
            rat.writeColumn(ratDS, "Selected", selected)
            
            rsgislib.rastergis.spatialLocation(thresImageClumpsFinal, "Eastings", "Northings")
            rsgislib.rastergis.selectClumpsOnGrid(thresImageClumpsFinal, "Selected", "PredictAOTFor", "Eastings", "Northings", "MinB7TOA", "min", 10, 10)
            
            MinB1TOA = rat.readColumn(ratDS, "MinB1TOA")
            MinB7TOA = rat.readColumn(ratDS, "MinB7TOA")
            MinB1RAD = rat.readColumn(ratDS, "MinB1RAD")
            PredictAOTFor = rat.readColumn(ratDS, "PredictAOTFor")
            
            PredB1Refl = (MinB7TOA/1000) * 0.33
            
            rat.writeColumn(ratDS, "PredB1Refl", PredB1Refl)
            
            numAOTValTests = int(math.ceil((aotValMax - aotValMin)/0.05))
            
            if not numAOTValTests >= 1:
                raise ARCSIException("min and max AOT range are too close together, they need to be at least 0.05 apart.")
            
            cAOT = aotValMin
            cDist = 0.0
            minAOT = 0.0
            minDist = 0.0
            
            aotVals = numpy.zeros_like(MinB7TOA, dtype=numpy.float)
            
            for i in range(len(PredB1Refl)):
                if PredictAOTFor[i] == 1:
                    print "Predicting AOD for Segment ", i
                    for j in range(numAOTValTests):
                        cAOT = aotValMin + (0.05 * j)
                        cDist = self.run6SToOptimiseAODValue(cAOT, MinB1RAD[i], PredB1Refl[i], aeroProfile, atmosProfile, grdRefl, surfaceAltitude)
                        if j == 0:
                            minAOT = cAOT
                            minDist = cDist
                        elif cDist < minDist:
                            minAOT = cAOT
                            minDist = cDist
                    predAOTArgs = list()
                    predAOTArgs.append(MinB1RAD[i])
                    predAOTArgs.append(PredB1Refl[i])
                    predAOTArgs.append(aeroProfile)
                    predAOTArgs.append(atmosProfile)
                    predAOTArgs.append(grdRefl)
                    predAOTArgs.append(surfaceAltitude)
                    res = minimize(self.run6SToOptimiseAODValue, minAOT, method='nelder-mead', options={'maxiter': 20, 'xtol': 0.001, 'disp': True}, args=predAOTArgs)
                    print "IDENTIFIED AOT: ", res.x[0]
                    aotVals[i] = res.x[0]
                else:
                    aotVals[i] = 0
            rat.writeColumn(ratDS, "AOT", aotVals)
            
            rsgislib.rastergis.interpolateClumpValues2Image(thresImageClumpsFinal, "PredictAOTFor", "Eastings", "Northings", "naturalnearestneighbour", "AOT", outputAOTImage, outFormat, rsgislib.TYPE_32FLOAT)
        
            gdalDriver = gdal.GetDriverByName(outFormat)
            gdalDriver.Delete(thresImage)
            gdalDriver.Delete(thresImageClumps)
            gdalDriver.Delete(thresImageClumpsRMSmall)
            gdalDriver.Delete(thresImageClumpsFinal)        
        
            return outputAOTImage
        except Exception as e:
            raise e


