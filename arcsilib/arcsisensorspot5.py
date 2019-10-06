"""
Module that contains the ARCSISPOT5Sensor class.
"""
############################################################################
#  arcsisensorspot5.py
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
# Purpose:  A class for read the SPOT5 sensor header file and applying
#           the pre-processing operations within ARCSI to the SPOT5
#           datasets.
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 07/10/2014
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

class ARCSISPOT5Sensor (ARCSIAbstractSensor):
    """
    A class which represents the RapidEye sensor to read
    header parameters and apply data processing operations.
    """
    def __init__(self, debugMode, inputImage):
        ARCSIAbstractSensor.__init__(self, debugMode, inputImage)
        self.sensor = "SPOT5"

        self.fileName = ""
        self.mission = ""
        self.missionIdx = ""
        self.instrument = ""
        self.instrumentIdx = ""
        self.acquIncidAngle = 0.0
        self.acquViewAngle = 0.0

        self.b1Bias = 0.0
        self.b2Bias = 0.0
        self.b3Bias = 0.0
        self.b4Bias = 0.0

        self.b1Gain = 0.0
        self.b2Gain = 0.0
        self.b3Gain = 0.0
        self.b4Gain = 0.0

        self.b1SolarIrradiance = 0.0
        self.b2SolarIrradiance = 0.0
        self.b3SolarIrradiance = 0.0
        self.b4SolarIrradiance = 0.0

        self.numXPxl = 0
        self.numYPxl = 0

        self.inImgHasGCPs = False

    def extractHeaderParameters(self, inputHeader, wktStr):
        """
        Understands and parses the RapidEye metadata.xml header file
        """
        try:
            self.headerFileName = os.path.split(inputHeader)[1]
            
            tree = ET.parse(inputHeader)
            root = tree.getroot()

            dimapVersion = root.find('Metadata_Id').find('METADATA_FORMAT').attrib['version']
            if dimapVersion != '1.1':
                print('DIMAP version is ' + dimapVersion)
                raise ARCSIException("Only DIMAP Version 1.1 is supported by this reader.")

            self.mission = root.find('Dataset_Sources').find('Source_Information').find('Scene_Source').find('MISSION').text.strip()
            self.missionIdx = root.find('Dataset_Sources').find('Source_Information').find('Scene_Source').find('MISSION_INDEX').text.strip()
            self.instrument = root.find('Dataset_Sources').find('Source_Information').find('Scene_Source').find('INSTRUMENT').text.strip()
            self.instrumentIdx = root.find('Dataset_Sources').find('Source_Information').find('Scene_Source').find('INSTRUMENT_INDEX').text.strip()

            if (self.mission != "SPOT") or (self.missionIdx != "5"):
                print("Mission = " + self.mission)
                print("Mission Index= " + self.missionIdx)
                raise ARCSIException("Only the SPOT 5 mission is supported")

            if (self.instrument != "HRG") or (self.instrumentIdx != "2"):
                print("Instrument = " + self.instrument)
                print("Instrument Index= " + self.instrumentIdx)
                raise ARCSIException("Only the HRG version 2 instrument is supported.")

            frameCorners = root.find('Dataset_Frame').findall('Vertex')
            count = 0
            for vertex in frameCorners:
                if count == 0:
                    self.latTL = float(vertex.find('FRAME_LAT').text.strip())
                    self.lonTL = float(vertex.find('FRAME_LON').text.strip())
                elif count == 1:
                    self.latTR = float(vertex.find('FRAME_LAT').text.strip())
                    self.lonTR = float(vertex.find('FRAME_LON').text.strip())
                elif count == 2:
                    self.latBR = float(vertex.find('FRAME_LAT').text.strip())
                    self.lonBR = float(vertex.find('FRAME_LON').text.strip())
                    self.numXPxl = int(vertex.find('FRAME_COL').text.strip())
                    self.numYPxl = int(vertex.find('FRAME_ROW').text.strip())
                else:
                    self.latBL = float(vertex.find('FRAME_LAT').text.strip())
                    self.lonBL = float(vertex.find('FRAME_LON').text.strip())
                count = count + 1

            if count != 4:
                print("Count == " + str(count))
                raise ARCSIException("The number of vertex's found was incorrect so the corners of the scene have not been correctly defined.")

            self.latCentre = float(root.find('Dataset_Frame').find('Scene_Center').find('FRAME_LAT').text.strip())
            self.lonCentre = float(root.find('Dataset_Frame').find('Scene_Center').find('FRAME_LON').text.strip())

            self.pxlXRes = (self.lonTR - self.lonTL)/float(self.numXPxl)
            self.pxlYRes = (self.latTL - self.latBL)/float(self.numYPxl)

            print("TL: " + str(self.lonTL) + "," + str(self.latTL))
            print("TR: " + str(self.lonTR) + "," + str(self.latTR))
            print("BR: " + str(self.lonBR) + "," + str(self.latBR))
            print("BL: " + str(self.lonBL) + "," + str(self.latBL))
            print("Centre Pt: " + str(self.lonCentre) + "," + str(self.latCentre))

            epsgCode = int(root.find('Coordinate_Reference_System').find('Horizontal_CS').find('HORIZONTAL_CS_CODE').text.strip().split(":")[1])
            sceneProj = osr.SpatialReference()
            sceneProj.ImportFromEPSG(epsgCode)
            self.inWKT = sceneProj.ExportToWkt()

            # Get date and time of the acquisition
            acDate = root.find('Dataset_Sources').find('Source_Information').find('Scene_Source').find('IMAGING_DATE').text.strip().split('-')
            acTime = root.find('Dataset_Sources').find('Source_Information').find('Scene_Source').find('IMAGING_TIME').text.strip().split(':')
            self.acquisitionTime = datetime.datetime(int(acDate[0]), int(acDate[1]), int(acDate[2]), int(acTime[0]), int(acTime[1]), int(acTime[2]))
            print("Aq Time = " + str(self.acquisitionTime))

            self.solarZenith = 90-float(root.find('Dataset_Sources').find('Source_Information').find('Scene_Source').find('SUN_ELEVATION').text.strip())
            self.solarAzimuth = float(root.find('Dataset_Sources').find('Source_Information').find('Scene_Source').find('SUN_AZIMUTH').text.strip())

            self.acquIncidAngle = float(root.find('Dataset_Sources').find('Source_Information').find('Scene_Source').find('INCIDENCE_ANGLE').text.strip())
            print("self.acquIncidAngle: ", self.acquIncidAngle)
            self.acquViewAngle = float(root.find('Dataset_Sources').find('Source_Information').find('Scene_Source').find('VIEWING_ANGLE').text.strip())
            print("self.acquViewAngle: ", self.acquViewAngle)

            self.sensorZenith = self.acquViewAngle
            print("self.sensorZenith: ", self.sensorZenith)
            self.sensorAzimuth = self.acquIncidAngle
            print("self.sensorAzimuth: ", self.sensorAzimuth)

            filesDIR = os.path.dirname(inputHeader)
            if not self.userSpInputImage is None:
                self.fileName = os.path.abspath(self.userSpInputImage)
            else:
                self.fileName = os.path.join(filesDIR, root.find('Data_Access').find('Data_File').find('DATA_FILE_PATH').attrib['href'])
            print("self.fileName = " + self.fileName)

            spectralBandInfo = root.find('Image_Interpretation').findall('Spectral_Band_Info')
            for bandInfo in spectralBandInfo:
                if bandInfo.find('BAND_DESCRIPTION').text.strip() == 'XS1':
                    self.b1Bias = float(bandInfo.find('PHYSICAL_BIAS').text.strip())
                    self.b1Gain = float(bandInfo.find('PHYSICAL_GAIN').text.strip())
                elif bandInfo.find('BAND_DESCRIPTION').text.strip() == 'XS2':
                    self.b2Bias = float(bandInfo.find('PHYSICAL_BIAS').text.strip())
                    self.b2Gain = float(bandInfo.find('PHYSICAL_GAIN').text.strip())
                elif bandInfo.find('BAND_DESCRIPTION').text.strip() == 'XS3':
                    self.b3Bias = float(bandInfo.find('PHYSICAL_BIAS').text.strip())
                    self.b3Gain = float(bandInfo.find('PHYSICAL_GAIN').text.strip())
                elif bandInfo.find('BAND_DESCRIPTION').text.strip() == 'SWIR':
                    self.b4Bias = float(bandInfo.find('PHYSICAL_BIAS').text.strip())
                    self.b4Gain = float(bandInfo.find('PHYSICAL_GAIN').text.strip())
                else:
                    print("Band Description = " + bandInfo.find('BAND_DESCRIPTION').text.strip())
                    raise ARCSIException("Image band description has not be recognised.")
            print("Gains and Bias' found.")

            #for child in root.find('Data_Strip').find('Sensor_Calibration').find('Solar_Irradiance'):
            #    print(child.tag, child.attrib)

            spectralIrrInfo = root.find('Data_Strip').find('Sensor_Calibration').find('Solar_Irradiance').findall('Band_Solar_Irradiance')
            for spectralIrr in spectralIrrInfo:
                if spectralIrr.find('BAND_INDEX').text.strip() == '1':
                    self.b1SolarIrradiance = float(spectralIrr.find('SOLAR_IRRADIANCE_VALUE').text.strip())
                elif spectralIrr.find('BAND_INDEX').text.strip() == '2':
                    self.b2SolarIrradiance = float(spectralIrr.find('SOLAR_IRRADIANCE_VALUE').text.strip())
                elif spectralIrr.find('BAND_INDEX').text.strip() == '3':
                    self.b3SolarIrradiance = float(spectralIrr.find('SOLAR_IRRADIANCE_VALUE').text.strip())
                elif spectralIrr.find('BAND_INDEX').text.strip() == '4':
                    self.b4SolarIrradiance = float(spectralIrr.find('SOLAR_IRRADIANCE_VALUE').text.strip())
                else:
                    print("Band Index = " + spectralIrr.find('BAND_INDEX').text.strip())
                    raise ARCSIException("Image band index has not be recognised.")
            print("Solar Irradiance values found...")

            dataset = gdal.Open(self.fileName, gdal.GA_ReadOnly)
            if not dataset is None:
                if dataset.GetGCPCount() > 0:
                    self.inImgHasGCPs = True
                else:
                    self.inImgHasGCPs = False
                dataset = None
            else:
                print("Could not open image to set band names: ", imageFile)



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
        Customises the generic name for the SPOT5 sensor
        """
        outname = self.defaultGenBaseOutFileName()
        return outname

    def expectedImageDataPresent(self):
        imageDataPresent = True

        if not os.path.exists(self.fileName):
            imageDataPresent = False

        return imageDataPresent

    def applyImageDataMask(self, inputHeader, outputPath, outputMaskName, outputImgName, outFormat, outWKTFile):
        raise ARCSIException("SPOT5 does not provide any image masks, do not use the MASK option.")

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

        spot5Band = collections.namedtuple('SPOT5Band', ['bandName', 'bandIndex', 'bias', 'gain'])
        bandDefnSeq.append(spot5Band(bandName="NIR", bandIndex=3, bias=self.b1Bias, gain=self.b1Gain))
        bandDefnSeq.append(spot5Band(bandName="Red", bandIndex=2, bias=self.b2Bias, gain=self.b2Gain))
        bandDefnSeq.append(spot5Band(bandName="Green", bandIndex=1, bias=self.b3Bias, gain=self.b3Gain))
        bandDefnSeq.append(spot5Band(bandName="SWIR", bandIndex=4, bias=self.b4Bias, gain=self.b4Gain))
        rsgislib.imagecalibration.spot5ToRadiance(self.fileName, outputImage, outFormat, bandDefnSeq)

        if self.inImgHasGCPs:
            arcsiUtils = ARCSIUtils()
            arcsiUtils.copyGCPs(self.fileName, outputImage)
            rsgislib.imageutils.assignSpatialInfo(outputImage, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0)

        return outputImage, None

    def generateImageSaturationMask(self, outputPath, outputName, outFormat):
        print("Generate Saturation Image")
        outputImage = os.path.join(outputPath, outputName)

        spot5Band = collections.namedtuple('SPOT5Band', ['bandName', 'fileName', 'bandIndex', 'satVal'])
        bandDefnSeq = list()
        bandDefnSeq.append(spot5Band(bandName="NIR", fileName=self.fileName, bandIndex=1, satVal=255.0))
        bandDefnSeq.append(spot5Band(bandName="Red", fileName=self.fileName, bandIndex=2, satVal=255.0))
        bandDefnSeq.append(spot5Band(bandName="Green", fileName=self.fileName, bandIndex=3, satVal=255.0))
        bandDefnSeq.append(spot5Band(bandName="SWIR", fileName=self.fileName, bandIndex=4, satVal=255.0))

        rsgislib.imagecalibration.saturatedPixelsMask(outputImage, outFormat, bandDefnSeq)

        if self.inImgHasGCPs:
            arcsiUtils = ARCSIUtils()
            arcsiUtils.copyGCPs(self.fileName, outputImage)
            rsgislib.imageutils.assignSpatialInfo(outputImage, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0)

        return outputImage

    def generateValidImageDataMask(self, outputPath, outputMaskName, viewAngleImg, outFormat):
        print("Generate valid image mask")
        outputImage = os.path.join(outputPath, outputMaskName)
        return outputImage

    def convertThermalToBrightness(self, inputRadImage, outputPath, outputName, outFormat, scaleFactor):
        raise ARCSIException("There are no thermal bands...")

    def convertImageToTOARefl(self, inputRadImage, outputPath, outputName, outFormat, scaleFactor):
        print("Converting to TOA")
        outputImage = os.path.join(outputPath, outputName)
        solarIrradianceVals = list()
        IrrVal = collections.namedtuple('SolarIrradiance', ['irradiance'])
        solarIrradianceVals.append(IrrVal(irradiance=self.b1SolarIrradiance))
        solarIrradianceVals.append(IrrVal(irradiance=self.b2SolarIrradiance))
        solarIrradianceVals.append(IrrVal(irradiance=self.b3SolarIrradiance))
        solarIrradianceVals.append(IrrVal(irradiance=self.b4SolarIrradiance))
        rsgislib.imagecalibration.radiance2TOARefl(inputRadImage, outputImage, outFormat, rsgislib.TYPE_16UINT, scaleFactor, self.acquisitionTime.year, self.acquisitionTime.month, self.acquisitionTime.day, self.solarZenith, solarIrradianceVals)

        if self.inImgHasGCPs:
            arcsiUtils = ARCSIUtils()
            arcsiUtils.copyGCPs(self.fileName, outputImage)
            rsgislib.imageutils.assignSpatialInfo(outputImage, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0)

        return outputImage

    def generateCloudMask(self, inputReflImage, inputSatImage, inputThermalImage, inputViewAngleImg, inputValidImg, outputPath, outputName, outFormat, tmpPath, scaleFactor, cloud_msk_methods=None):
        raise ARCSIException("SPOT5 does not have a cloud masking implementation in ARCSI.")

    def createCloudMaskDataArray(self, inImgDataArr):
        return inImgDataArr

    def defineDarkShadowImageBand(self):
        return 3

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
        s.wavelength = Py6S.Wavelength(0.450, 0.650, [0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.007200,0.007200,0.054000,0.054000,0.178800,0.178800,0.371900,0.371900,0.557100,0.557100,0.557100,0.696000,0.696000,0.795900,0.795900,0.865400,0.865400,0.920800,0.920800,0.962700,0.962700,0.983800,0.983800,1.000000,1.000000,0.996200,0.996200,0.978300,0.978300,0.949700,0.949700,0.904500,0.904500,0.846800,0.846800,0.789500,0.789500,0.722900,0.722900,0.661100,0.661100,0.593700,0.593700,0.525600,0.525600,0.456500,0.456500,0.377600,0.377600,0.294300,0.294300,0.205300,0.205300,0.132300,0.132300,0.076700,0.076700,0.041900,0.041900,0.022000,0.022000,0.011700,0.011700,0.006700,0.006700,0.003900,0.003900,0.000000,0.000000,0.000000])
        s.run()
        sixsCoeffs[0,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[0,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[0,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[0,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[0,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[0,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 2
        s.wavelength = Py6S.Wavelength(0.580, 0.7425, [0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.001700,0.001700,0.022000,0.022000,0.116500,0.116500,0.315300,0.315300,0.533600,0.533600,0.704900,0.704900,0.826400,0.826400,0.902100,0.902100,0.951500,0.951500,0.979600,0.979600,0.994200,0.994200,1.000000,1.000000,0.997300,0.997300,0.981800,0.981800,0.945000,0.945000,0.863200,0.863200,0.717100,0.717100,0.534000,0.534000,0.356000,0.356000,0.222900,0.222900,0.132700,0.132700,0.079200,0.079200,0.047600,0.047600,0.029100,0.029100,0.018000,0.018000,0.011100,0.011100,0.006800,0.006800,0.004200,0.004200,0.002500,0.002500])
        s.run()
        sixsCoeffs[1,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[1,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[1,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[1,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[1,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[1,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 3
        s.wavelength = Py6S.Wavelength(0.750, 0.945, [0.000000,0.000000,0.005400,0.005400,0.016100,0.016100,0.042900,0.042900,0.096600,0.096600,0.183500,0.183500,0.297400,0.297400,0.433700,0.433700,0.565000,0.565000,0.691000,0.691000,0.792200,0.792200,0.873300,0.873300,0.928600,0.928600,0.965000,0.965000,0.989700,0.989700,0.997900,0.997900,1.000000,1.000000,0.991200,0.991200,0.978700,0.978700,0.960700,0.960700,0.940400,0.940400,0.921500,0.921500,0.890800,0.890800,0.869500,0.869500,0.814700,0.814700,0.734400,0.734400,0.612700,0.612700,0.457600,0.457600,0.327200,0.327200,0.216800,0.216800,0.134300,0.134300,0.087400,0.087400,0.058200,0.058200,0.035800,0.035800,0.023700,0.023700,0.016700,0.016700,0.011400,0.011400,0.007800,0.007800,0.005300,0.005300,0.000000])
        s.run()
        sixsCoeffs[2,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[2,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[2,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[2,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[2,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[2,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 4
        s.wavelength = Py6S.Wavelength(1.500, 1.8025, [0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.011100,0.011100,0.015000,0.015000,0.019100,0.019100,0.024800,0.024800,0.030700,0.030700,0.043300,0.043300,0.056100,0.056100,0.087800,0.087800,0.119800,0.119800,0.182900,0.182900,0.246700,0.246700,0.331300,0.331300,0.415500,0.415500,0.513000,0.513000,0.610200,0.610200,0.722900,0.722900,0.835700,0.835700,0.918300,0.918300,1.000000,1.000000,0.994600,0.994600,0.987500,0.987500,0.932100,0.932100,0.878200,0.878200,0.857800,0.857800,0.834700,0.834700,0.859400,0.859400,0.882900,0.882900,0.894700,0.894700,0.905000,0.905000,0.848800,0.848800,0.792500,0.792500,0.685900,0.685900,0.578300,0.578300,0.468300,0.468300,0.359100,0.359100,0.286600,0.286600,0.214300,0.214300,0.171700,0.171700,0.128600,0.128600,0.101200,0.101200,0.073900,0.073900,0.058000,0.058000,0.042400,0.042400,0.034100,0.034100,0.025900,0.025900,0.021500,0.021500,0.017100,0.017100,0.014100,0.014100,0.011200,0.011200,0.009000,0.009000,0.006900,0.006900,0.005500,0.005500,0.004100,0.004100,0.003400,0.003400,0.002600,0.002600,0.002300,0.002300,0.001900,0.001900])
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

        imgBandCoeffs = list()

        sixsCoeffs = self.calc6SCoefficients(aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF)

        imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=1, aX=float(sixsCoeffs[0,0]), bX=float(sixsCoeffs[0,1]), cX=float(sixsCoeffs[0,2]), DirIrr=float(sixsCoeffs[0,3]), DifIrr=float(sixsCoeffs[0,4]), EnvIrr=float(sixsCoeffs[0,5])))
        imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=2, aX=float(sixsCoeffs[1,0]), bX=float(sixsCoeffs[1,1]), cX=float(sixsCoeffs[1,2]), DirIrr=float(sixsCoeffs[1,3]), DifIrr=float(sixsCoeffs[1,4]), EnvIrr=float(sixsCoeffs[1,5])))
        imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=3, aX=float(sixsCoeffs[2,0]), bX=float(sixsCoeffs[2,1]), cX=float(sixsCoeffs[2,2]), DirIrr=float(sixsCoeffs[2,3]), DifIrr=float(sixsCoeffs[2,4]), EnvIrr=float(sixsCoeffs[2,5])))
        imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=4, aX=float(sixsCoeffs[3,0]), bX=float(sixsCoeffs[3,1]), cX=float(sixsCoeffs[3,2]), DirIrr=float(sixsCoeffs[3,3]), DifIrr=float(sixsCoeffs[3,4]), EnvIrr=float(sixsCoeffs[3,5])))

        rsgislib.imagecalibration.apply6SCoeffSingleParam(inputRadImage, outputImage, outFormat, rsgislib.TYPE_16UINT, scaleFactor, 0, True, imgBandCoeffs)

        if self.inImgHasGCPs:
            arcsiUtils = ARCSIUtils()
            arcsiUtils.copyGCPs(self.fileName, outputImage)
            rsgislib.imageutils.assignSpatialInfo(outputImage, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0)

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

        # Band 2 (Red!)
        s.wavelength = Py6S.Wavelength(0.580, 0.7425, [0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.001700,0.001700,0.022000,0.022000,0.116500,0.116500,0.315300,0.315300,0.533600,0.533600,0.704900,0.704900,0.826400,0.826400,0.902100,0.902100,0.951500,0.951500,0.979600,0.979600,0.994200,0.994200,1.000000,1.000000,0.997300,0.997300,0.981800,0.981800,0.945000,0.945000,0.863200,0.863200,0.717100,0.717100,0.534000,0.534000,0.356000,0.356000,0.222900,0.222900,0.132700,0.132700,0.079200,0.079200,0.047600,0.047600,0.029100,0.029100,0.018000,0.018000,0.011100,0.011100,0.006800,0.006800,0.004200,0.004200,0.002500,0.002500])
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
        raise ARCSIException("SPOT5 does not provide an implement of a method to derive AOT from DDV.")

    def estimateImageToAODUsingDDV(self, inputRADImage, inputTOAImage, inputDEMFile, shadowMask, outputPath, outputName, outFormat, tmpPath, aeroProfile, atmosProfile, grdRefl, aotValMin, aotValMax):
        raise ARCSIException("SPOT5 does not provide an implement of a method to derive AOT from DDV.")

    def estimateImageToAODUsingDOS(self, inputRADImage, inputTOAImage, inputDEMFile, shadowMask, outputPath, outputName, outFormat, tmpPath, aeroProfile, atmosProfile, grdRefl, aotValMin, aotValMax, globalDOS, simpleDOS, dosOutRefl):
        try:
            print("Estimating AOD Using DOS")
            arcsiUtils = ARCSIUtils()
            outputAOTImage = os.path.join(outputPath, outputName)
            tmpBaseName = os.path.splitext(outputName)[0]
            imgExtension = arcsiUtils.getFileExtension(outFormat)

            dosGreenImage = ""
            minObjSize = 5
            darkPxlPercentile = 0.01
            blockSize = 1000
            if simpleDOS:
                outputDOSBlueName = tmpBaseName + "DOSGreen" + imgExtension
                dosGreenImage, bandOff = self.convertImageBandToReflectanceSimpleDarkSubtract(inputTOAImage, outputPath, outputDOSBlueName, outFormat, dosOutRefl, 1)
            elif globalDOS:
                dosGreenImage = self.performDOSOnSingleBand(inputTOAImage, 1, outputPath, tmpBaseName, "Green", "KEA", tmpPath, minObjSize, darkPxlPercentile, dosOutRefl)
            else:
                dosGreenImage = self.performLocalDOSOnSingleBand(inputTOAImage, 1, outputPath, tmpBaseName, "Green", "KEA", tmpPath, minObjSize, darkPxlPercentile, blockSize, dosOutRefl)

            thresImageClumpsFinal = os.path.join(tmpPath, tmpBaseName + "_clumps" + imgExtension)
            rsgislib.segmentation.segutils.runShepherdSegmentation(inputTOAImage, thresImageClumpsFinal, tmpath=tmpPath, gdalformat="KEA", numClusters=40, minPxls=10, bands=[5,4,1], processInMem=True)

            stats2CalcTOA = list()
            stats2CalcTOA.append(rsgislib.rastergis.BandAttStats(band=1, meanField="MeanElev"))
            rsgislib.rastergis.populateRATWithStats(inputDEMFile, thresImageClumpsFinal, stats2CalcTOA)

            stats2CalcTOA = list()
            stats2CalcTOA.append(rsgislib.rastergis.BandAttStats(band=1, meanField="MeanB1DOS"))
            rsgislib.rastergis.populateRATWithStats(dosGreenImage, thresImageClumpsFinal, stats2CalcTOA)

            stats2CalcRad = list()
            stats2CalcRad.append(rsgislib.rastergis.BandAttStats(band=1, meanField="MeanB1RAD"))
            stats2CalcRad.append(rsgislib.rastergis.BandAttStats(band=3, meanField="MeanB3RAD"))
            stats2CalcRad.append(rsgislib.rastergis.BandAttStats(band=2, meanField="MeanB2RAD"))
            rsgislib.rastergis.populateRATWithStats(inputRADImage, thresImageClumpsFinal, stats2CalcRad)

            ratDS = gdal.Open(thresImageClumpsFinal, gdal.GA_Update)
            Histogram = rat.readColumn(ratDS, "Histogram")
            MeanElev = rat.readColumn(ratDS, "MeanElev")

            MeanB5RAD = rat.readColumn(ratDS, "MeanB3RAD")
            MeanB3RAD = rat.readColumn(ratDS, "MeanB2RAD")

            radNDVI = (MeanB3RAD - MeanB2RAD)/(MeanB3RAD + MeanB2RAD)

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
                gdalDriver.Delete(dosGreenImage)

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
            dataset.GetRasterBand(1).SetDescription("Green")
            dataset.GetRasterBand(2).SetDescription("Red")
            dataset.GetRasterBand(3).SetDescription("NIR")
            dataset.GetRasterBand(4).SetDescription("SWIR")
            dataset = None
        else:
            print("Could not open image to set band names: ", imageFile)

    def cleanLocalFollowProcessing(self):
        print("")



