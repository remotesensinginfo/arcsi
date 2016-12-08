"""
Module that contains the ARCSISentinel2Sensor class.
"""
############################################################################
#  arcsisensorsentinel2.py
#
#  Copyright 2016 ARCSI.
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
# Date: 09/07/2016
# Version: 1.0
#
# History:
# Version 1.0 - Created.
#
############################################################################

# Import the future functionality (for Python 2)
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals
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
# Import the RIOS RAT module
from rios import rat

class ARCSISentinel2Sensor (ARCSIAbstractSensor):
    """
    A class which represents the Sentinel-2 sensor to read
    header parameters and apply data processing operations.
    """
    def __init__(self, debugMode, inputImage):
        ARCSIAbstractSensor.__init__(self, debugMode, inputImage)
        self.sensor = "SEN2"

    def extractHeaderParameters(self, inputHeader, wktStr):
        """
        Understands and parses the WorldView2 xml header file
        """
        try:
            self.headerFileName = os.path.split(inputHeader)[1]
            
            print("hello world...")
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
        return (0.0, 0.0)

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

    def imgNeedMosaicking(self):
        return True

    def applyImageDataMask(self, inputHeader, outputPath, outputMaskName, outputImgName, outFormat, outWKTFile):
        raise ARCSIException("Sentinel-2 does not provide any image masks, do not use the MASK option.")

    def mosaicImageTiles(self):
        raise ARCSIException("Image data does not need mosaicking")

    def generateValidImageDataMask(self, outputPath, outputMaskName, viewAngleImg, outFormat):
        raise ARCSIException("Don't know how to create a valid data mask for Sentinel-2")

    def convertImageToRadiance(self, outputPath, outputReflName, outputThermalName, outFormat):
        print("Converting to Radiance")
        raise ARCSIException("Not Implemented")

    def generateImageSaturationMask(self, outputPath, outputName, outFormat):
        print("Generate Saturation Image")
        raise ARCSIException("Not Implemented")

    def convertThermalToBrightness(self, inputRadImage, outputPath, outputName, outFormat, scaleFactor):
        raise ARCSIException("Not Implemented")

    def convertImageToTOARefl(self, inputRadImage, outputPath, outputName, outFormat, scaleFactor):
        print("Converting to TOA")
        raise ARCSIException("Not Implemented")

    def generateCloudMask(self, inputReflImage, inputSatImage, inputThermalImage, inputValidImg, outputPath, outputName, outFormat, tmpPath, scaleFactor):
        raise ARCSIException("Cloud Masking Not Implemented for Sentinel-2.")

    def calc6SCoefficients(self, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF):
        sixsCoeffs = numpy.zeros((12, 3), dtype=numpy.float32)
        # Set up 6S model
        s = Py6S.SixS()
        s.atmos_profile = atmosProfile
        s.aero_profile = aeroProfile
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
        s.wavelength = Py6S.Wavelength(0.430, 0.4575, [0.015297, 0.19593, 0.511598, 0.587385, 0.699199, 0.792047, 0.993878, 0.990178, 0.956726, 0.723933, 0.051814, 0.000303])
        s.run()
        sixsCoeffs[0,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[0,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[0,2] = float(s.outputs.values['coef_xc'])

        # Band 2
        s.wavelength = Py6S.Wavelength(0.440, 0.535, [0.001206, 0.002623, 0.002076, 0.002224, 0.002377, 0.002856, 0.009028, 0.038955, 0.292197, 0.382418, 0.400158, 0.424686, 0.505323, 0.529543, 0.534656, 0.543691, 0.601967, 0.621092, 0.575863, 0.546131, 0.571684, 0.633236, 0.738396, 0.768325, 0.788363, 0.809151, 0.844983, 0.840111, 0.78694, 0.761923, 0.810031, 0.901671, 1.0, 0.908308, 0.286992, 0.102833, 0.02508, 0.002585, 0.000441])
        s.run()
        sixsCoeffs[1,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[1,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[1,2] = float(s.outputs.values['coef_xc'])

        # Band 3
        s.wavelength = Py6S.Wavelength(0.5375, 0.5825, [0.00084, 0.080665, 0.341374, 0.828036, 0.888565, 0.860271, 0.834035, 0.867734, 0.933938, 1.0, 0.981107, 0.868656, 0.81291, 0.789606, 0.830458, 0.85799, 0.62498, 0.098293, 0.016512])
        s.run()
        sixsCoeffs[2,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[2,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[2,2] = float(s.outputs.values['coef_xc'])

        # Band 4
        s.wavelength = Py6S.Wavelength(0.6475, 0.6825, [0.034529, 0.817746, 0.983869, 0.995449, 0.977215, 0.814166, 0.764864, 0.830828, 0.883581, 0.955931, 0.973219, 0.965712, 0.944811, 0.422967, 0.063172])
        s.run()
        sixsCoeffs[3,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[3,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[3,2] = float(s.outputs.values['coef_xc'])

        # Band 5
        s.wavelength = Py6S.Wavelength(0.695, 0.7125, [0.04126, 0.478496, 1, 0.993239, 0.945953, 0.902399, 0.757197, 0.196706])
        s.run()
        sixsCoeffs[4,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[4,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[4,2] = float(s.outputs.values['coef_xc'])

        # Band 6
        s.wavelength = Py6S.Wavelength(0.7325, 0.7475, [0.085006, 0.920265, 0.934211, 0.981932, 0.993406, 0.962584, 0.506722])
        s.run()
        sixsCoeffs[5,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[5,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[5,2] = float(s.outputs.values['coef_xc'])

        # Band 7
        s.wavelength = Py6S.Wavelength(0.770, 0.7975, [0.014731, 0.199495, 0.898494, 0.994759, 0.964657, 0.846898, 0.777241, 0.800984, 0.757695, 0.536855, 0.077219, 0.003152])
        s.run()
        sixsCoeffs[6,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[6,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[6,2] = float(s.outputs.values['coef_xc'])

        # Band 8
        s.wavelength = Py6S.Wavelength(0.775, 0.9075, [0.019072, 0.056536, 0.203436, 0.450085, 0.81829, 0.960732, 0.985213, 0.93655, 0.941281, 0.962183, 0.959009, 0.945147, 0.945357, 0.937084, 0.900979, 0.86216, 0.801819, 0.755632, 0.708669, 0.690211, 0.682649, 0.67595, 0.660812, 0.65831, 0.685501, 0.720686, 0.776608, 0.78772, 0.776161, 0.759264, 0.720589, 0.69087, 0.649339, 0.627424, 0.604322, 0.591724, 0.581202, 0.580197, 0.589481, 0.596749, 0.605476, 0.613463, 0.637436, 0.659233, 0.659924, 0.615841, 0.526407, 0.49653, 0.529093, 0.537964, 0.326791, 0.14854, 0.033246, 0.007848])
        s.run()
        sixsCoeffs[7,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[7,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[7,2] = float(s.outputs.values['coef_xc'])

        # Band 9
        s.wavelength = Py6S.Wavelength(0.850, 0.880, [0.02471, 0.104944, 0.585731, 0.87843, 0.926043, 0.935962, 0.965458, 0.97988, 0.988474, 0.999626, 0.472189, 0.106955, 0.008819])
        s.run()
        sixsCoeffs[8,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[8,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[8,2] = float(s.outputs.values['coef_xc'])

        # Band 10
        s.wavelength = Py6S.Wavelength(0.9325, 0.9575, [0.018022, 0.408108, 0.873658, 0.983566, 0.996767, 0.998123, 1.0, 0.956408, 0.931094, 0.450443, 0.059807])
        s.run()
        sixsCoeffs[9,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[9,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[9,2] = float(s.outputs.values['coef_xc'])

        # Band 11
        s.wavelength = Py6S.Wavelength(1.540, 1.6825, [7.00E-06, 2.80E-05, 0.000147, 0.00048, 0.000911, 0.001684, 0.005345, 0.012628, 0.039584, 0.07493, 0.182597, 0.330736, 0.647173, 0.815215, 0.88703, 0.891417, 0.916528, 0.935322, 0.951416, 0.956429, 0.96348, 0.96818, 0.975915, 0.979878, 0.981412, 0.980705, 0.982736, 0.987807, 0.993288, 0.990405, 0.980023, 0.972568, 0.966371, 0.96605, 0.973463, 0.983472, 0.995476, 0.998568, 0.998804, 0.99973, 0.999814, 0.99162, 0.969903, 0.953287, 0.938586, 0.928114, 0.82498, 0.641891, 0.32371, 0.163972, 0.046194, 0.019359, 0.006523, 0.003409, 0.001423, 0.000498, 3.40E-05, 1.30E-05])
        s.run()
        sixsCoeffs[10,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[10,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[10,2] = float(s.outputs.values['coef_xc'])

        # Band 12
        s.wavelength = Py6S.Wavelength(2.080, 2.320, [0.002885, 0.006597, 0.00854, 0.010002, 0.013364, 0.017126, 0.027668, 0.040217, 0.073175, 0.11147, 0.203461, 0.284898, 0.408003, 0.476537, 0.543352, 0.568634, 0.598891, 0.621362, 0.663707, 0.696165, 0.741301, 0.772071, 0.809677, 0.828599, 0.851107, 0.854746, 0.859532, 0.863257, 0.869696, 0.878588, 0.889473, 0.896696, 0.904831, 0.905665, 0.904783, 0.903347, 0.901983, 0.904313, 0.908092, 0.91295, 0.921302, 0.927219, 0.934142, 0.937086, 0.937652, 0.942518, 0.942117, 0.938428, 0.933022, 0.921057, 0.908293, 0.908191, 0.922855, 0.919482, 0.924526, 0.931974, 0.946802, 0.954437, 0.962539, 0.966042, 0.96546, 0.963656, 0.957327, 0.953558, 0.951731, 0.952641, 0.960639, 0.968307, 0.982898, 0.990734, 0.998753, 0.999927, 0.993884, 0.983735, 0.958343, 0.938203, 0.905999, 0.881683, 0.84062, 0.809516, 0.749107, 0.688185, 0.566031, 0.474659, 0.342092, 0.263176, 0.16809, 0.124831, 0.082363, 0.062691, 0.042864, 0.034947, 0.027418, 0.023959, 0.016331, 0.007379, 0.002065])
        s.run()
        sixsCoeffs[11,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[11,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[11,2] = float(s.outputs.values['coef_xc'])

        return sixsCoeffs

    def convertImageToSurfaceReflSglParam(self, inputRadImage, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF, scaleFactor):
        print("Converting to Surface Reflectance")
        raise ARCSIException("Not Implemented")

    def convertImageToSurfaceReflDEMElevLUT(self, inputRadImage, inputDEMFile, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, aotVal, useBRDF, surfaceAltitudeMin, surfaceAltitudeMax, scaleFactor, elevCoeffs=None):
        print("Converting to Surface Reflectance")
        raise ARCSIException("Not Implemented")

    def convertImageToSurfaceReflAOTDEMElevLUT(self, inputRadImage, inputDEMFile, inputAOTImage, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, useBRDF, surfaceAltitudeMin, surfaceAltitudeMax, aotMin, aotMax, scaleFactor, elevAOTCoeffs=None):
        print("Converting to Surface Reflectance")
        raise ARCSIException("Not Implemented")

    def run6SToOptimiseAODValue(self, aotVal, radBlueVal, predBlueVal, aeroProfile, atmosProfile, grdRefl, surfaceAltitude):
        """Used as part of the optimastion for identifying values of AOD"""
        raise ARCSIException("Not Implemented")

    def findDDVTargets(self, inputTOAImage, outputPath, outputName, outFormat, tmpPath):
        print("Not implemented\n")
        sys.exit()

    def estimateImageToAODUsingDDV(self, inputRADImage, inputTOAImage, inputDEMFile, shadowMask, outputPath, outputName, outFormat, tmpPath, aeroProfile, atmosProfile, grdRefl, aotValMin, aotValMax):
        print("Not implemented\n")
        sys.exit()

    def estimateImageToAODUsingDOS(self, inputRADImage, inputTOAImage, inputDEMFile, shadowMask, outputPath, outputName, outFormat, tmpPath, aeroProfile, atmosProfile, grdRefl, aotValMin, aotValMax, globalDOS, simpleDOS, dosOutRefl):
        raise ARCSIException("Not Implemented")

    def estimateSingleAOTFromDOS(self, radianceImage, toaImage, inputDEMFile, tmpPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, minAOT, maxAOT, dosOutRefl):
        raise ARCSIException("Not Implemented")

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


