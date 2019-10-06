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
# Import the py6s module for running 6S from python.
import Py6S
# Import python math module
import math

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

        self.blueSpecRespFuncMin = 0.4
        self.blueSpecRespFuncMax = 0.99
        self.blueSpecRespFunc = [0.005715, 0.005715, 0.005715, 0.017144, 0.017144, 0.017144, 0.017144, 0.003742,
                                 0.003742, 0.003742, 0.003742, 0.005538, 0.005538, 0.005538, 0.005538, 0.034168,
                                 0.034168, 0.034168, 0.034168, 0.243317, 0.243317, 0.243317, 0.243317, 0.927788,
                                 0.927788, 0.927788, 0.927788, 0.998238, 0.998238, 0.998238, 0.998238, 0.981119,
                                 0.981119, 0.981119, 0.981119, 0.917197, 0.917197, 0.917197, 0.917197, 0.807515,
                                 0.807515, 0.807515, 0.807515, 0.679484, 0.679484, 0.679484, 0.679484, 0.502155,
                                 0.502155, 0.502155, 0.502155, 0.337688, 0.337688, 0.337688, 0.337688, 0.252910,
                                 0.252910, 0.252910, 0.252910, 0.168263, 0.168263, 0.168263, 0.168263, 0.102364,
                                 0.102364, 0.102364, 0.102364, 0.070991, 0.070991, 0.070991, 0.070991, 0.050316,
                                 0.050316, 0.050316, 0.050316, 0.043581, 0.043581, 0.043581, 0.043581, 0.037542,
                                 0.037542, 0.037542, 0.037542, 0.033158, 0.033158, 0.033158, 0.033158, 0.033498,
                                 0.033498, 0.033498, 0.033498, 0.036252, 0.036252, 0.036252, 0.036252, 0.045794,
                                 0.045794, 0.045794, 0.045794, 0.058726, 0.058726, 0.058726, 0.058726, 0.083455,
                                 0.083455, 0.083455, 0.083455, 0.093482, 0.093482, 0.093482, 0.093482, 0.068651,
                                 0.068651, 0.068651, 0.068651, 0.041995, 0.041995, 0.041995, 0.041995, 0.004385,
                                 0.004385, 0.004385, 0.004385, 0.000065, 0.000065, 0.000065, 0.000065, 0.000000,
                                 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000055,
                                 0.000055, 0.000055, 0.000055, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000390,
                                 0.000390, 0.000390, 0.000390, 0.001855, 0.001855, 0.001855, 0.001855, 0.001890,
                                 0.001890, 0.001890, 0.001890, 0.001268, 0.001268, 0.001268, 0.001268, 0.001457,
                                 0.001457, 0.001457, 0.001457, 0.002449, 0.002449, 0.002449, 0.002449, 0.004190,
                                 0.004190, 0.004190, 0.004190, 0.005073, 0.005073, 0.005073, 0.005073, 0.004262,
                                 0.004262, 0.004262, 0.004262, 0.003005, 0.003005, 0.003005, 0.003005, 0.002154,
                                 0.002154, 0.002154, 0.002154, 0.001753, 0.001753, 0.001753, 0.001753, 0.001567,
                                 0.001567, 0.001567, 0.001567, 0.001651, 0.001651, 0.001651, 0.001651, 0.001784,
                                 0.001784, 0.001784, 0.001784, 0.001802, 0.001802, 0.001802, 0.001802, 0.001647,
                                 0.001647, 0.001647, 0.001647, 0.001271, 0.001271, 0.001271, 0.001271, 0.001020,
                                 0.001020, 0.001020, 0.001020, 0.000747, 0.000747]

        self.greenSpecRespFuncMin = 0.4
        self.greenSpecRespFuncMax = 0.99
        self.greenSpecRespFunc = [0.002972, 0.002972, 0.002972, 0.003226, 0.003226, 0.003226, 0.003226, 0.002532,
                                  0.002532, 0.002532, 0.002532, 0.002815, 0.002815, 0.002815, 0.002815, 0.004382,
                                  0.004382, 0.004382, 0.004382, 0.055860, 0.055860, 0.055860, 0.055860, 0.198928,
                                  0.198928, 0.198928, 0.198928, 0.285744, 0.285744, 0.285744, 0.285744, 0.395509,
                                  0.395509, 0.395509, 0.395509, 0.505254, 0.505254, 0.505254, 0.505254, 0.639854,
                                  0.639854, 0.639854, 0.639854, 0.801640, 0.801640, 0.801640, 0.801640, 0.927988,
                                  0.927988, 0.927988, 0.927988, 0.983814, 0.983814, 0.983814, 0.983814, 1.000000,
                                  1.000000, 1.000000, 1.000000, 0.956976, 0.956976, 0.956976, 0.956976, 0.886456,
                                  0.886456, 0.886456, 0.886456, 0.796302, 0.796302, 0.796302, 0.796302, 0.681969,
                                  0.681969, 0.681969, 0.681969, 0.534918, 0.534918, 0.534918, 0.534918, 0.375797,
                                  0.375797, 0.375797, 0.375797, 0.224668, 0.224668, 0.224668, 0.224668, 0.150505,
                                  0.150505, 0.150505, 0.150505, 0.114317, 0.114317, 0.114317, 0.114317, 0.101420,
                                  0.101420, 0.101420, 0.101420, 0.081920, 0.081920, 0.081920, 0.081920, 0.079764,
                                  0.079764, 0.079764, 0.079764, 0.098114, 0.098114, 0.098114, 0.098114, 0.091845,
                                  0.091845, 0.091845, 0.091845, 0.073635, 0.073635, 0.073635, 0.073635, 0.008610,
                                  0.008610, 0.008610, 0.008610, 0.000563, 0.000563, 0.000563, 0.000563, 0.000000,
                                  0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                  0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000007,
                                  0.000007, 0.000007, 0.000007, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                  0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000357,
                                  0.000357, 0.000357, 0.000357, 0.001986, 0.001986, 0.001986, 0.001986, 0.001524,
                                  0.001524, 0.001524, 0.001524, 0.001071, 0.001071, 0.001071, 0.001071, 0.001259,
                                  0.001259, 0.001259, 0.001259, 0.002194, 0.002194, 0.002194, 0.002194, 0.003612,
                                  0.003612, 0.003612, 0.003612, 0.004444, 0.004444, 0.004444, 0.004444, 0.003703,
                                  0.003703, 0.003703, 0.003703, 0.002630, 0.002630, 0.002630, 0.002630, 0.001886,
                                  0.001886, 0.001886, 0.001886, 0.001527, 0.001527, 0.001527, 0.001527, 0.001404,
                                  0.001404, 0.001404, 0.001404, 0.001481, 0.001481, 0.001481, 0.001481, 0.001588,
                                  0.001588, 0.001588, 0.001588, 0.001600, 0.001600, 0.001600, 0.001600, 0.001449,
                                  0.001449, 0.001449, 0.001449, 0.001130, 0.001130, 0.001130, 0.001130, 0.000903,
                                  0.000903, 0.000903, 0.000903, 0.000694, 0.000694]

        self.redSpecRespFuncMin = 0.4
        self.redSpecRespFuncMax = 0.99
        self.redSpecRespFunc = [0.001561, 0.001561, 0.001561, 0.001745, 0.001745, 0.001745, 0.001745, 0.000704,
                                0.000704, 0.000704, 0.000704, 0.000833, 0.000833, 0.000833, 0.000833, 0.001107,
                                0.001107, 0.001107, 0.001107, 0.009164, 0.009164, 0.009164, 0.009164, 0.043702,
                                0.043702, 0.043702, 0.043702, 0.056346, 0.056346, 0.056346, 0.056346, 0.070534,
                                0.070534, 0.070534, 0.070534, 0.057347, 0.057347, 0.057347, 0.057347, 0.067886,
                                0.067886, 0.067886, 0.067886, 0.070524, 0.070524, 0.070524, 0.070524, 0.069002,
                                0.069002, 0.069002, 0.069002, 0.079968, 0.079968, 0.079968, 0.079968, 0.067132,
                                0.067132, 0.067132, 0.067132, 0.052190, 0.052190, 0.052190, 0.052190, 0.051645,
                                0.051645, 0.051645, 0.051645, 0.075141, 0.075141, 0.075141, 0.075141, 0.283971,
                                0.283971, 0.283971, 0.283971, 0.720800, 0.720800, 0.720800, 0.720800, 0.937774,
                                0.937774, 0.937774, 0.937774, 1.000000, 1.000000, 1.000000, 1.000000, 0.988703,
                                0.988703, 0.988703, 0.988703, 0.985597, 0.985597, 0.985597, 0.985597, 0.938495,
                                0.938495, 0.938495, 0.938495, 0.905978, 0.905978, 0.905978, 0.905978, 0.869700,
                                0.869700, 0.869700, 0.869700, 0.827064, 0.827064, 0.827064, 0.827064, 0.632247,
                                0.632247, 0.632247, 0.632247, 0.356836, 0.356836, 0.356836, 0.356836, 0.041825,
                                0.041825, 0.041825, 0.041825, 0.002260, 0.002260, 0.002260, 0.002260, 0.000473,
                                0.000473, 0.000473, 0.000473, 0.000165, 0.000165, 0.000165, 0.000165, 0.000002,
                                0.000002, 0.000002, 0.000002, 0.000000, 0.000000, 0.000000, 0.000000, 0.000045,
                                0.000045, 0.000045, 0.000045, 0.000041, 0.000041, 0.000041, 0.000041, 0.000023,
                                0.000023, 0.000023, 0.000023, 0.000061, 0.000061, 0.000061, 0.000061, 0.000563,
                                0.000563, 0.000563, 0.000563, 0.002845, 0.002845, 0.002845, 0.002845, 0.001958,
                                0.001958, 0.001958, 0.001958, 0.001308, 0.001308, 0.001308, 0.001308, 0.001519,
                                0.001519, 0.001519, 0.001519, 0.002449, 0.002449, 0.002449, 0.002449, 0.004086,
                                0.004086, 0.004086, 0.004086, 0.004939, 0.004939, 0.004939, 0.004939, 0.004081,
                                0.004081, 0.004081, 0.004081, 0.002846, 0.002846, 0.002846, 0.002846, 0.002036,
                                0.002036, 0.002036, 0.002036, 0.001629, 0.001629, 0.001629, 0.001629, 0.001523,
                                0.001523, 0.001523, 0.001523, 0.001623, 0.001623, 0.001623, 0.001623, 0.001732,
                                0.001732, 0.001732, 0.001732, 0.001694, 0.001694, 0.001694, 0.001694, 0.001549,
                                0.001549, 0.001549, 0.001549, 0.001246, 0.001246, 0.001246, 0.001246, 0.000969,
                                0.000969, 0.000969, 0.000969, 0.000722, 0.000722]

        self.nirSpecRespFuncMin = 0.4
        self.nirSpecRespFuncMax = 0.99
        self.nirSpecRespFunc = [0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                0.000000, 0.000000, 0.000000, 0.000523, 0.000523, 0.000523, 0.000523, 0.000000,
                                0.000000, 0.000000, 0.000000, 0.000861, 0.000861, 0.000861, 0.000861, 0.001066,
                                0.001066, 0.001066, 0.001066, 0.001109, 0.001109, 0.001109, 0.001109, 0.001710,
                                0.001710, 0.001710, 0.001710, 0.000288, 0.000288, 0.000288, 0.000288, 0.000497,
                                0.000497, 0.000497, 0.000497, 0.000824, 0.000824, 0.000824, 0.000824, 0.001211,
                                0.001211, 0.001211, 0.001211, 0.000619, 0.000619, 0.000619, 0.000619, 0.000934,
                                0.000934, 0.000934, 0.000934, 0.001321, 0.001321, 0.001321, 0.001321, 0.001184,
                                0.001184, 0.001184, 0.001184, 0.001484, 0.001484, 0.001484, 0.001484, 0.001756,
                                0.001756, 0.001756, 0.001756, 0.002109, 0.002109, 0.002109, 0.002109, 0.000683,
                                0.000683, 0.000683, 0.000683, 0.000617, 0.000617, 0.000617, 0.000617, 0.001981,
                                0.001981, 0.001981, 0.001981, 0.035315, 0.035315, 0.035315, 0.035315, 0.001457,
                                0.001457, 0.001457, 0.001457, 0.001801, 0.001801, 0.001801, 0.001801, 0.002020,
                                0.002020, 0.002020, 0.002020, 0.005506, 0.005506, 0.005506, 0.005506, 0.003080,
                                0.003080, 0.003080, 0.003080, 0.053163, 0.053163, 0.053163, 0.053163, 0.004693,
                                0.004693, 0.004693, 0.004693, 0.000518, 0.000518, 0.000518, 0.000518, 0.002206,
                                0.002206, 0.002206, 0.002206, 0.000360, 0.000360, 0.000360, 0.000360, 0.000555,
                                0.000555, 0.000555, 0.000555, 0.003616, 0.003616, 0.003616, 0.003616, 0.100389,
                                0.100389, 0.100389, 0.100389, 0.738535, 0.738535, 0.738535, 0.738535, 0.966667,
                                0.966667, 0.966667, 0.966667, 0.961212, 0.961212, 0.961212, 0.961212, 0.895448,
                                0.895448, 0.895448, 0.895448, 0.820819, 0.820819, 0.820819, 0.820819, 0.786388,
                                0.786388, 0.786388, 0.786388, 0.723819, 0.723819, 0.723819, 0.723819, 0.686119,
                                0.686119, 0.686119, 0.686119, 0.618513, 0.618513, 0.618513, 0.618513, 0.555248,
                                0.555248, 0.555248, 0.555248, 0.353697, 0.353697, 0.353697, 0.353697, 0.079050,
                                0.079050, 0.079050, 0.079050, 0.002591, 0.002591, 0.002591, 0.002591, 0.000043,
                                0.000043, 0.000043, 0.000043, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                0.000000, 0.000000, 0.000000, 0.000000, 0.000000]

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

                sensorNameID = ''
                if sensorNameID == "0f" or sensorNameID == "10":
                    self.blueSpecRespFuncMin = 0.4
                    self.blueSpecRespFuncMax = 0.99
                    self.blueSpecRespFunc = [0.004000, 0.004000, 0.004000, 0.004000, 0.004000, 0.004000, 0.004000,
                                             0.004000, 0.004000, 0.004000, 0.004000, 0.005000, 0.005000, 0.005000,
                                             0.005000, 0.010000, 0.010000, 0.010000, 0.010000, 0.473000, 0.473000,
                                             0.473000, 0.473000, 0.962000, 0.962000, 0.962000, 0.962000, 0.995000,
                                             0.995000, 0.995000, 0.995000, 1.000000, 1.000000, 1.000000, 1.000000,
                                             0.942000, 0.942000, 0.942000, 0.942000, 0.829000, 0.829000, 0.829000,
                                             0.829000, 0.712000, 0.712000, 0.712000, 0.712000, 0.549000, 0.549000,
                                             0.549000, 0.549000, 0.377000, 0.377000, 0.377000, 0.377000, 0.285000,
                                             0.285000, 0.285000, 0.285000, 0.215000, 0.215000, 0.215000, 0.215000,
                                             0.145000, 0.145000, 0.145000, 0.145000, 0.096000, 0.096000, 0.096000,
                                             0.096000, 0.077000, 0.077000, 0.077000, 0.077000, 0.065000, 0.065000,
                                             0.065000, 0.065000, 0.056000, 0.056000, 0.056000, 0.056000, 0.050000,
                                             0.050000, 0.050000, 0.050000, 0.049000, 0.049000, 0.049000, 0.049000,
                                             0.053000, 0.053000, 0.053000, 0.053000, 0.065000, 0.065000, 0.065000,
                                             0.065000, 0.084000, 0.084000, 0.084000, 0.084000, 0.107000, 0.107000,
                                             0.107000, 0.107000, 0.123000, 0.123000, 0.123000, 0.123000, 0.063000,
                                             0.063000, 0.063000, 0.063000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000]

                    self.greenSpecRespFuncMin = 0.4
                    self.greenSpecRespFuncMax = 0.99
                    self.greenSpecRespFunc = [0.004000, 0.004000, 0.004000, 0.004000, 0.004000, 0.004000, 0.004000,
                                              0.004000, 0.004000, 0.004000, 0.004000, 0.004000, 0.004000, 0.004000,
                                              0.004000, 0.005000, 0.005000, 0.005000, 0.005000, 0.100000, 0.100000,
                                              0.100000, 0.100000, 0.204000, 0.204000, 0.204000, 0.204000, 0.246000,
                                              0.246000, 0.246000, 0.246000, 0.368000, 0.368000, 0.368000, 0.368000,
                                              0.514000, 0.514000, 0.514000, 0.514000, 0.656000, 0.656000, 0.656000,
                                              0.656000, 0.814000, 0.814000, 0.814000, 0.814000, 0.938000, 0.938000,
                                              0.938000, 0.938000, 0.993000, 0.993000, 0.993000, 0.993000, 1.000000,
                                              1.000000, 1.000000, 1.000000, 0.974000, 0.974000, 0.974000, 0.974000,
                                              0.915000, 0.915000, 0.915000, 0.915000, 0.820000, 0.820000, 0.820000,
                                              0.820000, 0.711000, 0.711000, 0.711000, 0.711000, 0.567000, 0.567000,
                                              0.567000, 0.567000, 0.399000, 0.399000, 0.399000, 0.399000, 0.257000,
                                              0.257000, 0.257000, 0.257000, 0.171000, 0.171000, 0.171000, 0.171000,
                                              0.145000, 0.145000, 0.145000, 0.145000, 0.124000, 0.124000, 0.124000,
                                              0.124000, 0.110000, 0.110000, 0.110000, 0.110000, 0.104000, 0.104000,
                                              0.104000, 0.104000, 0.114000, 0.114000, 0.114000, 0.114000, 0.067000,
                                              0.067000, 0.067000, 0.067000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000]

                    self.redSpecRespFuncMin = 0.4
                    self.redSpecRespFuncMax = 0.99
                    self.redSpecRespFunc = [0.001000, 0.001000, 0.001000, 0.001000, 0.001000, 0.001000, 0.001000,
                                            0.001000, 0.001000, 0.001000, 0.001000, 0.001000, 0.001000, 0.001000,
                                            0.001000, 0.001000, 0.001000, 0.001000, 0.001000, 0.019000, 0.019000,
                                            0.019000, 0.019000, 0.030000, 0.030000, 0.030000, 0.030000, 0.029000,
                                            0.029000, 0.029000, 0.029000, 0.032000, 0.032000, 0.032000, 0.032000,
                                            0.031000, 0.031000, 0.031000, 0.031000, 0.028000, 0.028000, 0.028000,
                                            0.028000, 0.031000, 0.031000, 0.031000, 0.031000, 0.048000, 0.048000,
                                            0.048000, 0.048000, 0.067000, 0.067000, 0.067000, 0.067000, 0.049000,
                                            0.049000, 0.049000, 0.049000, 0.024000, 0.024000, 0.024000, 0.024000,
                                            0.021000, 0.021000, 0.021000, 0.021000, 0.048000, 0.048000, 0.048000,
                                            0.048000, 0.283000, 0.283000, 0.283000, 0.283000, 0.817000, 0.817000,
                                            0.817000, 0.817000, 0.982000, 0.982000, 0.982000, 0.982000, 1.000000,
                                            1.000000, 1.000000, 1.000000, 0.980000, 0.980000, 0.980000, 0.980000,
                                            0.979000, 0.979000, 0.979000, 0.979000, 0.945000, 0.945000, 0.945000,
                                            0.945000, 0.919000, 0.919000, 0.919000, 0.919000, 0.877000, 0.877000,
                                            0.877000, 0.877000, 0.833000, 0.833000, 0.833000, 0.833000, 0.397000,
                                            0.397000, 0.397000, 0.397000, 0.001000, 0.001000, 0.001000, 0.001000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000]

                    self.nirSpecRespFuncMin = 0.4
                    self.nirSpecRespFuncMax = 0.99
                    self.nirSpecRespFunc = [0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.002000, 0.002000,
                                            0.002000, 0.002000, 0.002000, 0.002000, 0.002000, 0.002000, 0.001000,
                                            0.001000, 0.001000, 0.001000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.001000, 0.001000, 0.001000, 0.001000, 0.002000, 0.002000,
                                            0.002000, 0.002000, 0.001000, 0.001000, 0.001000, 0.001000, 0.001000,
                                            0.001000, 0.001000, 0.001000, 0.002000, 0.002000, 0.002000, 0.002000,
                                            0.002000, 0.002000, 0.002000, 0.002000, 0.002000, 0.002000, 0.002000,
                                            0.002000, 0.001000, 0.001000, 0.001000, 0.001000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.001000, 0.001000, 0.001000, 0.001000, 0.002000, 0.002000, 0.002000,
                                            0.002000, 0.003000, 0.003000, 0.003000, 0.003000, 0.003000, 0.003000,
                                            0.003000, 0.003000, 0.003000, 0.003000, 0.003000, 0.003000, 0.002000,
                                            0.002000, 0.002000, 0.002000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.001000, 0.001000, 0.001000, 0.001000, 0.007000,
                                            0.007000, 0.007000, 0.007000, 0.162000, 0.162000, 0.162000, 0.162000,
                                            0.994000, 0.994000, 0.994000, 0.994000, 1.000000, 1.000000, 1.000000,
                                            1.000000, 0.924000, 0.924000, 0.924000, 0.924000, 0.869000, 0.869000,
                                            0.869000, 0.869000, 0.814000, 0.814000, 0.814000, 0.814000, 0.759000,
                                            0.759000, 0.759000, 0.759000, 0.711000, 0.711000, 0.711000, 0.711000,
                                            0.656000, 0.656000, 0.656000, 0.656000, 0.604000, 0.604000, 0.604000,
                                            0.604000, 0.557000, 0.557000, 0.557000, 0.557000, 0.239000, 0.239000,
                                            0.239000, 0.239000, 0.004000, 0.004000, 0.004000, 0.004000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000]
                elif sensorNameID == "0f" or sensorNameID == "10":
                    self.blueSpecRespFuncMin = 0.32
                    self.blueSpecRespFuncMax = 1.0
                    self.blueSpecRespFunc = [0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.027000, 0.027000, 0.027000, 0.027000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.054000, 0.054000,
                                             0.054000, 0.054000, 0.081000, 0.081000, 0.081000, 0.081000, 0.946000,
                                             0.946000, 0.946000, 0.946000, 1.000000, 1.000000, 1.000000, 1.000000,
                                             0.946000, 0.946000, 0.946000, 0.946000, 0.865000, 0.865000, 0.865000,
                                             0.865000, 0.757000, 0.757000, 0.757000, 0.757000, 0.622000, 0.622000,
                                             0.622000, 0.622000, 0.432000, 0.432000, 0.432000, 0.432000, 0.270000,
                                             0.270000, 0.270000, 0.270000, 0.189000, 0.189000, 0.189000, 0.189000,
                                             0.081000, 0.081000, 0.081000, 0.081000, 0.027000, 0.027000, 0.027000,
                                             0.027000, 0.027000, 0.027000, 0.027000, 0.027000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.027000, 0.027000, 0.027000, 0.027000, 0.027000,
                                             0.027000, 0.027000, 0.027000, 0.027000, 0.027000, 0.027000, 0.027000,
                                             0.027000, 0.027000, 0.027000, 0.027000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000]

                    self.greenSpecRespFuncMin = 0.32
                    self.greenSpecRespFuncMax = 1.0
                    self.greenSpecRespFunc = [0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.029000, 0.029000, 0.029000, 0.029000, 0.206000,
                                              0.206000, 0.206000, 0.206000, 0.353000, 0.353000, 0.353000, 0.353000,
                                              0.441000, 0.441000, 0.441000, 0.441000, 0.471000, 0.471000, 0.471000,
                                              0.471000, 0.588000, 0.588000, 0.588000, 0.588000, 0.765000, 0.765000,
                                              0.765000, 0.765000, 0.912000, 0.912000, 0.912000, 0.912000, 0.971000,
                                              0.971000, 0.971000, 0.971000, 1.000000, 1.000000, 1.000000, 1.000000,
                                              0.941000, 0.941000, 0.941000, 0.941000, 0.853000, 0.853000, 0.853000,
                                              0.853000, 0.765000, 0.765000, 0.765000, 0.765000, 0.647000, 0.647000,
                                              0.647000, 0.647000, 0.500000, 0.500000, 0.500000, 0.500000, 0.353000,
                                              0.353000, 0.353000, 0.353000, 0.176000, 0.176000, 0.176000, 0.176000,
                                              0.118000, 0.118000, 0.118000, 0.118000, 0.059000, 0.059000, 0.059000,
                                              0.059000, 0.059000, 0.059000, 0.059000, 0.059000, 0.029000, 0.029000,
                                              0.029000, 0.029000, 0.029000, 0.029000, 0.029000, 0.029000, 0.059000,
                                              0.059000, 0.059000, 0.059000, 0.059000, 0.059000, 0.059000, 0.059000,
                                              0.088000, 0.088000, 0.088000, 0.088000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000]

                    self.redSpecRespFuncMin = 0.32
                    self.redSpecRespFuncMax = 1.0
                    self.redSpecRespFunc = [0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.074000,
                                            0.074000, 0.074000, 0.074000, 0.111000, 0.111000, 0.111000, 0.111000,
                                            0.148000, 0.148000, 0.148000, 0.148000, 0.111000, 0.111000, 0.111000,
                                            0.111000, 0.148000, 0.148000, 0.148000, 0.148000, 0.148000, 0.148000,
                                            0.148000, 0.148000, 0.111000, 0.111000, 0.111000, 0.111000, 0.111000,
                                            0.111000, 0.111000, 0.111000, 0.111000, 0.111000, 0.111000, 0.111000,
                                            0.111000, 0.111000, 0.111000, 0.111000, 0.111000, 0.111000, 0.111000,
                                            0.111000, 0.111000, 0.111000, 0.111000, 0.111000, 0.185000, 0.185000,
                                            0.185000, 0.185000, 0.481000, 0.481000, 0.481000, 0.481000, 0.852000,
                                            0.852000, 0.852000, 0.852000, 1.000000, 1.000000, 1.000000, 1.000000,
                                            1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
                                            1.000000, 0.926000, 0.926000, 0.926000, 0.926000, 0.889000, 0.889000,
                                            0.889000, 0.889000, 0.852000, 0.852000, 0.852000, 0.852000, 0.815000,
                                            0.815000, 0.815000, 0.815000, 0.704000, 0.704000, 0.704000, 0.704000,
                                            0.519000, 0.519000, 0.519000, 0.519000, 0.037000, 0.037000, 0.037000,
                                            0.037000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000]

                    self.nirSpecRespFuncMin = 0.32
                    self.nirSpecRespFuncMax = 1.0
                    self.nirSpecRespFunc = [0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.100000, 0.100000, 0.100000,
                                            0.100000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.100000, 0.100000, 0.100000, 0.100000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.400000, 0.400000, 0.400000,
                                            0.400000, 0.900000, 0.900000, 0.900000, 0.900000, 1.000000, 1.000000,
                                            1.000000, 1.000000, 0.900000, 0.900000, 0.900000, 0.900000, 0.800000,
                                            0.800000, 0.800000, 0.800000, 0.800000, 0.800000, 0.800000, 0.800000,
                                            0.700000, 0.700000, 0.700000, 0.700000, 0.700000, 0.700000, 0.700000,
                                            0.700000, 0.600000, 0.600000, 0.600000, 0.600000, 0.500000, 0.500000,
                                            0.500000, 0.500000, 0.300000, 0.300000, 0.300000, 0.300000, 0.100000,
                                            0.100000, 0.100000, 0.100000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000]
                elif sensorNameID == "0e":
                    self.blueSpecRespFuncMin = 0.4
                    self.blueSpecRespFuncMax = 0.99
                    self.blueSpecRespFunc = [0.004000, 0.004000, 0.004000, 0.004000, 0.004000, 0.004000, 0.004000,
                                             0.004000, 0.004000, 0.004000, 0.004000, 0.005000, 0.005000, 0.005000,
                                             0.005000, 0.010000, 0.010000, 0.010000, 0.010000, 0.473000, 0.473000,
                                             0.473000, 0.473000, 0.962000, 0.962000, 0.962000, 0.962000, 0.995000,
                                             0.995000, 0.995000, 0.995000, 1.000000, 1.000000, 1.000000, 1.000000,
                                             0.942000, 0.942000, 0.942000, 0.942000, 0.829000, 0.829000, 0.829000,
                                             0.829000, 0.712000, 0.712000, 0.712000, 0.712000, 0.549000, 0.549000,
                                             0.549000, 0.549000, 0.377000, 0.377000, 0.377000, 0.377000, 0.285000,
                                             0.285000, 0.285000, 0.285000, 0.215000, 0.215000, 0.215000, 0.215000,
                                             0.145000, 0.145000, 0.145000, 0.145000, 0.096000, 0.096000, 0.096000,
                                             0.096000, 0.077000, 0.077000, 0.077000, 0.077000, 0.065000, 0.065000,
                                             0.065000, 0.065000, 0.056000, 0.056000, 0.056000, 0.056000, 0.050000,
                                             0.050000, 0.050000, 0.050000, 0.049000, 0.049000, 0.049000, 0.049000,
                                             0.053000, 0.053000, 0.053000, 0.053000, 0.065000, 0.065000, 0.065000,
                                             0.065000, 0.084000, 0.084000, 0.084000, 0.084000, 0.107000, 0.107000,
                                             0.107000, 0.107000, 0.123000, 0.123000, 0.123000, 0.123000, 0.063000,
                                             0.063000, 0.063000, 0.063000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                             0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000]

                    self.greenSpecRespFuncMin = 0.4
                    self.greenSpecRespFuncMax = 0.99
                    self.greenSpecRespFunc = [0.005000, 0.005000, 0.005000, 0.006000, 0.006000, 0.006000, 0.006000,
                                              0.004000, 0.004000, 0.004000, 0.004000, 0.005000, 0.005000, 0.005000,
                                              0.005000, 0.009000, 0.009000, 0.009000, 0.009000, 0.038000, 0.038000,
                                              0.038000, 0.038000, 0.187000, 0.187000, 0.187000, 0.187000, 0.258000,
                                              0.258000, 0.258000, 0.258000, 0.378000, 0.378000, 0.378000, 0.378000,
                                              0.531000, 0.531000, 0.531000, 0.531000, 0.675000, 0.675000, 0.675000,
                                              0.675000, 0.826000, 0.826000, 0.826000, 0.826000, 0.934000, 0.934000,
                                              0.934000, 0.934000, 0.988000, 0.988000, 0.988000, 0.988000, 1.000000,
                                              1.000000, 1.000000, 1.000000, 0.956000, 0.956000, 0.956000, 0.956000,
                                              0.892000, 0.892000, 0.892000, 0.892000, 0.804000, 0.804000, 0.804000,
                                              0.804000, 0.687000, 0.687000, 0.687000, 0.687000, 0.538000, 0.538000,
                                              0.538000, 0.538000, 0.376000, 0.376000, 0.376000, 0.376000, 0.241000,
                                              0.241000, 0.241000, 0.241000, 0.163000, 0.163000, 0.163000, 0.163000,
                                              0.139000, 0.139000, 0.139000, 0.139000, 0.121000, 0.121000, 0.121000,
                                              0.121000, 0.107000, 0.107000, 0.107000, 0.107000, 0.106000, 0.106000,
                                              0.106000, 0.106000, 0.122000, 0.122000, 0.122000, 0.122000, 0.150000,
                                              0.150000, 0.150000, 0.150000, 0.133000, 0.133000, 0.133000, 0.133000,
                                              0.026000, 0.026000, 0.026000, 0.026000, 0.002000, 0.002000, 0.002000,
                                              0.002000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.001000, 0.001000,
                                              0.001000, 0.001000, 0.006000, 0.006000, 0.006000, 0.006000, 0.005000,
                                              0.005000, 0.005000, 0.005000, 0.003000, 0.003000, 0.003000, 0.003000,
                                              0.004000, 0.004000, 0.004000, 0.004000, 0.007000, 0.007000, 0.007000,
                                              0.007000, 0.011000, 0.011000, 0.011000, 0.011000, 0.013000, 0.013000,
                                              0.013000, 0.013000, 0.011000, 0.011000, 0.011000, 0.011000, 0.008000,
                                              0.008000, 0.008000, 0.008000, 0.006000, 0.006000, 0.006000, 0.006000,
                                              0.005000, 0.005000, 0.005000, 0.005000, 0.004000, 0.004000, 0.004000,
                                              0.004000, 0.004000, 0.004000, 0.004000, 0.004000, 0.005000, 0.005000,
                                              0.005000, 0.005000, 0.005000, 0.005000, 0.005000, 0.005000, 0.004000,
                                              0.004000, 0.004000, 0.004000, 0.003000, 0.003000, 0.003000, 0.003000,
                                              0.003000, 0.003000, 0.003000, 0.003000, 0.002000, 0.002000]

                    self.redSpecRespFuncMin = 0.4
                    self.redSpecRespFuncMax = 0.99
                    self.redSpecRespFunc = [0.004000, 0.004000, 0.004000, 0.004000, 0.004000, 0.004000, 0.004000,
                                            0.001000, 0.001000, 0.001000, 0.001000, 0.002000, 0.002000, 0.002000,
                                            0.002000, 0.002000, 0.002000, 0.002000, 0.002000, 0.008000, 0.008000,
                                            0.008000, 0.008000, 0.027000, 0.027000, 0.027000, 0.027000, 0.029000,
                                            0.029000, 0.029000, 0.029000, 0.031000, 0.031000, 0.031000, 0.031000,
                                            0.030000, 0.030000, 0.030000, 0.030000, 0.027000, 0.027000, 0.027000,
                                            0.027000, 0.032000, 0.032000, 0.032000, 0.032000, 0.048000, 0.048000,
                                            0.048000, 0.048000, 0.062000, 0.062000, 0.062000, 0.062000, 0.042000,
                                            0.042000, 0.042000, 0.042000, 0.022000, 0.022000, 0.022000, 0.022000,
                                            0.023000, 0.023000, 0.023000, 0.023000, 0.066000, 0.066000, 0.066000,
                                            0.066000, 0.384000, 0.384000, 0.384000, 0.384000, 0.864000, 0.864000,
                                            0.864000, 0.864000, 0.980000, 0.980000, 0.980000, 0.980000, 1.000000,
                                            1.000000, 1.000000, 1.000000, 0.986000, 0.986000, 0.986000, 0.986000,
                                            0.978000, 0.978000, 0.978000, 0.978000, 0.944000, 0.944000, 0.944000,
                                            0.944000, 0.910000, 0.910000, 0.910000, 0.910000, 0.880000, 0.880000,
                                            0.880000, 0.880000, 0.833000, 0.833000, 0.833000, 0.833000, 0.796000,
                                            0.796000, 0.796000, 0.796000, 0.551000, 0.551000, 0.551000, 0.551000,
                                            0.088000, 0.088000, 0.088000, 0.088000, 0.007000, 0.007000, 0.007000,
                                            0.007000, 0.001000, 0.001000, 0.001000, 0.001000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.002000, 0.002000,
                                            0.002000, 0.002000, 0.009000, 0.009000, 0.009000, 0.009000, 0.006000,
                                            0.006000, 0.006000, 0.006000, 0.004000, 0.004000, 0.004000, 0.004000,
                                            0.005000, 0.005000, 0.005000, 0.005000, 0.007000, 0.007000, 0.007000,
                                            0.007000, 0.012000, 0.012000, 0.012000, 0.012000, 0.015000, 0.015000,
                                            0.015000, 0.015000, 0.012000, 0.012000, 0.012000, 0.012000, 0.009000,
                                            0.009000, 0.009000, 0.009000, 0.006000, 0.006000, 0.006000, 0.006000,
                                            0.005000, 0.005000, 0.005000, 0.005000, 0.005000, 0.005000, 0.005000,
                                            0.005000, 0.005000, 0.005000, 0.005000, 0.005000, 0.005000, 0.005000,
                                            0.005000, 0.005000, 0.005000, 0.005000, 0.005000, 0.005000, 0.005000,
                                            0.005000, 0.005000, 0.005000, 0.004000, 0.004000, 0.004000, 0.004000,
                                            0.003000, 0.003000, 0.003000, 0.003000, 0.002000, 0.002000]

                    self.nirSpecRespFuncMin = 0.4
                    self.nirSpecRespFuncMax = 0.99
                    self.nirSpecRespFunc = [0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.002000, 0.002000, 0.002000,
                                            0.002000, 0.000000, 0.000000, 0.000000, 0.000000, 0.001000, 0.001000,
                                            0.001000, 0.001000, 0.001000, 0.001000, 0.001000, 0.001000, 0.002000,
                                            0.002000, 0.002000, 0.002000, 0.005000, 0.005000, 0.005000, 0.005000,
                                            0.001000, 0.001000, 0.001000, 0.001000, 0.001000, 0.001000, 0.001000,
                                            0.001000, 0.002000, 0.002000, 0.002000, 0.002000, 0.002000, 0.002000,
                                            0.002000, 0.002000, 0.001000, 0.001000, 0.001000, 0.001000, 0.002000,
                                            0.002000, 0.002000, 0.002000, 0.002000, 0.002000, 0.002000, 0.002000,
                                            0.002000, 0.002000, 0.002000, 0.002000, 0.003000, 0.003000, 0.003000,
                                            0.003000, 0.004000, 0.004000, 0.004000, 0.004000, 0.006000, 0.006000,
                                            0.006000, 0.006000, 0.002000, 0.002000, 0.002000, 0.002000, 0.002000,
                                            0.002000, 0.002000, 0.002000, 0.005000, 0.005000, 0.005000, 0.005000,
                                            0.004000, 0.004000, 0.004000, 0.004000, 0.002000, 0.002000, 0.002000,
                                            0.002000, 0.003000, 0.003000, 0.003000, 0.003000, 0.003000, 0.003000,
                                            0.003000, 0.003000, 0.014000, 0.014000, 0.014000, 0.014000, 0.007000,
                                            0.007000, 0.007000, 0.007000, 0.059000, 0.059000, 0.059000, 0.059000,
                                            0.014000, 0.014000, 0.014000, 0.014000, 0.002000, 0.002000, 0.002000,
                                            0.002000, 0.007000, 0.007000, 0.007000, 0.007000, 0.001000, 0.001000,
                                            0.001000, 0.001000, 0.001000, 0.001000, 0.001000, 0.001000, 0.004000,
                                            0.004000, 0.004000, 0.004000, 0.139000, 0.139000, 0.139000, 0.139000,
                                            0.822000, 0.822000, 0.822000, 0.822000, 1.000000, 1.000000, 1.000000,
                                            1.000000, 0.960000, 0.960000, 0.960000, 0.960000, 0.918000, 0.918000,
                                            0.918000, 0.918000, 0.849000, 0.849000, 0.849000, 0.849000, 0.800000,
                                            0.800000, 0.800000, 0.800000, 0.761000, 0.761000, 0.761000, 0.761000,
                                            0.702000, 0.702000, 0.702000, 0.702000, 0.652000, 0.652000, 0.652000,
                                            0.652000, 0.609000, 0.609000, 0.609000, 0.609000, 0.522000, 0.522000,
                                            0.522000, 0.522000, 0.133000, 0.133000, 0.133000, 0.133000, 0.008000,
                                            0.008000, 0.008000, 0.008000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                                            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000]

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
        if nBands != self.numOfBands:
            raise ARCSIException('Input Image \'' + self.fileName + '\' does not match the number of image bands specified in header.')
        rasterDS = None

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
            bandDefnSeq = list()
            imgBand = collections.namedtuple('ImgBand', ['bandName', 'bandIndex', 'bias', 'gain'])
            bandDefnSeq.append(imgBand(bandName="Blue", bandIndex=1, bias=0.0, gain=self.radioCoefficentsDict[1].radCoef))
            bandDefnSeq.append(imgBand(bandName="Green", bandIndex=2, bias=0.0, gain=self.radioCoefficentsDict[2].radCoef))
            bandDefnSeq.append(imgBand(bandName="Red", bandIndex=3, bias=0.0, gain=self.radioCoefficentsDict[3].radCoef))
            if self.numOfBands == 4:
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

        reBand = collections.namedtuple('REBand', ['bandName', 'fileName', 'bandIndex', 'satVal'])
        bandDefnSeq = list()
        bandDefnSeq.append(reBand(bandName="Blue", fileName=self.fileName, bandIndex=1, satVal=65535.0))
        bandDefnSeq.append(reBand(bandName="Green", fileName=self.fileName, bandIndex=2, satVal=65535.0))
        bandDefnSeq.append(reBand(bandName="Red", fileName=self.fileName, bandIndex=3, satVal=65535.0))
        if self.numOfBands == 4:
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
        if self.numOfBands == 4:
            b4_rad = dn * self.radioCoefficentsDict[4].radCoef
            b4_toa = dn * self.radioCoefficentsDict[4].reflCoef
            b4_solar_irr = arcsiUtils.getESUNValue(b4_rad, b4_toa, self.acquisitionTime.day, self.acquisitionTime.month,
                                                   self.acquisitionTime.year, self.solarZenith)
            print("b4_solar_irr = {}".format(b4_solar_irr))
            solarIrradianceVals.append(IrrVal(irradiance=b4_solar_irr))

        rsgislib.imagecalibration.radiance2TOARefl(inputRadImage, outputImage, outFormat, rsgislib.TYPE_16UINT, scaleFactor, self.acquisitionTime.year, self.acquisitionTime.month, self.acquisitionTime.day, self.solarZenith, solarIrradianceVals)
        return outputImage

    def generateCloudMask(self, inputReflImage, inputSatImage, inputThermalImage, inputViewAngleImg, inputValidImg, outputPath, outputName, outFormat, tmpPath, scaleFactor, cloud_msk_methods=None):
        raise ARCSIException("Cloud Masking Not Implemented for PlanetScope.")

    def createCloudMaskDataArray(self, inImgDataArr):
        # Calc Whiteness
        meanArr = numpy.mean(inImgDataArr, axis=1)
        whitenessArr = numpy.absolute((inImgDataArr[...,0] - meanArr)/meanArr) + numpy.absolute((inImgDataArr[...,1] - meanArr)/meanArr) + numpy.absolute((inImgDataArr[...,2] - meanArr)/meanArr) + numpy.absolute((inImgDataArr[...,3] - meanArr)/meanArr)

        # Create and populate the output array.
        inShape = inImgDataArr.shape
        outShape = [inShape[0], inShape[1]+2]
        outArr = numpy.zeros(outShape, dtype=float)
        
        for i in range(inShape[1]):
            outArr[...,i] = inImgDataArr[...,i]
        
        idx = inShape[1]
        outArr[...,idx] = meanArr
        outArr[...,idx+1] = whitenessArr

        return outArr

    def defineDarkShadowImageBand(self):
        darkBand = 3
        if self.numOfBands == 4:
            darkBand = 4
        return darkBand

    def calc6SCoefficients(self, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF):
        sixsCoeffs = None
        if self.numOfBands == 3:
            sixsCoeffs = numpy.zeros((3, 6), dtype=numpy.float32)
        elif self.numOfBands == 4:
            sixsCoeffs = numpy.zeros((4, 6), dtype=numpy.float32)
        else:
            raise ARCSIException("calc6SCoefficients: Expecting either 3 or 4 image bands.")
        # Set up 6S model
        s = Py6S.SixS()
        s.atmos_profile = atmosProfile
        s.aero_profile = aeroProfile
        # s.ground_reflectance = Py6S.GroundReflectance.HomogeneousHapke(0.101, -0.263, 0.589, 0.046)
        s.ground_reflectance = grdRefl
        s.geometry = Py6S.Geometry.User()
        s.geometry.solar_z = self.solarZenith
        s.geometry.solar_a = self.solarAzimuth
        s.geometry.view_z = self.sensorZenith
        s.geometry.view_a = self.sensorAzimuth
        s.geometry.month = self.acquisitionTime.month
        s.geometry.day = self.acquisitionTime.day
        s.geometry.gmt_decimal_hour = float(self.acquisitionTime.hour) + float(self.acquisitionTime.minute) / 60.0
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
        s.wavelength = Py6S.Wavelength(self.blueSpecRespFuncMin, self.blueSpecRespFuncMax, self.blueSpecRespFunc)
        s.run()
        sixsCoeffs[0, 0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[0, 1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[0, 2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[0, 3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[0, 4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[0, 5] = float(s.outputs.values['environmental_irradiance'])

        # Band 2
        s.wavelength = Py6S.Wavelength(self.greenSpecRespFuncMin, self.greenSpecRespFuncMax, self.greenSpecRespFunc)
        s.run()
        sixsCoeffs[1, 0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[1, 1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[1, 2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[1, 3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[1, 4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[1, 5] = float(s.outputs.values['environmental_irradiance'])

        # Band 3
        s.wavelength = Py6S.Wavelength(self.redSpecRespFuncMin, self.redSpecRespFuncMax, self.redSpecRespFunc)
        s.run()
        sixsCoeffs[2, 0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[2, 1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[2, 2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[2, 3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[2, 4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[2, 5] = float(s.outputs.values['environmental_irradiance'])

        if self.numOfBands == 4:
            # Band 4
            s.wavelength = Py6S.Wavelength(self.nirSpecRespFuncMin, self.nirSpecRespFuncMax, self.nirSpecRespFunc)

            s.run()
            sixsCoeffs[3, 0] = float(s.outputs.values['coef_xa'])
            sixsCoeffs[3, 1] = float(s.outputs.values['coef_xb'])
            sixsCoeffs[3, 2] = float(s.outputs.values['coef_xc'])
            sixsCoeffs[3, 3] = float(s.outputs.values['direct_solar_irradiance'])
            sixsCoeffs[3, 4] = float(s.outputs.values['diffuse_solar_irradiance'])
            sixsCoeffs[3, 5] = float(s.outputs.values['environmental_irradiance'])

        return sixsCoeffs

    def convertImageToSurfaceReflSglParam(self, inputRadImage, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF, scaleFactor):
        print("Converting to Surface Reflectance")
        outputImage = os.path.join(outputPath, outputName)

        imgBandCoeffs = list()

        sixsCoeffs = self.calc6SCoefficients(aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF)

        imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=1, aX=float(sixsCoeffs[0,0]), bX=float(sixsCoeffs[0,1]), cX=float(sixsCoeffs[0,2]), DirIrr=float(sixsCoeffs[0,3]), DifIrr=float(sixsCoeffs[0,4]), EnvIrr=float(sixsCoeffs[0,5])))
        imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=2, aX=float(sixsCoeffs[1,0]), bX=float(sixsCoeffs[1,1]), cX=float(sixsCoeffs[1,2]), DirIrr=float(sixsCoeffs[1,3]), DifIrr=float(sixsCoeffs[1,4]), EnvIrr=float(sixsCoeffs[1,5])))
        imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=3, aX=float(sixsCoeffs[2,0]), bX=float(sixsCoeffs[2,1]), cX=float(sixsCoeffs[2,2]), DirIrr=float(sixsCoeffs[2,3]), DifIrr=float(sixsCoeffs[2,4]), EnvIrr=float(sixsCoeffs[2,5])))
        if self.numOfBands == 4:
            imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=4, aX=float(sixsCoeffs[3,0]), bX=float(sixsCoeffs[3,1]), cX=float(sixsCoeffs[3,2]), DirIrr=float(sixsCoeffs[3,3]), DifIrr=float(sixsCoeffs[3,4]), EnvIrr=float(sixsCoeffs[3,5])))

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
                if self.numOfBands == 4:
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
                    if self.numOfBands == 4:
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
        s.wavelength = Py6S.Wavelength(self.blueSpecRespFuncMin, self.blueSpecRespFuncMax, self.blueSpecRespFunc)
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
        raise ARCSIException("For the PlanetScope Sensor there is not implementation of findDDVTargets")

    def estimateImageToAODUsingDDV(self, inputRADImage, inputTOAImage, inputDEMFile, shadowMask, outputPath, outputName, outFormat, tmpPath, aeroProfile, atmosProfile, grdRefl, aotValMin, aotValMax):
        raise ARCSIException("For the PlanetScope Sensor there is not implementation of estimateImageToAODUsingDDV")

    def estimateImageToAODUsingDOS(self, inputRADImage, inputTOAImage, inputDEMFile, shadowMask, outputPath, outputName, outFormat, tmpPath, aeroProfile, atmosProfile, grdRefl, aotValMin, aotValMax, globalDOS, simpleDOS, dosOutRefl):
        raise ARCSIException("For the PlanetScope Sensor there is not implementation of estimateImageToAODUsingDOS")

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
            if self.numOfBands == 4:
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


