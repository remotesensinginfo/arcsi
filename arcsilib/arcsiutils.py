"""
Module that contains the ARCSIUtils class.
"""
############################################################################
#  arcsiutils.py
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
# Purpose:  A class with some useful utilites.
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


import numpy
import scipy.interpolate

from .arcsiexception import ARCSIException


def ARCSIEnum(*sequential, **named):
    """Handy way to fake an enumerated type in Python
    http://stackoverflow.com/questions/36932/how-can-i-represent-an-enum-in-python
    """
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type("ARCSIEnum", (), enums)


def readSpectralResponseFunc(inFile, seperator, ignoreLines, waveCol, respCol):
    specResp = list()
    try:
        specFile = open(inFile, "r")
        lineCount = 0
        for line in specFile:
            if lineCount >= ignoreLines:
                line = line.strip()
                if line:
                    lineVals = line.split(seperator)
                    if (len(lineVals) <= waveCol) or (len(lineVals) <= respCol):
                        raise ARCSIException("")
                    waveVal = float(lineVals[waveCol].strip())
                    respVal = float(lineVals[respCol].strip())
                    specResp.append([waveVal, respVal])
            lineCount += 1
        specFile.close()
    except ARCSIException as e:
        raise e
    except Exception as e:
        raise e
    return numpy.array(specResp)


def resampleSpectralResponseFunc(wvlens, respFunc, outSamp, sampleMethod):
    """
    Specifies the kind of interpolation as a string
    Options: 'linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic'
    Where, 'zero', 'slinear', 'quadratic' and 'cubic' refer to a spline
    interpolation of zeroth, first, second or third order) or as an integer
    specifying the order of the spline interpolator to use. Default is ‘linear’

    See scipy documentation for more information:
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html
    """
    resamFunc = scipy.interpolate.interp1d(
        wvlens,
        respFunc,
        kind=sampleMethod,
        axis=-1,
        copy=True,
        bounds_error=False,
        fill_value=0,
        assume_sorted=True,
    )
    oWVLens = numpy.arange(wvlens[0], wvlens[-1], outSamp)
    oSpecResp = resamFunc(oWVLens)
    return oWVLens, oSpecResp


def findMinimumElev(elev):
    elevVal = -550
    outElev = 0
    for i in range(100):
        if (elev > elevVal) & (elev < (elevVal + 100)):
            outElev = elevVal
            break
        elevVal = elevVal + 100
    return outElev


def findMaximumElev(elev):
    elevVal = -550
    outElev = 0
    for i in range(105):
        if (elev > elevVal) & (elev < (elevVal + 100)):
            outElev = elevVal + 100
            break
        elevVal = elevVal + 100
    return outElev


def findMinimumAOT(aot):
    aotVal = 0
    outAOT = 0
    for i in range(200):
        if (aot > aotVal) & (aot < (aotVal + 0.05)):
            outAOT = aotVal
            break
        aotVal = aotVal + 0.05
    return aotVal


def findMaximumAOT(aot):
    aotVal = 0
    outAOT = 0
    for i in range(200):
        if (aot > aotVal) & (aot < (aotVal + 0.05)):
            outAOT = aotVal + 0.05
            break
        aotVal = aotVal + 0.05
    return aotVal


class ARCSILandsatMetaUtils(object):
    """
    A class with common functions for parsing Landsat
    metadata
    """

    @staticmethod
    def getGeographicCorners(headerParams):
        """
        Function to get geographic coordinates of image from metatdata

        Returns array containing:

        * UL_LON
        * UL_LAT
        * UR_LON
        * UR_LAT
        * LL_LAT
        * LL_LON
        * LR_LAT
        * LR_LON

        """
        outCornerCoords = []
        geoVarList = [
            ["UL", "LAT"],
            ["UL", "LON"],
            ["UR", "LAT"],
            ["UR", "LON"],
            ["LL", "LAT"],
            ["LL", "LON"],
            ["LR", "LAT"],
            ["LR", "LON"],
        ]

        for geoItem in geoVarList:
            try:
                outCornerCoords.append(
                    float(
                        headerParams[
                            "CORNER_{0}_{1}_PRODUCT".format(geoItem[0], geoItem[1])
                        ]
                    )
                )
            except KeyError:
                outCornerCoords.append(
                    float(
                        headerParams[
                            "PRODUCT_{0}_CORNER_{1}".format(geoItem[0], geoItem[1])
                        ]
                    )
                )

        return outCornerCoords

    @staticmethod
    def getProjectedCorners(headerParams):
        """
        Function to get projected coordinates of image from metatdata

        Returns array containing:

        * UL_X
        * UL_Y
        * UR_X
        * UR_Y
        * LL_X
        * LL_Y
        * LR_X
        * LR_Y

        """
        outCornerCoords = []
        projectedVarList = [
            ["UL", "X"],
            ["UL", "Y"],
            ["UR", "X"],
            ["UR", "Y"],
            ["LL", "X"],
            ["LL", "Y"],
            ["LR", "X"],
            ["LR", "Y"],
        ]

        for projectedItem in projectedVarList:
            try:
                outCornerCoords.append(
                    float(
                        headerParams[
                            "CORNER_{0}_PROJECTION_{1}_PRODUCT".format(
                                projectedItem[0], projectedItem[1]
                            )
                        ]
                    )
                )
            except KeyError:
                outCornerCoords.append(
                    float(
                        headerParams[
                            "PRODUCT_{0}_CORNER_MAP{1}".format(
                                projectedItem[0], projectedItem[1]
                            )
                        ]
                    )
                )

        return outCornerCoords

    @staticmethod
    def getBandFilenames(headerParams, nBands):
        """
        Get filenames for individual bands

        Returns a list with a name for each band.
        """
        metaFilenames = []

        for i in range(1, nBands + 1):
            try:
                metaFilenames.append(headerParams["FILE_NAME_BAND_{}".format(i)])
            except KeyError:
                try:
                    metaFilenames.append(headerParams["BAND{}_FILE_NAME".format(i)])
                # For Landsat 7 ETM+ There are two band 6 files.
                # Just set to 'None' here and fetch separately.
                except Exception:
                    if i == 6:
                        metaFilenames.append(None)
                    else:
                        raise

        return metaFilenames


class ARCSISensorFactory(object):
    def getSensorClassFromName(self, sensor, debugMode, inputImage):
        sensorClass = None
        if sensor == "lsetm":
            from arcsilib.arcsisensorlandsat_etm import ARCSILandsatETMSensor

            sensorClass = ARCSILandsatETMSensor(debugMode, inputImage)
        elif sensor == "lstm":
            from arcsilib.arcsisensorlandsat_tm import ARCSILandsatTMSensor

            sensorClass = ARCSILandsatTMSensor(debugMode, inputImage)
        elif sensor == "lsmss":
            from arcsilib.arcsisensorlandsat_mss import ARCSILandsatMSSSensor

            sensorClass = ARCSILandsatMSSSensor(debugMode, inputImage)
        elif sensor == "lsoli":
            from arcsilib.arcsisensorlandsat_oli import ARCSILandsatOLISensor

            sensorClass = ARCSILandsatOLISensor(debugMode, inputImage)
        elif sensor == "sen2":
            from arcsilib.arcsisensorsentinel2 import ARCSISentinel2Sensor

            sensorClass = ARCSISentinel2Sensor(debugMode, inputImage)
        else:
            raise ARCSIException(
                "Could not get a class representing the sensor specified from the factory."
            )

        if sensorClass == None:
            raise ARCSIException("Something has gone wrong sensorClass is None!")

        return sensorClass
