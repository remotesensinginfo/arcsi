"""
Module that contains the ARCSIAbstractSensor class.
"""
############################################################################
#  arcsisensor.py
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
# Purpose:  An abstract class which provides the base class for each
#           sensor supported by the ARCSI system.
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

import collections
import datetime
import json
import math
import os
import shutil
import sys
from abc import ABCMeta, abstractmethod

import numpy
import rsgislib
import rsgislib.elevation
import rsgislib.imagecalc
import rsgislib.imagecalibration
import rsgislib.imagefilter
import rsgislib.imageutils
import rsgislib.rastergis
import rsgislib.segmentation
import rsgislib.tools.geometrytools
import rsgislib.tools.utils
import scipy.interpolate
import scipy.interpolate.rbf
from osgeo import gdal, ogr, osr
from rios import rat
from rios.imagereader import ImageReader
from rios.imagewriter import ImageWriter

from arcsilib import ARCSI_VERSION, ARCSI_WEBSITE

from .arcsiexception import ARCSIException


class ARCSIAbstractSensor(object):
    """
    An abstract class which represents a sensor and allows
    the various opperations required to be applied and standard
    variables (e.g., acqusiation date) stored and retrieved.
    """

    __metaclass__ = ABCMeta

    def __init__(self, debugMode, inputImage):
        self.sensor = "NA"
        self.headerFileName = ""
        self.userSpInputImage = inputImage
        self.debugMode = debugMode
        self.acquisitionTime = datetime.datetime.today()
        self.latTL = 0.0
        self.lonTL = 0.0
        self.latTR = 0.0
        self.lonTR = 0.0
        self.latBL = 0.0
        self.lonBL = 0.0
        self.latBR = 0.0
        self.lonBR = 0.0
        self.latCentre = 0.0
        self.lonCentre = 0.0
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
        self.inWKT = ""
        self.reprojectOutputs = False
        self.solarZenith = 0.0
        self.solarAzimuth = 0.0
        self.sensorZenith = 0.0
        self.sensorAzimuth = 0.0
        self.epsgCodes = dict()
        self.epsgCodes["WGS84UTM1N"] = 32601
        self.epsgCodes["WGS84UTM2N"] = 32602
        self.epsgCodes["WGS84UTM3N"] = 32603
        self.epsgCodes["WGS84UTM4N"] = 32604
        self.epsgCodes["WGS84UTM5N"] = 32605
        self.epsgCodes["WGS84UTM6N"] = 32606
        self.epsgCodes["WGS84UTM7N"] = 32607
        self.epsgCodes["WGS84UTM8N"] = 32608
        self.epsgCodes["WGS84UTM9N"] = 32609
        self.epsgCodes["WGS84UTM10N"] = 32610
        self.epsgCodes["WGS84UTM11N"] = 32611
        self.epsgCodes["WGS84UTM12N"] = 32612
        self.epsgCodes["WGS84UTM13N"] = 32613
        self.epsgCodes["WGS84UTM14N"] = 32614
        self.epsgCodes["WGS84UTM15N"] = 32615
        self.epsgCodes["WGS84UTM16N"] = 32616
        self.epsgCodes["WGS84UTM17N"] = 32617
        self.epsgCodes["WGS84UTM18N"] = 32618
        self.epsgCodes["WGS84UTM19N"] = 32619
        self.epsgCodes["WGS84UTM20N"] = 32620
        self.epsgCodes["WGS84UTM21N"] = 32621
        self.epsgCodes["WGS84UTM22N"] = 32622
        self.epsgCodes["WGS84UTM23N"] = 32623
        self.epsgCodes["WGS84UTM24N"] = 32624
        self.epsgCodes["WGS84UTM25N"] = 32625
        self.epsgCodes["WGS84UTM26N"] = 32626
        self.epsgCodes["WGS84UTM27N"] = 32627
        self.epsgCodes["WGS84UTM28N"] = 32628
        self.epsgCodes["WGS84UTM29N"] = 32629
        self.epsgCodes["WGS84UTM30N"] = 32630
        self.epsgCodes["WGS84UTM31N"] = 32631
        self.epsgCodes["WGS84UTM32N"] = 32632
        self.epsgCodes["WGS84UTM33N"] = 32633
        self.epsgCodes["WGS84UTM34N"] = 32634
        self.epsgCodes["WGS84UTM35N"] = 32635
        self.epsgCodes["WGS84UTM36N"] = 32636
        self.epsgCodes["WGS84UTM37N"] = 32637
        self.epsgCodes["WGS84UTM38N"] = 32638
        self.epsgCodes["WGS84UTM39N"] = 32639
        self.epsgCodes["WGS84UTM40N"] = 32640
        self.epsgCodes["WGS84UTM41N"] = 32641
        self.epsgCodes["WGS84UTM42N"] = 32642
        self.epsgCodes["WGS84UTM43N"] = 32643
        self.epsgCodes["WGS84UTM44N"] = 32644
        self.epsgCodes["WGS84UTM45N"] = 32645
        self.epsgCodes["WGS84UTM46N"] = 32646
        self.epsgCodes["WGS84UTM47N"] = 32647
        self.epsgCodes["WGS84UTM48N"] = 32648
        self.epsgCodes["WGS84UTM49N"] = 32649
        self.epsgCodes["WGS84UTM50N"] = 32650
        self.epsgCodes["WGS84UTM51N"] = 32651
        self.epsgCodes["WGS84UTM52N"] = 32652
        self.epsgCodes["WGS84UTM53N"] = 32653
        self.epsgCodes["WGS84UTM54N"] = 32654
        self.epsgCodes["WGS84UTM55N"] = 32655
        self.epsgCodes["WGS84UTM56N"] = 32656
        self.epsgCodes["WGS84UTM57N"] = 32657
        self.epsgCodes["WGS84UTM58N"] = 32658
        self.epsgCodes["WGS84UTM59N"] = 32659
        self.epsgCodes["WGS84UTM60N"] = 32660
        self.epsgCodes["WGS84UTM1S"] = 32701
        self.epsgCodes["WGS84UTM2S"] = 32702
        self.epsgCodes["WGS84UTM3S"] = 32703
        self.epsgCodes["WGS84UTM4S"] = 32704
        self.epsgCodes["WGS84UTM5S"] = 32705
        self.epsgCodes["WGS84UTM6S"] = 32706
        self.epsgCodes["WGS84UTM7S"] = 32707
        self.epsgCodes["WGS84UTM8S"] = 32708
        self.epsgCodes["WGS84UTM9S"] = 32709
        self.epsgCodes["WGS84UTM10S"] = 32710
        self.epsgCodes["WGS84UTM11S"] = 32711
        self.epsgCodes["WGS84UTM12S"] = 32712
        self.epsgCodes["WGS84UTM13S"] = 32713
        self.epsgCodes["WGS84UTM14S"] = 32714
        self.epsgCodes["WGS84UTM15S"] = 32715
        self.epsgCodes["WGS84UTM16S"] = 32716
        self.epsgCodes["WGS84UTM17S"] = 32717
        self.epsgCodes["WGS84UTM18S"] = 32718
        self.epsgCodes["WGS84UTM19S"] = 32719
        self.epsgCodes["WGS84UTM20S"] = 32720
        self.epsgCodes["WGS84UTM21S"] = 32721
        self.epsgCodes["WGS84UTM22S"] = 32722
        self.epsgCodes["WGS84UTM23S"] = 32723
        self.epsgCodes["WGS84UTM24S"] = 32724
        self.epsgCodes["WGS84UTM25S"] = 32725
        self.epsgCodes["WGS84UTM26S"] = 32726
        self.epsgCodes["WGS84UTM27S"] = 32727
        self.epsgCodes["WGS84UTM28S"] = 32728
        self.epsgCodes["WGS84UTM29S"] = 32729
        self.epsgCodes["WGS84UTM30S"] = 32730
        self.epsgCodes["WGS84UTM31S"] = 32731
        self.epsgCodes["WGS84UTM32S"] = 32732
        self.epsgCodes["WGS84UTM33S"] = 32733
        self.epsgCodes["WGS84UTM34S"] = 32734
        self.epsgCodes["WGS84UTM35S"] = 32735
        self.epsgCodes["WGS84UTM36S"] = 32736
        self.epsgCodes["WGS84UTM37S"] = 32737
        self.epsgCodes["WGS84UTM38S"] = 32738
        self.epsgCodes["WGS84UTM39S"] = 32739
        self.epsgCodes["WGS84UTM40S"] = 32740
        self.epsgCodes["WGS84UTM41S"] = 32741
        self.epsgCodes["WGS84UTM42S"] = 32742
        self.epsgCodes["WGS84UTM43S"] = 32743
        self.epsgCodes["WGS84UTM44S"] = 32744
        self.epsgCodes["WGS84UTM45S"] = 32745
        self.epsgCodes["WGS84UTM46S"] = 32746
        self.epsgCodes["WGS84UTM47S"] = 32747
        self.epsgCodes["WGS84UTM48S"] = 32748
        self.epsgCodes["WGS84UTM49S"] = 32749
        self.epsgCodes["WGS84UTM50S"] = 32750
        self.epsgCodes["WGS84UTM51S"] = 32751
        self.epsgCodes["WGS84UTM52S"] = 32752
        self.epsgCodes["WGS84UTM53S"] = 32753
        self.epsgCodes["WGS84UTM54S"] = 32754
        self.epsgCodes["WGS84UTM55S"] = 32755
        self.epsgCodes["WGS84UTM56S"] = 32756
        self.epsgCodes["WGS84UTM57S"] = 32757
        self.epsgCodes["WGS84UTM58S"] = 32758
        self.epsgCodes["WGS84UTM59S"] = 32759
        self.epsgCodes["WGS84UTM60S"] = 32760

    @abstractmethod
    def extractHeaderParameters(self, inputHeader, wktStr):
        pass

    def setReProjectOutputs(self, reproj=False):
        self.reprojectOutputs = reproj

    def getReProjectOutputs(self, reproj=False):
        return self.reprojectOutputs

    @abstractmethod
    def getSolarIrrStdSolarGeom(self):
        pass

    @abstractmethod
    def getSensorViewGeom(self):
        pass

    def checkInputImageValid(self):
        if not self.expectedImageDataPresent():
            raise ARCSIException("Not all of the image(s) are present.")

    def defaultGenBaseOutFileName(self):
        """
        A function to generate a generic standard file
        base name which will be sensible.

        It is expected that individual sensors may override this function.
        """
        date = self.acquisitionTime.strftime("%Y%m%d")

        east_west = "e"
        if self.lonCentre < 0:
            east_west = "w"
        north_south = "n"
        if self.latCentre < 0:
            north_south = "s"

        # pad lat and lons to a standardised format
        lat_pad = "{:04.1f}".format(round(abs(self.latCentre), 1))
        lon_pad = "{:05.1f}".format(round(abs(self.lonCentre), 1))

        pos = (
            "lat"
            + north_south
            + lat_pad.replace(".", "")
            + "lon"
            + east_west
            + lon_pad.replace(".", "")
        )

        outname = self.sensor + "_" + date + "_" + pos

        return outname

    def generateOutputBaseName(self):
        """
        Provides a default implementation for generating file name.
        """
        return self.defaultGenBaseOutFileName()

    def getJSONDictDefaultMetaData(
        self,
        productsStr,
        validMaskImage="",
        footprintCalc=False,
        calcdValuesDict=None,
        outFilesDict=None,
    ):
        """ """
        if outFilesDict is None:
            outFilesDict = dict()
        if calcdValuesDict is None:
            calcdValuesDict = dict()

        softwareDict = dict()
        softwareDict["Name"] = "ARCSI"
        softwareDict["URL"] = ARCSI_WEBSITE
        softwareDict["Version"] = ARCSI_VERSION

        productsDict = dict()
        productsDict["ARCSIProducts"] = productsStr
        nowDateTime = datetime.datetime.today()
        processDateDict = dict()
        processDateDict["Year"] = nowDateTime.year
        processDateDict["Month"] = nowDateTime.month
        processDateDict["Day"] = nowDateTime.day
        productsDict["ProcessDate"] = processDateDict
        processTimeDict = dict()
        processTimeDict["Hour"] = nowDateTime.hour
        processTimeDict["Minute"] = nowDateTime.minute
        processTimeDict["Second"] = nowDateTime.second
        productsDict["ProcessTime"] = processTimeDict
        for key in calcdValuesDict:
            productsDict[key] = calcdValuesDict[key]

        filesDict = dict()
        filesDict["FileBaseName"] = self.generateOutputBaseName()
        filesDict["ProviderMetadata"] = self.headerFileName

        for key in outFilesDict:
            filesDict[key] = os.path.basename(outFilesDict[key])

        sensorDict = dict()
        sensorDict["ARCSISensorName"] = self.sensor

        acqDict = dict()
        acqDateDict = dict()
        acqDateDict["Year"] = self.acquisitionTime.year
        acqDateDict["Month"] = self.acquisitionTime.month
        acqDateDict["Day"] = self.acquisitionTime.day
        acqDict["Date"] = acqDateDict
        acqTimeDict = dict()
        acqTimeDict["Hour"] = self.acquisitionTime.hour
        acqTimeDict["Minute"] = self.acquisitionTime.minute
        acqTimeDict["Second"] = self.acquisitionTime.second
        acqDict["Time"] = acqTimeDict
        acqDict["SolarZenith"] = self.solarZenith
        acqDict["SolarAzimuth"] = self.solarAzimuth
        acqDict["sensorZenith"] = self.sensorZenith
        acqDict["sensorAzimuth"] = self.sensorAzimuth

        locDict = dict()
        locGeogDict = dict()
        locGeogDict["CentreLon"] = self.lonCentre
        locGeogDict["CentreLat"] = self.latCentre
        locGeogBBOXDict = dict()
        locGeogBBOXDict["TLLat"] = self.latTL
        locGeogBBOXDict["TLLon"] = self.lonTL
        locGeogBBOXDict["TRLat"] = self.latTR
        locGeogBBOXDict["TRLon"] = self.lonTR
        locGeogBBOXDict["BLLat"] = self.latBL
        locGeogBBOXDict["BLLon"] = self.lonBL
        locGeogBBOXDict["BRLat"] = self.latBR
        locGeogBBOXDict["BRLon"] = self.lonBR
        locGeogDict["BBOX"] = locGeogBBOXDict
        locDict["Geographical"] = locGeogDict

        locProjDict = dict()
        if validMaskImage != "":
            if not footprintCalc:
                rsgislib.rastergis.clumps_spatial_extent(
                    clumps_img=validMaskImage,
                    min_xx="MinXX",
                    min_xy="MinXY",
                    max_xx="MaxXX",
                    max_xy="MaxXY",
                    min_yx="MinYX",
                    min_yy="MinYY",
                    max_yx="MaxYX",
                    max_yy="MaxYY",
                    rat_band=1,
                )

            ratDataset = gdal.Open(validMaskImage)
            minXXVals = rat.readColumn(ratDataset, "MinXX")
            minXYVals = rat.readColumn(ratDataset, "MinXY")
            maxXXVals = rat.readColumn(ratDataset, "MaxXX")
            maxXYVals = rat.readColumn(ratDataset, "MaxXY")
            minYXVals = rat.readColumn(ratDataset, "MinYX")
            minYYVals = rat.readColumn(ratDataset, "MinYY")
            maxYXVals = rat.readColumn(ratDataset, "MaxYX")
            maxYYVals = rat.readColumn(ratDataset, "MaxYY")
            ratDataset = None

            ## Remove First Row which is no data...
            dataMask = numpy.ones_like(minXXVals, dtype=numpy.int16)
            dataMask[0] = 0
            minXXVals = minXXVals[dataMask == 1]
            minXYVals = minXYVals[dataMask == 1]
            maxXXVals = maxXXVals[dataMask == 1]
            maxXYVals = maxXYVals[dataMask == 1]
            minYXVals = minYXVals[dataMask == 1]
            minYYVals = minYYVals[dataMask == 1]
            maxYXVals = maxYXVals[dataMask == 1]
            maxYYVals = maxYYVals[dataMask == 1]

            ## Remove any features which are all zero (i.e., polygon not present...)
            minXXValsSub = minXXVals[
                numpy.logical_not(
                    (minXXVals == 0)
                    & (minXYVals == 0)
                    & (maxXXVals == 0)
                    & (maxXYVals == 0)
                    & (minYXVals == 0)
                    & (minYYVals == 0)
                    & (maxYXVals == 0)
                    & (maxYYVals == 0)
                )
            ]
            minXYValsSub = minXYVals[
                numpy.logical_not(
                    (minXXVals == 0)
                    & (minXYVals == 0)
                    & (maxXXVals == 0)
                    & (maxXYVals == 0)
                    & (minYXVals == 0)
                    & (minYYVals == 0)
                    & (maxYXVals == 0)
                    & (maxYYVals == 0)
                )
            ]
            maxXXValsSub = maxXXVals[
                numpy.logical_not(
                    (minXXVals == 0)
                    & (minXYVals == 0)
                    & (maxXXVals == 0)
                    & (maxXYVals == 0)
                    & (minYXVals == 0)
                    & (minYYVals == 0)
                    & (maxYXVals == 0)
                    & (maxYYVals == 0)
                )
            ]
            maxXYValsSub = maxXYVals[
                numpy.logical_not(
                    (minXXVals == 0)
                    & (minXYVals == 0)
                    & (maxXXVals == 0)
                    & (maxXYVals == 0)
                    & (minYXVals == 0)
                    & (minYYVals == 0)
                    & (maxYXVals == 0)
                    & (maxYYVals == 0)
                )
            ]
            minYXValsSub = minYXVals[
                numpy.logical_not(
                    (minXXVals == 0)
                    & (minXYVals == 0)
                    & (maxXXVals == 0)
                    & (maxXYVals == 0)
                    & (minYXVals == 0)
                    & (minYYVals == 0)
                    & (maxYXVals == 0)
                    & (maxYYVals == 0)
                )
            ]
            minYYValsSub = minYYVals[
                numpy.logical_not(
                    (minXXVals == 0)
                    & (minXYVals == 0)
                    & (maxXXVals == 0)
                    & (maxXYVals == 0)
                    & (minYXVals == 0)
                    & (minYYVals == 0)
                    & (maxYXVals == 0)
                    & (maxYYVals == 0)
                )
            ]
            maxYXValsSub = maxYXVals[
                numpy.logical_not(
                    (minXXVals == 0)
                    & (minXYVals == 0)
                    & (maxXXVals == 0)
                    & (maxXYVals == 0)
                    & (minYXVals == 0)
                    & (minYYVals == 0)
                    & (maxYXVals == 0)
                    & (maxYYVals == 0)
                )
            ]
            maxYYValsSub = maxYYVals[
                numpy.logical_not(
                    (minXXVals == 0)
                    & (minXYVals == 0)
                    & (maxXXVals == 0)
                    & (maxXYVals == 0)
                    & (minYXVals == 0)
                    & (minYYVals == 0)
                    & (maxYXVals == 0)
                    & (maxYYVals == 0)
                )
            ]

            numFeats = minXXValsSub.shape[0]

            if numFeats == 1:
                locProjBBOXDict = dict()
                locProjBBOXDict["TLX"] = minXXValsSub[0]
                locProjBBOXDict["TLY"] = maxYYValsSub[0]
                locProjBBOXDict["TRX"] = maxXXValsSub[0]
                locProjBBOXDict["TRY"] = maxYYValsSub[0]
                locProjBBOXDict["BLX"] = minXXValsSub[0]
                locProjBBOXDict["BLY"] = minYYValsSub[0]
                locProjBBOXDict["BRX"] = maxXXValsSub[0]
                locProjBBOXDict["BRY"] = minYYValsSub[0]
                locProjDict["BBOX"] = locProjBBOXDict
                locProjDict["CentreX"] = minXXValsSub[0] + (
                    (maxXXValsSub[0] - minXXValsSub[0]) / 2
                )
                locProjDict["CentreY"] = minYYValsSub[0] + (
                    (maxYYValsSub[0] - minYYValsSub[0]) / 2
                )

                locProjVPolyDict = dict()
                locProjVPolyDict["MinXX"] = minXXValsSub[0]
                locProjVPolyDict["MinXY"] = minXYValsSub[0]
                locProjVPolyDict["MaxYX"] = maxYXValsSub[0]
                locProjVPolyDict["MaxYY"] = maxYYValsSub[0]
                locProjVPolyDict["MaxXX"] = maxXXValsSub[0]
                locProjVPolyDict["MaxXY"] = maxXYValsSub[0]
                locProjVPolyDict["MinYX"] = minYXValsSub[0]
                locProjVPolyDict["MinYY"] = minYYValsSub[0]
                locProjDict["VPOLY"] = locProjVPolyDict
        else:
            locProjDict["CentreX"] = self.xCentre
            locProjDict["CentreY"] = self.yCentre
            locProjBBOXDict = dict()
            locProjBBOXDict["TLX"] = self.xTL
            locProjBBOXDict["TLY"] = self.yTL
            locProjBBOXDict["TRX"] = self.xTR
            locProjBBOXDict["TRY"] = self.yTR
            locProjBBOXDict["BLX"] = self.xBL
            locProjBBOXDict["BLY"] = self.yBL
            locProjBBOXDict["BRX"] = self.xBR
            locProjBBOXDict["BRY"] = self.yBR
            locProjDict["BBOX"] = locProjBBOXDict
        locDict["Projected"] = locProjDict

        jsonBlock = dict()
        jsonBlock["FileInfo"] = filesDict
        jsonBlock["SensorInfo"] = sensorDict
        jsonBlock["AcquasitionInfo"] = acqDict
        jsonBlock["LocationInfo"] = locDict
        jsonBlock["ProductsInfo"] = productsDict
        jsonBlock["SoftwareInfo"] = softwareDict

        return jsonBlock

    def generateMetaDataFile(
        self,
        outputPath,
        outputFileName,
        productsStr,
        validMaskImage="",
        footprintCalc=False,
        calcdValuesDict=None,
        outFilesDict=None,
    ):
        """
        Provides a default implementation for generating file metadata.
        """
        if outFilesDict is None:
            outFilesDict = dict()
        if calcdValuesDict is None:
            calcdValuesDict = dict()
        outJSONFilePath = os.path.join(outputPath, outputFileName)
        jsonData = self.getJSONDictDefaultMetaData(
            productsStr, validMaskImage, footprintCalc, calcdValuesDict, outFilesDict
        )
        with open(outJSONFilePath, "w") as outfile:
            json.dump(
                jsonData,
                outfile,
                sort_keys=True,
                indent=4,
                separators=(",", ": "),
                ensure_ascii=False,
            )

    @abstractmethod
    def expectedImageDataPresent(self):
        pass

    def maskInputImages(self):
        return False

    def imgNeedMosaicking(self):
        return False

    def inImgsDiffRes(self):
        return False

    @abstractmethod
    def mosaicImageTiles(self, outputPath):
        pass

    @abstractmethod
    def resampleImgRes(
        self, outputPath, resampleToLowResImg, resampleMethod="cubic", multicore=False
    ):
        pass

    @abstractmethod
    def sharpenLowResRadImgBands(self, inputImg, outputImage, outFormat):
        pass

    def hasThermal(self):
        return False

    def getBBOXLatLon(self):
        """
        :return: (MinLon, MaxLon, MinLat, MaxLat)
        """
        minLon = self.lonTL
        if self.lonBL < minLon:
            minLon = self.lonBL

        maxLon = self.lonBR
        if self.lonTR > maxLon:
            maxLon = self.lonTR

        minLat = self.latBL
        if self.latBR < minLat:
            minLat = self.latBR

        maxLat = self.latTL
        if self.latTR > maxLat:
            maxLat = self.latTR

        return (minLon, maxLon, minLat, maxLat)

    def getBBOX(self):
        """
        :return: (MinX, MaxX, MinY, MaxY)
        """
        minX = self.xTL
        if self.xBL < minX:
            minX = self.xBL

        maxX = self.xBR
        if self.xTR > maxX:
            maxX = self.xTR

        minY = self.yBR
        if self.yBL < minY:
            minY = self.yBL

        maxY = self.yTL
        if self.yTR > maxY:
            maxY = self.yTR

        return (minX, maxX, minY, maxY)

    def getReProjBBOX(
        self, wktFile, proj4File, useWKT2Reproject, xPxlRes, yPxlRes, snap2Grid
    ):
        projImgBBOX = dict()
        projImgBBOX["MinX"] = 0.0
        projImgBBOX["MaxX"] = 0.0
        projImgBBOX["MinY"] = 0.0
        projImgBBOX["MaxY"] = 0.0

        srcProj = osr.SpatialReference()
        srcProj.ImportFromWkt(self.inWKT)

        tarProj = osr.SpatialReference()
        if useWKT2Reproject:
            wktStr = rsgislib.tools.utils.read_text_file_no_new_lines(wktFile)
            tarProj.ImportFromWkt(wktStr)
        else:
            proj4Str = rsgislib.tools.utils.read_text_file_no_new_lines(proj4File)
            tarProj.ImportFromProj4(proj4Str)

        in_bbox = self.getBBOX()
        out_bbox = rsgislib.tools.geometrytools.reproj_bbox(in_bbox, srcProj, tarProj)

        yPxlResAbs = math.fabs(yPxlRes)

        minXPt = math.floor(out_bbox[0] / xPxlRes) * xPxlRes
        maxYPt = math.ceil(out_bbox[3] / yPxlResAbs) * yPxlResAbs

        projImgBBOX["MinX"] = minXPt
        projImgBBOX["MaxY"] = maxYPt

        maxXPt = math.ceil(out_bbox[1] / xPxlRes) * xPxlRes
        minYPt = math.floor(out_bbox[2] / yPxlResAbs) * yPxlResAbs

        projImgBBOX["MaxX"] = maxXPt
        projImgBBOX["MinY"] = minYPt

        return projImgBBOX

    @abstractmethod
    def applyImageDataMask(
        self,
        inputHeader,
        inputImage,
        outputPath,
        outputMaskName,
        outputImgName,
        outFormat,
        outWKTFile,
    ):
        pass

    @abstractmethod
    def convertImageToRadiance(
        self, outputPath, outputReflName, outputThermalName, outFormat
    ):
        pass

    @abstractmethod
    def generateImageSaturationMask(self, outputPath, outputName, outFormat):
        pass

    @abstractmethod
    def generateValidImageDataMask(
        self, outputPath, outputMaskName, viewAngleImg, outFormat
    ):
        pass

    def generateImageFootprint(self, validMaskImage, outputPath, outputName):
        print("Creating Vector Footprint...")
        rsgislib.rastergis.clumps_spatial_extent(
            clumps_img=validMaskImage,
            min_xx="MinXX",
            min_xy="MinXY",
            max_xx="MaxXX",
            max_xy="MaxXY",
            min_yx="MinYX",
            min_yy="MinYY",
            max_yx="MaxYX",
            max_yy="MaxYY",
            rat_band=1,
        )

        ratDataset = gdal.Open(validMaskImage)

        minXXVals = rat.readColumn(ratDataset, "MinXX")
        minXYVals = rat.readColumn(ratDataset, "MinXY")
        maxXXVals = rat.readColumn(ratDataset, "MaxXX")
        maxXYVals = rat.readColumn(ratDataset, "MaxXY")
        minYXVals = rat.readColumn(ratDataset, "MinYX")
        minYYVals = rat.readColumn(ratDataset, "MinYY")
        maxYXVals = rat.readColumn(ratDataset, "MaxYX")
        maxYYVals = rat.readColumn(ratDataset, "MaxYY")

        ## Remove First Row which is no data...
        dataMask = numpy.ones_like(minXXVals, dtype=numpy.int16)
        dataMask[0] = 0
        minXXVals = minXXVals[dataMask == 1]
        minXYVals = minXYVals[dataMask == 1]
        maxXXVals = maxXXVals[dataMask == 1]
        maxXYVals = maxXYVals[dataMask == 1]
        minYXVals = minYXVals[dataMask == 1]
        minYYVals = minYYVals[dataMask == 1]
        maxYXVals = maxYXVals[dataMask == 1]
        maxYYVals = maxYYVals[dataMask == 1]

        ## Remove any features which are all zero (i.e., polygon not present...)
        minXXValsSub = minXXVals[
            numpy.logical_not(
                (minXXVals == 0)
                & (minXYVals == 0)
                & (maxXXVals == 0)
                & (maxXYVals == 0)
                & (minYXVals == 0)
                & (minYYVals == 0)
                & (maxYXVals == 0)
                & (maxYYVals == 0)
            )
        ]
        minXYValsSub = minXYVals[
            numpy.logical_not(
                (minXXVals == 0)
                & (minXYVals == 0)
                & (maxXXVals == 0)
                & (maxXYVals == 0)
                & (minYXVals == 0)
                & (minYYVals == 0)
                & (maxYXVals == 0)
                & (maxYYVals == 0)
            )
        ]
        maxXXValsSub = maxXXVals[
            numpy.logical_not(
                (minXXVals == 0)
                & (minXYVals == 0)
                & (maxXXVals == 0)
                & (maxXYVals == 0)
                & (minYXVals == 0)
                & (minYYVals == 0)
                & (maxYXVals == 0)
                & (maxYYVals == 0)
            )
        ]
        maxXYValsSub = maxXYVals[
            numpy.logical_not(
                (minXXVals == 0)
                & (minXYVals == 0)
                & (maxXXVals == 0)
                & (maxXYVals == 0)
                & (minYXVals == 0)
                & (minYYVals == 0)
                & (maxYXVals == 0)
                & (maxYYVals == 0)
            )
        ]
        minYXValsSub = minYXVals[
            numpy.logical_not(
                (minXXVals == 0)
                & (minXYVals == 0)
                & (maxXXVals == 0)
                & (maxXYVals == 0)
                & (minYXVals == 0)
                & (minYYVals == 0)
                & (maxYXVals == 0)
                & (maxYYVals == 0)
            )
        ]
        minYYValsSub = minYYVals[
            numpy.logical_not(
                (minXXVals == 0)
                & (minXYVals == 0)
                & (maxXXVals == 0)
                & (maxXYVals == 0)
                & (minYXVals == 0)
                & (minYYVals == 0)
                & (maxYXVals == 0)
                & (maxYYVals == 0)
            )
        ]
        maxYXValsSub = maxYXVals[
            numpy.logical_not(
                (minXXVals == 0)
                & (minXYVals == 0)
                & (maxXXVals == 0)
                & (maxXYVals == 0)
                & (minYXVals == 0)
                & (minYYVals == 0)
                & (maxYXVals == 0)
                & (maxYYVals == 0)
            )
        ]
        maxYYValsSub = maxYYVals[
            numpy.logical_not(
                (minXXVals == 0)
                & (minXYVals == 0)
                & (maxXXVals == 0)
                & (maxXYVals == 0)
                & (minYXVals == 0)
                & (minYYVals == 0)
                & (maxYXVals == 0)
                & (maxYYVals == 0)
            )
        ]

        numFeats = minXXValsSub.shape[0]

        outVecLayerNamePath = os.path.join(outputPath, outputName)
        driver = ogr.GetDriverByName("GeoJSON")
        if os.path.exists("{}.geojson".format(outVecLayerNamePath)):
            driver.DeleteDataSource("{}.geojson".format(outVecLayerNamePath))
        outDatasource = driver.CreateDataSource(
            "{}.geojson".format(outVecLayerNamePath)
        )
        raster_srs = osr.SpatialReference()
        raster_srs.ImportFromWkt(ratDataset.GetProjectionRef())
        outLayer = outDatasource.CreateLayer(outputName, srs=raster_srs)

        fieldYearDefn = ogr.FieldDefn("Year", ogr.OFTInteger)
        fieldYearDefn.SetWidth(6)
        outLayer.CreateField(fieldYearDefn)

        fieldMonthDefn = ogr.FieldDefn("Month", ogr.OFTInteger)
        fieldMonthDefn.SetWidth(6)
        outLayer.CreateField(fieldMonthDefn)

        fieldDayDefn = ogr.FieldDefn("Day", ogr.OFTInteger)
        fieldDayDefn.SetWidth(6)
        outLayer.CreateField(fieldDayDefn)

        fieldBaseNameDefn = ogr.FieldDefn("BaseName", ogr.OFTString)
        fieldBaseNameDefn.SetWidth(254)
        outLayer.CreateField(fieldBaseNameDefn)

        fieldSolZenDefn = ogr.FieldDefn("SolZen", ogr.OFTReal)
        fieldSolZenDefn.SetWidth(10)
        fieldSolZenDefn.SetPrecision(6)
        outLayer.CreateField(fieldSolZenDefn)

        fieldSolAziDefn = ogr.FieldDefn("SolAzi", ogr.OFTReal)
        fieldSolAziDefn.SetWidth(10)
        fieldSolAziDefn.SetPrecision(6)
        outLayer.CreateField(fieldSolAziDefn)

        fieldSenZenDefn = ogr.FieldDefn("SenZen", ogr.OFTReal)
        fieldSenZenDefn.SetWidth(10)
        fieldSenZenDefn.SetPrecision(6)
        outLayer.CreateField(fieldSenZenDefn)

        fieldSenAziDefn = ogr.FieldDefn("SenAzi", ogr.OFTReal)
        fieldSenAziDefn.SetWidth(10)
        fieldSenAziDefn.SetPrecision(6)
        outLayer.CreateField(fieldSenAziDefn)

        fieldCenLatDefn = ogr.FieldDefn("CenLat", ogr.OFTReal)
        fieldCenLatDefn.SetWidth(10)
        fieldCenLatDefn.SetPrecision(6)
        outLayer.CreateField(fieldCenLatDefn)

        fieldCenLonDefn = ogr.FieldDefn("CenLon", ogr.OFTReal)
        fieldCenLonDefn.SetWidth(10)
        fieldCenLonDefn.SetPrecision(6)
        outLayer.CreateField(fieldCenLonDefn)

        for i in range(numFeats):
            wktStr = (
                "POLYGON(("
                + str(minXXValsSub[i])
                + " "
                + str(minXYValsSub[i])
                + ", "
                + str(maxYXValsSub[i])
                + " "
                + str(maxYYValsSub[i])
                + ", "
                + str(maxXXValsSub[i])
                + " "
                + str(maxXYValsSub[i])
                + ", "
                + str(minYXValsSub[i])
                + " "
                + str(minYYValsSub[i])
                + ", "
                + str(minXXValsSub[i])
                + " "
                + str(minXYValsSub[i])
                + "))"
            )
            poly = ogr.CreateGeometryFromWkt(wktStr)
            feat = ogr.Feature(outLayer.GetLayerDefn())
            feat.SetGeometry(poly)
            feat.SetField("Year", self.acquisitionTime.year)
            feat.SetField("Month", self.acquisitionTime.month)
            feat.SetField("Day", self.acquisitionTime.day)
            feat.SetField("BaseName", self.generateOutputBaseName())
            feat.SetField("SolZen", self.solarZenith)
            feat.SetField("SolAzi", self.solarAzimuth)
            feat.SetField("SenZen", self.sensorZenith)
            feat.SetField("SenAzi", self.sensorAzimuth)
            feat.SetField("CenLat", self.latCentre)
            feat.SetField("CenLon", self.lonCentre)

            if outLayer.CreateFeature(feat) != 0:
                print(str(i) + ": " + wktStr)
                print("Failed to create feature in shapefile.\n")
                sys.exit(1)
            feat.Destroy()

        outDatasource.Destroy()
        ratDataset = None
        return "{}.geojson".format(outVecLayerNamePath)

    def generateTopoDirectShadowMask(
        self, inputDEMImage, outputPath, outputName, outFormat, tmpPath
    ):
        try:
            print("Calculating a direct topographic shadow mask.")
            solarAz, solarZen = self.getSolarIrrStdSolarGeom()
            print(
                "Solar Zenith = " + str(solarZen) + " Solar Azimuth = " + str(solarAz)
            )
            outputImage = os.path.join(outputPath, outputName)
            tmpBaseName = os.path.splitext(outputName)[0]
            imgExtension = rsgislib.imageutils.get_file_img_extension(outFormat)
            outputTmpFile = os.path.join(
                tmpPath, tmpBaseName + "_toposhadow_tmp." + imgExtension
            )

            demDataset = gdal.Open(inputDEMImage, gdal.GA_ReadOnly)
            if demDataset is None:
                raise ARCSIException("Could not open DEM dataset.")
            demBand = demDataset.GetRasterBand(1)
            demMax = demBand.GetMaximum()
            if demMax is None:
                (demMin, demMax) = demBand.ComputeRasterMinMax(0)
            demMax10p = demMax * 1.1  # Go 10% higher than max in DEM.
            rsgislib.elevation.shadow_mask(
                inputDEMImage, outputTmpFile, solarAz, solarZen, demMax10p, outFormat
            )
            rsgislib.imagefilter.apply_median_filter(
                outputTmpFile, outputImage, 3, outFormat, rsgislib.TYPE_8UINT
            )
            if not self.debugMode:
                gdalDriver = gdal.GetDriverByName(outFormat)
                gdalDriver.Delete(outputTmpFile)

            return outputImage
        except Exception as e:
            raise e

    @abstractmethod
    def convertThermalToBrightness(
        self, inputRadImage, outputPath, outputName, outFormat, scaleFactor
    ):
        pass

    @abstractmethod
    def convertImageToTOARefl(
        self, inputRadImage, outputPath, outputName, outFormat, scaleFactor
    ):
        pass

    @abstractmethod
    def generateCloudMask(
        self,
        inputReflImage,
        inputSatImage,
        inputThermalImage,
        inputViewAngleImg,
        inputValidImg,
        outputPath,
        outputName,
        outFormat,
        tmpPath,
        scaleFactor,
        cloud_msk_methods=None,
    ):
        pass

    @abstractmethod
    def createCloudMaskDataArray(self, inImgDataArr):
        pass

    @abstractmethod
    def defineDarkShadowImageBand(self):
        pass

    def generateClearSkyMask(
        self,
        cloudsImg,
        inputValidImg,
        outputPath,
        outputName,
        outFormat,
        tmpPath,
        initClearSkyRegionDist=3000,
        initClearSkyRegionMinSize=3000,
        finalClearSkyRegionDist=1000,
        morphSize=21,
    ):
        try:
            outputImage = os.path.join(outputPath, outputName)
            tmpBaseName = os.path.splitext(outputName)[0]
            tmpBaseDIR = os.path.join(tmpPath, tmpBaseName)

            tmpDIRExisted = True
            if not os.path.exists(tmpBaseDIR):
                os.makedirs(tmpBaseDIR)
                tmpDIRExisted = False

            deleteTmpFiles = not self.debugMode

            rsgislib.imagecalibration.calc_clear_sky_regions(
                cloudsImg,
                inputValidImg,
                outputImage,
                outFormat,
                tmpBaseDIR,
                deleteTmpFiles,
                initClearSkyRegionDist,
                initClearSkyRegionMinSize,
                finalClearSkyRegionDist,
                morphSize,
            )

            if not self.debugMode:
                if not tmpDIRExisted:
                    shutil.rmtree(tmpBaseDIR, ignore_errors=True)

            return outputImage
        except Exception as e:
            raise e

    @abstractmethod
    def convertImageToSurfaceReflSglParam(
        self,
        inputRadImage,
        outputPath,
        outputName,
        outFormat,
        aeroProfile,
        atmosProfile,
        grdRefl,
        surfaceAltitude,
        aotVal,
        useBRDF,
        scaleFactor,
    ):
        pass

    @abstractmethod
    def calc6SCoefficients(
        self, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF
    ):
        pass

    def buildElevation6SCoeffLUT(
        self,
        aeroProfile,
        atmosProfile,
        grdRefl,
        aotVal,
        useBRDF,
        surfaceAltitudeMin,
        surfaceAltitudeMax,
    ):
        lut = list()
        elevRange = (surfaceAltitudeMax - surfaceAltitudeMin) / 100
        numElevSteps = int(math.ceil(elevRange) + 1)
        elevVal = surfaceAltitudeMin
        for i in range(numElevSteps):
            print("Building LUT Elevation ", elevVal)
            lut.append(
                rsgislib.imagecalibration.ElevLUTFeat(
                    Elev=elevVal,
                    Coeffs=self.calc6SCoefficients(
                        aeroProfile,
                        atmosProfile,
                        grdRefl,
                        (float(elevVal) / 1000),
                        aotVal,
                        useBRDF,
                    ),
                )
            )
            elevVal = elevVal + 100
        return lut

    @abstractmethod
    def convertImageToSurfaceReflDEMElevLUT(
        self,
        inputRadImage,
        inputDEMFile,
        outputPath,
        outputName,
        outFormat,
        aeroProfile,
        atmosProfile,
        grdRefl,
        aotVal,
        useBRDF,
        surfaceAltitudeMin,
        surfaceAltitudeMax,
        scaleFactor,
        elevCoeffs=None,
    ):
        pass

    def buildElevationAOT6SCoeffLUT(
        self,
        aeroProfile,
        atmosProfile,
        grdRefl,
        useBRDF,
        surfaceAltitudeMin,
        surfaceAltitudeMax,
        aotMin,
        aotMax,
    ):
        lut = list()
        elevRange = (surfaceAltitudeMax - surfaceAltitudeMin) / 100
        numElevSteps = int(math.ceil(elevRange) + 1)
        elevVal = surfaceAltitudeMin

        aotRange = (aotMax - aotMin) / 0.05
        numAOTSteps = int(math.ceil(aotRange) + 1) + 1
        aotVal = aotMin

        for i in range(numElevSteps):
            print("Building LUT Elevation ", elevVal)
            aotVal = aotMin
            aotCoeffLUT = list()
            for j in range(numAOTSteps):
                aotCoeffLUT.append(
                    rsgislib.imagecalibration.AOTLUTFeat(
                        AOT=aotVal,
                        Coeffs=self.calc6SCoefficients(
                            aeroProfile,
                            atmosProfile,
                            grdRefl,
                            (float(elevVal) / 1000),
                            aotVal,
                            useBRDF,
                        ),
                    )
                )
                aotVal = aotVal + 0.05
            lut.append(
                rsgislib.imagecalibration.ElevLUTFeat(Elev=elevVal, Coeffs=aotCoeffLUT)
            )
            elevVal = elevVal + 100
        return lut

    @abstractmethod
    def convertImageToSurfaceReflAOTDEMElevLUT(
        self,
        inputRadImage,
        inputDEMFile,
        inputAOTImage,
        outputPath,
        outputName,
        outFormat,
        aeroProfile,
        atmosProfile,
        grdRefl,
        useBRDF,
        surfaceAltitudeMin,
        surfaceAltitudeMax,
        aotMin,
        aotMax,
        scaleFactor,
        elevAOTCoeffs=None,
    ):
        pass

    def convertSREF2StdisedSREF(
        self,
        inputSREFImage,
        inputSREFWholeImage,
        inputDEMFile,
        inputTopoShadowMask,
        outputPath,
        outputName,
        outputWholeName,
        outFormat,
        tmpPath,
        sixsLUTCoeffs,
        aotLUT,
        scaleFactor,
        brdfBeta=1,
        outIncidenceAngle=0,
        outExitanceAngle=0,
    ):
        print("Converting to Standardised Reflectance")
        try:
            outputStdSREFImage = os.path.join(outputPath, outputName)
            outputStdSREFWholeImage = os.path.join(outputPath, outputWholeName)
            tmpBaseName = os.path.splitext(outputName)[0]
            tmpWholeBaseName = os.path.splitext(outputWholeName)[0]
            imgExtension = rsgislib.imageutils.get_file_img_extension(outFormat)
            tmpBaseDIR = os.path.join(tmpPath, tmpBaseName)

            tmpDIRExisted = True
            if not os.path.exists(tmpBaseDIR):
                os.makedirs(tmpBaseDIR)
                tmpDIRExisted = False

            # Define Angles with the correct origin.
            solarAz, solarZen = self.getSolarIrrStdSolarGeom()
            viewAz = 0.0
            viewZen = 0.0

            # Derive layers from DEM
            incidAngleImg = os.path.join(
                tmpBaseDIR, tmpBaseName + "_incidangle." + imgExtension
            )
            rsgislib.elevation.local_incidence_angle(
                inputDEMFile, incidAngleImg, solarAz, solarZen, outFormat
            )

            existAngleImg = os.path.join(
                tmpBaseDIR, tmpBaseName + "_existangle." + imgExtension
            )
            rsgislib.elevation.local_existance_angle(
                inputDEMFile, existAngleImg, viewAz, viewZen, outFormat
            )

            slopeImg = os.path.join(tmpBaseDIR, tmpBaseName + "_slope." + imgExtension)
            rsgislib.elevation.slope(inputDEMFile, slopeImg, "degrees", outFormat)

            # Derive the valid area mask
            validMaskSREF = os.path.join(
                tmpBaseDIR, tmpBaseName + "_validMask." + imgExtension
            )
            validMaskSREFWhole = os.path.join(
                tmpBaseDIR, tmpWholeBaseName + "_validMask." + imgExtension
            )
            if inputSREFWholeImage is not None:
                rsgislib.imageutils.gen_valid_mask(
                    inputSREFWholeImage, validMaskSREFWhole, outFormat, 0.0
                )
            rsgislib.imageutils.gen_valid_mask(
                inputSREFImage, validMaskSREF, outFormat, 0.0
            )

            # Calculate the solar irradiance
            solarIrradianceImg = os.path.join(
                tmpBaseDIR, tmpBaseName + "_solarirr." + imgExtension
            )
            solarIrradianceWholeImg = os.path.join(
                tmpBaseDIR, tmpWholeBaseName + "_solarirr." + imgExtension
            )
            if aotLUT:
                raise ARCSIException(
                    "Doh! Currently don't have an implementation of rsgislib.imagecalibration.calc_irradiance_img_elev_lut for using an Elev and AOT LUT..."
                )
            else:
                if inputSREFWholeImage is not None:
                    rsgislib.imagecalibration.calc_irradiance_img_elev_lut(
                        validMaskSREFWhole,
                        inputDEMFile,
                        incidAngleImg,
                        slopeImg,
                        inputSREFWholeImage,
                        inputTopoShadowMask,
                        solarIrradianceWholeImg,
                        outFormat,
                        solarZen,
                        scaleFactor,
                        sixsLUTCoeffs,
                    )
                rsgislib.imagecalibration.calc_irradiance_img_elev_lut(
                    validMaskSREF,
                    inputDEMFile,
                    incidAngleImg,
                    slopeImg,
                    inputSREFImage,
                    inputTopoShadowMask,
                    solarIrradianceImg,
                    outFormat,
                    solarZen,
                    scaleFactor,
                    sixsLUTCoeffs,
                )

            rsgislib.imagecalibration.calc_standardised_reflectance_sd2010(
                validMaskSREF,
                inputSREFImage,
                solarIrradianceImg,
                incidAngleImg,
                existAngleImg,
                outputStdSREFImage,
                outFormat,
                scaleFactor,
                brdfBeta,
                outIncidenceAngle,
                outExitanceAngle,
            )
            if inputSREFWholeImage is not None:
                rsgislib.imagecalibration.calc_standardised_reflectance_sd2010(
                    validMaskSREFWhole,
                    inputSREFWholeImage,
                    solarIrradianceWholeImg,
                    incidAngleImg,
                    existAngleImg,
                    outputStdSREFWholeImage,
                    outFormat,
                    scaleFactor,
                    brdfBeta,
                    outIncidenceAngle,
                    outExitanceAngle,
                )
            else:
                outputStdSREFWholeImage = ""

            if not self.debugMode:
                gdalDriver = gdal.GetDriverByName(outFormat)
                gdalDriver.Delete(incidAngleImg)
                gdalDriver.Delete(existAngleImg)
                gdalDriver.Delete(slopeImg)
                gdalDriver.Delete(validMaskSREF)
                gdalDriver.Delete(solarIrradianceImg)
                if not tmpDIRExisted:
                    shutil.rmtree(tmpBaseDIR, ignore_errors=True)

            return (outputStdSREFImage, outputStdSREFWholeImage)
        except Exception as e:
            raise e

    def calcDarkTargetOffsetsForBand(
        self,
        inputTOAImage,
        offsetImage,
        band,
        outFormat,
        histBinWidth,
        minObjSize,
        darkPxlPercentile,
        tmpDarkPxlsImg,
        tmpDarkPxlsClumpsImg,
        tmpDarkPxlsClumpsRMSmallImg,
        tmpDarkObjsImg,
    ):
        print("Band: ", band)
        bandHist = rsgislib.imagecalc.get_histogram(
            inputTOAImage, band, histBinWidth, False, 1, 10000
        )
        sumPxls = numpy.sum(bandHist[0])
        findTargets = True

        while findTargets:
            findTargets = False
            numPxlThreshold = sumPxls * darkPxlPercentile  # 0.001 # 1th percentile

            pxlThreshold = 0
            pxlCount = 0
            for bin in bandHist[0]:
                pxlCount = pxlCount + bin
                if pxlCount < numPxlThreshold:
                    pxlThreshold = pxlThreshold + histBinWidth
                else:
                    break
            print("Image Band Threshold (For Dark Pixels) = ", pxlThreshold)

            dataType = rsgislib.TYPE_8UINT
            expression = str("(b1!=0)&&(b1<=") + str(pxlThreshold) + str(")?1:0")
            bandDefns = []
            bandDefns.append(rsgislib.imagecalc.BandDefn("b1", inputTOAImage, band))
            rsgislib.imagecalc.band_math(
                tmpDarkPxlsImg, expression, outFormat, dataType, bandDefns
            )
            rsgislib.segmentation.clump(
                tmpDarkPxlsImg, tmpDarkPxlsClumpsImg, outFormat, False, 0.0
            )
            rsgislib.rastergis.pop_rat_img_stats(tmpDarkPxlsClumpsImg, True, False)
            rsgislib.segmentation.rm_small_clumps(
                tmpDarkPxlsClumpsImg, tmpDarkPxlsClumpsRMSmallImg, minObjSize, outFormat
            )
            rsgislib.segmentation.relabel_clumps(
                tmpDarkPxlsClumpsRMSmallImg, tmpDarkObjsImg, outFormat, False
            )
            rsgislib.rastergis.pop_rat_img_stats(tmpDarkObjsImg, True, False)
            stats2CalcTOA = list()
            stats2CalcTOA.append(
                rsgislib.rastergis.BandAttStats(
                    band=(band), min_field="MinTOARefl", mean_field="MeanTOARefl"
                )
            )
            rsgislib.rastergis.populate_rat_with_stats(
                inputTOAImage, tmpDarkObjsImg, stats2CalcTOA
            )

            ratDS = gdal.Open(tmpDarkObjsImg, gdal.GA_Update)
            Histogram = rat.readColumn(ratDS, "Histogram")
            if Histogram.size > 10:
                findTargets = False
                ratDS = None
                break
            else:
                findTargets = True
                darkPxlPercentile = darkPxlPercentile * 2
                print(
                    "Trying Dark Pixel Percentile Threshold: " + str(darkPxlPercentile)
                )
                ratDS = None

        ratDS = gdal.Open(tmpDarkObjsImg, gdal.GA_Update)
        Histogram = rat.readColumn(ratDS, "Histogram")
        selected = Histogram * 2
        selected[...] = 1
        selected[0] = 0
        rat.writeColumn(ratDS, "Selected", selected)
        ratDS = None

        rsgislib.rastergis.clumps_spatial_location(
            tmpDarkObjsImg, "Eastings", "Northings"
        )
        rsgislib.rastergis.select_clumps_on_grid(
            tmpDarkObjsImg,
            "Selected",
            "SelectedGrid",
            "Eastings",
            "Northings",
            "MeanTOARefl",
            "min",
            20,
            20,
        )

        ratDS = gdal.Open(tmpDarkObjsImg, gdal.GA_Update)
        Eastings = rat.readColumn(ratDS, "Eastings")
        Northings = rat.readColumn(ratDS, "Northings")
        MinTOARefl = rat.readColumn(ratDS, "MinTOARefl")
        SelectedGrid = rat.readColumn(ratDS, "SelectedGrid")
        ratDS = None

        Eastings = Eastings[SelectedGrid != 0]
        Northings = Northings[SelectedGrid != 0]
        MinTOARefl = MinTOARefl[SelectedGrid != 0]

        interpSmoothing = 10.0
        self.interpolateImageFromPointData(
            inputTOAImage,
            Eastings,
            Northings,
            MinTOARefl,
            offsetImage,
            outFormat,
            interpSmoothing,
            True,
            0.0,
        )

    def findPerBandDarkTargetsOffsets(
        self,
        inputTOAImage,
        numBands,
        outputPath,
        outputName,
        outFormat,
        tmpPath,
        minObjSize,
        darkPxlPercentile,
    ):
        try:
            tmpBaseName = os.path.splitext(outputName)[0]
            binWidth = 1

            bandDarkTargetOffsetImages = list()
            imgExtension = rsgislib.imageutils.get_file_img_extension(outFormat)
            tmpDarkPxlsImg = os.path.join(
                tmpPath, tmpBaseName + "_darkpxls." + imgExtension
            )
            tmpDarkPxlsClumpsImg = os.path.join(
                tmpPath, tmpBaseName + "_darkclumps." + imgExtension
            )
            tmpDarkPxlsClumpsRMSmallImg = os.path.join(
                tmpPath, tmpBaseName + "_darkclumpsrmsmall." + imgExtension
            )
            tmpDarkObjsImg = os.path.join(
                tmpPath, tmpBaseName + "_darkobjs." + imgExtension
            )

            for band in range(numBands):
                offsetImage = os.path.join(
                    tmpPath, f"{tmpBaseName}_darktargetoffs_b{band + 1}.{imgExtension}",
                )
                self.calcDarkTargetOffsetsForBand(
                    inputTOAImage,
                    offsetImage,
                    band + 1,
                    outFormat,
                    binWidth,
                    minObjSize,
                    darkPxlPercentile,
                    tmpDarkPxlsImg,
                    tmpDarkPxlsClumpsImg,
                    tmpDarkPxlsClumpsRMSmallImg,
                    tmpDarkObjsImg,
                )
                bandDarkTargetOffsetImages.append(offsetImage)

            outputImage = os.path.join(
                outputPath, tmpBaseName + "_dosuboffs." + imgExtension
            )
            print(outputImage)
            rsgislib.imageutils.stack_img_bands(
                bandDarkTargetOffsetImages,
                None,
                outputImage,
                None,
                0,
                outFormat,
                rsgislib.TYPE_32FLOAT,
            )

            if not self.debugMode:
                gdalDriver = gdal.GetDriverByName(outFormat)
                gdalDriver.Delete(tmpDarkPxlsImg)
                gdalDriver.Delete(tmpDarkPxlsClumpsImg)
                gdalDriver.Delete(tmpDarkPxlsClumpsRMSmallImg)
                gdalDriver.Delete(tmpDarkObjsImg)
                for image in bandDarkTargetOffsetImages:
                    gdalDriver.Delete(image)

            return outputImage

        except Exception as e:
            raise e

    def performDOSOnSingleBand(
        self,
        inputTOAImage,
        band,
        outputPath,
        tmpBaseName,
        bandName,
        outFormat,
        tmpPath,
        minObjSize,
        darkPxlPercentile,
        dosOutRefl,
    ):
        try:
            binWidth = 1

            imgExtension = rsgislib.imageutils.get_file_img_extension(outFormat)
            tmpDarkPxlsImg = os.path.join(
                tmpPath, tmpBaseName + "_darkpxls." + imgExtension
            )
            tmpDarkPxlsClumpsImg = os.path.join(
                tmpPath, tmpBaseName + "_darkclumps." + imgExtension
            )
            tmpDarkPxlsClumpsRMSmallImg = os.path.join(
                tmpPath, tmpBaseName + "_darkclumpsrmsmall." + imgExtension
            )
            tmpDarkObjsImg = os.path.join(
                tmpPath, tmpBaseName + "_darkobjs." + imgExtension
            )

            offsetImage = os.path.join(
                outputPath, f"{tmpBaseName}_darktargetoffs_b{band}.{imgExtension}"
            )
            self.calcDarkTargetOffsetsForBand(
                inputTOAImage,
                offsetImage,
                band,
                outFormat,
                binWidth,
                minObjSize,
                darkPxlPercentile,
                tmpDarkPxlsImg,
                tmpDarkPxlsClumpsImg,
                tmpDarkPxlsClumpsRMSmallImg,
                tmpDarkObjsImg,
            )

            outputName = f"{tmpBaseName}DOS{bandName}.{imgExtension}"
            outputImage = os.path.join(outputPath, outputName)
            print(outputImage)

            bandDefns = []
            bandDefns.append(rsgislib.imagecalc.BandDefn("Off", offsetImage, 1))
            bandDefns.append(rsgislib.imagecalc.BandDefn("TOA", inputTOAImage, band))
            expression = (
                "(TOA==0)?0:((TOA-Off)+"
                + str(dosOutRefl)
                + ")<=0?1.0:(TOA-Off)+"
                + str(dosOutRefl)
            )
            rsgislib.imagecalc.band_math(
                outputImage, expression, outFormat, rsgislib.TYPE_16UINT, bandDefns
            )

            if not self.debugMode:
                gdalDriver = gdal.GetDriverByName(outFormat)
                gdalDriver.Delete(tmpDarkPxlsImg)
                gdalDriver.Delete(tmpDarkPxlsClumpsImg)
                gdalDriver.Delete(tmpDarkPxlsClumpsRMSmallImg)
                gdalDriver.Delete(tmpDarkObjsImg)

            return outputImage

        except Exception as e:
            raise e

    def performLocalDOSOnSingleBand(
        self,
        inputTOAImage,
        band,
        outputPath,
        tmpBaseName,
        bandName,
        outFormat,
        tmpPath,
        minObjSize,
        darkPxlPercentile,
        blockSize,
        dosOutRefl,
    ):
        try:
            bandDarkTargetOffsetImages = list()
            imgExtension = rsgislib.imageutils.get_file_img_extension(outFormat)
            tmpDarkTargetAllImage = os.path.join(
                tmpPath, tmpBaseName + "_darkpxls_allbands." + imgExtension
            )
            tmpDarkPxlsImg = os.path.join(
                tmpPath, tmpBaseName + "_darkpxls." + imgExtension
            )
            tmpDarkPxlsClumpsImg = os.path.join(
                tmpPath, tmpBaseName + "_darkclumps." + imgExtension
            )
            tmpDarkPxlsClumpsRMSmallImg = os.path.join(
                tmpPath, tmpBaseName + "_darkclumpsrmsmall." + imgExtension
            )
            tmpDarkObjsImg = os.path.join(
                tmpPath, tmpBaseName + "_darkobjs." + imgExtension
            )

            binWidth = 1
            self.findDOSLocalDarkTargets(
                inputTOAImage,
                tmpDarkTargetAllImage,
                blockSize,
                "KEA",
                binWidth,
                darkPxlPercentile,
            )

            print("Band ", band)
            rsgislib.imageutils.select_img_bands(
                tmpDarkTargetAllImage,
                tmpDarkPxlsImg,
                "KEA",
                rsgislib.TYPE_8UINT,
                [band],
            )
            rsgislib.segmentation.clump(
                tmpDarkPxlsImg, tmpDarkPxlsClumpsImg, "KEA", False, 0.0
            )
            rsgislib.rastergis.pop_rat_img_stats(tmpDarkPxlsClumpsImg, True, False)
            rsgislib.segmentation.rm_small_clumps(
                tmpDarkPxlsClumpsImg, tmpDarkPxlsClumpsRMSmallImg, minObjSize, "KEA"
            )
            rsgislib.segmentation.relabel_clumps(
                tmpDarkPxlsClumpsRMSmallImg, tmpDarkObjsImg, "KEA", False
            )
            rsgislib.rastergis.pop_rat_img_stats(tmpDarkObjsImg, True, False)
            stats2CalcTOA = list()
            stats2CalcTOA.append(
                rsgislib.rastergis.BandAttStats(
                    band=(band + 1), min_field="MinTOARefl", mean_field="MeanTOARefl"
                )
            )
            rsgislib.rastergis.populate_rat_with_stats(
                inputTOAImage, tmpDarkObjsImg, stats2CalcTOA
            )

            ratDS = gdal.Open(tmpDarkObjsImg, gdal.GA_Update)
            Histogram = rat.readColumn(ratDS, "Histogram")
            selected = Histogram * 2
            selected[...] = 1
            selected[0] = 0
            rat.writeColumn(ratDS, "Selected", selected)
            ratDS = None

            rsgislib.rastergis.clumps_spatial_location(
                tmpDarkObjsImg, "Eastings", "Northings"
            )
            rsgislib.rastergis.select_clumps_on_grid(
                tmpDarkObjsImg,
                "Selected",
                "SelectedGrid",
                "Eastings",
                "Northings",
                "MeanTOARefl",
                "min",
                20,
                20,
            )

            print("Interpolating the offset image...")

            offsetImage = os.path.join(
                tmpPath, f"{tmpBaseName}_darktargetoffs_b{band}.{imgExtension}"
            )

            ratDS = gdal.Open(tmpDarkObjsImg, gdal.GA_Update)
            Eastings = rat.readColumn(ratDS, "Eastings")
            Northings = rat.readColumn(ratDS, "Northings")
            MinTOARefl = rat.readColumn(ratDS, "MinTOARefl")
            SelectedGrid = rat.readColumn(ratDS, "SelectedGrid")
            ratDS = None

            Eastings = Eastings[SelectedGrid != 0]
            Northings = Northings[SelectedGrid != 0]
            MinTOARefl = MinTOARefl[SelectedGrid != 0]

            interpSmoothing = 10.0
            self.interpolateImageFromPointData(
                inputTOAImage,
                Eastings,
                Northings,
                MinTOARefl,
                offsetImage,
                "KEA",
                interpSmoothing,
                True,
                0.0,
            )

            outputName = f"{tmpBaseName}DOS{bandName}.{imgExtension}"
            outputImage = os.path.join(outputPath, outputName)
            print(outputImage)

            bandDefns = []
            bandDefns.append(rsgislib.imagecalc.BandDefn("Off", offsetImage, 1))
            bandDefns.append(rsgislib.imagecalc.BandDefn("TOA", inputTOAImage, band))
            expression = (
                "(TOA==0)?0:((TOA-Off)+"
                + str(dosOutRefl)
                + ")<=0?1.0:(TOA-Off)+"
                + str(dosOutRefl)
            )
            rsgislib.imagecalc.band_math(
                outputImage, expression, outFormat, rsgislib.TYPE_16UINT, bandDefns
            )

            if not self.debugMode:
                gdalDriver = gdal.GetDriverByName("KEA")
                gdalDriver.Delete(tmpDarkPxlsImg)
                gdalDriver.Delete(tmpDarkPxlsClumpsImg)
                gdalDriver.Delete(tmpDarkPxlsClumpsRMSmallImg)
                gdalDriver.Delete(tmpDarkObjsImg)
                gdalDriver.Delete(tmpDarkTargetAllImage)
                gdalDriver.Delete(offsetImage)

            return outputImage

        except Exception as e:
            raise e

    def findDOSLocalDarkTargets(
        self,
        inputTOAImage,
        darkTargetImage,
        blockSize,
        outFormat,
        histBinWidth,
        darkPxlPercentile,
    ):
        reader = ImageReader(
            inputTOAImage, windowxsize=blockSize, windowysize=blockSize
        )
        writer = None
        for (info, block) in reader:
            out = numpy.zeros_like(block)

            # Iterate through the image bands
            for i in range(len(out)):
                minVal = numpy.min(block[i])
                maxVal = numpy.max(block[i])
                if (maxVal - minVal) > 5:
                    data = block[i].flatten()
                    data = data[data != 0]

                    minVal = numpy.min(data)
                    maxVal = numpy.max(data)

                    numBins = (math.ceil(maxVal - minVal) / histBinWidth) + 1
                    if data.shape[0] > ((blockSize * blockSize) * 0.1):
                        # out[i,...] = 1
                        histo, histoBins = numpy.histogram(
                            data, bins=numBins, range=(float(minVal), float(maxVal))
                        )
                        numValues = numpy.sum(histo)
                        numValsPercentile = math.floor(numValues * darkPxlPercentile)
                        binValCount = 0
                        threshold = 0.0
                        for n in range(histo.shape[0]):
                            if (binValCount + histo[n]) > numValsPercentile:
                                break
                            else:
                                binValCount = binValCount + histo[n]
                                threshold = histoBins[n + 1]
                        out[i, ((block[i] <= threshold) & (block[i] > 0))] = 1
                    else:
                        out[i, ...] = 0

            if writer is None:
                writer = ImageWriter(
                    darkTargetImage,
                    info=info,
                    firstblock=out,
                    drivername=outFormat,
                    creationoptions=[],
                )
            else:
                writer.write(out)
        writer.close(calcStats=False)

    def findPerBandLocalDarkTargetsOffsets(
        self,
        inputTOAImage,
        numBands,
        outputPath,
        outputName,
        outFormat,
        tmpPath,
        blockSize,
        minObjSize,
        darkPxlPercentile,
    ):
        try:
            tmpBaseName = os.path.splitext(outputName)[0]
            binWidth = 1

            bandDarkTargetOffsetImages = list()
            imgExtension = rsgislib.imageutils.get_file_img_extension(outFormat)
            tmpDarkTargetAllImage = os.path.join(
                tmpPath, tmpBaseName + "_darkpxls_allbands." + imgExtension
            )
            tmpDarkPxlsImg = os.path.join(
                tmpPath, tmpBaseName + "_darkpxls." + imgExtension
            )
            tmpDarkPxlsClumpsImg = os.path.join(
                tmpPath, tmpBaseName + "_darkclumps." + imgExtension
            )
            tmpDarkPxlsClumpsRMSmallImg = os.path.join(
                tmpPath, tmpBaseName + "_darkclumpsrmsmall." + imgExtension
            )
            tmpDarkObjsImg = os.path.join(
                tmpPath, tmpBaseName + "_darkobjs." + imgExtension
            )

            self.findDOSLocalDarkTargets(
                inputTOAImage,
                tmpDarkTargetAllImage,
                blockSize,
                outFormat,
                binWidth,
                darkPxlPercentile,
            )

            for band in range(numBands):
                print("Band ", band + 1)
                rsgislib.imageutils.select_img_bands(
                    tmpDarkTargetAllImage,
                    tmpDarkPxlsImg,
                    outFormat,
                    rsgislib.TYPE_8UINT,
                    [band + 1],
                )
                rsgislib.segmentation.clump(
                    tmpDarkPxlsImg, tmpDarkPxlsClumpsImg, outFormat, False, 0.0
                )
                rsgislib.rastergis.pop_rat_img_stats(tmpDarkPxlsClumpsImg, True, False)
                rsgislib.segmentation.rm_small_clumps(
                    tmpDarkPxlsClumpsImg,
                    tmpDarkPxlsClumpsRMSmallImg,
                    minObjSize,
                    outFormat,
                )
                rsgislib.segmentation.relabel_clumps(
                    tmpDarkPxlsClumpsRMSmallImg, tmpDarkObjsImg, outFormat, False
                )
                rsgislib.rastergis.pop_rat_img_stats(tmpDarkObjsImg, True, False)
                stats2CalcTOA = list()
                stats2CalcTOA.append(
                    rsgislib.rastergis.BandAttStats(
                        band=(band + 1),
                        min_field="MinTOARefl",
                        mean_field="MeanTOARefl",
                    )
                )
                rsgislib.rastergis.populate_rat_with_stats(
                    inputTOAImage, tmpDarkObjsImg, stats2CalcTOA
                )

                ratDS = gdal.Open(tmpDarkObjsImg, gdal.GA_Update)
                Histogram = rat.readColumn(ratDS, "Histogram")
                selected = Histogram * 2
                selected[...] = 1
                selected[0] = 0
                rat.writeColumn(ratDS, "Selected", selected)
                ratDS = None

                rsgislib.rastergis.clumps_spatial_location(
                    tmpDarkObjsImg, "Eastings", "Northings"
                )
                rsgislib.rastergis.select_clumps_on_grid(
                    tmpDarkObjsImg,
                    "Selected",
                    "SelectedGrid",
                    "Eastings",
                    "Northings",
                    "MeanTOARefl",
                    "min",
                    20,
                    20,
                )

                print("Interpolating the offset image...")

                offsetImage = os.path.join(
                    tmpPath, f"{tmpBaseName}_darktargetoffs_b{band + 1}.{imgExtension}",
                )

                ratDS = gdal.Open(tmpDarkObjsImg, gdal.GA_Update)
                Eastings = rat.readColumn(ratDS, "Eastings")
                Northings = rat.readColumn(ratDS, "Northings")
                MinTOARefl = rat.readColumn(ratDS, "MinTOARefl")
                SelectedGrid = rat.readColumn(ratDS, "SelectedGrid")
                ratDS = None

                Eastings = Eastings[SelectedGrid != 0]
                Northings = Northings[SelectedGrid != 0]
                MinTOARefl = MinTOARefl[SelectedGrid != 0]

                interpSmoothing = 10.0
                self.interpolateImageFromPointData(
                    inputTOAImage,
                    Eastings,
                    Northings,
                    MinTOARefl,
                    offsetImage,
                    outFormat,
                    interpSmoothing,
                    True,
                    0.0,
                )

                bandDarkTargetOffsetImages.append(offsetImage)

            outputImage = os.path.join(
                outputPath, tmpBaseName + "_dosuboffs." + imgExtension
            )
            print(outputImage)
            rsgislib.imageutils.stack_img_bands(
                bandDarkTargetOffsetImages,
                None,
                outputImage,
                None,
                0,
                outFormat,
                rsgislib.TYPE_32FLOAT,
            )

            if not self.debugMode:
                gdalDriver = gdal.GetDriverByName(outFormat)
                gdalDriver.Delete(tmpDarkPxlsImg)
                gdalDriver.Delete(tmpDarkPxlsClumpsImg)
                gdalDriver.Delete(tmpDarkPxlsClumpsRMSmallImg)
                gdalDriver.Delete(tmpDarkObjsImg)
                gdalDriver.Delete(tmpDarkTargetAllImage)
                for image in bandDarkTargetOffsetImages:
                    gdalDriver.Delete(image)

            return outputImage

        except Exception as e:
            raise e

    def convertImageToReflectanceSimpleDarkSubtract(
        self,
        inputTOAImage,
        outputPath,
        outputName,
        outFormat,
        dosOutRefl,
        offsetsList=None,
    ):
        try:
            print("Perform Simple Dark Object Subtraction")
            outputImage = os.path.join(outputPath, outputName)

            if offsetsList is None:
                OffValDOS = collections.namedtuple("DOSOffset", ["offset"])
                percentiles = rsgislib.imagecalc.calc_band_percentile(
                    inputTOAImage, 0.01, 0
                )
                offsetsList = list()
                for val in percentiles:
                    offsetsList.append(OffValDOS(offset=val))

            rsgislib.imagecalibration.apply_subtract_single_offsets(
                inputTOAImage,
                outputImage,
                outFormat,
                rsgislib.TYPE_16UINT,
                True,
                True,
                0.0,
                dosOutRefl,
                offsetsList,
            )

            return outputImage, offsetsList
        except Exception as e:
            raise e

    def convertImageBandToReflectanceSimpleDarkSubtract(
        self, inputTOAImage, outputPath, outputName, outFormat, dosOutRefl, imgBand
    ):
        try:
            print("Perform Simple Dark Object Subtraction on image band")
            outputImage = os.path.join(outputPath, outputName)

            percentiles = rsgislib.imagecalc.calc_band_percentile(
                inputTOAImage, 0.01, 0
            )

            print("Band offset = " + str(percentiles[imgBand - 1]))
            expression = (
                "(b1 == 0.0)?0.0:((b1-"
                + str(percentiles[imgBand - 1])
                + ") + "
                + str(dosOutRefl)
                + ") < 0?1.0: (b1-"
                + str(percentiles[imgBand - 1])
                + ") + "
                + str(dosOutRefl)
            )
            bandDefns = []
            bandDefns.append(rsgislib.imagecalc.BandDefn("b1", inputTOAImage, imgBand))
            rsgislib.imagecalc.band_math(
                outputImage, expression, outFormat, rsgislib.TYPE_16UINT, bandDefns
            )

            return outputImage, percentiles[imgBand - 1]
        except Exception as e:
            raise e

    @abstractmethod
    def findDDVTargets(self, inputTOAImage, outputPath, outputName, outFormat, tmpPath):
        pass

    @abstractmethod
    def estimateImageToAODUsingDDV(
        self,
        inputRADImage,
        inputTOAImage,
        inputDEMFile,
        shadowMask,
        outputPath,
        outputName,
        outFormat,
        tmpPath,
        aeroProfile,
        atmosProfile,
        grdRefl,
        aotValMin,
        aotValMax,
    ):
        pass

    @abstractmethod
    def estimateImageToAODUsingDOS(
        self,
        inputRADImage,
        inputTOAImage,
        inputDEMFile,
        shadowMask,
        outputPath,
        outputName,
        outFormat,
        tmpPath,
        aeroProfile,
        atmosProfile,
        grdRefl,
        aotValMin,
        aotValMax,
        globalDOS,
        simpleDOS,
        dosOutRefl,
    ):
        pass

    @abstractmethod
    def estimateSingleAOTFromDOS(
        self,
        radianceImage,
        toaImage,
        inputDEMFile,
        tmpPath,
        outputName,
        outFormat,
        aeroProfile,
        atmosProfile,
        grdRefl,
        minAOT,
        maxAOT,
        dosOutRefl,
    ):
        pass

    def estimateSingleAOTFromDOSBandImpl(
        self,
        radianceImage,
        toaImage,
        inputDEMFile,
        tmpPath,
        outputName,
        outFormat,
        aeroProfile,
        atmosProfile,
        grdRefl,
        aotValMin,
        aotValMax,
        dosOutRefl,
        imgBand,
    ):
        try:
            print("Estimating a single AOD value Using DOS")
            # Using a simple DOS find RAD and SREF a value from within the image.
            radVal = 0.0
            srefVal = 0.0
            imgExtension = rsgislib.imageutils.get_file_img_extension(outFormat)

            dosBandFileName = outputName + "_simbanddos." + imgExtension
            dosBandFile, bandOff = self.convertImageBandToReflectanceSimpleDarkSubtract(
                toaImage, tmpPath, dosBandFileName, outFormat, dosOutRefl, imgBand
            )

            darkROIMask = os.path.join(
                tmpPath, outputName + "_darkROIMask." + imgExtension,
            )
            expression = "((b1 != 0) && (b1 < " + str(dosOutRefl + 5) + "))?1.0:0.0"
            bandDefns = []
            bandDefns.append(rsgislib.imagecalc.BandDefn("b1", dosBandFile, 1))
            rsgislib.imagecalc.band_math(
                darkROIMask, expression, outFormat, rsgislib.TYPE_8UINT, bandDefns
            )

            darkROIMaskClumps = os.path.join(
                tmpPath, outputName + "_darkROIMaskClumps." + imgExtension,
            )
            rsgislib.segmentation.clump(
                darkROIMask, darkROIMaskClumps, outFormat, False, 0.0
            )
            rsgislib.rastergis.pop_rat_img_stats(darkROIMaskClumps, True, False)

            darkROIMaskClumpsRMSmall = os.path.join(
                tmpPath, outputName + "_darkROIMaskClumpsRMSmall." + imgExtension,
            )
            rsgislib.segmentation.rm_small_clumps(
                darkROIMaskClumps, darkROIMaskClumpsRMSmall, 5, outFormat
            )
            darkROIMaskClumpsFinal = os.path.join(
                tmpPath, outputName + "_darkROIMaskClumpsFinal." + imgExtension,
            )
            rsgislib.segmentation.relabel_clumps(
                darkROIMaskClumpsRMSmall, darkROIMaskClumpsFinal, outFormat, False
            )
            rsgislib.rastergis.pop_rat_img_stats(darkROIMaskClumpsFinal, True, False)

            stats2CalcDOS = list()
            stats2CalcDOS.append(
                rsgislib.rastergis.BandAttStats(band=imgBand, mean_field="MeanTOARefl")
            )
            rsgislib.rastergis.populate_rat_with_stats(
                toaImage, darkROIMaskClumpsFinal, stats2CalcDOS
            )
            stats2CalcRAD = list()
            stats2CalcRAD.append(
                rsgislib.rastergis.BandAttStats(band=imgBand, mean_field="MeanRad")
            )
            rsgislib.rastergis.populate_rat_with_stats(
                radianceImage, darkROIMaskClumpsFinal, stats2CalcRAD
            )
            stats2CalcElev = list()
            stats2CalcElev.append(
                rsgislib.rastergis.BandAttStats(band=1, mean_field="MeanElev")
            )
            rsgislib.rastergis.populate_rat_with_stats(
                inputDEMFile, darkROIMaskClumpsFinal, stats2CalcElev
            )

            ratDS = gdal.Open(darkROIMaskClumpsFinal, gdal.GA_ReadOnly)
            Histogram = rat.readColumn(ratDS, "Histogram")
            MeanTOARefl = rat.readColumn(ratDS, "MeanTOARefl")
            MeanRad = rat.readColumn(ratDS, "MeanRad")
            MeanElev = rat.readColumn(ratDS, "MeanElev")
            ratDS = None

            maxObjSizeArrIdx = numpy.where(Histogram == numpy.max(Histogram))

            reflTOA = MeanTOARefl[maxObjSizeArrIdx][0]
            reflDOS = reflTOA - bandOff
            if reflDOS < dosOutRefl:
                reflDOS = dosOutRefl
            reflDOS = reflDOS / 1000
            radVal = MeanRad[maxObjSizeArrIdx][0]
            elevVal = MeanElev[maxObjSizeArrIdx][0]

            if not self.debugMode:
                gdalDriver = gdal.GetDriverByName(outFormat)
                gdalDriver.Delete(dosBandFile)
                gdalDriver.Delete(darkROIMask)
                gdalDriver.Delete(darkROIMaskClumps)
                gdalDriver.Delete(darkROIMaskClumpsRMSmall)
                gdalDriver.Delete(darkROIMaskClumpsFinal)

            # Second Step - estimate AOT to get RAD value to SREF.
            numAOTValTests = int(math.ceil((aotValMax - aotValMin) / 0.05)) + 1
            if not numAOTValTests >= 1:
                raise ARCSIException(
                    "min and max AOT range are too close together, they need to be at least 0.05 apart."
                )

            cAOT = aotValMin
            cDist = 0.0
            minAOT = 0.0
            minDist = 0.0

            for j in range(numAOTValTests):
                cAOT = aotValMin + (0.05 * j)
                cDist = self.run6SToOptimiseAODValue(
                    cAOT,
                    radVal,
                    reflDOS,
                    aeroProfile,
                    atmosProfile,
                    grdRefl,
                    elevVal / 1000,
                )
                if j == 0:
                    minAOT = cAOT
                    minDist = cDist
                elif cDist < minDist:
                    minAOT = cAOT
                    minDist = cDist
            aotVal = minAOT
            print("IDENTIFIED AOT: ", aotVal)

            return aotVal
        except Exception as e:
            raise

    @abstractmethod
    def setBandNames(self, imageFile):
        pass

    def interpolateImageFromPointData(
        self,
        templateInImage,
        xVals,
        yVals,
        zVals,
        outputImage,
        outFormat,
        smoothingParam,
        notNegOut,
        notNegMinVal,
    ):
        print("Interpolating Image: Number of Features = ", xVals.shape[0])

        reader = ImageReader(templateInImage, windowxsize=200, windowysize=200)
        writer = None
        for (info, block) in reader:
            pxlCoords = info.getBlockCoordArrays()
            interZnn = scipy.interpolate.griddata(
                (xVals, yVals),
                zVals,
                (pxlCoords[0].flatten(), pxlCoords[1].flatten()),
                method="nearest",
            )
            interZcub = scipy.interpolate.griddata(
                (xVals, yVals),
                zVals,
                (pxlCoords[0].flatten(), pxlCoords[1].flatten()),
                method="cubic",
            )
            interZ = numpy.where(numpy.isnan(interZcub), interZnn, interZcub)
            if notNegOut:
                interZ = numpy.where(interZ < 0, notNegMinVal, interZ)
            out = numpy.reshape(interZ, block[0].shape)
            out = numpy.expand_dims(out, axis=0)

            if writer is None:
                writer = ImageWriter(
                    outputImage,
                    info=info,
                    firstblock=out,
                    drivername=outFormat,
                    creationoptions=[],
                )
            else:
                writer.write(out)
        writer.close(calcStats=True)
        print("Interpolating Image - Complete")

    @abstractmethod
    def cleanLocalFollowProcessing(self):
        pass

    def cleanFollowProcessing(self, outputDIR=None, fileEnding2Keep=None):
        self.cleanLocalFollowProcessing()
        if (outputDIR is not None) and (fileEnding2Keep is not None):
            baseName = self.generateOutputBaseName()
            import glob

            allFiles = glob.glob(os.path.join(outputDIR, baseName + "*"))
            for file in allFiles:
                foundEnd = False
                for fileEnd in fileEnding2Keep:
                    if fileEnd in file:
                        foundEnd = True
                        break
                if not foundEnd:
                    print("REMOVE: " + file)
                    os.remove(file)
