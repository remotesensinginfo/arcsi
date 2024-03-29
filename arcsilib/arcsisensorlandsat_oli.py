"""
Module that contains the ARCSILandsatOLISensor class.
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

import collections
import datetime
import json
import math
import os
import shutil

import fmask.config
import fmask.fmask
import fmask.landsatangles
import numpy
import Py6S
import rios.fileinfo
import rsgislib
import rsgislib.imagecalc
import rsgislib.imagecalibration
import rsgislib.imagecalibration.solarangles
import rsgislib.imageutils
import rsgislib.rastergis
import rsgislib.segmentation
import rsgislib.segmentation.shepherdseg
import rsgislib.tools.geometrytools
import rsgislib.tools.utils
from osgeo import gdal, osr
from rios import rat

from .arcsiexception import ARCSIException
from .arcsisensor import ARCSIAbstractSensor


class ARCSILandsatOLISensor(ARCSIAbstractSensor):
    """
    A class which represents the landsat 8 sensor to read
    header parameters and apply data processing operations.
    """

    def __init__(self, debugMode, inputImage):
        ARCSIAbstractSensor.__init__(self, debugMode, inputImage)
        self.sensor = "LS_OLI"
        self.collection_num = 0

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

        self.b1CalMin = 0.0
        self.b1CalMax = 0.0
        self.b2CalMin = 0.0
        self.b2CalMax = 0.0
        self.b3CalMin = 0.0
        self.b3CalMax = 0.0
        self.b4CalMin = 0.0
        self.b4CalMax = 0.0
        self.b5CalMin = 0.0
        self.b5CalMax = 0.0
        self.b6CalMin = 0.0
        self.b6CalMax = 0.0
        self.b7CalMin = 0.0
        self.b7CalMax = 0.0
        self.b8CalMin = 0.0
        self.b8CalMax = 0.0
        self.b9CalMin = 0.0
        self.b9CalMax = 0.0
        self.b10CalMin = 0.0
        self.b10CalMax = 0.0
        self.b11CalMin = 0.0
        self.b11CalMax = 0.0

        self.k1ConstB10 = 0.0
        self.k1ConstB11 = 0.0
        self.k2ConstB10 = 0.0
        self.k2ConstB11 = 0.0

        self.sensorID = ""
        self.spacecraftID = ""
        self.cloudCover = 0.0
        self.cloudCoverLand = 0.0
        self.earthSunDistance = 0.0
        self.gridCellSizePan = 0.0
        self.gridCellSizeRefl = 0.0
        self.gridCellSizeTherm = 0.0

    def extractHeaderParameters(self, inputHeader, wktStr):
        """
        Understands and parses the Landsat MTL header files
        """
        try:
            if self.userSpInputImage is not None:
                raise ARCSIException(
                    "Landsat sensor cannot accept a user specified image file - only the images in the header file will be used."
                )
            self.headerFileName = os.path.split(inputHeader)[1]

            print("Reading header file")
            hFile = open(inputHeader, "r")
            headerParams = dict()
            for line in hFile:
                line = line.strip()
                if line:
                    lineVals = line.split("=")
                    if len(lineVals) == 2:
                        if (lineVals[0].strip() != "GROUP") or (
                            lineVals[0].strip() != "END_GROUP"
                        ):
                            headerParams[lineVals[0].strip()] = (
                                lineVals[1].strip().replace('"', "")
                            )
            hFile.close()
            print("Extracting Header Values")
            # Get the sensor info.
            if (
                (headerParams["SPACECRAFT_ID"].upper() == "LANDSAT_8")
                or (headerParams["SPACECRAFT_ID"].upper() == "LANDSAT8")
            ) and (headerParams["SENSOR_ID"].upper() == "OLI_TIRS"):
                self.sensor = "LS8"
            elif (
                (headerParams["SPACECRAFT_ID"].upper() == "LANDSAT_9")
                or (headerParams["SPACECRAFT_ID"].upper() == "LANDSAT9")
            ) and (headerParams["SENSOR_ID"].upper() == "OLI_TIRS"):
                self.sensor = "LS9"
            else:
                raise ARCSIException(
                    "Do no recognise the spacecraft and sensor or combination."
                )

            self.sensorID = headerParams["SENSOR_ID"]
            self.spacecraftID = headerParams["SPACECRAFT_ID"]

            if headerParams["COLLECTION_NUMBER"] == "01":
                self.collection_num = 1
            elif headerParams["COLLECTION_NUMBER"] == "02":
                self.collection_num = 2
            else:
                raise ARCSIException(
                    "Can only process collection 1 and 2 data: {}".format(
                        headerParams["COLLECTION_NUMBER"]
                    )
                )

            # Get row/path
            self.row = int(headerParams["WRS_ROW"])
            self.path = int(headerParams["WRS_PATH"])

            # Get date and time of the acquisition
            acData = headerParams["DATE_ACQUIRED"].split("-")
            acTime = headerParams["SCENE_CENTER_TIME"].split(":")
            secsTime = acTime[2].split(".")
            self.acquisitionTime = datetime.datetime(
                int(acData[0]),
                int(acData[1]),
                int(acData[2]),
                int(acTime[0]),
                int(acTime[1]),
                int(secsTime[0]),
            )

            self.solarZenith = 90 - rsgislib.tools.utils.str_to_float(
                headerParams["SUN_ELEVATION"]
            )
            self.solarAzimuth = rsgislib.tools.utils.str_to_float(
                headerParams["SUN_AZIMUTH"]
            )

            # Get the geographic lat/long corners of the image.
            self.latTL = rsgislib.tools.utils.str_to_float(
                headerParams["CORNER_UL_LAT_PRODUCT"]
            )
            self.lonTL = rsgislib.tools.utils.str_to_float(
                headerParams["CORNER_UL_LON_PRODUCT"]
            )
            self.latTR = rsgislib.tools.utils.str_to_float(
                headerParams["CORNER_UR_LAT_PRODUCT"]
            )
            self.lonTR = rsgislib.tools.utils.str_to_float(
                headerParams["CORNER_UR_LON_PRODUCT"]
            )
            self.latBL = rsgislib.tools.utils.str_to_float(
                headerParams["CORNER_LL_LAT_PRODUCT"]
            )
            self.lonBL = rsgislib.tools.utils.str_to_float(
                headerParams["CORNER_LL_LON_PRODUCT"]
            )
            self.latBR = rsgislib.tools.utils.str_to_float(
                headerParams["CORNER_LR_LAT_PRODUCT"]
            )
            self.lonBR = rsgislib.tools.utils.str_to_float(
                headerParams["CORNER_LR_LON_PRODUCT"]
            )

            # Get the projected X/Y corners of the image
            self.xTL = rsgislib.tools.utils.str_to_float(
                headerParams["CORNER_UL_PROJECTION_X_PRODUCT"]
            )
            self.yTL = rsgislib.tools.utils.str_to_float(
                headerParams["CORNER_UL_PROJECTION_Y_PRODUCT"]
            )
            self.xTR = rsgislib.tools.utils.str_to_float(
                headerParams["CORNER_UR_PROJECTION_X_PRODUCT"]
            )
            self.yTR = rsgislib.tools.utils.str_to_float(
                headerParams["CORNER_UR_PROJECTION_Y_PRODUCT"]
            )
            self.xBL = rsgislib.tools.utils.str_to_float(
                headerParams["CORNER_LL_PROJECTION_X_PRODUCT"]
            )
            self.yBL = rsgislib.tools.utils.str_to_float(
                headerParams["CORNER_LL_PROJECTION_Y_PRODUCT"]
            )
            self.xBR = rsgislib.tools.utils.str_to_float(
                headerParams["CORNER_LR_PROJECTION_X_PRODUCT"]
            )
            self.yBR = rsgislib.tools.utils.str_to_float(
                headerParams["CORNER_LR_PROJECTION_Y_PRODUCT"]
            )

            # Get projection
            inProj = osr.SpatialReference()
            if (
                (headerParams["MAP_PROJECTION"] == "UTM")
                and (headerParams["DATUM"] == "WGS84")
                and (headerParams["ELLIPSOID"] == "WGS84")
            ):
                utmZone = int(headerParams["UTM_ZONE"])
                utmCode = "WGS84UTM" + str(utmZone) + str("N")
                # print("UTM: ", utmCode)
                inProj.ImportFromEPSG(self.epsgCodes[utmCode])
            elif (
                (headerParams["MAP_PROJECTION"] == "PS")
                and (headerParams["DATUM"] == "WGS84")
                and (headerParams["ELLIPSOID"] == "WGS84")
            ):
                inProj.ImportFromWkt(
                    'PROJCS["PS WGS84", GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563, AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433],AUTHORITY["EPSG","4326"]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",-71],PARAMETER["central_meridian",0],PARAMETER["scale_factor",1],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]]]'
                )
            else:
                raise ARCSIException(
                    "Expecting Landsat to be projected in UTM or PolarStereographic (PS) with datum=WGS84 and ellipsoid=WGS84."
                )

            if self.inWKT == "":
                self.inWKT = inProj.ExportToWkt()

            # Check image is square!
            if not (
                (self.xTL == self.xBL)
                and (self.yTL == self.yTR)
                and (self.xTR == self.xBR)
                and (self.yBL == self.yBR)
            ):
                raise ARCSIException("Image is not square in projected coordinates.")

            self.xCentre = self.xTL + ((self.xTR - self.xTL) / 2)
            self.yCentre = self.yBR + ((self.yTL - self.yBR) / 2)

            (
                self.lonCentre,
                self.latCentre,
            ) = rsgislib.tools.geometrytools.reproj_point_to_wgs84(
                inProj, self.xCentre, self.yCentre
            )

            # print("Lat: " + str(self.latCentre) + " Long: " + str(self.lonCentre))

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

            if "FILE_NAME_BAND_QUALITY" in headerParams:
                self.bandQAFile = os.path.join(
                    filesDIR, headerParams["FILE_NAME_BAND_QUALITY"]
                )
            elif "FILE_NAME_QUALITY_L1_PIXEL" in headerParams:
                self.bandQAFile = os.path.join(
                    filesDIR, headerParams["FILE_NAME_QUALITY_L1_PIXEL"]
                )
            else:
                print(
                    "Warning - the quality band is not available. Are you using collection 1 or 2 data?"
                )
                self.bandQAFile = ""

            self.b1RadMulti = rsgislib.tools.utils.str_to_float(
                headerParams["RADIANCE_MULT_BAND_1"]
            )
            self.b2RadMulti = rsgislib.tools.utils.str_to_float(
                headerParams["RADIANCE_MULT_BAND_2"]
            )
            self.b3RadMulti = rsgislib.tools.utils.str_to_float(
                headerParams["RADIANCE_MULT_BAND_3"]
            )
            self.b4RadMulti = rsgislib.tools.utils.str_to_float(
                headerParams["RADIANCE_MULT_BAND_4"]
            )
            self.b5RadMulti = rsgislib.tools.utils.str_to_float(
                headerParams["RADIANCE_MULT_BAND_5"]
            )
            self.b6RadMulti = rsgislib.tools.utils.str_to_float(
                headerParams["RADIANCE_MULT_BAND_6"]
            )
            self.b7RadMulti = rsgislib.tools.utils.str_to_float(
                headerParams["RADIANCE_MULT_BAND_7"]
            )
            self.b8RadMulti = rsgislib.tools.utils.str_to_float(
                headerParams["RADIANCE_MULT_BAND_8"]
            )
            self.b9RadMulti = rsgislib.tools.utils.str_to_float(
                headerParams["RADIANCE_MULT_BAND_9"]
            )
            self.b10RadMulti = rsgislib.tools.utils.str_to_float(
                headerParams["RADIANCE_MULT_BAND_10"]
            )
            self.b11RadMulti = rsgislib.tools.utils.str_to_float(
                headerParams["RADIANCE_MULT_BAND_11"]
            )

            self.b1RadAdd = rsgislib.tools.utils.str_to_float(
                headerParams["RADIANCE_ADD_BAND_1"]
            )
            self.b2RadAdd = rsgislib.tools.utils.str_to_float(
                headerParams["RADIANCE_ADD_BAND_2"]
            )
            self.b3RadAdd = rsgislib.tools.utils.str_to_float(
                headerParams["RADIANCE_ADD_BAND_3"]
            )
            self.b4RadAdd = rsgislib.tools.utils.str_to_float(
                headerParams["RADIANCE_ADD_BAND_4"]
            )
            self.b5RadAdd = rsgislib.tools.utils.str_to_float(
                headerParams["RADIANCE_ADD_BAND_5"]
            )
            self.b6RadAdd = rsgislib.tools.utils.str_to_float(
                headerParams["RADIANCE_ADD_BAND_6"]
            )
            self.b7RadAdd = rsgislib.tools.utils.str_to_float(
                headerParams["RADIANCE_ADD_BAND_7"]
            )
            self.b8RadAdd = rsgislib.tools.utils.str_to_float(
                headerParams["RADIANCE_ADD_BAND_8"]
            )
            self.b9RadAdd = rsgislib.tools.utils.str_to_float(
                headerParams["RADIANCE_ADD_BAND_9"]
            )
            self.b10RadAdd = rsgislib.tools.utils.str_to_float(
                headerParams["RADIANCE_ADD_BAND_10"]
            )
            self.b11RadAdd = rsgislib.tools.utils.str_to_float(
                headerParams["RADIANCE_ADD_BAND_11"]
            )

            self.k1ConstB10 = rsgislib.tools.utils.str_to_float(
                headerParams["K1_CONSTANT_BAND_10"]
            )
            self.k1ConstB11 = rsgislib.tools.utils.str_to_float(
                headerParams["K1_CONSTANT_BAND_11"]
            )
            self.k2ConstB10 = rsgislib.tools.utils.str_to_float(
                headerParams["K2_CONSTANT_BAND_10"]
            )
            self.k2ConstB11 = rsgislib.tools.utils.str_to_float(
                headerParams["K2_CONSTANT_BAND_11"]
            )

            self.b1ReflMulti = rsgislib.tools.utils.str_to_float(
                headerParams["REFLECTANCE_MULT_BAND_1"]
            )
            self.b2ReflMulti = rsgislib.tools.utils.str_to_float(
                headerParams["REFLECTANCE_MULT_BAND_2"]
            )
            self.b3ReflMulti = rsgislib.tools.utils.str_to_float(
                headerParams["REFLECTANCE_MULT_BAND_3"]
            )
            self.b4ReflMulti = rsgislib.tools.utils.str_to_float(
                headerParams["REFLECTANCE_MULT_BAND_4"]
            )
            self.b5ReflMulti = rsgislib.tools.utils.str_to_float(
                headerParams["REFLECTANCE_MULT_BAND_5"]
            )
            self.b6ReflMulti = rsgislib.tools.utils.str_to_float(
                headerParams["REFLECTANCE_MULT_BAND_6"]
            )
            self.b7ReflMulti = rsgislib.tools.utils.str_to_float(
                headerParams["REFLECTANCE_MULT_BAND_7"]
            )
            self.b8ReflMulti = rsgislib.tools.utils.str_to_float(
                headerParams["REFLECTANCE_MULT_BAND_8"]
            )
            self.b9ReflMulti = rsgislib.tools.utils.str_to_float(
                headerParams["REFLECTANCE_MULT_BAND_9"]
            )

            self.b1ReflAdd = rsgislib.tools.utils.str_to_float(
                headerParams["REFLECTANCE_ADD_BAND_1"]
            )
            self.b2ReflAdd = rsgislib.tools.utils.str_to_float(
                headerParams["REFLECTANCE_ADD_BAND_2"]
            )
            self.b3ReflAdd = rsgislib.tools.utils.str_to_float(
                headerParams["REFLECTANCE_ADD_BAND_3"]
            )
            self.b4ReflAdd = rsgislib.tools.utils.str_to_float(
                headerParams["REFLECTANCE_ADD_BAND_4"]
            )
            self.b5ReflAdd = rsgislib.tools.utils.str_to_float(
                headerParams["REFLECTANCE_ADD_BAND_5"]
            )
            self.b6ReflAdd = rsgislib.tools.utils.str_to_float(
                headerParams["REFLECTANCE_ADD_BAND_6"]
            )
            self.b7ReflAdd = rsgislib.tools.utils.str_to_float(
                headerParams["REFLECTANCE_ADD_BAND_7"]
            )
            self.b8ReflAdd = rsgislib.tools.utils.str_to_float(
                headerParams["REFLECTANCE_ADD_BAND_8"]
            )
            self.b9ReflAdd = rsgislib.tools.utils.str_to_float(
                headerParams["REFLECTANCE_ADD_BAND_9"]
            )

            self.b1CalMin = rsgislib.tools.utils.str_to_float(
                headerParams["QUANTIZE_CAL_MIN_BAND_1"]
            )
            self.b1CalMax = rsgislib.tools.utils.str_to_float(
                headerParams["QUANTIZE_CAL_MAX_BAND_1"]
            )
            self.b2CalMin = rsgislib.tools.utils.str_to_float(
                headerParams["QUANTIZE_CAL_MIN_BAND_2"]
            )
            self.b2CalMax = rsgislib.tools.utils.str_to_float(
                headerParams["QUANTIZE_CAL_MAX_BAND_2"]
            )
            self.b3CalMin = rsgislib.tools.utils.str_to_float(
                headerParams["QUANTIZE_CAL_MIN_BAND_3"]
            )
            self.b3CalMax = rsgislib.tools.utils.str_to_float(
                headerParams["QUANTIZE_CAL_MAX_BAND_3"]
            )
            self.b4CalMin = rsgislib.tools.utils.str_to_float(
                headerParams["QUANTIZE_CAL_MIN_BAND_4"]
            )
            self.b4CalMax = rsgislib.tools.utils.str_to_float(
                headerParams["QUANTIZE_CAL_MAX_BAND_4"]
            )
            self.b5CalMin = rsgislib.tools.utils.str_to_float(
                headerParams["QUANTIZE_CAL_MIN_BAND_5"]
            )
            self.b5CalMax = rsgislib.tools.utils.str_to_float(
                headerParams["QUANTIZE_CAL_MAX_BAND_5"]
            )
            self.b6CalMin = rsgislib.tools.utils.str_to_float(
                headerParams["QUANTIZE_CAL_MIN_BAND_6"]
            )
            self.b6CalMax = rsgislib.tools.utils.str_to_float(
                headerParams["QUANTIZE_CAL_MAX_BAND_6"]
            )
            self.b7CalMin = rsgislib.tools.utils.str_to_float(
                headerParams["QUANTIZE_CAL_MIN_BAND_7"]
            )
            self.b7CalMax = rsgislib.tools.utils.str_to_float(
                headerParams["QUANTIZE_CAL_MAX_BAND_7"]
            )
            self.b8CalMin = rsgislib.tools.utils.str_to_float(
                headerParams["QUANTIZE_CAL_MIN_BAND_8"]
            )
            self.b8CalMax = rsgislib.tools.utils.str_to_float(
                headerParams["QUANTIZE_CAL_MAX_BAND_8"]
            )
            self.b9CalMin = rsgislib.tools.utils.str_to_float(
                headerParams["QUANTIZE_CAL_MIN_BAND_9"]
            )
            self.b9CalMax = rsgislib.tools.utils.str_to_float(
                headerParams["QUANTIZE_CAL_MAX_BAND_9"]
            )
            self.b10CalMin = rsgislib.tools.utils.str_to_float(
                headerParams["QUANTIZE_CAL_MIN_BAND_10"]
            )
            self.b10CalMax = rsgislib.tools.utils.str_to_float(
                headerParams["QUANTIZE_CAL_MAX_BAND_10"]
            )
            self.b11CalMin = rsgislib.tools.utils.str_to_float(
                headerParams["QUANTIZE_CAL_MIN_BAND_11"]
            )
            self.b11CalMax = rsgislib.tools.utils.str_to_float(
                headerParams["QUANTIZE_CAL_MAX_BAND_11"]
            )

            if "CLOUD_COVER" in headerParams:
                self.cloudCover = rsgislib.tools.utils.str_to_float(
                    headerParams["CLOUD_COVER"], 0.0
                )
            if "CLOUD_COVER_LAND" in headerParams:
                self.cloudCoverLand = rsgislib.tools.utils.str_to_float(
                    headerParams["CLOUD_COVER_LAND"], 0.0
                )
            if "EARTH_SUN_DISTANCE" in headerParams:
                self.earthSunDistance = rsgislib.tools.utils.str_to_float(
                    headerParams["EARTH_SUN_DISTANCE"], 0.0
                )
            if "GRID_CELL_SIZE_REFLECTIVE" in headerParams:
                self.gridCellSizeRefl = rsgislib.tools.utils.str_to_float(
                    headerParams["GRID_CELL_SIZE_REFLECTIVE"], 60.0
                )
            if "GRID_CELL_SIZE_THERMAL" in headerParams:
                self.gridCellSizeTherm = rsgislib.tools.utils.str_to_float(
                    headerParams["GRID_CELL_SIZE_THERMAL"], 30.0
                )
            if "GRID_CELL_SIZE_PANCHROMATIC" in headerParams:
                self.gridCellSizePan = rsgislib.tools.utils.str_to_float(
                    headerParams["GRID_CELL_SIZE_PANCHROMATIC"], 15.0
                )

            # Read MTL header into python dict for python-fmask
            self.fmaskMTLInfo = fmask.config.readMTLFile(inputHeader)

            if "FILE_DATE" in headerParams:
                fileDateStr = headerParams["FILE_DATE"].strip()
            else:
                fileDateStr = headerParams["DATE_PRODUCT_GENERATED"].strip()
            fileDateStr = fileDateStr.replace("Z", "")
            self.fileDateObj = datetime.datetime.strptime(
                fileDateStr, "%Y-%m-%dT%H:%M:%S"
            )

        except Exception as e:
            raise e

    def getSolarIrrStdSolarGeom(self):
        """
        Get Solar Azimuth and Zenith as standard geometry.
        Azimuth: N=0, E=90, S=180, W=270.
        """
        solarAz = rsgislib.imagecalibration.solarangles.get_solar_irr_convention_solar_azimuth_from_usgs(
            self.solarAzimuth
        )
        return (solarAz, self.solarZenith)

    def getSensorViewGeom(self):
        """
        Get sensor viewing angles
        returns (viewAzimuth, viewZenith)
        """
        return (0.0, 0.0)

    def generateOutputBaseName(self):
        """
        Provides an implementation for the landsat sensor
        """
        rowpath = "r" + str(self.row) + "p" + str(self.path)
        outname = self.defaultGenBaseOutFileName()
        outname = outname + str("_") + rowpath
        return outname

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
        Generate file metadata.
        """
        if outFilesDict is None:
            outFilesDict = dict()
        if calcdValuesDict is None:
            calcdValuesDict = dict()
        outJSONFilePath = os.path.join(outputPath, outputFileName)
        jsonData = self.getJSONDictDefaultMetaData(
            productsStr, validMaskImage, footprintCalc, calcdValuesDict, outFilesDict
        )
        sensorInfo = jsonData["SensorInfo"]
        sensorInfo["Row"] = self.row
        sensorInfo["Path"] = self.path
        sensorInfo["SensorID"] = self.sensorID
        sensorInfo["SpacecraftID"] = self.spacecraftID
        acqDict = jsonData["AcquasitionInfo"]
        acqDict["EarthSunDistance"] = self.earthSunDistance
        imgInfo = dict()
        imgInfo["CloudCover"] = self.cloudCover
        imgInfo["CloudCoverLand"] = self.cloudCoverLand
        imgInfo["CellSizePan"] = self.gridCellSizePan
        imgInfo["CellSizeRefl"] = self.gridCellSizeRefl
        imgInfo["CellSizeTherm"] = self.gridCellSizeTherm
        jsonData["ImageInfo"] = imgInfo

        with open(outJSONFilePath, "w") as outfile:
            json.dump(
                jsonData,
                outfile,
                sort_keys=True,
                indent=4,
                separators=(",", ": "),
                ensure_ascii=False,
            )

    def expectedImageDataPresent(self):
        imageDataPresent = True

        if not os.path.exists(self.band1File):
            imageDataPresent = False
        if not os.path.exists(self.band2File):
            imageDataPresent = False
        if not os.path.exists(self.band3File):
            imageDataPresent = False
        if not os.path.exists(self.band4File):
            imageDataPresent = False
        if not os.path.exists(self.band5File):
            imageDataPresent = False
        if not os.path.exists(self.band6File):
            imageDataPresent = False
        if not os.path.exists(self.band7File):
            imageDataPresent = False
        # if not os.path.exists(self.band8File):
        #    imageDataPresent = False
        if not os.path.exists(self.band9File):
            imageDataPresent = False
        if not os.path.exists(self.band10File):
            imageDataPresent = False
        if not os.path.exists(self.band11File):
            imageDataPresent = False
        # if not os.path.exists(self.bandQAFile):
        #    imageDataPresent = False

        return imageDataPresent

    def hasThermal(self):
        return True

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
        raise ARCSIException(
            "Landsat 8 does not provide any image masks, do not use the MASK option."
        )

    def mosaicImageTiles(self, outputPath):
        raise ARCSIException("Image data does not need mosaicking")

    def resampleImgRes(
        self,
        outputPath,
        resampleToLowResImg,
        resampleMethod=rsgislib.INTERP_CUBIC,
        multicore=False,
    ):
        raise ARCSIException("Image data does not need resampling")

    def sharpenLowResRadImgBands(self, inputImg, outputImage, outFormat):
        raise ARCSIException("Image sharpening is not available for this sensor.")

    def generateValidImageDataMask(
        self, outputPath, outputMaskName, viewAngleImg, outFormat
    ):
        print("Create the valid data mask")
        tmpBaseName = os.path.splitext(outputMaskName)[0]
        tmpValidPxlMsk = os.path.join(outputPath, tmpBaseName + "vldpxlmsk.kea")
        outputImage = os.path.join(outputPath, outputMaskName)
        inImages = [
            self.band1File,
            self.band2File,
            self.band3File,
            self.band4File,
            self.band5File,
            self.band6File,
            self.band7File,
            self.band10File,
            self.band11File,
        ]
        rsgislib.imageutils.gen_valid_mask(
            input_imgs=inImages,
            output_img=tmpValidPxlMsk,
            gdalformat="KEA",
            no_data_val=0.0,
        )
        rsgislib.rastergis.pop_rat_img_stats(tmpValidPxlMsk, True, False, True)
        # Check there is valid data
        ratDS = gdal.Open(tmpValidPxlMsk, gdal.GA_ReadOnly)
        Histogram = rat.readColumn(ratDS, "Histogram")
        ratDS = None
        if Histogram.shape[0] < 2:
            raise ARCSIException("There is no valid data in this image.")
        if not os.path.exists(viewAngleImg):
            print("Calculate Image Angles.")
            imgInfo = rios.fileinfo.ImageInfo(tmpValidPxlMsk)
            corners = fmask.landsatangles.findImgCorners(tmpValidPxlMsk, imgInfo)
            nadirLine = fmask.landsatangles.findNadirLine(corners)
            extentSunAngles = fmask.landsatangles.sunAnglesForExtent(
                imgInfo, self.fmaskMTLInfo
            )
            satAzimuth = fmask.landsatangles.satAzLeftRight(nadirLine)
            fmask.landsatangles.makeAnglesImage(
                tmpValidPxlMsk,
                viewAngleImg,
                nadirLine,
                extentSunAngles,
                satAzimuth,
                imgInfo,
            )
            dataset = gdal.Open(viewAngleImg, gdal.GA_Update)
            if not dataset is None:
                dataset.GetRasterBand(1).SetDescription("SatelliteAzimuth")
                dataset.GetRasterBand(2).SetDescription("SatelliteZenith")
                dataset.GetRasterBand(3).SetDescription("SolorAzimuth")
                dataset.GetRasterBand(4).SetDescription("SolorZenith")
            dataset = None
        rsgislib.imagecalc.band_math(
            outputImage,
            "(VA<14)&&(VM==1)?1:0",
            outFormat,
            rsgislib.TYPE_8UINT,
            [
                rsgislib.imagecalc.BandDefn("VA", viewAngleImg, 2),
                rsgislib.imagecalc.BandDefn("VM", tmpValidPxlMsk, 1),
            ],
        )

        rsgislib.imageutils.delete_gdal_layer(tmpValidPxlMsk)
        return outputImage

    def convertImageToRadiance(
        self, outputPath, outputReflName, outputThermalName, outFormat
    ):
        print("Converting to Radiance")
        outputReflImage = os.path.join(outputPath, outputReflName)
        outputThermalImage = None
        bandDefnSeq = list()

        lsBand = collections.namedtuple(
            "LSBand", ["band_name", "input_img", "img_band", "add_val", "multi_val"]
        )
        bandDefnSeq.append(
            lsBand(
                band_name="Coastal",
                input_img=self.band1File,
                img_band=1,
                add_val=self.b1RadAdd,
                multi_val=self.b1RadMulti,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="Blue",
                input_img=self.band2File,
                img_band=1,
                add_val=self.b2RadAdd,
                multi_val=self.b2RadMulti,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="Green",
                input_img=self.band3File,
                img_band=1,
                add_val=self.b3RadAdd,
                multi_val=self.b3RadMulti,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="Red",
                input_img=self.band4File,
                img_band=1,
                add_val=self.b4RadAdd,
                multi_val=self.b4RadMulti,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="NIR",
                input_img=self.band5File,
                img_band=1,
                add_val=self.b5RadAdd,
                multi_val=self.b5RadMulti,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="SWIR1",
                input_img=self.band6File,
                img_band=1,
                add_val=self.b6RadAdd,
                multi_val=self.b6RadMulti,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="SWIR2",
                input_img=self.band7File,
                img_band=1,
                add_val=self.b7RadAdd,
                multi_val=self.b7RadMulti,
            )
        )
        rsgislib.imagecalibration.landsat_to_radiance_multi_add(
            outputReflImage, outFormat, bandDefnSeq
        )

        if not outputThermalName == None:
            outputThermalImage = os.path.join(outputPath, outputThermalName)
            bandDefnSeq = list()
            lsBand = collections.namedtuple(
                "LSBand", ["band_name", "input_img", "img_band", "add_val", "multi_val"]
            )
            bandDefnSeq.append(
                lsBand(
                    band_name="ThermalB10",
                    input_img=self.band10File,
                    img_band=1,
                    add_val=self.b10RadAdd,
                    multi_val=self.b10RadMulti,
                )
            )
            bandDefnSeq.append(
                lsBand(
                    band_name="ThermalB11",
                    input_img=self.band11File,
                    img_band=1,
                    add_val=self.b11RadAdd,
                    multi_val=self.b11RadMulti,
                )
            )
            rsgislib.imagecalibration.landsat_to_radiance_multi_add(
                outputThermalImage, outFormat, bandDefnSeq
            )

        return outputReflImage, outputThermalImage

    def generateImageSaturationMask(self, outputPath, outputName, outFormat):
        print("Generate Saturation Image")
        outputImage = os.path.join(outputPath, outputName)

        lsBand = collections.namedtuple(
            "LSBand", ["band_name", "input_img", "img_band", "sat_val"]
        )
        bandDefnSeq = list()
        bandDefnSeq.append(
            lsBand(
                band_name="Coastal",
                input_img=self.band1File,
                img_band=1,
                sat_val=self.b1CalMax,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="Blue",
                input_img=self.band2File,
                img_band=1,
                sat_val=self.b2CalMax,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="Green",
                input_img=self.band3File,
                img_band=1,
                sat_val=self.b3CalMax,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="Red",
                input_img=self.band4File,
                img_band=1,
                sat_val=self.b4CalMax,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="NIR",
                input_img=self.band5File,
                img_band=1,
                sat_val=self.b5CalMax,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="SWIR1",
                input_img=self.band6File,
                img_band=1,
                sat_val=self.b6CalMax,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="SWIR2",
                input_img=self.band7File,
                img_band=1,
                sat_val=self.b7CalMax,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="ThermalB10",
                input_img=self.band10File,
                img_band=1,
                sat_val=self.b10CalMax,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="ThermalB11",
                input_img=self.band11File,
                img_band=1,
                sat_val=self.b11CalMax,
            )
        )

        rsgislib.imagecalibration.saturated_pixels_mask(
            outputImage, outFormat, bandDefnSeq
        )

        return outputImage

    def convertThermalToBrightness(
        self, inputRadImage, outputPath, outputName, outFormat, scaleFactor
    ):
        print("Converting to Thermal Brightness")
        outputThermalImage = os.path.join(outputPath, outputName)
        bandDefnSeq = list()

        lsBand = collections.namedtuple("LSBand", ["band_name", "img_band", "k1", "k2"])
        bandDefnSeq.append(
            lsBand(
                band_name="ThermalB10",
                img_band=1,
                k1=self.k1ConstB10,
                k2=self.k2ConstB10,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="ThermalB11",
                img_band=2,
                k1=self.k1ConstB11,
                k2=self.k2ConstB11,
            )
        )
        rsgislib.imagecalibration.landsat_thermal_rad_to_brightness(
            inputRadImage,
            outputThermalImage,
            outFormat,
            rsgislib.TYPE_32INT,
            scaleFactor,
            bandDefnSeq,
        )
        return outputThermalImage

    def convertImageToTOARefl(
        self, inputRadImage, outputPath, outputName, outFormat, scaleFactor
    ):
        print("Converting to TOA")
        outputImage = os.path.join(outputPath, outputName)
        solarIrradianceVals = list()
        IrrVal = collections.namedtuple("SolarIrradiance", ["irradiance"])
        solarIrradianceVals.append(IrrVal(irradiance=1876.61))
        solarIrradianceVals.append(IrrVal(irradiance=1970.03))
        solarIrradianceVals.append(IrrVal(irradiance=1848.9))
        solarIrradianceVals.append(IrrVal(irradiance=1571.3))
        solarIrradianceVals.append(IrrVal(irradiance=967.66))
        solarIrradianceVals.append(IrrVal(irradiance=245.73))
        solarIrradianceVals.append(IrrVal(irradiance=82.03))
        rsgislib.imagecalibration.radiance_to_toa_refl(
            inputRadImage,
            outputImage,
            outFormat,
            rsgislib.TYPE_16UINT,
            scaleFactor,
            self.acquisitionTime.year,
            self.acquisitionTime.month,
            self.acquisitionTime.day,
            self.solarZenith,
            solarIrradianceVals,
        )
        return outputImage

    def generateCloudMask(
        self,
        inputReflImage,
        inputSatImage,
        inputThermalImage,
        inputViewAngleImg,
        inputValidImg,
        outputPath,
        outputCloudName,
        outputCloudProb,
        outFormat,
        tmpPath,
        scaleFactor,
        cloud_msk_methods=None,
    ):
        import rsgislib.imageutils

        try:
            outputImage = os.path.join(outputPath, outputCloudName)
            tmpBaseName = os.path.splitext(outputCloudName)[0]
            tmpBaseDIR = os.path.join(tmpPath, tmpBaseName)

            tmpDIRExisted = True
            if not os.path.exists(tmpBaseDIR):
                os.makedirs(tmpBaseDIR)
                tmpDIRExisted = False

            if (cloud_msk_methods is None) or (cloud_msk_methods == "FMASK"):
                tmpFMaskOut = os.path.join(tmpBaseDIR, tmpBaseName + "_pyfmaskout.kea")

                # Create tmp TOA stack with Band 9 (Cirrus)
                tmpTOAB9Out = os.path.join(tmpBaseDIR, tmpBaseName + "_TOAB9.kea")
                toaExp = (
                    "((b1*"
                    + str(self.b9ReflMulti)
                    + ")+"
                    + str(self.b9ReflAdd)
                    + ")*"
                    + str(scaleFactor)
                )
                rsgislib.imagecalc.image_math(
                    self.band9File,
                    tmpTOAB9Out,
                    toaExp,
                    "KEA",
                    rsgislib.TYPE_16UINT,
                    False,
                )
                if not rsgislib.imageutils.do_gdal_layers_have_same_proj(
                    inputReflImage, tmpTOAB9Out
                ):
                    tmpTOAB9OutNotProj = tmpTOAB9Out
                    tmpTOAB9Out = os.path.join(
                        tmpBaseDIR, tmpBaseName + "_TOAB9_reproj.kea"
                    )
                    rsgislib.imageutils.resample_img_to_match(
                        inputReflImage,
                        tmpTOAB9OutNotProj,
                        tmpTOAB9Out,
                        "KEA",
                        rsgislib.INTERP_CUBIC,
                        rsgislib.TYPE_16UINT,
                    )

                tmpReflStackOut = os.path.join(
                    tmpBaseDIR, tmpBaseName + "_TOAreflStack.kea"
                )
                rsgislib.imageutils.stack_img_bands(
                    [inputReflImage, tmpTOAB9Out],
                    None,
                    tmpReflStackOut,
                    None,
                    0,
                    "KEA",
                    rsgislib.TYPE_16UINT,
                )

                tmpThermStack = os.path.join(
                    tmpBaseDIR, tmpBaseName + "_ThermB10B11Stack.kea"
                )
                rsgislib.imageutils.stack_img_bands(
                    [self.band10File, self.band11File],
                    None,
                    tmpThermStack,
                    None,
                    0,
                    "KEA",
                    rsgislib.TYPE_16UINT,
                )

                tmpThermalLayer = tmpThermStack
                if not rsgislib.imageutils.do_gdal_layers_have_same_proj(
                    inputThermalImage, tmpThermStack
                ):
                    tmpThermalLayer = os.path.join(
                        tmpBaseDIR, tmpBaseName + "_thermalresample.kea"
                    )
                    rsgislib.imageutils.resample_img_to_match(
                        inputThermalImage,
                        tmpThermStack,
                        tmpThermalLayer,
                        "KEA",
                        rsgislib.INTERP_CUBIC,
                        rsgislib.TYPE_32FLOAT,
                    )

                minCloudSize = 0
                cloudBufferDistance = 150
                shadowBufferDistance = 300

                fmaskFilenames = fmask.config.FmaskFilenames()
                fmaskFilenames.setTOAReflectanceFile(tmpReflStackOut)
                fmaskFilenames.setThermalFile(tmpThermalLayer)
                fmaskFilenames.setSaturationMask(inputSatImage)
                fmaskFilenames.setOutputCloudMaskFile(tmpFMaskOut)

                thermalGain1040um = self.b10RadMulti
                thermalOffset1040um = self.b10RadAdd
                thermalBand1040um = 0
                thermalInfo = fmask.config.ThermalFileInfo(
                    thermalBand1040um,
                    thermalGain1040um,
                    thermalOffset1040um,
                    self.k1ConstB10,
                    self.k2ConstB10,
                )

                anglesInfo = fmask.config.AnglesFileInfo(
                    inputViewAngleImg,
                    3,
                    inputViewAngleImg,
                    2,
                    inputViewAngleImg,
                    1,
                    inputViewAngleImg,
                    0,
                )

                fmaskConfig = fmask.config.FmaskConfig(fmask.config.FMASK_LANDSAT8)
                fmaskConfig.setTOARefScaling(float(scaleFactor))
                fmaskConfig.setThermalInfo(thermalInfo)
                fmaskConfig.setAnglesInfo(anglesInfo)
                fmaskConfig.setKeepIntermediates(False)
                fmaskConfig.setVerbose(True)
                fmaskConfig.setTempDir(tmpBaseDIR)
                fmaskConfig.setMinCloudSize(minCloudSize)
                fmaskConfig.setEqn17CloudProbThresh(
                    fmask.config.FmaskConfig.Eqn17CloudProbThresh
                )
                fmaskConfig.setEqn20NirSnowThresh(
                    fmask.config.FmaskConfig.Eqn20NirSnowThresh
                )
                fmaskConfig.setEqn20GreenSnowThresh(
                    fmask.config.FmaskConfig.Eqn20GreenSnowThresh
                )

                # Work out a suitable buffer size, in pixels, dependent on the resolution of the input TOA image
                toaImgInfo = rios.fileinfo.ImageInfo(inputReflImage)
                fmaskConfig.setCloudBufferSize(
                    int(cloudBufferDistance / toaImgInfo.xRes)
                )
                fmaskConfig.setShadowBufferSize(
                    int(shadowBufferDistance / toaImgInfo.xRes)
                )

                fmask.fmask.doFmask(fmaskFilenames, fmaskConfig)

                rsgislib.imagecalc.image_math(
                    tmpFMaskOut,
                    outputImage,
                    "(b1==2)?1:(b1==3)?2:0",
                    outFormat,
                    rsgislib.TYPE_8UINT,
                )

            elif cloud_msk_methods == "LSMSK":
                if (self.bandQAFile == "") or (not os.path.exists(self.bandQAFile)):
                    raise ARCSIException(
                        "The QA band is not present - cannot use this for cloud masking."
                    )

                bqa_img_file = self.bandQAFile
                if not rsgislib.imageutils.do_gdal_layers_have_same_proj(
                    bqa_img_file, inputReflImage
                ):
                    bqa_img_file = os.path.join(tmpBaseDIR, tmpBaseName + "_BQA.kea")
                    rsgislib.imageutils.resample_img_to_match(
                        inputReflImage,
                        self.bandQAFile,
                        bqa_img_file,
                        "KEA",
                        rsgislib.INTERP_NEAREST_NEIGHBOUR,
                        rsgislib.TYPE_16UINT,
                        no_data_val=0,
                        multicore=False,
                    )

                if self.collection_num == 1:
                    exp = (
                        "(b1==2800)||(b1==2804)||(b1==2808)||(b1==2812)||(b1==6896)||(b1==6900)||(b1==6904)||(b1==6908)?1:"
                        "(b1==2976)||(b1==2980)||(b1==2984)||(b1==2988)||(b1==3008)||(b1==3012)||(b1==3016)||(b1==3020)||"
                        "(b1==7072)||(b1==7076)||(b1==7080)||(b1==7084)||(b1==7104)||(b1==7108)||(b1==7112)||(b1==7116)?2:0"
                    )
                    rsgislib.imagecalc.image_math(
                        bqa_img_file, outputImage, exp, outFormat, rsgislib.TYPE_8UINT
                    )
                elif self.collection_num == 2:
                    import rsgislib.imagecalibration.sensorlvl2data

                    c2_bqa_ind_img_file = os.path.join(
                        tmpBaseDIR, tmpBaseName + "c2_qa_ind_bands.kea"
                    )
                    rsgislib.imagecalibration.sensorlvl2data.parse_landsat_c2_qa_pixel_img(
                        bqa_img_file, c2_bqa_ind_img_file, gdalformat="KEA"
                    )
                    band_defns = list()
                    band_defns.append(
                        rsgislib.imagecalc.BandDefn(
                            "DilatedCloud", c2_bqa_ind_img_file, 2
                        )
                    )
                    band_defns.append(
                        rsgislib.imagecalc.BandDefn("Cloud", c2_bqa_ind_img_file, 4)
                    )
                    band_defns.append(
                        rsgislib.imagecalc.BandDefn(
                            "CloudShadow", c2_bqa_ind_img_file, 5
                        )
                    )
                    rsgislib.imagecalc.band_math(
                        outputImage,
                        "(DilatedCloud == 1)||(Cloud == 1)?1:(CloudShadow == 1)?2:0",
                        "KEA",
                        rsgislib.TYPE_8UINT,
                        band_defns,
                    )
                else:
                    raise ARCSIException(
                        "Can only read Collection 1 and 2 cloud masks."
                    )

            else:
                raise ARCSIException(
                    "Landsat only has FMASK and LSMSK cloud masking options; option provided is unknown."
                )

            if outFormat == "KEA":
                rsgislib.rastergis.pop_rat_img_stats(outputImage, True, True)
                ratDataset = gdal.Open(outputImage, gdal.GA_Update)
                red = rat.readColumn(ratDataset, "Red")
                green = rat.readColumn(ratDataset, "Green")
                blue = rat.readColumn(ratDataset, "Blue")
                ClassName = numpy.empty_like(red, dtype=numpy.dtype("a255"))

                red[0] = 0
                green[0] = 0
                blue[0] = 0

                if (red.shape[0] == 2) or (red.shape[0] == 3):
                    red[1] = 0
                    green[1] = 0
                    blue[1] = 255
                    ClassName[1] = "Clouds"

                    if red.shape[0] == 3:
                        red[2] = 0
                        green[2] = 255
                        blue[2] = 255
                        ClassName[2] = "Shadows"

                rat.writeColumn(ratDataset, "Red", red)
                rat.writeColumn(ratDataset, "Green", green)
                rat.writeColumn(ratDataset, "Blue", blue)
                rat.writeColumn(ratDataset, "ClassName", ClassName)
                ratDataset = None
            rsgislib.imageutils.copy_proj_from_img(outputImage, inputReflImage)

            if not self.debugMode:
                if not tmpDIRExisted:
                    shutil.rmtree(tmpBaseDIR, ignore_errors=True)

            return outputImage, None
        except Exception as e:
            raise e

    def createCloudMaskDataArray(self, inImgDataArr):
        return inImgDataArr

    def defineDarkShadowImageBand(self):
        return 5

    def calc6SCoefficients(
        self, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF
    ):
        sixsCoeffs = numpy.zeros((7, 6), dtype=numpy.float32)
        # Set up 6S model
        s = Py6S.SixS()
        s.atmos_profile = atmosProfile
        s.aero_profile = aeroProfile
        # s.ground_reflectance = Py6S.GroundReflectance.HomogeneousHapke(0.101, -0.263, 0.589, 0.046)
        s.ground_reflectance = grdRefl
        s.geometry = Py6S.Geometry.Landsat_TM()
        s.geometry.month = self.acquisitionTime.month
        s.geometry.day = self.acquisitionTime.day
        s.geometry.gmt_decimal_hour = (
            float(self.acquisitionTime.hour) + float(self.acquisitionTime.minute) / 60.0
        )
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
        s.wavelength = Py6S.Wavelength(
            0.427,
            0.4595,
            [
                0.000073,
                0.001628,
                0.024767,
                0.254149,
                0.908749,
                0.977393,
                0.986713,
                0.993137,
                0.982780,
                0.905808,
                0.226412,
                0.036603,
                0.002414,
                0.000255,
            ],
        )
        s.run()
        sixsCoeffs[0, 0] = float(s.outputs.values["coef_xa"])
        sixsCoeffs[0, 1] = float(s.outputs.values["coef_xb"])
        sixsCoeffs[0, 2] = float(s.outputs.values["coef_xc"])
        sixsCoeffs[0, 3] = float(s.outputs.values["direct_solar_irradiance"])
        sixsCoeffs[0, 4] = float(s.outputs.values["diffuse_solar_irradiance"])
        sixsCoeffs[0, 5] = float(s.outputs.values["environmental_irradiance"])

        # Band 2
        s.wavelength = Py6S.Wavelength(
            0.436,
            0.5285,
            [
                0.000010,
                0.000117,
                0.000455,
                0.001197,
                0.006869,
                0.027170,
                0.271370,
                0.723971,
                0.903034,
                0.909880,
                0.889667,
                0.877453,
                0.879688,
                0.891913,
                0.848533,
                0.828339,
                0.868497,
                0.912538,
                0.931726,
                0.954248,
                0.956424,
                0.978564,
                0.989469,
                0.968801,
                0.988729,
                0.967361,
                0.966125,
                0.981834,
                0.963135,
                0.996498,
                0.844893,
                0.190738,
                0.005328,
                0.001557,
                0.000516,
                0.000162,
                0.000023,
                -0.000016,
            ],
        )
        s.run()
        sixsCoeffs[1, 0] = float(s.outputs.values["coef_xa"])
        sixsCoeffs[1, 1] = float(s.outputs.values["coef_xb"])
        sixsCoeffs[1, 2] = float(s.outputs.values["coef_xc"])
        sixsCoeffs[1, 3] = float(s.outputs.values["direct_solar_irradiance"])
        sixsCoeffs[1, 4] = float(s.outputs.values["diffuse_solar_irradiance"])
        sixsCoeffs[1, 5] = float(s.outputs.values["environmental_irradiance"])

        # Band 3
        s.wavelength = Py6S.Wavelength(
            0.512,
            0.6095,
            [
                -0.000046,
                0.00011,
                0.000648,
                0.001332,
                0.003446,
                0.007024,
                0.025513,
                0.070551,
                0.353885,
                0.741205,
                0.954627,
                0.959215,
                0.969873,
                0.961397,
                0.977001,
                0.990784,
                0.982642,
                0.977765,
                0.946245,
                0.959038,
                0.966447,
                0.958314,
                0.983397,
                0.974522,
                0.978208,
                0.974392,
                0.969181,
                0.982956,
                0.968886,
                0.986657,
                0.904478,
                0.684974,
                0.190467,
                0.035393,
                0.002574,
                0.000394,
                -0.000194,
                -0.000292,
                -0.000348,
                -0.000351,
            ],
        )
        s.run()
        sixsCoeffs[2, 0] = float(s.outputs.values["coef_xa"])
        sixsCoeffs[2, 1] = float(s.outputs.values["coef_xb"])
        sixsCoeffs[2, 2] = float(s.outputs.values["coef_xc"])
        sixsCoeffs[2, 3] = float(s.outputs.values["direct_solar_irradiance"])
        sixsCoeffs[2, 4] = float(s.outputs.values["diffuse_solar_irradiance"])
        sixsCoeffs[2, 5] = float(s.outputs.values["environmental_irradiance"])

        # Band 4
        s.wavelength = Py6S.Wavelength(
            0.625,
            0.690,
            [
                -0.000342,
                0.000895,
                0.007197,
                0.030432,
                0.299778,
                0.764443,
                0.950823,
                0.951831,
                0.984173,
                0.983434,
                0.959441,
                0.955548,
                0.981688,
                0.992388,
                0.97696,
                0.98108,
                0.980678,
                0.962154,
                0.966928,
                0.848855,
                0.123946,
                0.017702,
                0.001402,
                0.000117,
                -0.000376,
                -0.000458,
                -0.000429,
            ],
        )
        s.run()
        sixsCoeffs[3, 0] = float(s.outputs.values["coef_xa"])
        sixsCoeffs[3, 1] = float(s.outputs.values["coef_xb"])
        sixsCoeffs[3, 2] = float(s.outputs.values["coef_xc"])
        sixsCoeffs[3, 3] = float(s.outputs.values["direct_solar_irradiance"])
        sixsCoeffs[3, 4] = float(s.outputs.values["diffuse_solar_irradiance"])
        sixsCoeffs[3, 5] = float(s.outputs.values["environmental_irradiance"])

        # Band 5
        s.wavelength = Py6S.Wavelength(
            0.829,
            0.899,
            [
                -0.000034,
                0.000050,
                0.000314,
                0.000719,
                0.002107,
                0.004744,
                0.017346,
                0.048191,
                0.249733,
                0.582623,
                0.960215,
                0.973133,
                1.000000,
                0.980733,
                0.957357,
                0.947044,
                0.948450,
                0.950632,
                0.969821,
                0.891066,
                0.448364,
                0.174619,
                0.034532,
                0.012440,
                0.002944,
                0.001192,
                0.000241,
                0.000044,
                -0.000084,
            ],
        )
        s.run()
        sixsCoeffs[4, 0] = float(s.outputs.values["coef_xa"])
        sixsCoeffs[4, 1] = float(s.outputs.values["coef_xb"])
        sixsCoeffs[4, 2] = float(s.outputs.values["coef_xc"])
        sixsCoeffs[4, 3] = float(s.outputs.values["direct_solar_irradiance"])
        sixsCoeffs[4, 4] = float(s.outputs.values["diffuse_solar_irradiance"])
        sixsCoeffs[4, 5] = float(s.outputs.values["environmental_irradiance"])

        # Band 6
        s.wavelength = Py6S.Wavelength(
            1.515,
            1.6975,
            [
                -0.00002,
                0.00015,
                0.00047,
                0.00076,
                0.00137,
                0.00186,
                0.00288,
                0.00377,
                0.00553,
                0.00732,
                0.01099,
                0.01430,
                0.02183,
                0.02995,
                0.04786,
                0.06573,
                0.10189,
                0.13864,
                0.22026,
                0.29136,
                0.42147,
                0.52568,
                0.67668,
                0.75477,
                0.85407,
                0.89183,
                0.91301,
                0.92295,
                0.92641,
                0.92368,
                0.92283,
                0.92206,
                0.92661,
                0.94253,
                0.94618,
                0.94701,
                0.95286,
                0.94967,
                0.95905,
                0.96005,
                0.96147,
                0.96018,
                0.96470,
                0.96931,
                0.97691,
                0.98126,
                0.98861,
                0.99802,
                0.99964,
                0.99344,
                0.96713,
                0.93620,
                0.84097,
                0.75189,
                0.57323,
                0.45197,
                0.29175,
                0.21115,
                0.12846,
                0.09074,
                0.05275,
                0.03731,
                0.02250,
                0.01605,
                0.00959,
                0.00688,
                0.00426,
                0.00306,
                0.00178,
                0.00124,
                0.00068,
                0.00041,
                0.00011,
                -0.00003,
            ],
        )
        s.run()
        sixsCoeffs[5, 0] = float(s.outputs.values["coef_xa"])
        sixsCoeffs[5, 1] = float(s.outputs.values["coef_xb"])
        sixsCoeffs[5, 2] = float(s.outputs.values["coef_xc"])
        sixsCoeffs[5, 3] = float(s.outputs.values["direct_solar_irradiance"])
        sixsCoeffs[5, 4] = float(s.outputs.values["diffuse_solar_irradiance"])
        sixsCoeffs[5, 5] = float(s.outputs.values["environmental_irradiance"])

        # Band 7
        s.wavelength = Py6S.Wavelength(
            2.037,
            2.3545,
            [
                -0.000010,
                0.000083,
                0.000240,
                0.000368,
                0.000599,
                0.000814,
                0.001222,
                0.001546,
                0.002187,
                0.002696,
                0.003733,
                0.004627,
                0.006337,
                0.007996,
                0.011005,
                0.013610,
                0.018899,
                0.023121,
                0.032071,
                0.040206,
                0.056429,
                0.070409,
                0.100640,
                0.128292,
                0.179714,
                0.227234,
                0.311347,
                0.377044,
                0.488816,
                0.554715,
                0.663067,
                0.722284,
                0.792667,
                0.836001,
                0.867845,
                0.886411,
                0.906527,
                0.911091,
                0.929693,
                0.936544,
                0.942952,
                0.943194,
                0.948776,
                0.949643,
                0.956635,
                0.947423,
                0.950874,
                0.947014,
                0.957717,
                0.946412,
                0.951641,
                0.948644,
                0.940311,
                0.947923,
                0.938737,
                0.941859,
                0.944482,
                0.951661,
                0.939939,
                0.935493,
                0.938955,
                0.929162,
                0.930508,
                0.933908,
                0.936472,
                0.933523,
                0.946217,
                0.955661,
                0.963135,
                0.964365,
                0.962905,
                0.962473,
                0.957814,
                0.958041,
                0.951706,
                0.960212,
                0.947696,
                0.959060,
                0.955750,
                0.953245,
                0.966786,
                0.960173,
                0.977637,
                0.982760,
                0.985056,
                0.999600,
                0.992469,
                0.995894,
                0.997261,
                0.991127,
                0.986037,
                0.984536,
                0.972794,
                0.976540,
                0.974409,
                0.967502,
                0.955095,
                0.955588,
                0.922405,
                0.894940,
                0.823876,
                0.744025,
                0.602539,
                0.502693,
                0.355569,
                0.278260,
                0.186151,
                0.141435,
                0.092029,
                0.069276,
                0.046332,
                0.035634,
                0.024000,
                0.018688,
                0.012930,
                0.010155,
                0.007088,
                0.005643,
                0.003903,
                0.003025,
                0.002047,
                0.001554,
                0.000974,
                0.000680,
                0.000320,
                0.000119,
                -0.000134,
                -0.000263,
            ],
        )
        s.run()
        sixsCoeffs[6, 0] = float(s.outputs.values["coef_xa"])
        sixsCoeffs[6, 1] = float(s.outputs.values["coef_xb"])
        sixsCoeffs[6, 2] = float(s.outputs.values["coef_xc"])
        sixsCoeffs[6, 3] = float(s.outputs.values["direct_solar_irradiance"])
        sixsCoeffs[6, 4] = float(s.outputs.values["diffuse_solar_irradiance"])
        sixsCoeffs[6, 5] = float(s.outputs.values["environmental_irradiance"])

        return sixsCoeffs

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
        print("Converting to Surface Reflectance")
        outputImage = os.path.join(outputPath, outputName)

        imgBandCoeffs = list()

        sixsCoeffs = self.calc6SCoefficients(
            aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF
        )

        imgBandCoeffs.append(
            rsgislib.imagecalibration.Band6SCoeff(
                band=1,
                aX=float(sixsCoeffs[0, 0]),
                bX=float(sixsCoeffs[0, 1]),
                cX=float(sixsCoeffs[0, 2]),
                DirIrr=float(sixsCoeffs[0, 3]),
                DifIrr=float(sixsCoeffs[0, 4]),
                EnvIrr=float(sixsCoeffs[0, 5]),
            )
        )
        imgBandCoeffs.append(
            rsgislib.imagecalibration.Band6SCoeff(
                band=2,
                aX=float(sixsCoeffs[1, 0]),
                bX=float(sixsCoeffs[1, 1]),
                cX=float(sixsCoeffs[1, 2]),
                DirIrr=float(sixsCoeffs[1, 3]),
                DifIrr=float(sixsCoeffs[1, 4]),
                EnvIrr=float(sixsCoeffs[1, 5]),
            )
        )
        imgBandCoeffs.append(
            rsgislib.imagecalibration.Band6SCoeff(
                band=3,
                aX=float(sixsCoeffs[2, 0]),
                bX=float(sixsCoeffs[2, 1]),
                cX=float(sixsCoeffs[2, 2]),
                DirIrr=float(sixsCoeffs[2, 3]),
                DifIrr=float(sixsCoeffs[2, 4]),
                EnvIrr=float(sixsCoeffs[2, 5]),
            )
        )
        imgBandCoeffs.append(
            rsgislib.imagecalibration.Band6SCoeff(
                band=4,
                aX=float(sixsCoeffs[3, 0]),
                bX=float(sixsCoeffs[3, 1]),
                cX=float(sixsCoeffs[3, 2]),
                DirIrr=float(sixsCoeffs[3, 3]),
                DifIrr=float(sixsCoeffs[3, 4]),
                EnvIrr=float(sixsCoeffs[3, 5]),
            )
        )
        imgBandCoeffs.append(
            rsgislib.imagecalibration.Band6SCoeff(
                band=5,
                aX=float(sixsCoeffs[4, 0]),
                bX=float(sixsCoeffs[4, 1]),
                cX=float(sixsCoeffs[4, 2]),
                DirIrr=float(sixsCoeffs[4, 3]),
                DifIrr=float(sixsCoeffs[4, 4]),
                EnvIrr=float(sixsCoeffs[4, 5]),
            )
        )
        imgBandCoeffs.append(
            rsgislib.imagecalibration.Band6SCoeff(
                band=6,
                aX=float(sixsCoeffs[5, 0]),
                bX=float(sixsCoeffs[5, 1]),
                cX=float(sixsCoeffs[5, 2]),
                DirIrr=float(sixsCoeffs[5, 3]),
                DifIrr=float(sixsCoeffs[5, 4]),
                EnvIrr=float(sixsCoeffs[5, 5]),
            )
        )
        imgBandCoeffs.append(
            rsgislib.imagecalibration.Band6SCoeff(
                band=7,
                aX=float(sixsCoeffs[6, 0]),
                bX=float(sixsCoeffs[6, 1]),
                cX=float(sixsCoeffs[6, 2]),
                DirIrr=float(sixsCoeffs[6, 3]),
                DifIrr=float(sixsCoeffs[6, 4]),
                EnvIrr=float(sixsCoeffs[6, 5]),
            )
        )

        rsgislib.imagecalibration.apply_6s_coeff_single_param(
            inputRadImage,
            outputImage,
            outFormat,
            rsgislib.TYPE_16UINT,
            scaleFactor,
            0,
            True,
            imgBandCoeffs,
        )
        return outputImage

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
        print("Converting to Surface Reflectance")
        outputImage = os.path.join(outputPath, outputName)

        if elevCoeffs is None:
            print("Build an LUT for elevation values.")
            elev6SCoeffsLUT = self.buildElevation6SCoeffLUT(
                aeroProfile,
                atmosProfile,
                grdRefl,
                aotVal,
                useBRDF,
                surfaceAltitudeMin,
                surfaceAltitudeMax,
            )
            print("LUT has been built.")

            elevCoeffs = list()
            for elevLUT in elev6SCoeffsLUT:
                imgBandCoeffs = list()
                sixsCoeffs = elevLUT.Coeffs
                elevVal = elevLUT.Elev
                imgBandCoeffs.append(
                    rsgislib.imagecalibration.Band6SCoeff(
                        band=1,
                        aX=float(sixsCoeffs[0, 0]),
                        bX=float(sixsCoeffs[0, 1]),
                        cX=float(sixsCoeffs[0, 2]),
                        DirIrr=float(sixsCoeffs[0, 3]),
                        DifIrr=float(sixsCoeffs[0, 4]),
                        EnvIrr=float(sixsCoeffs[0, 5]),
                    )
                )
                imgBandCoeffs.append(
                    rsgislib.imagecalibration.Band6SCoeff(
                        band=2,
                        aX=float(sixsCoeffs[1, 0]),
                        bX=float(sixsCoeffs[1, 1]),
                        cX=float(sixsCoeffs[1, 2]),
                        DirIrr=float(sixsCoeffs[1, 3]),
                        DifIrr=float(sixsCoeffs[1, 4]),
                        EnvIrr=float(sixsCoeffs[1, 5]),
                    )
                )
                imgBandCoeffs.append(
                    rsgislib.imagecalibration.Band6SCoeff(
                        band=3,
                        aX=float(sixsCoeffs[2, 0]),
                        bX=float(sixsCoeffs[2, 1]),
                        cX=float(sixsCoeffs[2, 2]),
                        DirIrr=float(sixsCoeffs[2, 3]),
                        DifIrr=float(sixsCoeffs[2, 4]),
                        EnvIrr=float(sixsCoeffs[2, 5]),
                    )
                )
                imgBandCoeffs.append(
                    rsgislib.imagecalibration.Band6SCoeff(
                        band=4,
                        aX=float(sixsCoeffs[3, 0]),
                        bX=float(sixsCoeffs[3, 1]),
                        cX=float(sixsCoeffs[3, 2]),
                        DirIrr=float(sixsCoeffs[3, 3]),
                        DifIrr=float(sixsCoeffs[3, 4]),
                        EnvIrr=float(sixsCoeffs[3, 5]),
                    )
                )
                imgBandCoeffs.append(
                    rsgislib.imagecalibration.Band6SCoeff(
                        band=5,
                        aX=float(sixsCoeffs[4, 0]),
                        bX=float(sixsCoeffs[4, 1]),
                        cX=float(sixsCoeffs[4, 2]),
                        DirIrr=float(sixsCoeffs[4, 3]),
                        DifIrr=float(sixsCoeffs[4, 4]),
                        EnvIrr=float(sixsCoeffs[4, 5]),
                    )
                )
                imgBandCoeffs.append(
                    rsgislib.imagecalibration.Band6SCoeff(
                        band=6,
                        aX=float(sixsCoeffs[5, 0]),
                        bX=float(sixsCoeffs[5, 1]),
                        cX=float(sixsCoeffs[5, 2]),
                        DirIrr=float(sixsCoeffs[5, 3]),
                        DifIrr=float(sixsCoeffs[5, 4]),
                        EnvIrr=float(sixsCoeffs[5, 5]),
                    )
                )
                imgBandCoeffs.append(
                    rsgislib.imagecalibration.Band6SCoeff(
                        band=7,
                        aX=float(sixsCoeffs[6, 0]),
                        bX=float(sixsCoeffs[6, 1]),
                        cX=float(sixsCoeffs[6, 2]),
                        DirIrr=float(sixsCoeffs[6, 3]),
                        DifIrr=float(sixsCoeffs[6, 4]),
                        EnvIrr=float(sixsCoeffs[6, 5]),
                    )
                )

                elevCoeffs.append(
                    rsgislib.imagecalibration.ElevLUTFeat(
                        Elev=float(elevVal), Coeffs=imgBandCoeffs
                    )
                )

        rsgislib.imagecalibration.apply_6s_coeff_elev_lut_param(
            inputRadImage,
            inputDEMFile,
            outputImage,
            outFormat,
            rsgislib.TYPE_16UINT,
            scaleFactor,
            0,
            True,
            elevCoeffs,
        )
        return outputImage, elevCoeffs

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
        print("Converting to Surface Reflectance")
        outputImage = os.path.join(outputPath, outputName)

        if elevAOTCoeffs is None:
            print("Build an LUT for elevation and AOT values.")
            elevAOT6SCoeffsLUT = self.buildElevationAOT6SCoeffLUT(
                aeroProfile,
                atmosProfile,
                grdRefl,
                useBRDF,
                surfaceAltitudeMin,
                surfaceAltitudeMax,
                aotMin,
                aotMax,
            )

            elevAOTCoeffs = list()
            for elevLUT in elevAOT6SCoeffsLUT:
                elevVal = elevLUT.Elev
                aotLUT = elevLUT.Coeffs
                aot6SCoeffsOut = list()
                for aotFeat in aotLUT:
                    sixsCoeffs = aotFeat.Coeffs
                    aotVal = aotFeat.AOT
                    imgBandCoeffs = list()
                    imgBandCoeffs.append(
                        rsgislib.imagecalibration.Band6SCoeff(
                            band=1,
                            aX=float(sixsCoeffs[0, 0]),
                            bX=float(sixsCoeffs[0, 1]),
                            cX=float(sixsCoeffs[0, 2]),
                            DirIrr=float(sixsCoeffs[0, 3]),
                            DifIrr=float(sixsCoeffs[0, 4]),
                            EnvIrr=float(sixsCoeffs[0, 5]),
                        )
                    )
                    imgBandCoeffs.append(
                        rsgislib.imagecalibration.Band6SCoeff(
                            band=2,
                            aX=float(sixsCoeffs[1, 0]),
                            bX=float(sixsCoeffs[1, 1]),
                            cX=float(sixsCoeffs[1, 2]),
                            DirIrr=float(sixsCoeffs[1, 3]),
                            DifIrr=float(sixsCoeffs[1, 4]),
                            EnvIrr=float(sixsCoeffs[1, 5]),
                        )
                    )
                    imgBandCoeffs.append(
                        rsgislib.imagecalibration.Band6SCoeff(
                            band=3,
                            aX=float(sixsCoeffs[2, 0]),
                            bX=float(sixsCoeffs[2, 1]),
                            cX=float(sixsCoeffs[2, 2]),
                            DirIrr=float(sixsCoeffs[2, 3]),
                            DifIrr=float(sixsCoeffs[2, 4]),
                            EnvIrr=float(sixsCoeffs[2, 5]),
                        )
                    )
                    imgBandCoeffs.append(
                        rsgislib.imagecalibration.Band6SCoeff(
                            band=4,
                            aX=float(sixsCoeffs[3, 0]),
                            bX=float(sixsCoeffs[3, 1]),
                            cX=float(sixsCoeffs[3, 2]),
                            DirIrr=float(sixsCoeffs[3, 3]),
                            DifIrr=float(sixsCoeffs[3, 4]),
                            EnvIrr=float(sixsCoeffs[3, 5]),
                        )
                    )
                    imgBandCoeffs.append(
                        rsgislib.imagecalibration.Band6SCoeff(
                            band=5,
                            aX=float(sixsCoeffs[4, 0]),
                            bX=float(sixsCoeffs[4, 1]),
                            cX=float(sixsCoeffs[4, 2]),
                            DirIrr=float(sixsCoeffs[4, 3]),
                            DifIrr=float(sixsCoeffs[4, 4]),
                            EnvIrr=float(sixsCoeffs[4, 5]),
                        )
                    )
                    imgBandCoeffs.append(
                        rsgislib.imagecalibration.Band6SCoeff(
                            band=6,
                            aX=float(sixsCoeffs[5, 0]),
                            bX=float(sixsCoeffs[5, 1]),
                            cX=float(sixsCoeffs[5, 2]),
                            DirIrr=float(sixsCoeffs[5, 3]),
                            DifIrr=float(sixsCoeffs[5, 4]),
                            EnvIrr=float(sixsCoeffs[5, 5]),
                        )
                    )
                    imgBandCoeffs.append(
                        rsgislib.imagecalibration.Band6SCoeff(
                            band=7,
                            aX=float(sixsCoeffs[6, 0]),
                            bX=float(sixsCoeffs[6, 1]),
                            cX=float(sixsCoeffs[6, 2]),
                            DirIrr=float(sixsCoeffs[6, 3]),
                            DifIrr=float(sixsCoeffs[6, 4]),
                            EnvIrr=float(sixsCoeffs[6, 5]),
                        )
                    )
                    aot6SCoeffsOut.append(
                        rsgislib.imagecalibration.AOTLUTFeat(
                            AOT=float(aotVal), Coeffs=imgBandCoeffs
                        )
                    )
                elevAOTCoeffs.append(
                    rsgislib.imagecalibration.ElevLUTFeat(
                        Elev=float(elevVal), Coeffs=aot6SCoeffsOut
                    )
                )

        rsgislib.imagecalibration.apply_6s_coeff_elev_aot_lut_param(
            inputRadImage,
            inputDEMFile,
            inputAOTImage,
            outputImage,
            outFormat,
            rsgislib.TYPE_16UINT,
            scaleFactor,
            0,
            True,
            elevAOTCoeffs,
        )

        return outputImage, elevAOTCoeffs

    def run6SToOptimiseAODValue(
        self,
        aotVal,
        radBlueVal,
        predBlueVal,
        aeroProfile,
        atmosProfile,
        grdRefl,
        surfaceAltitude,
    ):
        """Used as part of the optimastion for identifying values of AOD"""
        print(
            "Testing AOD Val: ",
            aotVal,
        )
        s = Py6S.SixS()
        s.atmos_profile = atmosProfile
        s.aero_profile = aeroProfile
        s.ground_reflectance = grdRefl
        s.geometry = Py6S.Geometry.Landsat_TM()
        s.geometry.month = self.acquisitionTime.month
        s.geometry.day = self.acquisitionTime.day
        s.geometry.gmt_decimal_hour = (
            float(self.acquisitionTime.hour) + float(self.acquisitionTime.minute) / 60.0
        )
        s.geometry.latitude = self.latCentre
        s.geometry.longitude = self.lonCentre
        s.altitudes = Py6S.Altitudes()
        s.altitudes.set_target_custom_altitude(surfaceAltitude)
        s.altitudes.set_sensor_satellite_level()
        s.atmos_corr = Py6S.AtmosCorr.AtmosCorrLambertianFromRadiance(200)
        s.aot550 = aotVal

        # Band 2 (Blue!)
        s.wavelength = Py6S.Wavelength(
            0.436,
            0.5285,
            [
                0.000010,
                0.000117,
                0.000455,
                0.001197,
                0.006869,
                0.027170,
                0.271370,
                0.723971,
                0.903034,
                0.909880,
                0.889667,
                0.877453,
                0.879688,
                0.891913,
                0.848533,
                0.828339,
                0.868497,
                0.912538,
                0.931726,
                0.954248,
                0.956424,
                0.978564,
                0.989469,
                0.968801,
                0.988729,
                0.967361,
                0.966125,
                0.981834,
                0.963135,
                0.996498,
                0.844893,
                0.190738,
                0.005328,
                0.001557,
                0.000516,
                0.000162,
                0.000023,
                -0.000016,
            ],
        )
        s.run()
        aX = float(s.outputs.values["coef_xa"])
        bX = float(s.outputs.values["coef_xb"])
        cX = float(s.outputs.values["coef_xc"])

        tmpVal = (aX * radBlueVal) - bX
        reflBlueVal = tmpVal / (1.0 + cX * tmpVal)
        outDist = math.sqrt(math.pow((reflBlueVal - predBlueVal), 2))
        print("\taX: ", aX, " bX: ", bX, " cX: ", cX, "     Dist = ", outDist)
        return outDist

    def findDDVTargets(self, inputTOAImage, outputPath, outputName, outFormat, tmpPath):
        raise ARCSIException("Not Implemented")

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
        raise ARCSIException("Not Implemented")

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
        raise ARCSIException("Not Implemented")

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
        try:
            return self.estimateSingleAOTFromDOSBandImpl(
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
                2,
            )
        except Exception as e:
            raise

    def setBandNames(self, imageFile):
        dataset = gdal.Open(imageFile, gdal.GA_Update)
        dataset.GetRasterBand(1).SetDescription("Coastal")
        dataset.GetRasterBand(2).SetDescription("Blue")
        dataset.GetRasterBand(3).SetDescription("Green")
        dataset.GetRasterBand(4).SetDescription("Red")
        dataset.GetRasterBand(5).SetDescription("NIR")
        dataset.GetRasterBand(6).SetDescription("SWIR1")
        dataset.GetRasterBand(7).SetDescription("SWIR2")
        dataset = None

    def cleanLocalFollowProcessing(self):
        print("")
