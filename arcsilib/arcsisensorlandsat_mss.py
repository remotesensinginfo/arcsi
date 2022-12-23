"""
Module that contains the ARCSILandsatMSSSensor class.
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
#           the pre-processing operations within ARCSI to the landsat 5 MSS
#           datasets.
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 13/07/2013
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

import numpy
import Py6S
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


class ARCSILandsatMSSSensor(ARCSIAbstractSensor):
    """
    A class which represents the landsat 5 MSS sensor to read
    header parameters and apply data processing operations.
    """

    def __init__(self, debugMode, inputImage):
        ARCSIAbstractSensor.__init__(self, debugMode, inputImage)
        self.sensor = "LS_MSS"
        self.collection_num = 0

        self.band1File = ""
        self.band2File = ""
        self.band3File = ""
        self.band4File = ""
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

        self.b1MinRad = 0.0
        self.b1MaxRad = 0.0
        self.b2MinRad = 0.0
        self.b2MaxRad = 0.0
        self.b3MinRad = 0.0
        self.b3MaxRad = 0.0
        self.b4MinRad = 0.0
        self.b4MaxRad = 0.0

        self.sensorID = ""
        self.spacecraftID = ""
        self.cloudCover = 0.0
        self.cloudCoverLand = 0.0
        self.earthSunDistance = 0.0
        self.gridCellSizeRefl = 0.0

    def extractHeaderParameters(self, inputHeader, wktStr):
        """
        Understands and parses the Landsat MTL header files
        """
        try:
            if not self.userSpInputImage is None:
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
                (headerParams["SPACECRAFT_ID"].upper() == "LANDSAT_5")
                or (headerParams["SPACECRAFT_ID"].upper() == "LANDSAT5")
            ) and (headerParams["SENSOR_ID"].upper() == "MSS"):
                self.sensor = "LS5MSS"
            elif (
                (headerParams["SPACECRAFT_ID"].upper() == "LANDSAT_4")
                or (headerParams["SPACECRAFT_ID"].upper() == "LANDSAT4")
            ) and (headerParams["SENSOR_ID"].upper() == "MSS"):
                self.sensor = "LS4MSS"
            elif (
                (headerParams["SPACECRAFT_ID"].upper() == "LANDSAT_3")
                or (headerParams["SPACECRAFT_ID"].upper() == "LANDSAT3")
            ) and (headerParams["SENSOR_ID"].upper() == "MSS"):
                self.sensor = "LS3MSS"
            elif (
                (headerParams["SPACECRAFT_ID"].upper() == "LANDSAT_2")
                or (headerParams["SPACECRAFT_ID"].upper() == "LANDSAT2")
            ) and (headerParams["SENSOR_ID"].upper() == "MSS"):
                self.sensor = "LS2MSS"
            elif (
                (headerParams["SPACECRAFT_ID"].upper() == "LANDSAT_1")
                or (headerParams["SPACECRAFT_ID"].upper() == "LANDSAT1")
            ) and (headerParams["SENSOR_ID"].upper() == "MSS"):
                self.sensor = "LS1MSS"
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
            try:
                self.row = int(headerParams["WRS_ROW"])
            except Exception as e:
                self.row = int(headerParams["STARTING_ROW"])
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

            filesDIR = os.path.dirname(inputHeader)

            if "FILE_NAME_BAND_7" in headerParams:
                self.band1File = os.path.join(
                    filesDIR, headerParams["FILE_NAME_BAND_4"]
                )
                self.band2File = os.path.join(
                    filesDIR, headerParams["FILE_NAME_BAND_5"]
                )
                self.band3File = os.path.join(
                    filesDIR, headerParams["FILE_NAME_BAND_6"]
                )
                self.band4File = os.path.join(
                    filesDIR, headerParams["FILE_NAME_BAND_7"]
                )

                self.b1CalMin = rsgislib.tools.utils.str_to_float(
                    headerParams["QUANTIZE_CAL_MIN_BAND_4"], 1.0
                )
                self.b1CalMax = rsgislib.tools.utils.str_to_float(
                    headerParams["QUANTIZE_CAL_MAX_BAND_4"], 255.0
                )
                self.b2CalMin = rsgislib.tools.utils.str_to_float(
                    headerParams["QUANTIZE_CAL_MIN_BAND_5"], 1.0
                )
                self.b2CalMax = rsgislib.tools.utils.str_to_float(
                    headerParams["QUANTIZE_CAL_MAX_BAND_5"], 255.0
                )
                self.b3CalMin = rsgislib.tools.utils.str_to_float(
                    headerParams["QUANTIZE_CAL_MIN_BAND_6"], 1.0
                )
                self.b3CalMax = rsgislib.tools.utils.str_to_float(
                    headerParams["QUANTIZE_CAL_MAX_BAND_6"], 255.0
                )
                self.b4CalMin = rsgislib.tools.utils.str_to_float(
                    headerParams["QUANTIZE_CAL_MIN_BAND_7"], 1.0
                )
                self.b4CalMax = rsgislib.tools.utils.str_to_float(
                    headerParams["QUANTIZE_CAL_MAX_BAND_7"], 255.0
                )

                self.b1MinRad = rsgislib.tools.utils.str_to_float(
                    headerParams["RADIANCE_MINIMUM_BAND_4"], 2.500
                )
                self.b1MaxRad = rsgislib.tools.utils.str_to_float(
                    headerParams["RADIANCE_MAXIMUM_BAND_4"], 220.800
                )
                self.b2MinRad = rsgislib.tools.utils.str_to_float(
                    headerParams["RADIANCE_MINIMUM_BAND_5"], 2.700
                )
                self.b2MaxRad = rsgislib.tools.utils.str_to_float(
                    headerParams["RADIANCE_MAXIMUM_BAND_5"], 163.600
                )
                self.b3MinRad = rsgislib.tools.utils.str_to_float(
                    headerParams["RADIANCE_MINIMUM_BAND_6"], 4.700
                )
                self.b3MaxRad = rsgislib.tools.utils.str_to_float(
                    headerParams["RADIANCE_MAXIMUM_BAND_6"], 140.300
                )
                self.b4MinRad = rsgislib.tools.utils.str_to_float(
                    headerParams["RADIANCE_MINIMUM_BAND_7"], 2.900
                )
                self.b4MaxRad = rsgislib.tools.utils.str_to_float(
                    headerParams["RADIANCE_MAXIMUM_BAND_7"], 117.500
                )
            else:
                self.band1File = os.path.join(
                    filesDIR, headerParams["FILE_NAME_BAND_1"]
                )
                self.band2File = os.path.join(
                    filesDIR, headerParams["FILE_NAME_BAND_2"]
                )
                self.band3File = os.path.join(
                    filesDIR, headerParams["FILE_NAME_BAND_3"]
                )
                self.band4File = os.path.join(
                    filesDIR, headerParams["FILE_NAME_BAND_4"]
                )

                self.b1CalMin = rsgislib.tools.utils.str_to_float(
                    headerParams["QUANTIZE_CAL_MIN_BAND_1"], 1.0
                )
                self.b1CalMax = rsgislib.tools.utils.str_to_float(
                    headerParams["QUANTIZE_CAL_MAX_BAND_1"], 255.0
                )
                self.b2CalMin = rsgislib.tools.utils.str_to_float(
                    headerParams["QUANTIZE_CAL_MIN_BAND_2"], 1.0
                )
                self.b2CalMax = rsgislib.tools.utils.str_to_float(
                    headerParams["QUANTIZE_CAL_MAX_BAND_2"], 255.0
                )
                self.b3CalMin = rsgislib.tools.utils.str_to_float(
                    headerParams["QUANTIZE_CAL_MIN_BAND_3"], 1.0
                )
                self.b3CalMax = rsgislib.tools.utils.str_to_float(
                    headerParams["QUANTIZE_CAL_MAX_BAND_3"], 255.0
                )
                self.b4CalMin = rsgislib.tools.utils.str_to_float(
                    headerParams["QUANTIZE_CAL_MIN_BAND_4"], 1.0
                )
                self.b4CalMax = rsgislib.tools.utils.str_to_float(
                    headerParams["QUANTIZE_CAL_MAX_BAND_4"], 255.0
                )

                self.b1MinRad = rsgislib.tools.utils.str_to_float(
                    headerParams["RADIANCE_MINIMUM_BAND_1"], 2.500
                )
                self.b1MaxRad = rsgislib.tools.utils.str_to_float(
                    headerParams["RADIANCE_MAXIMUM_BAND_1"], 220.800
                )
                self.b2MinRad = rsgislib.tools.utils.str_to_float(
                    headerParams["RADIANCE_MINIMUM_BAND_2"], 2.700
                )
                self.b2MaxRad = rsgislib.tools.utils.str_to_float(
                    headerParams["RADIANCE_MAXIMUM_BAND_2"], 163.600
                )
                self.b3MinRad = rsgislib.tools.utils.str_to_float(
                    headerParams["RADIANCE_MINIMUM_BAND_3"], 4.700
                )
                self.b3MaxRad = rsgislib.tools.utils.str_to_float(
                    headerParams["RADIANCE_MAXIMUM_BAND_3"], 140.300
                )
                self.b4MinRad = rsgislib.tools.utils.str_to_float(
                    headerParams["RADIANCE_MINIMUM_BAND_4"], 2.900
                )
                self.b4MaxRad = rsgislib.tools.utils.str_to_float(
                    headerParams["RADIANCE_MAXIMUM_BAND_4"], 117.500
                )

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
        imgInfo["CellSizeRefl"] = self.gridCellSizeRefl
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

        return imageDataPresent

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
            "Landsat 5 MSS does not provide any image masks, do not use the MASK option."
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
        inImages = [self.band1File, self.band2File, self.band3File, self.band4File]
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
            rsgislib.rastergis.clumps_spatial_extent(
                clumps_img=tmpValidPxlMsk,
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
            rsgislib.imagecalibration.calc_nadir_img_view_angle(
                tmpValidPxlMsk,
                viewAngleImg,
                "KEA",
                705000.0,
                "MinXX",
                "MinXY",
                "MaxXX",
                "MaxXY",
                "MinYX",
                "MinYY",
                "MaxYX",
                "MaxYY",
            )
        rsgislib.imagecalc.image_math(
            viewAngleImg, outputImage, "b1<7.65?1:0", outFormat, rsgislib.TYPE_8UINT
        )

        rsgislib.imageutils.delete_gdal_layer(tmpValidPxlMsk)
        return outputImage

    def convertImageToRadiance(
        self, outputPath, outputReflName, outputThermalName, outFormat
    ):
        print("Converting to Radiance")
        outputImage = os.path.join(outputPath, outputReflName)
        bandDefnSeq = list()
        lsBand = collections.namedtuple(
            "LSBand",
            [
                "band_name",
                "input_img",
                "img_band",
                "l_min",
                "l_max",
                "q_cal_min",
                "q_cal_max",
            ],
        )
        bandDefnSeq.append(
            lsBand(
                band_name="Green",
                input_img=self.band1File,
                img_band=1,
                l_min=self.b1MinRad,
                l_max=self.b1MaxRad,
                q_cal_min=self.b1CalMin,
                q_cal_max=self.b1CalMax,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="Red",
                input_img=self.band2File,
                img_band=1,
                l_min=self.b2MinRad,
                l_max=self.b2MaxRad,
                q_cal_min=self.b2CalMin,
                q_cal_max=self.b2CalMax,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="NIR1",
                input_img=self.band3File,
                img_band=1,
                l_min=self.b3MinRad,
                l_max=self.b3MaxRad,
                q_cal_min=self.b3CalMin,
                q_cal_max=self.b3CalMax,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="NIR2",
                input_img=self.band4File,
                img_band=1,
                l_min=self.b4MinRad,
                l_max=self.b4MaxRad,
                q_cal_min=self.b4CalMin,
                q_cal_max=self.b4CalMax,
            )
        )
        rsgislib.imagecalibration.landsat_to_radiance(
            outputImage, outFormat, bandDefnSeq
        )
        return outputImage, None

    def generateImageSaturationMask(self, outputPath, outputName, outFormat):
        print("Generate Saturation Image")
        outputImage = os.path.join(outputPath, outputName)

        lsBand = collections.namedtuple(
            "LSBand", ["band_name", "input_img", "img_band", "sat_val"]
        )
        bandDefnSeq = list()
        bandDefnSeq.append(
            lsBand(
                band_name="Green",
                input_img=self.band1File,
                img_band=1,
                sat_val=self.b1CalMax,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="Red",
                input_img=self.band2File,
                img_band=1,
                sat_val=self.b2CalMax,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="NIR1",
                input_img=self.band3File,
                img_band=1,
                sat_val=self.b3CalMax,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="NIR2",
                input_img=self.band4File,
                img_band=1,
                sat_val=self.b4CalMax,
            )
        )

        rsgislib.imagecalibration.saturated_pixels_mask(
            outputImage, outFormat, bandDefnSeq
        )

        return outputImage

    def convertThermalToBrightness(
        self, inputRadImage, outputPath, outputName, outFormat, scaleFactor
    ):
        raise ARCSIException("There are no thermal bands...")

    def convertImageToTOARefl(
        self, inputRadImage, outputPath, outputName, outFormat, scaleFactor
    ):
        print("Converting to TOA")
        outputImage = os.path.join(outputPath, outputName)
        solarIrradianceVals = list()
        IrrVal = collections.namedtuple("SolarIrradiance", ["irradiance"])
        solarIrradianceVals.append(IrrVal(irradiance=1824.0))
        solarIrradianceVals.append(IrrVal(irradiance=1570.0))
        solarIrradianceVals.append(IrrVal(irradiance=1249.0))
        solarIrradianceVals.append(IrrVal(irradiance=853.4))
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
        raise ARCSIException("Cloud Masking Not Implemented for LS MSS.")

    def createCloudMaskDataArray(self, inImgDataArr):
        # Calc Whiteness
        meanArr = numpy.mean(inImgDataArr, axis=1)
        whitenessArr = (
            numpy.absolute((inImgDataArr[..., 0] - meanArr) / meanArr)
            + numpy.absolute((inImgDataArr[..., 1] - meanArr) / meanArr)
            + numpy.absolute((inImgDataArr[..., 2] - meanArr) / meanArr)
            + numpy.absolute((inImgDataArr[..., 3] - meanArr) / meanArr)
        )
        # Calc NDVI
        ndvi = (inImgDataArr[..., 3] - inImgDataArr[..., 1]) / (
            inImgDataArr[..., 3] + inImgDataArr[..., 1]
        )

        # Create and populate the output array.
        inShape = inImgDataArr.shape
        outShape = [inShape[0], inShape[1] + 3]
        outArr = numpy.zeros(outShape, dtype=float)

        for i in range(inShape[1]):
            outArr[..., i] = inImgDataArr[..., i]

        idx = inShape[1]
        outArr[..., idx] = meanArr
        outArr[..., idx + 1] = whitenessArr
        outArr[..., idx + 2] = ndvi

        return outArr

    def defineDarkShadowImageBand(self):
        return 4

    def calc6SCoefficients(
        self, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF
    ):
        sixsCoeffs = numpy.zeros((4, 6), dtype=numpy.float32)
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
        s.wavelength = Py6S.Wavelength(Py6S.PredefinedWavelengths.LANDSAT_MSS_B1)
        s.run()
        sixsCoeffs[0, 0] = float(s.outputs.values["coef_xa"])
        sixsCoeffs[0, 1] = float(s.outputs.values["coef_xb"])
        sixsCoeffs[0, 2] = float(s.outputs.values["coef_xc"])
        sixsCoeffs[0, 3] = float(s.outputs.values["direct_solar_irradiance"])
        sixsCoeffs[0, 4] = float(s.outputs.values["diffuse_solar_irradiance"])
        sixsCoeffs[0, 5] = float(s.outputs.values["environmental_irradiance"])

        # Band 2
        s.wavelength = Py6S.Wavelength(Py6S.PredefinedWavelengths.LANDSAT_MSS_B2)
        s.run()
        sixsCoeffs[1, 0] = float(s.outputs.values["coef_xa"])
        sixsCoeffs[1, 1] = float(s.outputs.values["coef_xb"])
        sixsCoeffs[1, 2] = float(s.outputs.values["coef_xc"])
        sixsCoeffs[1, 3] = float(s.outputs.values["direct_solar_irradiance"])
        sixsCoeffs[1, 4] = float(s.outputs.values["diffuse_solar_irradiance"])
        sixsCoeffs[1, 5] = float(s.outputs.values["environmental_irradiance"])

        # Band 3
        s.wavelength = Py6S.Wavelength(Py6S.PredefinedWavelengths.LANDSAT_MSS_B3)
        s.run()
        sixsCoeffs[2, 0] = float(s.outputs.values["coef_xa"])
        sixsCoeffs[2, 1] = float(s.outputs.values["coef_xb"])
        sixsCoeffs[2, 2] = float(s.outputs.values["coef_xc"])
        sixsCoeffs[2, 3] = float(s.outputs.values["direct_solar_irradiance"])
        sixsCoeffs[2, 4] = float(s.outputs.values["diffuse_solar_irradiance"])
        sixsCoeffs[2, 5] = float(s.outputs.values["environmental_irradiance"])

        # Band 4
        s.wavelength = Py6S.Wavelength(Py6S.PredefinedWavelengths.LANDSAT_MSS_B4)
        s.run()
        sixsCoeffs[3, 0] = float(s.outputs.values["coef_xa"])
        sixsCoeffs[3, 1] = float(s.outputs.values["coef_xb"])
        sixsCoeffs[3, 2] = float(s.outputs.values["coef_xc"])
        sixsCoeffs[3, 3] = float(s.outputs.values["direct_solar_irradiance"])
        sixsCoeffs[3, 4] = float(s.outputs.values["diffuse_solar_irradiance"])
        sixsCoeffs[3, 5] = float(s.outputs.values["environmental_irradiance"])

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

        # Band 1 (Blue!)
        s.wavelength = Py6S.Wavelength(Py6S.PredefinedWavelengths.LANDSAT_MSS_B1)
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
        raise ARCSIException("findDDVTargets Not implemented\n")

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
        raise ARCSIException("estimateImageToAODUsingDDV Not implemented\n")

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
                1,
            )
        except Exception as e:
            raise

    def setBandNames(self, imageFile):
        """
        Set band names for Landsat 5 MSS

        Green: 500 - 600 nm
        Red: 600 - 700 nm
        NIR1: 700 - 800 nm
        NIR2: 800 - 1100 nm

        http://landsat.usgs.gov/about_landsat5.php

        """
        dataset = gdal.Open(imageFile, gdal.GA_Update)
        dataset.GetRasterBand(1).SetDescription("Green")
        dataset.GetRasterBand(2).SetDescription("Red")
        dataset.GetRasterBand(3).SetDescription("NIR1")
        dataset.GetRasterBand(4).SetDescription("NIR2")
        dataset = None

    def cleanLocalFollowProcessing(self):
        print("")
