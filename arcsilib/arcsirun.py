"""
Module that contains the ARCSIRun class.
"""
############################################################################
#  arcsirun.py
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
# Purpose: Class which does the actural hard work in performing the
#          atmospheric correction using the input parameters. 
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 08/12/2016
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
# Import the system library
import sys
# Import the subprocess module
import subprocess
# Import the OS python module
import os
# Import the os.path module
import os.path
# Import the time module
import time
# Import the copy module
import copy
# Import the glob module
import glob
# Import ARCSI library
import arcsilib
# Import the ARCSI exception class
from arcsilib.arcsiexception import ARCSIException
# Import the ARCSI sensor factory class
from arcsilib.arcsiutils import ARCSISensorFactory
# Import the ARCSI utilities class
from arcsilib.arcsiutils import ARCSIUtils
# Import the sensor classes
from arcsilib.arcsisensor import ARCSIAbstractSensor
# Import the image utilities module from rsgislib
import rsgislib.imageutils
# Import the image calculations module from rsgislib
import rsgislib.imagecalc
# Import the raster GIS module from rsgislib
import rsgislib.rastergis
# Import the base rsgislib module
import rsgislib
# Import the py6s module for running 6S from python.
import Py6S
# Import the osgeo gdal library
import osgeo.gdal as gdal
# Import python math library
import math
# Import the arcsi version number
from arcsilib import ARCSI_VERSION
# Import the arcsi copyright year
from arcsilib import ARCSI_COPYRIGHT_YEAR
# Import the arcsi support email
from arcsilib import ARCSI_SUPPORT_EMAIL
# Import the arcsi website
from arcsilib import ARCSI_WEBSITE
# Import the arcsi copyright names
from arcsilib import ARCSI_COPYRIGHT_NAMES
# Import the list of sensors arcsi supports
from arcsilib import ARCSI_SENSORS_LIST
# Import the list of products arcsi supports
from arcsilib import ARCSI_PRODUCTS_LIST

class ARCSIRun (object):
    """
    The \'main\' class which executes the whole ARCSI package.
    """

    def convertVisabilityToAOD(self, vis):
        return (3.9449/vis)+0.08498

    def findMinimumElev(self, elev):
        elevVal = -500
        outElev = 0
        for i in range(50):
            if (elev > elevVal) & (elev < (elevVal+100)):
                outElev = elevVal
                break
            elevVal = elevVal + 100
        return outElev

    def findMaximumElev(self, elev):
        elevVal = -500
        outElev = 0
        for i in range(90):
            if (elev > elevVal) & (elev < (elevVal+100)):
                outElev = elevVal + 100
                break
            elevVal = elevVal + 100
        return outElev


    def findMinimumAOT(self, aot):
        aotVal = 0
        outAOT = 0
        for i in range(200):
            if (aot > aotVal) & (aot < (aotVal+0.05)):
                outAOT = aotVal
                break
            aotVal = aotVal + 0.05
        return aotVal

    def findMaximumAOT(self, aot):
        aotVal = 0
        outAOT = 0
        for i in range(200):
            if (aot > aotVal) & (aot < (aotVal+0.05)):
                outAOT = aotVal+ 0.05
                break
            aotVal = aotVal + 0.05
        return aotVal


    def runARCSI(self, inputHeader, inputImage, cloudMaskUsrImg, sensorStr, inWKTFile, outFormat, outFilePath, outBaseName,
            outWKTFile, outProj4File, projAbbv, xPxlResUsr, yPxlResUsr, productsStr, calcStatsPy, aeroProfileOption,
            atmosProfileOption, aeroProfileOptionImg, atmosProfileOptionImg,  grdReflOption, surfaceAltitude,
            atmosOZoneVal,atmosWaterVal, atmosOZoneWaterSpecified, aeroWaterVal, aeroDustVal, aeroOceanicVal,
            aeroSootVal, aeroComponentsSpecified, aotVal, visVal, tmpPath, minAOT, maxAOT, lowAOT, upAOT,
            demFile, aotFile, globalDOS, dosOutRefl, simpleDOS, debugMode, scaleFactor, interpAlgor,
            initClearSkyRegionDist, initClearSkyRegionMinSize, finalClearSkyRegionDist, clearSkyMorphSize, fullImgOuts,
            checkOutputs, classmlclouds, cloudtrainclouds, cloudtrainother):
        """
        A function contains the main flow of the software
        """

        
        arcsiUtils = ARCSIUtils()
        rsgisUtils = rsgislib.RSGISPyUtils()
        # Create list to store products to be calculated and those actually calculated.
        prodsToCalc = dict()
        prodsCalculated = dict()
        try:
            # Read WKT file if provided.
            wktStr = None
            if inWKTFile != None:
                wktStr = arcsiUtils.readTextFile(inWKTFile)

            # Step 1: Get the Sensor specific class from factory
            sensorFact = ARCSISensorFactory()
            sensorClass = sensorFact.getSensorClassFromName(sensorStr, debugMode, inputImage)

            # Step 2: Read header parameters
            sensorClass.extractHeaderParameters(inputHeader, wktStr)
            print("")

            if not sensorClass.expectedImageDataPresent():
                raise Exception("Not all the expected input images are present as listed in the header file.")
            else:
                print("Input imagery as listed in header file is present.\n")

            # Step 3: If aerosol and atmosphere images are specified then sample them to find
            #         the aerosol and atmosphere generic model to use for conversion to SREF
            if not aeroProfileOptionImg == None:
                print("Get aero profile from image...")
                aeroProfileMode = int(rsgislib.imagecalc.getImageBandModeInEnv(aeroProfileOptionImg, 1, 1, None, sensorClass.latTL, sensorClass.latBR, sensorClass.lonBR, sensorClass.lonTL)[0])

                if aeroProfileMode == 1:
                    aeroProfileOption = "Maritime"
                elif aeroProfileMode == 2:
                    aeroProfileOption = "Continental"
                else:
                    raise Exception("The aerosol profile from the input image was not recognised.")
                print("Aerosol Profile = ", aeroProfileOption)
                print("")
            if not atmosProfileOptionImg == None:
                print("Get atmos profile from image...")
                atmosProfileMode = int(rsgislib.imagecalc.getImageBandModeInEnv(atmosProfileOptionImg, 1, 1, None, sensorClass.latTL, sensorClass.latBR, sensorClass.lonBR, sensorClass.lonTL)[0])
                summerWinter = arcsiUtils.isSummerOrWinter(sensorClass.latCentre, sensorClass.lonCentre, sensorClass.acquisitionTime )
                if atmosProfileMode == 1:
                    atmosProfileOption = "Tropical"
                elif atmosProfileMode == 2:
                    if summerWinter == 1:
                        atmosProfileOption = "MidlatitudeSummer"
                    elif summerWinter == 2:
                        atmosProfileOption = "MidlatitudeWinter"
                    else:
                        raise Exception("Not recognised as being summer or winter.")
                elif atmosProfileMode == 3:
                    if summerWinter == 1:
                        atmosProfileOption = "SubarcticSummer"
                    elif summerWinter == 2:
                        atmosProfileOption = "SubarcticWinter"
                    else:
                        raise Exception("Not recognised as being summer or winter.")
                else:
                    raise Exception("The atmosphere profile from the input image was not recognised.")
                print("Atmosphere Profile = ", atmosProfileOption)
                print("")

            # Step 3: Get Output Image Base Name.
            if (outBaseName is None) or (outBaseName is ""):
                outBaseName = sensorClass.generateOutputBaseName()
                if not projAbbv is None:
                    outBaseNameProj = outBaseName + "_" + str(projAbbv)
            print("Image Base Name: " + outBaseName + "\n")

            # Check whether output files for this input already exist - if checkOutputs is True.
            if checkOutputs:
                prevOutFiles = glob.glob(os.path.join(outFilePath, outBaseName+'*'))
                if len(prevOutFiles) > 0:
                    sys.stderr.write('Error outputs already exist: \'' + inputHeader + '\'\n')
                    for tmpFile in prevOutFiles:
                        sys.stderr.write('\tFile: \'' + tmpFile + '\'\n')
                    raise Exception("Output files already exist and the \'check outputs\' option (--checkouts) was specified so can't continue.")

            # Step 4: Find the products which are to be generated.
            prodsToCalc["RAD"] = False
            prodsToCalc["TOA"] = False
            prodsToCalc["CLOUDS"] = False
            prodsToCalc["CLEARSKY"] = False
            prodsToCalc["DDVAOT"] = False
            prodsToCalc["DOSAOT"] = False
            prodsToCalc["DOSAOTSGL"] = False
            prodsToCalc["SREF"] = False
            prodsToCalc["DOS"] = False
            prodsToCalc["THERMAL"] = False
            prodsToCalc["SATURATE"] = False
            prodsToCalc["TOPOSHADOW"] = False
            prodsToCalc["FOOTPRINT"] = False
            prodsToCalc["METADATA"] = False
            prodsToCalc["STDSREF"] = False

            # Make a copy of the dictionary to store calculated products.
            prodsCalculated = copy.copy(prodsToCalc)
            needAtmModel = False
            for prod in productsStr:
                if prod == 'RAD':
                    prodsToCalc["RAD"] = True
                elif prod == 'SATURATE':
                    prodsToCalc["SATURATE"] = True
                elif prod == 'TOA':
                    prodsToCalc["RAD"] = True
                    prodsToCalc["TOA"] = True
                elif prod == 'CLOUDS':
                    prodsToCalc["RAD"] = True
                    prodsToCalc["TOA"] = True
                    prodsToCalc["CLOUDS"] = True
                    prodsToCalc["SATURATE"] = True
                    if sensorClass.hasThermal():
                        prodsToCalc["THERMAL"] = True
                elif prod == 'CLEARSKY':
                    prodsToCalc["RAD"] = True
                    prodsToCalc["TOA"] = True
                    prodsToCalc["CLOUDS"] = True
                    prodsToCalc["SATURATE"] = True
                    if sensorClass.hasThermal():
                        prodsToCalc["THERMAL"] = True
                    prodsToCalc["CLEARSKY"] = True
                elif prod == 'DDVAOT':
                    prodsToCalc["RAD"] = True
                    prodsToCalc["TOA"] = True
                    prodsToCalc["DDVAOT"] = True
                    needAtmModel = True
                elif prod == 'DOSAOTSGL':
                    prodsToCalc["RAD"] = True
                    prodsToCalc["TOA"] = True
                    prodsToCalc["DOSAOTSGL"] = True
                    needAtmModel = True
                elif prod == 'SREF':
                    prodsToCalc["RAD"] = True
                    prodsToCalc["SREF"] = True
                    needAtmModel = True
                elif prod == 'STDSREF':
                    prodsToCalc["RAD"] = True
                    prodsToCalc["SREF"] = True
                    prodsToCalc["STDSREF"] = True
                    prodsToCalc["TOPOSHADOW"] = True
                    needAtmModel = True
                    if (demFile is None) or (demFile is ""):
                        raise ARCSIException("STDSREF requires a DEM file to be provided.")
                elif prod == 'DOS':
                    prodsToCalc["RAD"] = True
                    prodsToCalc["TOA"] = True
                    prodsToCalc["DOS"] = True
                elif prod == 'DOSAOT':
                    prodsToCalc["RAD"] = True
                    prodsToCalc["TOA"] = True
                    prodsToCalc["DOSAOT"] = True
                    needAtmModel = True
                elif prod == 'THERMAL':
                    if sensorClass.hasThermal():
                        prodsToCalc["THERMAL"] = True
                    else:
                        raise ARCSIException("The sensor does not have thermal bands. Check you inputs.")
                elif prod == 'TOPOSHADOW':
                    prodsToCalc["RAD"] = True
                    prodsToCalc["TOPOSHADOW"] = True
                    if (demFile is None) or (demFile is ""):
                        raise ARCSIException("STDSREF requires a DEM file to be provided.")
                elif prod == 'FOOTPRINT':
                    prodsToCalc["FOOTPRINT"] = True
                elif prod == 'METADATA':
                    prodsToCalc["METADATA"] = True
                
            if prodsToCalc["DOSAOT"] and prodsToCalc["DDVAOT"]:
                raise ARCSIException("You cannot specify both the DOSAOT and DDVAOT products, you must choose one or the other.")

            aeroProfile = None
            atmosProfile = None
            grdRefl = None
            useBRDF = False

            if needAtmModel:
                if aeroProfileOption == None:
                    raise ARCSIException("An aersol profile has not been specified.")
                elif aeroProfileOption == "NoAerosols":
                    aeroProfile = Py6S.AeroProfile.PredefinedType(Py6S.AeroProfile.NoAerosols)
                elif aeroProfileOption == "Continental":
                    aeroProfile = Py6S.AeroProfile.PredefinedType(Py6S.AeroProfile.Continental)
                elif aeroProfileOption == "Maritime":
                    aeroProfile = Py6S.AeroProfile.PredefinedType(Py6S.AeroProfile.Maritime)
                elif aeroProfileOption == "Urban":
                    aeroProfile = Py6S.AeroProfile.PredefinedType(Py6S.AeroProfile.Urban)
                elif aeroProfileOption == "Desert":
                    aeroProfile = Py6S.AeroProfile.PredefinedType(Py6S.AeroProfile.Desert)
                elif aeroProfileOption == "BiomassBurning":
                    aeroProfile = Py6S.AeroProfile.PredefinedType(Py6S.AeroProfile.BiomassBurning)
                elif aeroProfileOption == "Stratospheric":
                    aeroProfile = Py6S.AeroProfile.PredefinedType(Py6S.AeroProfile.Stratospheric)
                else:
                    raise ARCSIException("The specified aersol profile is unknown.")

                if atmosProfileOption == None:
                    raise ARCSIException("An atmospheric profile has not been specified.")
                elif atmosProfileOption == "NoGaseousAbsorption":
                    atmosProfile = Py6S.AtmosProfile.PredefinedType(Py6S.AtmosProfile.NoGaseousAbsorption)
                elif atmosProfileOption == "Tropical":
                    atmosProfile = Py6S.AtmosProfile.PredefinedType(Py6S.AtmosProfile.Tropical)
                elif atmosProfileOption == "MidlatitudeSummer":
                    atmosProfile = Py6S.AtmosProfile.PredefinedType(Py6S.AtmosProfile.MidlatitudeSummer)
                elif atmosProfileOption == "MidlatitudeWinter":
                    atmosProfile = Py6S.AtmosProfile.PredefinedType(Py6S.AtmosProfile.MidlatitudeWinter)
                elif atmosProfileOption == "SubarcticSummer":
                    atmosProfile = Py6S.AtmosProfile.PredefinedType(Py6S.AtmosProfile.SubarcticSummer)
                elif atmosProfileOption == "SubarcticWinter":
                    atmosProfile = Py6S.AtmosProfile.PredefinedType(Py6S.AtmosProfile.SubarcticWinter)
                elif atmosProfileOption == "USStandard1962":
                    atmosProfile = Py6S.AtmosProfile.PredefinedType(Py6S.AtmosProfile.USStandard1962)
                else:
                    raise ARCSIException("The specified atmospheric profile is unknown.")

                if atmosOZoneWaterSpecified:
                    atmosProfile = Py6S.AtmosProfile.UserWaterAndOzone(atmosWaterVal, atmosOZoneVal)

                if grdReflOption == None:
                    raise ARCSIException("A ground reflectance has not been specified.")
                elif grdReflOption == "GreenVegetation":
                      grdRefl = Py6S.GroundReflectance.HomogeneousLambertian(Py6S.GroundReflectance.GreenVegetation)
                elif grdReflOption == "ClearWater":
                    grdRefl = Py6S.GroundReflectance.HomogeneousLambertian(Py6S.GroundReflectance.ClearWater)
                elif grdReflOption == "Sand":
                    grdRefl = Py6S.GroundReflectance.HomogeneousLambertian(Py6S.GroundReflectance.Sand)
                elif grdReflOption == "LakeWater":
                    grdRefl = Py6S.GroundReflectance.HomogeneousLambertian(Py6S.GroundReflectance.LakeWater)
                elif grdReflOption == "BRDFHapke":
                    grdRefl = Py6S.GroundReflectance.HomogeneousHapke(0.101, -0.263, 0.589, 0.046)
                    useBRDF = True
                else:
                    raise ARCSIException("The specified ground reflectance is unknown.")

            # Decide is the image needs to be reprojected, if so define parameters.
            reproject = False
            useWKT2Reproject = True
            xPxlRes = 0.0
            yPxlRes = 0.0
            pxlResDefd = False
            reProjStr = ''
            if not outWKTFile is None:
                reproject = True
                useWKT2Reproject = True
                reProjStr = outWKTFile
            elif not outProj4File is None:
                reproject = True
                useWKT2Reproject = False
                reProjStr = '"' + arcsiUtils.readTextFile(outProj4File) + '"'

            if (xPxlResUsr is None) or (yPxlResUsr is None):
                pxlResDefd = False
            else:
                pxlResDefd = True
                xPxlRes = xPxlResUsr
                yPxlRes = yPxlResUsr


            validMaskImage=None
            validMaskImageProj=""
            viewAngleImg=""
            viewAngleImgProj=""
            radianceImage=""
            radianceImageWhole=""
            saturateImage=""
            saturateImageProj=""
            thermalRadImage=""
            thermalBrightImage=""
            thermalBrightImageWhole=""
            maskImage=""
            toaImage=""
            toaImageWhole=""
            srefImage=""
            srefDOSWholeImage=""
            sref6SWholeImage=""
            aotFile=""
            cloudsImage=""
            clearskyImage=""
            outDEMName=""
            topoShadowImage=""
            footprintShpFile=""
            metaDataFile=""
            propOfCloud = 0.0
            propOfClearSky = 0.0
            projImgBBOX = dict()
            projImgBBOX['MinX'] = 0.0
            projImgBBOX['MaxX'] = 0.0
            projImgBBOX['MinY'] = 0.0
            projImgBBOX['MaxY'] = 0.0
            processStageStr = ""
            processStageWholeImgStr = ""
            processSREFStr = ""
            finalOutFiles = dict()
            calcdOutVals = dict()
            sixsLUTCoeffs = None
            aotLUT = False

            if reproject:
                if useWKT2Reproject:
                    calcdOutVals["REPROJECT"] = arcsiUtils.readTextFile(outWKTFile)
                else:
                    calcdOutVals["REPROJECT"] = arcsiUtils.readTextFile(outProj4File)

            # Check Input image(s) is valid before proceeding.
            print('Checking Input Images are valid')
            sensorClass.checkInputImageValid()

            if sensorClass.imgNeedMosaicking():
                print("Mosacking Input Image Tiles.")
                sensorClass.mosaicImageTiles(outFilePath)

            # Get the valid image data maskImage
            outName = outBaseName + "_valid" + arcsiUtils.getFileExtension(outFormat)
            viewAngleImg = os.path.join(outFilePath, outBaseName + "_viewangle" + arcsiUtils.getFileExtension(outFormat))
            validMaskImage = sensorClass.generateValidImageDataMask(outFilePath, outName, viewAngleImg, outFormat)
            if not os.path.exists(viewAngleImg):
                viewAngleImg = None
            if not validMaskImage is None:
                rsgislib.rastergis.populateStats(validMaskImage, True, True)
            if not viewAngleImg is None:
                rsgislib.imageutils.popImageStats(viewAngleImg, usenodataval=True, nodataval=99999, calcpyramids=True)
            print("")
            if not validMaskImage is None:
                finalOutFiles["VALID_MASK"] = validMaskImage
            if not viewAngleImg is None:
                finalOutFiles["VIEW_ANGLE"] = viewAngleImg

            if reproject and (not validMaskImage is None):
                if not pxlResDefd:
                    validImgDS = gdal.Open(validMaskImage, gdal.GA_ReadOnly)
                    if validImgDS is None:
                        raise ARCSIException('Could not open the Valid Image Mask ' + validMaskImage)
                    geoTransform = validImgDS.GetGeoTransform()
                    if geoTransform is None:
                        raise ARCSIException('Could read the geotransform from the Valid Image Mask ' + validMaskImage)
                    xPxlRes = geoTransform[1]
                    yPxlRes = geoTransform[5]
                    validImgDS = None
                print("Define re-projected image BBOX.")
                projImgBBOX = sensorClass.getReProjBBOX(outWKTFile, outProj4File, useWKT2Reproject, xPxlRes, yPxlRes, True)
                print("Output Image BBOX (TL, BR): [({0:10.2f}, {1:10.2f}), ({2:10.2f}, {3:10.2f})]".format(projImgBBOX['MinX'], projImgBBOX['MaxY'], projImgBBOX['MaxX'], projImgBBOX['MinY']))
                outName = outBaseNameProj + "_valid" + arcsiUtils.getFileExtension(outFormat)
                validMaskImageProj = os.path.join(outFilePath, outName)

                cmd = 'gdalwarp -t_srs ' + reProjStr + ' -tr ' + str(xPxlRes) + ' ' + str(yPxlRes) + ' -ot Byte -wt Float32 ' \
                    + '-te ' + str(projImgBBOX['MinX']) + ' ' + str(projImgBBOX['MinY']) + ' ' + str(projImgBBOX['MaxX']) + ' ' + str(projImgBBOX['MaxY']) \
                    + ' -r near -tap -srcnodata 0 -dstnodata 0 -of ' + outFormat + ' -overwrite ' \
                    + validMaskImage + ' ' + validMaskImageProj
                print(cmd)
                try:
                    subprocess.call(cmd, shell=True)
                except OSError as e:
                    raise ARCSIException('Could not re-projection valid image mask: ' + cmd)
                if not os.path.exists(validMaskImageProj):
                    raise ARCSIException('Reprojected valid image mask is not present: ' + validMaskImageProj)
                else:
                    rsgislib.rastergis.populateStats(validMaskImageProj, True, True)
                finalOutFiles["VALID_MASK"] = validMaskImageProj
                print("")

            if reproject and ((not viewAngleImg is "") and (not viewAngleImg is None)):
                if not pxlResDefd:
                    viewAngleImgDS = gdal.Open(viewAngleImg, gdal.GA_ReadOnly)
                    if viewAngleImgDS is None:
                        raise ARCSIException('Could not open the Valid Image Mask ' + viewAngleImg)
                    geoTransform = viewAngleImgDS.GetGeoTransform()
                    if geoTransform is None:
                        raise ARCSIException('Could read the geotransform from the Valid Image Mask ' + viewAngleImg)
                    xPxlRes = geoTransform[1]
                    yPxlRes = geoTransform[5]
                    viewAngleImgDS = None
                print("Define re-projected image BBOX.")
                projImgBBOX = sensorClass.getReProjBBOX(outWKTFile, outProj4File, useWKT2Reproject, xPxlRes, yPxlRes, True)
                print("Output Image BBOX (TL, BR): [({0:10.2f}, {1:10.2f}), ({2:10.2f}, {3:10.2f})]".format(projImgBBOX['MinX'], projImgBBOX['MaxY'], projImgBBOX['MaxX'], projImgBBOX['MinY']))
                outName = outBaseNameProj + "_viewangle" + arcsiUtils.getFileExtension(outFormat)
                viewAngleImgProj = os.path.join(outFilePath, outName)

                cmd = 'gdalwarp -t_srs ' + reProjStr + ' -tr ' + str(xPxlRes) + ' ' + str(yPxlRes) + ' -ot Float32 -wt Float32 ' \
                    + '-te ' + str(projImgBBOX['MinX']) + ' ' + str(projImgBBOX['MinY']) + ' ' + str(projImgBBOX['MaxX']) + ' ' + str(projImgBBOX['MaxY']) \
                    + ' -r cubicspline -tap -srcnodata 99999 -dstnodata 99999 -of ' + outFormat + ' -overwrite ' \
                    + viewAngleImg + ' ' + viewAngleImgProj
                print(cmd)
                try:
                    subprocess.call(cmd, shell=True)
                except OSError as e:
                    raise ARCSIException('Could not re-projection valid image mask: ' + cmd)
                if not os.path.exists(viewAngleImgProj):
                    raise ARCSIException('Reprojected valid image mask is not present: ' + viewAngleImgProj)
                else:
                    rsgislib.imageutils.popImageStats(viewAngleImg, usenodataval=True, nodataval=99999, calcpyramids=True)
                rsgisUtils.deleteFileWithBasename(viewAngleImg)
                viewAngleImg = viewAngleImgProj
                finalOutFiles["VIEW_ANGLE"] = viewAngleImgProj
                print("")

            if prodsToCalc["FOOTPRINT"]:
                if validMaskImage is None:
                    raise ARCSIException("To generate a footprint a valid image mask is required - not supported by this sensor?")
                outFootprintLyrName = ''
                if reproject:
                    outFootprintLyrName = outBaseNameProj + "_footprint"
                    footprintShpFile = sensorClass.generateImageFootprint(validMaskImageProj, outFilePath, outFootprintLyrName)
                else:
                    outFootprintLyrName = outBaseName + "_footprint"
                    footprintShpFile = sensorClass.generateImageFootprint(validMaskImage, outFilePath, outFootprintLyrName)
                finalOutFiles["FOOTPRINT"] = validMaskImageProj
                prodsCalculated["FOOTPRINT"] = True
                print("")

            if prodsToCalc["SATURATE"]:
                # Execute generation of the saturation image
                outName = outBaseName + "_sat" + arcsiUtils.getFileExtension(outFormat)
                saturateImage = sensorClass.generateImageSaturationMask(outFilePath, outName, outFormat)

                if reproject:
                    outNameProj = outBaseNameProj + "_sat" + arcsiUtils.getFileExtension(outFormat)
                    saturateImageProj = os.path.join(outFilePath, outNameProj)
                    cmd = 'gdalwarp -t_srs ' + reProjStr + ' -tr ' + str(xPxlRes) + ' ' + str(yPxlRes) + ' -ot Byte ' \
                    + '-te ' + str(projImgBBOX['MinX']) + ' ' + str(projImgBBOX['MinY']) + ' ' + str(projImgBBOX['MaxX']) + ' ' + str(projImgBBOX['MaxY']) \
                    + ' -r near -tap -of ' + outFormat + ' -overwrite ' \
                    + saturateImage + ' ' + saturateImageProj
                    print(cmd)
                    try:
                        subprocess.call(cmd, shell=True)
                    except OSError as e:
                        raise ARCSIException('Could not re-projection saturated image mask: ' + cmd)
                    if not os.path.exists(saturateImageProj):
                        raise ARCSIException('Reprojected saturated image mask is not present: ' + saturateImageProj)
                    else:
                        rsgisUtils.deleteFileWithBasename(saturateImage)
                        saturateImage = saturateImageProj
                if calcStatsPy:
                    print("Calculating Statistics...")
                    rsgislib.imageutils.popImageStats(saturateImage, usenodataval=False, nodataval=0, calcpyramids=True)
                prodsCalculated["SATURATE"] = True
                print("")

            # Step 5: Convert to Radiance
            if prodsToCalc["RAD"]:
                # Execute conversion to radiance
                outName = outBaseName + "_rad" + arcsiUtils.getFileExtension(outFormat)
                outThermName = None
                if prodsToCalc["THERMAL"]:
                    outThermName = outBaseName + "_therm_rad" + arcsiUtils.getFileExtension(outFormat)
                radianceImage, thermalRadImage = sensorClass.convertImageToRadiance(outFilePath, outName, outThermName, outFormat)

                if sensorClass.maskInputImages():
                    processStageStr = processStageStr + "_msk"
                    outImgName = outBaseName + processStageStr + "_rad" + arcsiUtils.getFileExtension(outFormat)
                    outMaskName = outBaseName + "_mask" + arcsiUtils.getFileExtension(outFormat)
                    radianceImageTmp, maskImage = sensorClass.applyImageDataMask(inputHeader, radianceImage, outFilePath, outMaskName, outImgName, outFormat, None)
                    if not radianceImageTmp is radianceImage:
                        rsgisUtils.deleteFileWithBasename(radianceImage)
                        radianceImage = radianceImageTmp
                    if not maskImage == None:
                        print("Setting Band Names...")
                        sensorClass.setBandNames(radianceImage)
                        if calcStatsPy:
                            print("Calculating Statistics...")
                            rsgislib.imageutils.popImageStats(radianceImage, True, 0.0, True)
                            rsgislib.rastergis.populateStats(maskImage, True, True)

                if not validMaskImage is None:
                    print("Masking to valid data area.")
                    processStageStr = processStageStr + "_vmsk"
                    outRadPathName = os.path.join(outFilePath, outBaseName + processStageStr + "_rad" + arcsiUtils.getFileExtension(outFormat))
                    rsgislib.imageutils.maskImage(radianceImage, validMaskImage, outRadPathName, outFormat, rsgisUtils.getRSGISLibDataTypeFromImg(radianceImage), 0.0, 0.0)
                    rsgisUtils.deleteFileWithBasename(radianceImage)
                    radianceImage = outRadPathName
                    if not thermalRadImage == None:
                        outThermPathName = os.path.join(outFilePath, outBaseName + processStageStr + "_thermrad" + arcsiUtils.getFileExtension(outFormat))
                        rsgislib.imageutils.maskImage(thermalRadImage, validMaskImage, outThermPathName, outFormat, rsgisUtils.getRSGISLibDataTypeFromImg(thermalRadImage), 0.0, 0.0)
                        rsgisUtils.deleteFileWithBasename(thermalRadImage)
                        thermalRadImage = outThermPathName
                    if reproject:
                        rsgisUtils.deleteFileWithBasename(validMaskImage)
                        validMaskImage = validMaskImageProj

                if reproject:
                    if not radianceImage is None:
                        outName = outBaseNameProj + processStageStr + "_rad" + arcsiUtils.getFileExtension(outFormat)

                        outRadImagePath = os.path.join(outFilePath, outName)
                        cmd = 'gdalwarp -t_srs ' + reProjStr + ' -tr ' + str(xPxlRes) + ' ' + str(yPxlRes) + ' -ot Float32 -wt Float32 ' \
                        + '-te ' + str(projImgBBOX['MinX']) + ' ' + str(projImgBBOX['MinY']) + ' ' + str(projImgBBOX['MaxX']) + ' ' + str(projImgBBOX['MaxY']) \
                        + ' -r ' + interpAlgor + ' -tap -srcnodata 0 -dstnodata 0 -of ' + outFormat + ' -overwrite ' \
                        + radianceImage + ' ' + outRadImagePath
                        print(cmd)
                        try:
                            subprocess.call(cmd, shell=True)
                        except OSError as e:
                            raise ARCSIException('Could not re-projection radiance image: ' + cmd)
                        if not os.path.exists(outRadImagePath):
                            raise ARCSIException('Reprojected radiance image is not present: ' + outRadImagePath)
                        else:
                            rsgisUtils.deleteFileWithBasename(radianceImage)
                            radianceImage = outRadImagePath
                    if not thermalRadImage is None:
                        outName = outBaseNameProj + processStageStr + "_thrad" + arcsiUtils.getFileExtension(outFormat)
                        outThermRadImagePath = os.path.join(outFilePath, outName)
                        cmd = 'gdalwarp -t_srs ' + reProjStr + ' -tr ' + str(xPxlRes) + ' ' + str(yPxlRes) + ' -ot Float32 -wt Float32 ' \
                        + '-te ' + str(projImgBBOX['MinX']) + ' ' + str(projImgBBOX['MinY']) + ' ' + str(projImgBBOX['MaxX']) + ' ' + str(projImgBBOX['MaxY']) \
                        + ' -r ' + interpAlgor + ' -tap -srcnodata 0 -dstnodata 0 -of ' + outFormat + ' -overwrite ' \
                        + thermalRadImage + ' ' + outThermRadImagePath
                        try:
                            subprocess.call(cmd, shell=True)
                        except OSError as e:
                            raise ARCSIException('Could not re-projection thermal radiance image: ' + cmd)
                        if not os.path.exists(outThermRadImagePath):
                            raise ARCSIException('Reprojected thermal radiance image is not present: ' + outThermRadImagePath)
                        else:
                            rsgisUtils.deleteFileWithBasename(thermalRadImage)
                            thermalRadImage = outThermRadImagePath
                    outBaseName = outBaseNameProj

                print("Setting Band Names...")
                sensorClass.setBandNames(radianceImage)
                if calcStatsPy:
                    print("Calculating Statistics...")
                    rsgislib.imageutils.popImageStats(radianceImage, True, 0.0, True)
                    if not thermalRadImage == None:
                        rsgislib.imageutils.popImageStats(thermalRadImage, True, 0.0, True)
                
                finalOutFiles["RADIANCE_WHOLE"] = radianceImage
                finalOutFiles["RADIANCE"] = radianceImage
                if not thermalRadImage == None:
                    finalOutFiles["THERM_RADIANCE_WHOLE"] = thermalRadImage

                radianceImageWhole = radianceImage
                prodsCalculated["RAD"] = True
                print("")

            # Execute calibrate thermal to brightness
            if prodsToCalc["THERMAL"]:
                outName = outBaseName + processStageStr + "_thrad_thermbright" + arcsiUtils.getFileExtension(outFormat)
                thermalBrightImage = sensorClass.convertThermalToBrightness(thermalRadImage, outFilePath, outName, outFormat, scaleFactor)
                if calcStatsPy:
                    print("Calculating Statistics...")
                    rsgislib.imageutils.popImageStats(thermalBrightImage, True, 0.0, True)
                finalOutFiles["THERMAL_BRIGHT_WHOLE"] = thermalBrightImage
                finalOutFiles["THERMAL_BRIGHT"] = thermalBrightImage
                thermalBrightImageWhole = thermalBrightImage
                prodsCalculated["THERMAL"] = True
                print("")

            # Execute conversion to top of atmosphere reflectance
            if prodsToCalc["TOA"]:
                outName = outBaseName + processStageStr +"_rad_toa" + arcsiUtils.getFileExtension(outFormat)
                toaImage = sensorClass.convertImageToTOARefl(radianceImage, outFilePath, outName, outFormat, scaleFactor)
                print("Setting Band Names...")
                sensorClass.setBandNames(toaImage)
                if calcStatsPy:
                    print("Calculating Statistics...")
                    rsgislib.imageutils.popImageStats(toaImage, True, 0.0, True)
                finalOutFiles["TOA_WHOLE"] = toaImage
                finalOutFiles["TOA"] = toaImage
                toaImageWhole = toaImage
                prodsCalculated["TOA"] = True
                print("")

            # Save the process stage string for using with whole image outputs.
            processStageWholeImgStr = processStageStr

            # Step 7: Generate Cloud Masks
            if prodsToCalc["CLOUDS"]:
                outName = outBaseName + "_clouds" + arcsiUtils.getFileExtension(outFormat)
                if cloudMaskUsrImg == None:
                    if classmlclouds:
                        cloudsImage = sensorClass.generateCloudMaskML(toaImage, validMaskImage, outFilePath, outName, outFormat, tmpPath, cloudtrainclouds, cloudtrainother, numCores=1)
                    else:
                        cloudsImage = sensorClass.generateCloudMask(toaImage, saturateImage, thermalBrightImage, validMaskImage, outFilePath, outName, outFormat, tmpPath, scaleFactor)
                    if calcStatsPy:
                        print("Calculating Statistics...")
                        rsgislib.rastergis.populateStats(cloudsImage, False, True)
                else:
                    cloudsImage = cloudMaskUsrImg
                finalOutFiles["CLOUD_MASK"] = cloudsImage

                # Calculate the proportion of the scene cover by cloud.
                propOfCloud = rsgislib.imagecalc.calcPropTrueExp('b1==1?1:b1==2?1:0', [rsgislib.imagecalc.BandDefn('b1', cloudsImage, 1)], validMaskImage)
                print("The scene is " + str(propOfCloud*100) + "% cloud.")
                calcdOutVals['ARCSI_CLOUD_COVER'] = propOfCloud

                if propOfCloud < 0.98: # Less than 98% cloud cover then process.
                    print("Applying cloud masks to images...")
                    processStageStr = processStageStr + "_mclds"
                    outputRADImage = os.path.join(outFilePath, outBaseName + processStageStr + "_rad" + arcsiUtils.getFileExtension(outFormat))
                    rsgislib.imageutils.maskImage(radianceImage, cloudsImage, outputRADImage, outFormat, rsgislib.TYPE_32FLOAT, 0, [1,2])
                    radianceImage = outputRADImage
                    sensorClass.setBandNames(radianceImage)
                    finalOutFiles["RADIANCE"] = radianceImage
                    if calcStatsPy:
                        print("Calculating Statistics...")
                        rsgislib.imageutils.popImageStats(radianceImage, True, 0.0, True)
                    outputTOAImage = os.path.join(outFilePath, outBaseName + processStageStr + "_rad_toa" + arcsiUtils.getFileExtension(outFormat))
                    rsgislib.imageutils.maskImage(toaImage, cloudsImage, outputTOAImage, outFormat, rsgislib.TYPE_16UINT, 0, [1,2])
                    toaImage = outputTOAImage
                    sensorClass.setBandNames(toaImage)
                    finalOutFiles["TOA"] = toaImage
                    if calcStatsPy:
                        print("Calculating Statistics...")
                        rsgislib.imageutils.popImageStats(toaImage, True, 0.0, True)
                prodsCalculated["CLOUDS"] = True
                print("")

            # Don't continue further if there is more than 95% cloud cover in the scene.
            if  (not prodsToCalc["CLOUDS"]) or (prodsToCalc["CLOUDS"] and propOfCloud < 0.95):
                if prodsToCalc["CLEARSKY"]:
                    outName = outBaseName + "_clearsky" + arcsiUtils.getFileExtension(outFormat)
                    clearskyImage = sensorClass.generateClearSkyMask(cloudsImage, validMaskImage, outFilePath, outName, outFormat, tmpPath, initClearSkyRegionDist, initClearSkyRegionMinSize, finalClearSkyRegionDist, clearSkyMorphSize)
                    if calcStatsPy:
                        print("Calculating Statistics...")
                        rsgislib.rastergis.populateStats(clearskyImage, True, True)
                    finalOutFiles["CLEARSKY_MASK"] = clearskyImage

                    # Calculate the proportion of the scene which is clear sky.
                    propOfClearSky = rsgislib.imagecalc.calcPropTrueExp('b1==1?1:0', [rsgislib.imagecalc.BandDefn('b1', clearskyImage, 1)], validMaskImage)
                    print("The scene is " + str(propOfClearSky*100) + "% clear-sky.")
                    calcdOutVals['ARCSI_CLEARSKY_COVER'] = propOfClearSky

                    if propOfClearSky > 0.05: # Keep going if at least 5% of the scene is clear sky
                        print("Applying clear-sky masks to images...")
                        processStageStr = processStageStr + "_clearsky"
                        outputRADImage = os.path.join(outFilePath, outBaseName + processStageStr + "_rad" + arcsiUtils.getFileExtension(outFormat))
                        rsgislib.imageutils.maskImage(radianceImage, clearskyImage, outputRADImage, outFormat, rsgislib.TYPE_32FLOAT, 0, 0)
                        radianceImage = outputRADImage
                        sensorClass.setBandNames(radianceImage)
                        finalOutFiles["RADIANCE"] = radianceImage
                        if calcStatsPy:
                            print("Calculating Statistics...")
                            rsgislib.imageutils.popImageStats(radianceImage, True, 0.0, True)
                        outputTOAImage = os.path.join(outFilePath, outBaseName + processStageStr + "_rad_toa" + arcsiUtils.getFileExtension(outFormat))
                        rsgislib.imageutils.maskImage(toaImage, clearskyImage, outputTOAImage, outFormat, rsgislib.TYPE_16UINT, 0, 0)
                        toaImage = outputTOAImage
                        sensorClass.setBandNames(toaImage)
                        finalOutFiles["TOA"] = toaImage
                        if calcStatsPy:
                            print("Calculating Statistics...")
                            rsgislib.imageutils.popImageStats(toaImage, True, 0.0, True)
                    prodsCalculated["CLEARSKY"] = True
                    print("")

                # Don't continue further if there is less than 5% clear sky in the scene.
                if  (not prodsToCalc["CLEARSKY"]) or (prodsToCalc["CLEARSKY"] and propOfClearSky > 0.05):
                    # Interpolate DEM image to output refl image resolution and convert projection to the same as the output image.
                    if (not (demFile == None)) and (outDEMName == ""):
                        outDEMNameTmp = os.path.join(outFilePath, (outBaseName + "_demtmp" + arcsiUtils.getFileExtension(outFormat)))
                        rsgislib.imageutils.createCopyImage(radianceImage, outDEMNameTmp, 1, -32768.0, outFormat, rsgislib.TYPE_32FLOAT)

                        inDEMDS = gdal.Open(demFile, gdal.GA_ReadOnly)
                        outDEMDS = gdal.Open(outDEMNameTmp, gdal.GA_Update)

                        print("Subset and reproject DEM...")
                        gdal.ReprojectImage(inDEMDS, outDEMDS, None, None, gdal.GRA_CubicSpline)
                        inDEMDS = None
                        outDEMDS = None

                        outDEMName = os.path.join(outFilePath, (outBaseName + "_dem" + arcsiUtils.getFileExtension(outFormat)))
                        print("Output DEM: ", outDEMName)
                        rsgislib.imageutils.maskImage(outDEMNameTmp, validMaskImage, outDEMName, outFormat, rsgislib.TYPE_32FLOAT, -32768.0, 0)

                        outDEMNameMsk = os.path.join(outFilePath, (outBaseName + "_demmsk" + arcsiUtils.getFileExtension(outFormat)))
                        if prodsToCalc["CLEARSKY"]:
                            rsgislib.imageutils.maskImage(outDEMName, clearskyImage, outDEMNameMsk, outFormat, rsgislib.TYPE_32FLOAT, -32768.0, 0)
                        elif prodsToCalc["CLOUDS"]:
                            rsgislib.imageutils.maskImage(outDEMName, cloudsImage, outDEMNameMsk, outFormat, rsgislib.TYPE_32FLOAT, -32768.0, [1,2,3])
                        else:
                            outDEMNameMsk = outDEMName

                        # Calculate DEM statistics and set no data value.
                        if prodsToCalc["CLEARSKY"] or prodsToCalc["CLOUDS"]:
                            rsgislib.imageutils.popImageStats(outDEMNameMsk, True, -32768.0, True)
                        rsgislib.imageutils.popImageStats(outDEMName, True, -32768.0, True)

                        # Remove tmp DEM file.
                        rsgisUtils.deleteFileWithBasename(outDEMNameTmp)
                        finalOutFiles["IMAGE_DEM"] = outDEMName

                    # Execute generation of the topographic shadow image
                    if prodsToCalc["TOPOSHADOW"]:
                        outName = outBaseName + "_toposhad" + arcsiUtils.getFileExtension(outFormat)
                        topoShadowImage = sensorClass.generateTopoDirectShadowMask(outDEMNameMsk, outFilePath, outName, outFormat, tmpPath)
                        if calcStatsPy:
                            print("Calculating Statistics...")
                            rsgislib.rastergis.populateStats(topoShadowImage, True, True)
                        finalOutFiles["TOPO_SHADOW_MASK"] = topoShadowImage

                        processStageStr = processStageStr + "_topshad"
                        if prodsToCalc["RAD"]:
                            outputRADImage = os.path.join(outFilePath, outBaseName + processStageStr + "_rad"  + arcsiUtils.getFileExtension(outFormat))
                            rsgislib.imageutils.maskImage(radianceImage, topoShadowImage, outputRADImage, outFormat, rsgislib.TYPE_32FLOAT, 0, 1)
                            radianceImage = outputRADImage
                            sensorClass.setBandNames(radianceImage)
                            if calcStatsPy:
                                print("Calculating Statistics...")
                                rsgislib.imageutils.popImageStats(radianceImage, True, 0.0, True)
                        if prodsToCalc["TOA"]:
                            outputTOAImage = os.path.join(outFilePath, outBaseName + processStageStr + "_rad_toa" + arcsiUtils.getFileExtension(outFormat))
                            rsgislib.imageutils.maskImage(toaImage, topoShadowImage, outputTOAImage, outFormat, rsgislib.TYPE_16UINT, 0, 1)
                            toaImage = outputTOAImage
                            sensorClass.setBandNames(toaImage)
                            if calcStatsPy:
                                print("Calculating Statistics...")
                                rsgislib.imageutils.popImageStats(toaImage, True, 0.0, True)

                        prodsCalculated["TOPOSHADOW"] = True
                        print("")

                    # Step 8: Convert to an approximation of Surface Reflectance using a dark object subtraction
                    if prodsToCalc["DOS"]:
                        print("Convert to reflectance using dark object subtraction.")
                        outName = outBaseName + processStageStr + "_rad_toa_dos" + arcsiUtils.getFileExtension(outFormat)
                        outWholeName = outBaseName + processStageWholeImgStr + "_rad_toa_dos" + arcsiUtils.getFileExtension(outFormat)
                        srefImage, offVals = sensorClass.convertImageToReflectanceSimpleDarkSubtract(toaImage, outFilePath, outName, outFormat, dosOutRefl)
                        calcdOutVals['ARCSI_DOS_OFFSETS'] = offVals
                        if fullImgOuts:
                            srefDOSWholeImage, offVals = sensorClass.convertImageToReflectanceSimpleDarkSubtract(toaImageWhole, outFilePath, outWholeName, outFormat, dosOutRefl, offVals)

                        print("Setting Band Names...")
                        sensorClass.setBandNames(srefImage)
                        finalOutFiles["SREF_DOS_IMG"] = srefImage
                        if fullImgOuts:
                            sensorClass.setBandNames(srefDOSWholeImage)
                            finalOutFiles["SREF_DOS_IMG_WHOLE"] = srefDOSWholeImage

                        if calcStatsPy:
                            print("Calculating Statistics...")
                            rsgislib.imageutils.popImageStats(srefImage, True, 0.0, True)
                            if fullImgOuts:
                                rsgislib.imageutils.popImageStats(srefDOSWholeImage, True, 0.0, True)
                            print("")
                        prodsCalculated["DOS"] = True

                    # Step 9: Use image to estimate AOD values
                    if prodsToCalc["DOSAOTSGL"]:
                        calcdOutVals['ARCSI_AOT_RANGE_MIN'] = minAOT
                        calcdOutVals['ARCSI_AOT_RANGE_MAX'] = maxAOT
                        aotVal = sensorClass.estimateSingleAOTFromDOS(radianceImage, toaImage, outDEMName, tmpPath, outBaseName, outFormat, aeroProfile, atmosProfile, grdRefl, minAOT, maxAOT, dosOutRefl)
                        minAOT = aotVal - lowAOT
                        if minAOT < 0.01:
                            minAOT = 0.05
                        maxAOT = aotVal + upAOT
                        print("AOT Search Range = [" + str(minAOT) + ", " + str(maxAOT) + "]")
                        calcdOutVals['ARCSI_AOT_VALUE'] = aotVal
                        prodsCalculated["DOSAOTSGL"] = True

                    if prodsToCalc["DDVAOT"]:
                        outName = outBaseName + "_ddvaod" + arcsiUtils.getFileExtension(outFormat)
                        aotFile = sensorClass.estimateImageToAODUsingDDV(radianceImage, toaImage, outDEMName, topoShadowImage, outFilePath, outName, outFormat, tmpPath, aeroProfile, atmosProfile, grdRefl, minAOT, maxAOT)
                        dataset = gdal.Open(aotFile, gdal.GA_Update)
                        dataset.GetRasterBand(1).SetDescription("AOT")
                        dataset = None
                        if calcStatsPy:
                            print("Calculating Statistics...")
                            rsgislib.imageutils.popImageStats(aotFile, True, 0.0, True)
                        finalOutFiles["AOTIMG_DDV"] = aotFile
                        calcdOutVals['ARCSI_AOT_RANGE_MIN'] = minAOT
                        calcdOutVals['ARCSI_AOT_RANGE_MAX'] = maxAOT
                        prodsCalculated["DDVAOT"] = True
                        print("")

                    if prodsToCalc["DOSAOT"]:
                        outName = outBaseName + "_dosaod" + arcsiUtils.getFileExtension(outFormat)
                        aotFile = sensorClass.estimateImageToAODUsingDOS(radianceImage, toaImage, outDEMName, topoShadowImage, outFilePath, outName, outFormat, tmpPath, aeroProfile, atmosProfile, grdRefl, minAOT, maxAOT, globalDOS, simpleDOS, dosOutRefl)
                        dataset = gdal.Open(aotFile, gdal.GA_Update)
                        dataset.GetRasterBand(1).SetDescription("AOT")
                        dataset = None
                        if calcStatsPy:
                            print("Calculating Statistics...")
                            rsgislib.imageutils.popImageStats(aotFile, True, 0.0, True)
                        finalOutFiles["AOTIMG_DOS"] = aotFile
                        calcdOutVals['ARCSI_AOT_RANGE_MIN'] = minAOT
                        calcdOutVals['ARCSI_AOT_RANGE_MAX'] = maxAOT
                        prodsCalculated["DOSAOT"] = True
                        print("")

                    # Step 10: Convert to Surface Reflectance using 6S Standard Models
                    if prodsToCalc["SREF"]:
                        if (prodsToCalc["DDVAOT"] or prodsToCalc["DOSAOT"]) and (not aotFile == ""):
                            imgDS = gdal.Open(aotFile, gdal.GA_ReadOnly )
                            imgBand = imgDS.GetRasterBand(1)
                            (min,max,mean,stddev) = imgBand.ComputeStatistics(False)
                            print("AOT Mean (Std Dev) = " + str(mean) + " (" + str(stddev) + ")")
                            print("AOT [Min, Max] = [" + str(min) + "," + str(max) + "]")
                            aotVal = mean
                            if aotVal < 0.01:
                                print("WARNING: Something has gone wrong as AOT value is 0 or below. Setting to 0.05")
                                aotVal = 0.05
                            imgDS = None

                        if (aotVal == None) and (visVal == None) and (aotFile == ""):
                            raise ARCSIException("Either the AOT or the visability need to specified.")
                        elif (aotVal == None) and (aotFile == ""):
                            print("Convert to vis to aot...")
                            aotVal = self.convertVisabilityToAOD(visVal)

                        if not (aotVal == None):
                            print("AOT Value: {}".format(aotVal))
                            calcdOutVals['ARCSI_AOT_VALUE'] = aotVal

                        if (demFile == None):
                            processSREFStr = '_rad_sref'
                            outName = outBaseName + processStageStr + processSREFStr + arcsiUtils.getFileExtension(outFormat)
                            srefImage = sensorClass.convertImageToSurfaceReflSglParam(radianceImage, outFilePath, outName, outFormat, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF, scaleFactor)
                            if fullImgOuts:
                                outName = outBaseName + processStageWholeImgStr + processSREFStr + arcsiUtils.getFileExtension(outFormat)
                                sref6SWholeImage = sensorClass.convertImageToSurfaceReflSglParam(radianceImageWhole, outFilePath, outName, outFormat, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF, scaleFactor)
                            calcdOutVals['ARCSI_ELEVATION_VALUE'] = surfaceAltitude
                        else:
                            # Calc Min, Max Elevation for region intersecting with the image.
                            statsElev = rsgislib.imagecalc.getImageStatsInEnv(outDEMName, 1, -32768.0, sensorClass.latTL, sensorClass.latBR, sensorClass.lonBR, sensorClass.lonTL)

                            print("Minimum Elevation = ", statsElev[0])
                            print("Maximum Elevation = ", statsElev[1])

                            minElev = self.findMinimumElev(statsElev[0])
                            maxElev = self.findMaximumElev(statsElev[1])

                            calcdOutVals['ARCSI_LUT_ELEVATION_MIN'] = minElev
                            calcdOutVals['ARCSI_LUT_ELEVATION_MAX'] = maxElev

                            elevRange = (maxElev - minElev) / 100
                            numElevSteps = math.ceil(elevRange) + 1
                            print("Elevation Ranges from ", minElev, " to ", maxElev, " an LUT with ", numElevSteps, " will be created.")

                            if (aotFile == None) or (aotFile == ""):
                                print("Build an DEM LUT with AOT == " + str(aotVal) + "...")
                                processSREFStr = '_rad_srefdem'
                                outName = outBaseName + processStageStr + processSREFStr + arcsiUtils.getFileExtension(outFormat)
                                srefImage, sixsLUTCoeffs = sensorClass.convertImageToSurfaceReflDEMElevLUT(radianceImage, outDEMName, outFilePath, outName, outFormat, aeroProfile, atmosProfile, grdRefl, aotVal, useBRDF, minElev, maxElev, scaleFactor)
                                if fullImgOuts:
                                    outName = outBaseName + processStageWholeImgStr + processSREFStr + arcsiUtils.getFileExtension(outFormat)
                                    sref6SWholeImage, sixsLUTCoeffs = sensorClass.convertImageToSurfaceReflDEMElevLUT(radianceImageWhole, outDEMName, outFilePath, outName, outFormat, aeroProfile, atmosProfile, grdRefl, aotVal, useBRDF, minElev, maxElev, scaleFactor, sixsLUTCoeffs)
                                calcdOutVals['ARCSI_6S_COEFFICENTS'] = sixsLUTCoeffs
                                aotLUT = False
                            else:
                                print("Build an AOT and DEM LUT...")
                                statsAOT = rsgislib.imagecalc.getImageStatsInEnv(aotFile, 1, -9999, sensorClass.latTL, sensorClass.latBR, sensorClass.lonBR, sensorClass.lonTL)

                                minAOT = self.findMinimumAOT(statsAOT[0])
                                if minAOT < 0.01:
                                    minAOT = 0.05
                                maxAOT = self.findMaximumAOT(statsAOT[1])

                                calcdOutVals['ARCSI_LUT_AOT_MIN'] = minAOT
                                calcdOutVals['ARCSI_LUT_AOT_MAX'] = maxAOT

                                aotRange = (maxAOT - minAOT) / 0.05
                                numAOTSteps = math.ceil(aotRange) + 1
                                print("AOT Ranges from ", minAOT, " to ", maxAOT, " an LUT with ", numAOTSteps, " will be created.")
                                processSREFStr = '_rad_srefdemaot'
                                outName = outBaseName + processStageStr + processSREFStr + arcsiUtils.getFileExtension(outFormat)
                                srefImage, sixsLUTCoeffs = sensorClass.convertImageToSurfaceReflAOTDEMElevLUT(radianceImage, outDEMName, aotFile, outFilePath, outName, outFormat, aeroProfile, atmosProfile, grdRefl, useBRDF, minElev, maxElev, minAOT, maxAOT, scaleFactor)
                                if fullImgOuts:
                                    outName = outBaseName + processStageWholeImgStr + processSREFStr + arcsiUtils.getFileExtension(outFormat)
                                    sref6SWholeImage, sixsLUTCoeffs = sensorClass.convertImageToSurfaceReflAOTDEMElevLUT(radianceImageWhole, outDEMName, aotFile, outFilePath, outName, outFormat, aeroProfile, atmosProfile, grdRefl, useBRDF, minElev, maxElev, minAOT, maxAOT, scaleFactor, sixsLUTCoeffs)
                                calcdOutVals['ARCSI_6S_COEFFICENTS'] = sixsLUTCoeffs
                                aotLUT = True

                        print("Setting Band Names...")
                        sensorClass.setBandNames(srefImage)
                        if fullImgOuts:
                            sensorClass.setBandNames(sref6SWholeImage)
                            finalOutFiles["SREF_DOS_IMG_WHOLE"] = sref6SWholeImage

                        if calcStatsPy:
                            print("Calculating Statistics...")
                            rsgislib.imageutils.popImageStats(srefImage, True, 0.0, True)
                            if fullImgOuts:
                                rsgislib.imageutils.popImageStats(sref6SWholeImage, True, 0.0, True)
                        finalOutFiles["SREF_6S_IMG"] = srefImage
                        if fullImgOuts:
                            finalOutFiles["SREF_6S_WHOLE_IMG"] = sref6SWholeImage
                        prodsCalculated["SREF"] = True
                        print("")

                    if prodsToCalc["STDSREF"]:
                        if sixsLUTCoeffs is None:
                            raise ARCSIException("To calculate standardised reflectance \'STDSREF\' you need to have a LUT within 6S.\nTherefore, a DEM is required and optional an AOT surface should be calculated.")
                        if not fullImgOuts:
                            sref6SWholeImage = None
                        processSREFStr = processSREFStr + '_stdsref'
                        outName = outBaseName + processStageStr + processSREFStr + arcsiUtils.getFileExtension(outFormat)
                        outNameWhole = outBaseName + processStageWholeImgStr + processSREFStr + arcsiUtils.getFileExtension(outFormat)
                        stdSREFImg, stdSREFWholeImg = sensorClass.convertSREF2StdisedSREF(srefImage, sref6SWholeImage, outDEMName, topoShadowImage, outFilePath, outName, outNameWhole, outFormat, tmpPath, sixsLUTCoeffs, aotLUT, scaleFactor, brdfBeta=1, outIncidenceAngle=0, outExitanceAngle=0)

                        print("Setting Band Names...")
                        sensorClass.setBandNames(stdSREFImg)
                        if fullImgOuts:
                            sensorClass.setBandNames(stdSREFWholeImg)

                        if calcStatsPy:
                            print("Calculating Statistics...")
                            rsgislib.imageutils.popImageStats(stdSREFImg, True, 0.0, True)
                            if fullImgOuts:
                                rsgislib.imageutils.popImageStats(stdSREFWholeImg, True, 0.0, True)
                        finalOutFiles["STD_SREF_IMG"] = stdSREFImg
                        if fullImgOuts:
                            finalOutFiles["STD_SREF_WHOLE_IMG"] = stdSREFWholeImg
                        prodsCalculated["STDSREF"] = True
                        print("")
                        
                else:
                    keys2Del = []
                    for key in prodsToCalc.keys():
                        if prodsToCalc[key] is not prodsCalculated[key]:
                            if key in ['DDVAOT', 'DOSAOT', 'DOSAOTSGL', 'SREF', 'DOS', 'STDSREF', 'TOPOSHADOW']:
                                keys2Del.append(key)
                    for key in keys2Del:
                        del prodsToCalc[key]
            else:
                keys2Del = []
                for key in prodsToCalc.keys():
                    if prodsToCalc[key] is not prodsCalculated[key]:
                        if key in ['CLEARSKY', 'DDVAOT', 'DOSAOT', 'DOSAOTSGL', 'SREF', 'DOS', 'STDSREF', 'TOPOSHADOW']:
                            keys2Del.append(key)
                for key in keys2Del:
                    del prodsToCalc[key]

            if prodsToCalc["METADATA"]:
                print("Exporting Meta-data file")
                outName = outBaseName + "_meta.json"
                finalOutFiles["METADATA"] = outName

                validMaskImagePath = ""
                if not validMaskImage is None:
                    if reproject:
                        validMaskImagePath = validMaskImageProj
                    else:
                        validMaskImagePath = validMaskImage

                sensorClass.generateMetaDataFile(outFilePath, outName, productsStr, validMaskImagePath, prodsToCalc["FOOTPRINT"], calcdOutVals, finalOutFiles)
                prodsCalculated["METADATA"] = True
                print("")

            print('Clean up anything left over...')
            sensorClass.cleanFollowProcessing()

        except ARCSIException as e:
            print('Input Header: \'' + inputHeader + '\'', file=sys.stderr)
            if outBaseName is not None:
                print('Output Basename: \'' + outBaseName + '\'', file=sys.stderr)
            print("Error: {}".format(e), file=sys.stderr)
        except Exception as e:
            print('Input Header: \'' + inputHeader + '\'', file=sys.stderr)
            if outBaseName is not None:
                print('Output Basename: \'' + outBaseName + '\'', file=sys.stderr)
            print("Error: {}".format(e), file=sys.stderr)
        finally:
            failedProdsList = []
            # Check all requested products have been created
            for key in prodsToCalc.keys():
                if prodsToCalc[key] is not prodsCalculated[key]:
                    failedProdsList.append(key)
            if len(failedProdsList) > 0:
                print("Error: The following products were not generated:", file=sys.stderr)
                print(" ".join(failedProdsList), file=sys.stderr)
                print('Input Header: \'' + inputHeader + '\'', file=sys.stderr)
                if outBaseName is not None:
                    print('Output Basename: \'' + outBaseName + '\'', file=sys.stderr)
                print('\n\n', file=sys.stderr) # Put gap into output log to easier to see where one ends and the next starts.


    def print2ConsoleListSensors(self):
        """
        A function which lists the currently supported sensors
        and the names by which they should be specified to the
        ARCSI command line argument.
        """
        print("Supported Sensors are:")
        print("\t----------------------------------------------------------------------------------------------------------------------------------")
        print("\tSensor        | Shorthand     | Functions")
        print("\t----------------------------------------------------------------------------------------------------------------------------------")
        print("\tLandsat 1 MSS | \'ls1\'       | RAD, TOA, DOSAOT, DOSAOTSGL, SREF, STDSREF, DOS, TOPOSHADOW, FOOTPRINT, METADATA")
        print("\tLandsat 2 MSS | \'ls2\'       | RAD, TOA, DOSAOT, DOSAOTSGL, SREF, STDSREF, DOS, TOPOSHADOW, FOOTPRINT, METADATA")
        print("\tLandsat 3 MSS | \'ls3\'       | RAD, TOA, DOSAOT, DOSAOTSGL, SREF, STDSREF, DOS, TOPOSHADOW, FOOTPRINT, METADATA")
        print("\tLandsat 4 MSS | \'ls4mss\'    | RAD, TOA, DOSAOT, DOSAOTSGL, SREF, STDSREF, DOS, TOPOSHADOW, FOOTPRINT, METADATA")
        print("\tLandsat 4 TM  | \'ls4tm\'     | RAD, TOA, DOSAOT, DDVAOT, DOSAOTSGL, STDSREF, SREF, DOS, THERMAL, TOPOSHADOW, FOOTPRINT, METADATA")
        print("\tLandsat 5 MSS | \'ls5mss\'    | RAD, TOA, DOSAOT, DOSAOTSGL, SREF, STDSREF, DOS, TOPOSHADOW, FOOTPRINT, METADATA")
        print("\tLandsat 5 TM  | \'ls5tm\'     | RAD, TOA, DOSAOT, DDVAOT, DOSAOTSGL, STDSREF, SREF, DOS, THERMAL, TOPOSHADOW, FOOTPRINT, METADATA")
        print("\tLandsat 7 ETM | \'ls7\'       | RAD, TOA, DOSAOT, DDVAOT, DOSAOTSGL, STDSREF, SREF, DOS, THERMAL, TOPOSHADOW, FOOTPRINT, METADATA")
        print("\tLandsat 8     | \'ls8\'       | RAD, TOA, DOSAOT, DDVAOT, DOSAOTSGL, STDSREF, SREF, DOS, THERMAL, TOPOSHADOW, FOOTPRINT, METADATA")
        print("\tRapideye      | \'rapideye\'  | RAD, TOA, DOSAOT, DOSAOTSGL, SREF, STDSREF, DOS, TOPOSHADOW, METADATA")
        print("\tWorldView2    | \'wv2\'       | RAD, TOA, DOSAOT, DOSAOTSGL, SREF, STDSREF, DOS, TOPOSHADOW, METADATA")
        print("\tSPOT5         | \'spot5\'     | RAD, TOA, DOSAOT, DOSAOTSGL, SREF, STDSREF, DOS, TOPOSHADOW, METADATA")
        print("\t----------------------------------------------------------------------------------------------------------------------------------")

    def print2ConsoleListProductDescription(self, product):
        """
        A function which lists the currently supported products
        and describes what that are and the parameters they require.
        """
        print("Hello World. this has not be written yet!")

    def print2ConsoleListEnvVars(self):
        """
        A function which lists the available environmental variables for ARCSI.
        """
        print("ARCSI_OUT_FORMAT       in place of the -f, --format option")
        print("ARCSI_OUTPUT_PATH      in place of the -o, --outpath option")
        print("ARCSI_TMP_PATH         in place of the --tmpath option")
        print("ARCSI_DEM_PATH         in place of the -d, --dem option")
        print("ARCSI_AEROIMG_PATH     in place of the --aeroimg option")
        print("ARCSI_ATMOSIMG_PATH    in place of the --atmosimg option")
        print("ARCSI_MIN_AOT          in place of the --minaot option")
        print("ARCSI_MAX_AOT          in place of the --maxaot option")
        print("ARCSI_LOW_AOT          in place of the --lowaot option")
        print("ARCSI_UP_AOT           in place of the --upaot option")
        print("ARCSI_OUTDOS_REFL      in place of the --dosout option")
        print("                       Note reflectance values are multiplied")
        print("                       by 1000 so a value of 20 is 2 %")
        print("ARCSI_USE_LOCALDOS     in place of the --localdos (variable ")
        print("                       values can be either `TRUE' or `FALSE') option")
        print("ARCSI_USE_SIMPLEDOS    in place of the --simpledos (variable ")
        print("                       values can be either `TRUE' or `FALSE') option")
        print("")
