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

# Import updated print function into python 2.7
from __future__ import print_function
# Import updated division operator into python 2.7
from __future__ import division
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
# Import the multiprocessing Pool module
from multiprocessing import Pool


class ARCSIParamsObj (object):

    def __init__(self):
        inputHeader = ""
        inputImage = None
        sensorClass = None
        sensorStr = ""
        productsStr = ""
        tmpPath = ""
        outFilePath = ""
        outFormat = 'KEA'
        outFormatExt = '.kea'
        calcStatsPy = True
        outBaseName = None
        cloudMaskUsrImg = None
        wktStr = None
        inWKTFile = None
        outWKTFile = None
        outProj4File = None
        projAbbv = None
        outBaseNameProj = None
        xPxlResUsr = None
        yPxlResUsr = None        
        aeroProfileOption = None
        atmosProfileOption = None
        aeroProfileOptionImg = None
        atmosProfileOptionImg = None
        grdReflOption = 'GreenVegetation'
        surfaceAltitude = 0.0
        atmosOZoneVal = None
        atmosWaterVal = None
        atmosOZoneWaterSpecified = None
        aeroWaterVal = None
        aeroDustVal = None
        aeroOceanicVal = None
        aeroSootVal = None
        aeroComponentsSpecified = False
        aotVal = None
        visVal = None
        minAOT = 0.05
        maxAOT = 0.5
        lowAOT = 0.1
        upAOT = 0.4
        demFile = None
        aotFile = None
        globalDOS = True
        dosOutRefl = 20
        simpleDOS = False
        debugMode = False
        scaleFactor = 1000
        interpAlgor = 'cubic'
        interpAlgorResample = 'near'
        initClearSkyRegionDist = 3000
        initClearSkyRegionMinSize = 3000
        finalClearSkyRegionDist = 1000
        clearSkyMorphSize = 21
        fullImgOuts = False
        checkOutputs = False
        classmlclouds = False
        cloudtrainclouds = None
        cloudtrainother = None
        resample2LowResImg = False
        needAtmModel = False
        prodsToCalc = dict()
        prodsCalculated = dict()
        aeroProfile = None
        atmosProfile = None
        grdRefl = None
        useBRDF = False
        reproject = False
        useWKT2Reproject = True
        xPxlRes = 0.0
        yPxlRes = 0.0
        pxlResDefd = False
        reProjStr = ''
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
        srefDOSImage=""
        srefImage=""
        srefDOSWholeImage=""
        sref6SWholeImage=""
        aotFile=""
        cloudsImage=""
        clearskyImage=""
        outDEMName=""
        outDEMNameMsk=""
        demNoDataVal = -32768.0
        topoShadowImage=""
        footprintVecFile=""
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
        fileEnding2Keep = None
        cloud_methods = None

def prepParametersObj(inputHeader, inputImage, cloudMaskUsrImg, sensorStr, inWKTFile, outFormat, outFilePath, outBaseName, outWKTFile, outProj4File, projAbbv, xPxlResUsr, yPxlResUsr, productsStr, calcStatsPy, aeroProfileOption, atmosProfileOption, aeroProfileOptionImg, atmosProfileOptionImg,  grdReflOption, surfaceAltitude, atmosOZoneVal, atmosWaterVal, atmosOZoneWaterSpecified, aeroWaterVal, aeroDustVal, aeroOceanicVal, aeroSootVal, aeroComponentsSpecified, aotVal, visVal, tmpPath, minAOT, maxAOT, lowAOT, upAOT, demFile, demNoDataUsrVal, aotFile, globalDOS, dosOutRefl, simpleDOS, debugMode, scaleFactor, interpAlgor, interpAlgorResample, initClearSkyRegionDist, initClearSkyRegionMinSize, finalClearSkyRegionDist, clearSkyMorphSize, fullImgOuts, checkOutputs, classmlclouds, cloudtrainclouds, cloudtrainother, resample2LowResImg, fileEnding2Keep, cloud_methods):
    """
    """
    arcsiUtils = ARCSIUtils()
    rsgisUtils = rsgislib.RSGISPyUtils()

    paramsObj = ARCSIParamsObj()

    paramsObj.inputHeader = inputHeader
    paramsObj.inputImage = inputImage
    paramsObj.sensorStr = sensorStr
    paramsObj.productsStr = productsStr
    paramsObj.tmpPath = tmpPath
    paramsObj.outFilePath = outFilePath
    paramsObj.outFormat = outFormat
    paramsObj.outFormatExt = arcsiUtils.getFileExtension(outFormat)
    paramsObj.calcStatsPy = calcStatsPy
    paramsObj.outBaseName = outBaseName
    paramsObj.cloudMaskUsrImg = cloudMaskUsrImg
    paramsObj.inWKTFile = inWKTFile
    paramsObj.outWKTFile = outWKTFile
    paramsObj.outProj4File = outProj4File
    paramsObj.projAbbv = projAbbv
    paramsObj.xPxlResUsr = xPxlResUsr
    paramsObj.yPxlResUsr = yPxlResUsr        
    paramsObj.aeroProfileOption = aeroProfileOption
    paramsObj.atmosProfileOption = atmosProfileOption
    paramsObj.aeroProfileOptionImg = aeroProfileOptionImg
    paramsObj.atmosProfileOptionImg = atmosProfileOptionImg
    paramsObj.grdReflOption = grdReflOption
    paramsObj.surfaceAltitude = surfaceAltitude
    paramsObj.atmosOZoneVal = atmosOZoneVal
    paramsObj.atmosWaterVal = atmosWaterVal
    paramsObj.atmosOZoneWaterSpecified = atmosOZoneWaterSpecified
    paramsObj.aeroWaterVal = aeroWaterVal
    paramsObj.aeroDustVal = aeroDustVal
    paramsObj.aeroOceanicVal = aeroOceanicVal
    paramsObj.aeroSootVal = aeroSootVal
    paramsObj.aeroComponentsSpecified = aeroComponentsSpecified
    paramsObj.aotVal = aotVal
    paramsObj.visVal = visVal
    paramsObj.minAOT = minAOT
    paramsObj.maxAOT = maxAOT
    paramsObj.lowAOT = lowAOT
    paramsObj.upAOT = upAOT
    paramsObj.demFile = demFile
    paramsObj.aotFile = aotFile
    paramsObj.globalDOS = globalDOS
    paramsObj.dosOutRefl = dosOutRefl
    paramsObj.simpleDOS = simpleDOS
    paramsObj.debugMode = debugMode
    paramsObj.scaleFactor = scaleFactor
    paramsObj.interpAlgor = interpAlgor
    paramsObj.interpAlgorResample = interpAlgorResample
    paramsObj.initClearSkyRegionDist = initClearSkyRegionDist
    paramsObj.initClearSkyRegionMinSize = initClearSkyRegionMinSize
    paramsObj.finalClearSkyRegionDist = finalClearSkyRegionDist
    paramsObj.clearSkyMorphSize = clearSkyMorphSize
    paramsObj.fullImgOuts = fullImgOuts
    paramsObj.checkOutputs = checkOutputs
    paramsObj.classmlclouds = classmlclouds
    paramsObj.cloudtrainclouds = cloudtrainclouds
    paramsObj.cloudtrainother = cloudtrainother
    paramsObj.resample2LowResImg = resample2LowResImg
    paramsObj.fileEnding2Keep = fileEnding2Keep
    paramsObj.cloud_methods = cloud_methods

    # Read WKT file if provided.
    paramsObj.wktStr = None
    if paramsObj.inWKTFile is not None:
        paramsObj.wktStr = arcsiUtils.readTextFile(paramsObj.inWKTFile)

    # Step 1: Get the Sensor specific class from factory
    sensorFact = ARCSISensorFactory()
    paramsObj.sensorClass = sensorFact.getSensorClassFromName(paramsObj.sensorStr, paramsObj.debugMode, paramsObj.inputImage)
    # Step 2: Read header parameters
    paramsObj.sensorClass.extractHeaderParameters(paramsObj.inputHeader, paramsObj.wktStr)
    print("")

    if not paramsObj.sensorClass.expectedImageDataPresent():
        raise ARCSIException("Not all the expected input images are present as listed in the header file.")
    else:
        print("Input imagery as listed in header file is present.\n")

    # Step 3: If aerosol and atmosphere images are specified then sample them to find
    #         the aerosol and atmosphere generic model to use for conversion to SREF
    if paramsObj.aeroProfileOptionImg is not None:
        print("Get aero profile from image...")
        aeroProfileMode = int(rsgislib.imagecalc.getImageBandModeInEnv(paramsObj.aeroProfileOptionImg, 1, 1, None, paramsObj.sensorClass.lonTL, paramsObj.sensorClass.lonBR, paramsObj.sensorClass.latBR, paramsObj.sensorClass.latTL)[0])

        if aeroProfileMode is 1:
            paramsObj.aeroProfileOption = "Maritime"
        elif aeroProfileMode is 2:
            paramsObj.aeroProfileOption = "Continental"
        else:
            raise ARCSIException("The aerosol profile from the input image was not recognised.")
        print("Aerosol Profile = ", paramsObj.aeroProfileOption)
        print("")
    if paramsObj.atmosProfileOptionImg is not None:
        print("Get atmos profile from image...")
        atmosProfileMode = int(rsgislib.imagecalc.getImageBandModeInEnv(paramsObj.atmosProfileOptionImg, 1, 1, None, paramsObj.sensorClass.lonTL, paramsObj.sensorClass.lonBR, paramsObj.sensorClass.latBR, paramsObj.sensorClass.latTL)[0])
        summerWinter = arcsiUtils.isSummerOrWinter(paramsObj.sensorClass.latCentre, paramsObj.sensorClass.lonCentre, paramsObj.sensorClass.acquisitionTime )
        if atmosProfileMode is 1:
            paramsObj.atmosProfileOption = "Tropical"
        elif atmosProfileMode is 2:
            if summerWinter is 1:
                paramsObj.atmosProfileOption = "MidlatitudeSummer"
            elif summerWinter is 2:
                paramsObj.atmosProfileOption = "MidlatitudeWinter"
            else:
                raise ARCSIException("Not recognised as being summer or winter.")
        elif atmosProfileMode is 3:
            if summerWinter is 1:
                paramsObj.atmosProfileOption = "SubarcticSummer"
            elif summerWinter is 2:
                paramsObj.atmosProfileOption = "SubarcticWinter"
            else:
                raise ARCSIException("Not recognised as being summer or winter.")
        else:
            raise ARCSIException("The atmosphere profile from the input image was not recognised.")
        print("Atmosphere Profile = ", paramsObj.atmosProfileOption)
        print("")

    # Step 3: Get Output Image Base Name.
    if (paramsObj.outBaseName is None) or (paramsObj.outBaseName == ""):
        paramsObj.outBaseName = paramsObj.sensorClass.generateOutputBaseName()
        if not paramsObj.projAbbv is None:
            paramsObj.outBaseNameProj = paramsObj.outBaseName + "_" + str(paramsObj.projAbbv)
    print("Image Base Name: " + paramsObj.outBaseName + "\n")

    # Check whether output files for this input already exist - if checkOutputs is True.
    if paramsObj.checkOutputs:
        prevOutFiles = glob.glob(os.path.join(outFilePath, paramsObj.outBaseName+'*'))
        if len(prevOutFiles) > 0:
            sys.stderr.write('Error outputs already exist: \'' + inputHeader + '\'\n')
            for tmpFile in prevOutFiles:
                sys.stderr.write('\tFile: \'' + tmpFile + '\'\n')
            raise ARCSIException("Output files already exist and the \'check outputs\' option (--checkouts) was specified so can't continue.")

    # Step 4: Find the products which are to be generated.
    # Make a copy of the dictionary to store calculated products.
    paramsObj.prodsToCalc = dict()
    paramsObj.prodsToCalc["RAD"] = False
    paramsObj.prodsToCalc["TOA"] = False
    paramsObj.prodsToCalc["CLOUDS"] = False
    paramsObj.prodsToCalc["CLEARSKY"] = False
    paramsObj.prodsToCalc["DDVAOT"] = False
    paramsObj.prodsToCalc["DOSAOT"] = False
    paramsObj.prodsToCalc["DOSAOTSGL"] = False
    paramsObj.prodsToCalc["SREF"] = False
    paramsObj.prodsToCalc["DOS"] = False
    paramsObj.prodsToCalc["THERMAL"] = False
    paramsObj.prodsToCalc["SATURATE"] = False
    paramsObj.prodsToCalc["TOPOSHADOW"] = False
    paramsObj.prodsToCalc["FOOTPRINT"] = False
    paramsObj.prodsToCalc["METADATA"] = False
    paramsObj.prodsToCalc["STDSREF"] = False
    paramsObj.prodsToCalc["SHARP"] = False

    # Make a copy of the dictionary to store calculated products.
    paramsObj.prodsCalculated = copy.copy(paramsObj.prodsToCalc)
    paramsObj.needAtmModel = False
    for prod in paramsObj.productsStr:
        if prod == 'RAD':
            paramsObj.prodsToCalc["RAD"] = True
        elif prod == 'SATURATE':
            paramsObj.prodsToCalc["SATURATE"] = True
        elif prod == 'TOA':
            paramsObj.prodsToCalc["RAD"] = True
            paramsObj.prodsToCalc["TOA"] = True
        elif prod == 'CLOUDS':
            paramsObj.prodsToCalc["RAD"] = True
            paramsObj.prodsToCalc["TOA"] = True
            paramsObj.prodsToCalc["CLOUDS"] = True
            paramsObj.prodsToCalc["SATURATE"] = True
            if paramsObj.sensorClass.hasThermal():
                paramsObj.prodsToCalc["THERMAL"] = True
        elif prod == 'CLEARSKY':
            paramsObj.prodsToCalc["RAD"] = True
            paramsObj.prodsToCalc["TOA"] = True
            paramsObj.prodsToCalc["CLOUDS"] = True
            paramsObj.prodsToCalc["SATURATE"] = True
            if paramsObj.sensorClass.hasThermal():
                paramsObj.prodsToCalc["THERMAL"] = True
            paramsObj.prodsToCalc["CLEARSKY"] = True
        elif prod == 'DDVAOT':
            paramsObj.prodsToCalc["RAD"] = True
            paramsObj.prodsToCalc["TOA"] = True
            paramsObj.prodsToCalc["DDVAOT"] = True
            paramsObj.needAtmModel = True
        elif prod == 'DOSAOTSGL':
            paramsObj.prodsToCalc["RAD"] = True
            paramsObj.prodsToCalc["TOA"] = True
            paramsObj.prodsToCalc["DOSAOTSGL"] = True
            paramsObj.needAtmModel = True
        elif prod == 'SREF':
            paramsObj.prodsToCalc["RAD"] = True
            paramsObj.prodsToCalc["SREF"] = True
            paramsObj.needAtmModel = True
        elif prod == 'STDSREF':
            paramsObj.prodsToCalc["RAD"] = True
            paramsObj.prodsToCalc["SREF"] = True
            paramsObj.prodsToCalc["STDSREF"] = True
            paramsObj.prodsToCalc["TOPOSHADOW"] = True
            paramsObj.needAtmModel = True
            if (paramsObj.demFile is None) or (paramsObj.demFile == ""):
                raise ARCSIException("STDSREF requires a DEM file to be provided.")
        elif prod == 'DOS':
            paramsObj.prodsToCalc["RAD"] = True
            paramsObj.prodsToCalc["TOA"] = True
            paramsObj.prodsToCalc["DOS"] = True
        elif prod == 'DOSAOT':
            paramsObj.prodsToCalc["RAD"] = True
            paramsObj.prodsToCalc["TOA"] = True
            paramsObj.prodsToCalc["DOSAOT"] = True
            paramsObj.needAtmModel = True
        elif prod == 'THERMAL':
            if paramsObj.sensorClass.hasThermal():
                paramsObj.prodsToCalc["THERMAL"] = True
            else:
                raise ARCSIException("The sensor does not have thermal bands. Check you inputs.")
        elif prod == 'TOPOSHADOW':
            paramsObj.prodsToCalc["RAD"] = True
            paramsObj.prodsToCalc["TOPOSHADOW"] = True
            if (paramsObj.demFile is None) or (paramsObj.demFile == ""):
                raise ARCSIException("STDSREF requires a DEM file to be provided.")
        elif prod == 'FOOTPRINT':
            paramsObj.prodsToCalc["FOOTPRINT"] = True
        elif prod == 'METADATA':
            paramsObj.prodsToCalc["METADATA"] = True
        elif prod == 'SHARP':
            paramsObj.prodsToCalc["RAD"] = True
            paramsObj.prodsToCalc["SHARP"] = True
            if paramsObj.resample2LowResImg:
                raise ARCSIException("'SHARP' product is only appropriate if lower resolution images bands are being sharpened to a higher resolution.")

    if paramsObj.prodsToCalc["DOSAOT"] and paramsObj.prodsToCalc["DDVAOT"]:
        raise ARCSIException("You cannot specify both the DOSAOT and DDVAOT products, you must choose one or the other.")
    
    paramsObj.aeroProfile = None
    paramsObj.atmosProfile = None
    paramsObj.grdRefl = None
    paramsObj.useBRDF = False

    if paramsObj.needAtmModel:
        if paramsObj.aeroProfileOption is None:
            raise ARCSIException("An aersol profile has not been specified.")
        elif paramsObj.aeroProfileOption == "NoAerosols":
            paramsObj.aeroProfile = Py6S.AeroProfile.PredefinedType(Py6S.AeroProfile.NoAerosols)
        elif paramsObj.aeroProfileOption == "Continental":
            paramsObj.aeroProfile = Py6S.AeroProfile.PredefinedType(Py6S.AeroProfile.Continental)
        elif paramsObj.aeroProfileOption == "Maritime":
            paramsObj.aeroProfile = Py6S.AeroProfile.PredefinedType(Py6S.AeroProfile.Maritime)
        elif paramsObj.aeroProfileOption == "Urban":
            paramsObj.aeroProfile = Py6S.AeroProfile.PredefinedType(Py6S.AeroProfile.Urban)
        elif paramsObj.aeroProfileOption == "Desert":
            paramsObj.aeroProfile = Py6S.AeroProfile.PredefinedType(Py6S.AeroProfile.Desert)
        elif paramsObj.aeroProfileOption == "BiomassBurning":
            paramsObj.aeroProfile = Py6S.AeroProfile.PredefinedType(Py6S.AeroProfile.BiomassBurning)
        elif paramsObj.aeroProfileOption == "Stratospheric":
            paramsObj.aeroProfile = Py6S.AeroProfile.PredefinedType(Py6S.AeroProfile.Stratospheric)
        else:
            raise ARCSIException("The specified aersol profile is unknown.")

        if paramsObj.atmosProfileOption is None:
            raise ARCSIException("An atmospheric profile has not been specified.")
        elif paramsObj.atmosProfileOption == "NoGaseousAbsorption":
            paramsObj.atmosProfile = Py6S.AtmosProfile.PredefinedType(Py6S.AtmosProfile.NoGaseousAbsorption)
        elif paramsObj.atmosProfileOption == "Tropical":
            paramsObj.atmosProfile = Py6S.AtmosProfile.PredefinedType(Py6S.AtmosProfile.Tropical)
        elif paramsObj.atmosProfileOption == "MidlatitudeSummer":
            paramsObj.atmosProfile = Py6S.AtmosProfile.PredefinedType(Py6S.AtmosProfile.MidlatitudeSummer)
        elif paramsObj.atmosProfileOption == "MidlatitudeWinter":
            paramsObj.atmosProfile = Py6S.AtmosProfile.PredefinedType(Py6S.AtmosProfile.MidlatitudeWinter)
        elif paramsObj.atmosProfileOption == "SubarcticSummer":
            paramsObj.atmosProfile = Py6S.AtmosProfile.PredefinedType(Py6S.AtmosProfile.SubarcticSummer)
        elif paramsObj.atmosProfileOption == "SubarcticWinter":
            paramsObj.atmosProfile = Py6S.AtmosProfile.PredefinedType(Py6S.AtmosProfile.SubarcticWinter)
        elif paramsObj.atmosProfileOption == "USStandard1962":
            paramsObj.atmosProfile = Py6S.AtmosProfile.PredefinedType(Py6S.AtmosProfile.USStandard1962)
        else:
            raise ARCSIException("The specified atmospheric profile is unknown.")

        if paramsObj.atmosOZoneWaterSpecified:
            paramsObj.atmosProfile = Py6S.AtmosProfile.UserWaterAndOzone(atmosWaterVal, atmosOZoneVal)

        if paramsObj.grdReflOption is None:
            raise ARCSIException("A ground reflectance has not been specified.")
        elif paramsObj.grdReflOption == "GreenVegetation":
              paramsObj.grdRefl = Py6S.GroundReflectance.HomogeneousLambertian(Py6S.GroundReflectance.GreenVegetation)
        elif paramsObj.grdReflOption == "ClearWater":
            paramsObj.grdRefl = Py6S.GroundReflectance.HomogeneousLambertian(Py6S.GroundReflectance.ClearWater)
        elif paramsObj.grdReflOption == "Sand":
            paramsObj.grdRefl = Py6S.GroundReflectance.HomogeneousLambertian(Py6S.GroundReflectance.Sand)
        elif paramsObj.grdReflOption == "LakeWater":
            paramsObj.grdRefl = Py6S.GroundReflectance.HomogeneousLambertian(Py6S.GroundReflectance.LakeWater)
        elif paramsObj.grdReflOption == "BRDFHapke":
            paramsObj.grdRefl = Py6S.GroundReflectance.HomogeneousHapke(0.101, -0.263, 0.589, 0.046)
            paramsObj.useBRDF = True
        else:
            raise ARCSIException("The specified ground reflectance is unknown.")

    # Decide is the image needs to be reprojected, if so define parameters.
    paramsObj.reproject = False
    paramsObj.useWKT2Reproject = True
    paramsObj.xPxlRes = 0.0
    paramsObj.yPxlRes = 0.0
    paramsObj.pxlResDefd = False
    paramsObj.reProjStr = ''
    if paramsObj.outWKTFile is not None:
        paramsObj.reproject = True
        paramsObj.useWKT2Reproject = True
        paramsObj.reProjStr = paramsObj.outWKTFile
    elif outProj4File is not None:
        paramsObj.reproject = True
        paramsObj.useWKT2Reproject = False
        paramsObj.reProjStr = '"' + arcsiUtils.readTextFile(paramsObj.outProj4File) + '"'

    if (paramsObj.xPxlResUsr is None) or (paramsObj.yPxlResUsr is None):
        paramsObj.pxlResDefd = False
    else:
        paramsObj.pxlResDefd = True
        paramsObj.xPxlRes = paramsObj.xPxlResUsr
        paramsObj.yPxlRes = paramsObj.yPxlResUsr

    paramsObj.sensorClass.setReProjectOutputs(paramsObj.reproject)

    paramsObj.validMaskImage=None
    paramsObj.validMaskImageProj=""
    paramsObj.viewAngleImg=""
    paramsObj.viewAngleImgProj=""
    paramsObj.radianceImage=""
    paramsObj.radianceImageWhole=""
    paramsObj.saturateImage=""
    paramsObj.saturateImageProj=""
    paramsObj.thermalRadImage=""
    paramsObj.thermalBrightImage=""
    paramsObj.thermalBrightImageWhole=""
    paramsObj.maskImage=""
    paramsObj.toaImage=""
    paramsObj.toaImageWhole=""
    paramsObj.srefDOSImage = ""
    paramsObj.srefImage=""
    paramsObj.srefDOSWholeImage=""
    paramsObj.sref6SWholeImage=""
    paramsObj.aotFile=""
    paramsObj.cloudsImage=""
    paramsObj.clearskyImage=""
    paramsObj.outDEMName=""
    paramsObj.outDEMNameMsk=""
    paramsObj.demNoDataVal = -32768.0
    if demNoDataUsrVal is not None:
        paramsObj.demNoDataVal = demNoDataUsrVal
    elif paramsObj.demFile is not None:
        paramsObj.demNoDataVal = rsgisUtils.getImageNoDataValue(paramsObj.demFile)
        if paramsObj.demNoDataVal is None:
            raise ARCSIException("A no data value for the inputted DEM has not been defined - cannot continue without a no data value. A no data value can be define using the --demnodata option.")
    paramsObj.topoShadowImage=""
    paramsObj.footprintVecFile=""
    paramsObj.metaDataFile=""
    paramsObj.propOfCloud = 0.0
    paramsObj.propOfClearSky = 0.0
    paramsObj.projImgBBOX = dict()
    paramsObj.projImgBBOX['MinX'] = 0.0
    paramsObj.projImgBBOX['MaxX'] = 0.0
    paramsObj.projImgBBOX['MinY'] = 0.0
    paramsObj.projImgBBOX['MaxY'] = 0.0
    paramsObj.processStageStr = ""
    paramsObj.processStageWholeImgStr = ""
    paramsObj.processSREFStr = ""
    paramsObj.finalOutFiles = dict()
    paramsObj.calcdOutVals = dict()
    paramsObj.sixsLUTCoeffs = None
    paramsObj.aotLUT = False

    if paramsObj.reproject:
        if paramsObj.useWKT2Reproject:
            paramsObj.calcdOutVals["REPROJECT"] = arcsiUtils.readTextFile(paramsObj.outWKTFile)
        else:
            paramsObj.calcdOutVals["REPROJECT"] = arcsiUtils.readTextFile(paramsObj.outProj4File)
    return paramsObj

def checkForValidInput(paramsObj):
    print('Checking Input Images are valid')
    paramsObj.sensorClass.checkInputImageValid()

def resampleBands(paramsObj):
    # Check if bands need resampling
    if paramsObj.sensorClass.inImgsDiffRes():
        print('Resampling image bands to match one another.')
        if paramsObj.interpAlgorResample == 'near':
            paramsObj.interpAlgorResample = 'nearestneighbour'
        paramsObj.sensorClass.resampleImgRes(paramsObj.outFilePath, paramsObj.resample2LowResImg, paramsObj.interpAlgorResample, False)

def mosaicInputImages(paramsObj):
    if paramsObj.sensorClass.imgNeedMosaicking():
        print("Mosacking Input Image Tiles.")
        paramsObj.sensorClass.mosaicImageTiles(paramsObj.outFilePath)

def createValidMaskViewAngle(paramsObj):
    # Get the valid image data maskImage
    rsgisUtils = rsgislib.RSGISPyUtils()
    outName = paramsObj.outBaseName + "_valid" + paramsObj.outFormatExt
    paramsObj.viewAngleImg = os.path.join(paramsObj.outFilePath, paramsObj.outBaseName + "_viewangle" + paramsObj.outFormatExt)
    paramsObj.validMaskImage = paramsObj.sensorClass.generateValidImageDataMask(paramsObj.outFilePath, outName, paramsObj.viewAngleImg, paramsObj.outFormat)
    if not os.path.exists(paramsObj.viewAngleImg):
        paramsObj.viewAngleImg = None
    if paramsObj.validMaskImage is not None:
        rsgislib.rastergis.populateStats(paramsObj.validMaskImage, True, True)
    if paramsObj.viewAngleImg is not None:
        rsgislib.imageutils.popImageStats(paramsObj.viewAngleImg, usenodataval=True, nodataval=1000, calcpyramids=True)
    print("")
    if paramsObj.validMaskImage is not None:
        paramsObj.finalOutFiles["VALID_MASK"] = paramsObj.validMaskImage
    if paramsObj.viewAngleImg is not None:
        paramsObj.finalOutFiles["VIEW_ANGLE"] = paramsObj.viewAngleImg

    if paramsObj.reproject and (paramsObj.validMaskImage is not None):
        if not paramsObj.pxlResDefd:
            validImgDS = gdal.Open(paramsObj.validMaskImage, gdal.GA_ReadOnly)
            if validImgDS is None:
                raise ARCSIException('Could not open the Valid Image Mask ' + paramsObj.validMaskImage)
            geoTransform = validImgDS.GetGeoTransform()
            if geoTransform is None:
                raise ARCSIException('Could read the geotransform from the Valid Image Mask ' + paramsObj.validMaskImage)
            paramsObj.xPxlRes = geoTransform[1]
            paramsObj.yPxlRes = geoTransform[5]
            paramsObj.pxlResDefd = True
            paramsObj.validImgDS = None
        print("Define re-projected image BBOX.")
        paramsObj.projImgBBOX = paramsObj.sensorClass.getReProjBBOX(paramsObj.outWKTFile, paramsObj.outProj4File, paramsObj.useWKT2Reproject, paramsObj.xPxlRes, paramsObj.yPxlRes, True)
        print("Output Image BBOX (TL, BR): [({0:10.2f}, {1:10.2f}), ({2:10.2f}, {3:10.2f})]".format(paramsObj.projImgBBOX['MinX'], paramsObj.projImgBBOX['MaxY'], paramsObj.projImgBBOX['MaxX'], paramsObj.projImgBBOX['MinY']))
        outName = paramsObj.outBaseNameProj + "_valid" + paramsObj.outFormatExt
        paramsObj.validMaskImageProj = os.path.join(paramsObj.outFilePath, outName)

        cmd = 'gdalwarp -t_srs ' + paramsObj.reProjStr + ' -tr ' + str(paramsObj.xPxlRes) + ' ' + str(paramsObj.yPxlRes) + ' -ot Byte -wt Float32 ' \
            + '-te ' + str(paramsObj.projImgBBOX['MinX']) + ' ' + str(paramsObj.projImgBBOX['MinY']) + ' ' + str(paramsObj.projImgBBOX['MaxX']) + ' ' + str(paramsObj.projImgBBOX['MaxY']) \
            + ' -r near -tap -srcnodata 0 -dstnodata 0 -of ' + paramsObj.outFormat + ' -overwrite ' \
            + paramsObj.validMaskImage + ' ' + paramsObj.validMaskImageProj
        print(cmd)
        try:
            subprocess.call(cmd, shell=True)
        except OSError as e:
            raise ARCSIException('Could not re-projection valid image mask: ' + cmd)
        if not os.path.exists(paramsObj.validMaskImageProj):
            raise ARCSIException('Reprojected valid image mask is not present: ' + paramsObj.validMaskImageProj)
        else:
            rsgislib.rastergis.populateStats(paramsObj.validMaskImageProj, True, True)
        paramsObj.finalOutFiles["VALID_MASK"] = paramsObj.validMaskImageProj
        print("")

    if paramsObj.reproject and ((paramsObj.viewAngleImg != "") and (paramsObj.viewAngleImg is not None)):
        if not paramsObj.pxlResDefd:
            viewAngleImgDS = gdal.Open(paramsObj.viewAngleImg, gdal.GA_ReadOnly)
            if viewAngleImgDS is None:
                raise ARCSIException('Could not open the Valid Image Mask ' + paramsObj.viewAngleImg)
            geoTransform = viewAngleImgDS.GetGeoTransform()
            if geoTransform is None:
                raise ARCSIException('Could read the geotransform from the Valid Image Mask ' + paramsObj.viewAngleImg)
            paramsObj.xPxlRes = geoTransform[1]
            paramsObj.yPxlRes = geoTransform[5]
            paramsObj.pxlResDefd = True
            viewAngleImgDS = None
        print("Define re-projected image BBOX.")
        paramsObj.projImgBBOX = paramsObj.sensorClass.getReProjBBOX(paramsObj.outWKTFile, paramsObj.outProj4File, paramsObj.useWKT2Reproject, paramsObj.xPxlRes, paramsObj.yPxlRes, True)
        print("Output Image BBOX (TL, BR): [({0:10.2f}, {1:10.2f}), ({2:10.2f}, {3:10.2f})]".format(paramsObj.projImgBBOX['MinX'], paramsObj.projImgBBOX['MaxY'], paramsObj.projImgBBOX['MaxX'], paramsObj.projImgBBOX['MinY']))
        outName = paramsObj.outBaseNameProj + "_viewangle" + paramsObj.outFormatExt
        paramsObj.viewAngleImgProj = os.path.join(paramsObj.outFilePath, outName)

        cmd = 'gdalwarp -t_srs ' + paramsObj.reProjStr + ' -tr ' + str(paramsObj.xPxlRes) + ' ' + str(paramsObj.yPxlRes) + ' -ot Float32 -wt Float32 ' \
            + '-te ' + str(paramsObj.projImgBBOX['MinX']) + ' ' + str(paramsObj.projImgBBOX['MinY']) + ' ' + str(paramsObj.projImgBBOX['MaxX']) + ' ' + str(paramsObj.projImgBBOX['MaxY']) \
            + ' -r cubicspline -tap -srcnodata 99999 -dstnodata 99999 -of ' + paramsObj.outFormat + ' -overwrite ' \
            + paramsObj.viewAngleImg + ' ' + paramsObj.viewAngleImgProj
        print(cmd)
        try:
            subprocess.call(cmd, shell=True)
        except OSError as e:
            raise ARCSIException('Could not re-projection valid image mask: ' + cmd)
        if not os.path.exists(paramsObj.viewAngleImgProj):
            raise ARCSIException('Reprojected valid image mask is not present: ' + paramsObj.viewAngleImgProj)
        else:
            rsgislib.imageutils.popImageStats(paramsObj.viewAngleImgProj, usenodataval=True, nodataval=99999, calcpyramids=True)
        rsgisUtils.deleteFileWithBasename(paramsObj.viewAngleImg)
        paramsObj.viewAngleImg = paramsObj.viewAngleImgProj
        paramsObj.finalOutFiles["VIEW_ANGLE"] = paramsObj.viewAngleImgProj
        print("")

def createFootprint(paramsObj):
    if paramsObj.prodsToCalc["FOOTPRINT"]:
        if paramsObj.validMaskImage is None:
            raise ARCSIException("To generate a footprint a valid image mask is required - not supported by this sensor?")
        outFootprintLyrName = ''
        if paramsObj.reproject:
            outFootprintLyrName = paramsObj.outBaseNameProj + "_footprint"
            paramsObj.footprintVecFile = paramsObj.sensorClass.generateImageFootprint(paramsObj.validMaskImageProj, paramsObj.outFilePath, outFootprintLyrName)
        else:
            outFootprintLyrName = paramsObj.outBaseName + "_footprint"
            paramsObj.footprintVecFile = paramsObj.sensorClass.generateImageFootprint(paramsObj.validMaskImage, paramsObj.outFilePath, outFootprintLyrName)
        paramsObj.finalOutFiles["FOOTPRINT"] = paramsObj.footprintVecFile
        paramsObj.prodsCalculated["FOOTPRINT"] = True
        print("")

def createSaturatedImage(paramsObj):
    if paramsObj.prodsToCalc["SATURATE"]:
        # Execute generation of the saturation image
        rsgisUtils = rsgislib.RSGISPyUtils()
        outName = paramsObj.outBaseName + "_sat" + paramsObj.outFormatExt
        paramsObj.saturateImage = paramsObj.sensorClass.generateImageSaturationMask(paramsObj.outFilePath, outName, paramsObj.outFormat)

        if paramsObj.reproject:
            outNameProj = paramsObj.outBaseNameProj + "_sat" + paramsObj.outFormatExt
            paramsObj.saturateImageProj = os.path.join(paramsObj.outFilePath, outNameProj)
            cmd = 'gdalwarp -t_srs ' + paramsObj.reProjStr + ' -tr ' + str(paramsObj.xPxlRes) + ' ' + str(paramsObj.yPxlRes) + ' -ot Byte ' \
            + '-te ' + str(paramsObj.projImgBBOX['MinX']) + ' ' + str(paramsObj.projImgBBOX['MinY']) + ' ' + str(paramsObj.projImgBBOX['MaxX']) + ' ' + str(paramsObj.projImgBBOX['MaxY']) \
            + ' -r near -tap -of ' + paramsObj.outFormat + ' -overwrite ' + paramsObj.saturateImage + ' ' + paramsObj.saturateImageProj
            print(cmd)
            try:
                subprocess.call(cmd, shell=True)
            except OSError as e:
                raise ARCSIException('Could not re-projection saturated image mask: ' + cmd)
            if not os.path.exists(paramsObj.saturateImageProj):
                raise ARCSIException('Reprojected saturated image mask is not present: ' + paramsObj.saturateImageProj)
            else:
                rsgisUtils.deleteFileWithBasename(paramsObj.saturateImage)
                paramsObj.saturateImage = paramsObj.saturateImageProj
        if paramsObj.calcStatsPy:
            print("Calculating Statistics...")
            rsgislib.rastergis.populateStats(paramsObj.saturateImage, True, True)
        paramsObj.prodsCalculated["SATURATE"] = True
        print("")

def convertInputImageToRadiance(paramsObj):
    # Convert to Radiance
    if paramsObj.prodsToCalc["RAD"]:
        rsgisUtils = rsgislib.RSGISPyUtils()
        # Execute conversion to radiance
        outName = paramsObj.outBaseName + "_rad" + paramsObj.outFormatExt
        outThermName = None
        if paramsObj.prodsToCalc["THERMAL"]:
            outThermName = paramsObj.outBaseName + "_therm_rad" + paramsObj.outFormatExt
        paramsObj.radianceImage, paramsObj.thermalRadImage = paramsObj.sensorClass.convertImageToRadiance(paramsObj.outFilePath, outName, outThermName, paramsObj.outFormat)

        if paramsObj.sensorClass.maskInputImages():
            paramsObj.processStageStr = paramsObj.processStageStr + "_msk"
            outImgName = paramsObj.outBaseName + paramsObj.processStageStr + "_rad" + paramsObj.outFormatExt
            outMaskName = paramsObj.outBaseName + "_mask" + paramsObj.outFormatExt
            radianceImageTmp, paramsObj.maskImage = paramsObj.sensorClass.applyImageDataMask(paramsObj.inputHeader, paramsObj.radianceImage, paramsObj.outFilePath, outMaskName, outImgName, paramsObj.outFormat, None)
            if radianceImageTmp is not paramsObj.radianceImage:
                rsgisUtils.deleteFileWithBasename(paramsObj.radianceImage)
                paramsObj.radianceImage = radianceImageTmp
            if paramsObj.maskImage is not None:
                print("Setting Band Names...")
                paramsObj.sensorClass.setBandNames(paramsObj.radianceImage)
                if paramsObj.calcStatsPy:
                    print("Calculating Statistics...")
                    rsgislib.imageutils.popImageStats(paramsObj.radianceImage, True, 0.0, True)
                    rsgislib.rastergis.populateStats(paramsObj.maskImage, True, True)

        if paramsObj.validMaskImage is not None:
            print("Masking to valid data area.")
            paramsObj.processStageStr = paramsObj.processStageStr + "_vmsk"
            outRadPathName = os.path.join(paramsObj.outFilePath, paramsObj.outBaseName + paramsObj.processStageStr + "_rad" + paramsObj.outFormatExt)
            rsgislib.imageutils.maskImage(paramsObj.radianceImage, paramsObj.validMaskImage, outRadPathName, paramsObj.outFormat, rsgisUtils.getRSGISLibDataTypeFromImg(paramsObj.radianceImage), 0.0, 0.0)
            rsgisUtils.deleteFileWithBasename(paramsObj.radianceImage)
            paramsObj.radianceImage = outRadPathName
            if paramsObj.thermalRadImage is not None:
                outThermPathName = os.path.join(paramsObj.outFilePath, paramsObj.outBaseName + paramsObj.processStageStr + "_thermrad" + paramsObj.outFormatExt)
                rsgislib.imageutils.maskImage(paramsObj.thermalRadImage, paramsObj.validMaskImage, outThermPathName, paramsObj.outFormat, rsgisUtils.getRSGISLibDataTypeFromImg(paramsObj.thermalRadImage), 0.0, 0.0)
                rsgisUtils.deleteFileWithBasename(paramsObj.thermalRadImage)
                paramsObj.thermalRadImage = outThermPathName
            if paramsObj.reproject:
                rsgisUtils.deleteFileWithBasename(paramsObj.validMaskImage)
                paramsObj.validMaskImage = paramsObj.validMaskImageProj
        
        if paramsObj.prodsToCalc["SHARP"]:
            print("Sharpen radiance image...")
            paramsObj.processStageStr = paramsObj.processStageStr + "_sharp"
            outRadSharpImageName = paramsObj.outBaseName + paramsObj.processStageStr + "_rad" + paramsObj.outFormatExt
            outRadSharpImage = os.path.join(paramsObj.outFilePath, outRadSharpImageName)
            paramsObj.sensorClass.sharpenLowResRadImgBands(paramsObj.radianceImage, outRadSharpImage, paramsObj.outFormat)
            rsgisUtils.deleteFileWithBasename(paramsObj.radianceImage)
            paramsObj.radianceImage = outRadSharpImage
            paramsObj.prodsCalculated['SHARP'] = True

        if paramsObj.reproject:
            if paramsObj.radianceImage is not None:
                outName = paramsObj.outBaseNameProj + paramsObj.processStageStr + "_rad" + paramsObj.outFormatExt

                outRadImagePath = os.path.join(paramsObj.outFilePath, outName)
                cmd = 'gdalwarp -t_srs ' + paramsObj.reProjStr + ' -tr ' + str(paramsObj.xPxlRes) + ' ' + str(paramsObj.yPxlRes) + ' -ot Float32 -wt Float32 ' \
                + '-te ' + str(paramsObj.projImgBBOX['MinX']) + ' ' + str(paramsObj.projImgBBOX['MinY']) + ' ' + str(paramsObj.projImgBBOX['MaxX']) + ' ' + str(paramsObj.projImgBBOX['MaxY']) \
                + ' -r ' + paramsObj.interpAlgor + ' -tap -srcnodata 0 -dstnodata 0 -of ' + paramsObj.outFormat + ' -overwrite ' \
                + paramsObj.radianceImage + ' ' + outRadImagePath
                print(cmd)
                try:
                    subprocess.call(cmd, shell=True)
                except OSError as e:
                    raise ARCSIException('Could not re-projection radiance image: ' + cmd)
                if not os.path.exists(outRadImagePath):
                    raise ARCSIException('Reprojected radiance image is not present: ' + outRadImagePath)
                else:
                    rsgisUtils.deleteFileWithBasename(paramsObj.radianceImage)
                    paramsObj.radianceImage = outRadImagePath
            if paramsObj.thermalRadImage is not None:
                outName = paramsObj.outBaseNameProj + paramsObj.processStageStr + "_thrad" + paramsObj.outFormatExt
                outThermRadImagePath = os.path.join(paramsObj.outFilePath, outName)
                cmd = 'gdalwarp -t_srs ' + paramsObj.reProjStr + ' -tr ' + str(paramsObj.xPxlRes) + ' ' + str(paramsObj.yPxlRes) + ' -ot Float32 -wt Float32 ' \
                + '-te ' + str(paramsObj.projImgBBOX['MinX']) + ' ' + str(paramsObj.projImgBBOX['MinY']) + ' ' + str(paramsObj.projImgBBOX['MaxX']) + ' ' + str(paramsObj.projImgBBOX['MaxY']) \
                + ' -r ' + paramsObj.interpAlgor + ' -tap -srcnodata 0 -dstnodata 0 -of ' + paramsObj.outFormat + ' -overwrite ' \
                + paramsObj.thermalRadImage + ' ' + outThermRadImagePath
                try:
                    subprocess.call(cmd, shell=True)
                except OSError as e:
                    raise ARCSIException('Could not re-projection thermal radiance image: ' + cmd)
                if not os.path.exists(outThermRadImagePath):
                    raise ARCSIException('Reprojected thermal radiance image is not present: ' + outThermRadImagePath)
                else:
                    rsgisUtils.deleteFileWithBasename(paramsObj.thermalRadImage)
                    paramsObj.thermalRadImage = outThermRadImagePath
            paramsObj.outBaseName = paramsObj.outBaseNameProj

        print("Setting Band Names...")
        paramsObj.sensorClass.setBandNames(paramsObj.radianceImage)
        if paramsObj.calcStatsPy:
            print("Calculating Statistics...")
            rsgislib.imageutils.popImageStats(paramsObj.radianceImage, True, 0.0, True)
            if paramsObj.thermalRadImage is not None:
                rsgislib.imageutils.popImageStats(paramsObj.thermalRadImage, True, 0.0, True)
        
        paramsObj.finalOutFiles["RADIANCE_WHOLE"] = paramsObj.radianceImage
        paramsObj.finalOutFiles["RADIANCE"] = paramsObj.radianceImage
        if paramsObj.thermalRadImage is not None:
            paramsObj.finalOutFiles["THERM_RADIANCE_WHOLE"] = paramsObj.thermalRadImage

        paramsObj.radianceImageWhole = paramsObj.radianceImage
        paramsObj.prodsCalculated["RAD"] = True
        print("")

def calcThermalBrightness(paramsObj):
    # Execute calibrate thermal to brightness
    if paramsObj.prodsToCalc["THERMAL"]:
        outName = paramsObj.outBaseName + paramsObj.processStageStr + "_thrad_thermbright" + paramsObj.outFormatExt
        paramsObj.thermalBrightImage = paramsObj.sensorClass.convertThermalToBrightness(paramsObj.thermalRadImage, paramsObj.outFilePath, outName, paramsObj.outFormat, paramsObj.scaleFactor)
        if paramsObj.calcStatsPy:
            print("Calculating Statistics...")
            rsgislib.imageutils.popImageStats(paramsObj.thermalBrightImage, True, 0.0, True)
        paramsObj.finalOutFiles["THERMAL_BRIGHT_WHOLE"] = paramsObj.thermalBrightImage
        paramsObj.finalOutFiles["THERMAL_BRIGHT"] = paramsObj.thermalBrightImage
        paramsObj.thermalBrightImageWhole = paramsObj.thermalBrightImage
        paramsObj.prodsCalculated["THERMAL"] = True
        print("")

def calcTOAReflectance(paramsObj):
    # Execute conversion to top of atmosphere reflectance
    if paramsObj.prodsToCalc["TOA"]:
        outName = paramsObj.outBaseName + paramsObj.processStageStr +"_rad_toa" + paramsObj.outFormatExt
        paramsObj.toaImage = paramsObj.sensorClass.convertImageToTOARefl(paramsObj.radianceImage, paramsObj.outFilePath, outName, paramsObj.outFormat, paramsObj.scaleFactor)

        rsgisUtils = rsgislib.RSGISPyUtils()
        if not rsgisUtils.doGDALLayersHaveSameProj(paramsObj.radianceImage, paramsObj.toaImage):
            tmpImg = os.path.join(paramsObj.tmpPath, paramsObj.outBaseName+rsgisUtils.uidGenerator()+'_toaimg'+paramsObj.outFormatExt)
            rsgisUtils.renameGDALLayer(paramsObj.toaImage, tmpImg)
            interpAlgorOpt = paramsObj.interpAlgor
            if paramsObj.interpAlgor == 'near':
                interpAlgorOpt = 'nearestneighbour'
            rsgislib.imageutils.resampleImage2Match(paramsObj.radianceImage, tmpImg, paramsObj.toaImage, 'KEA', interpAlgorOpt, rsgislib.TYPE_16UINT, multicore=False)
            rsgisUtils.deleteFileWithBasename(tmpImg)

        print("Setting Band Names...")
        paramsObj.sensorClass.setBandNames(paramsObj.toaImage)
        if paramsObj.calcStatsPy:
            print("Calculating Statistics...")
            rsgislib.imageutils.popImageStats(paramsObj.toaImage, True, 0.0, True)
        paramsObj.finalOutFiles["TOA_WHOLE"] = paramsObj.toaImage
        paramsObj.finalOutFiles["TOA"] = paramsObj.toaImage
        paramsObj.toaImageWhole = paramsObj.toaImage
        paramsObj.prodsCalculated["TOA"] = True
        print("")

def performCloudMasking(paramsObj):
    # Generate Cloud and Apply Cloud Masks
    if paramsObj.prodsToCalc["CLOUDS"]:
        print("Perform Cloud Classification...")
        rsgisUtils = rsgislib.RSGISPyUtils()
        outName = paramsObj.outBaseName + "_clouds" + paramsObj.outFormatExt
        if paramsObj.cloudMaskUsrImg is None:
            if paramsObj.classmlclouds:
                paramsObj.cloudsImage = paramsObj.sensorClass.generateCloudMaskML(paramsObj.toaImage, paramsObj.validMaskImage, paramsObj.outFilePath, outName, paramsObj.outFormat, paramsObj.tmpPath, paramsObj.cloudtrainclouds, paramsObj.cloudtrainother, paramsObj.scaleFactor, numCores=1)
            else:
                paramsObj.cloudsImage = paramsObj.sensorClass.generateCloudMask(paramsObj.toaImage, paramsObj.saturateImage, paramsObj.thermalBrightImage, paramsObj.viewAngleImg, paramsObj.validMaskImage, paramsObj.outFilePath, outName, paramsObj.outFormat, paramsObj.tmpPath, paramsObj.scaleFactor, paramsObj.cloud_methods)
            if paramsObj.calcStatsPy:
                print("Calculating Statistics...")
                rsgislib.rastergis.populateStats(paramsObj.cloudsImage, False, True)
        else:
            paramsObj.cloudsImage = paramsObj.cloudMaskUsrImg
        paramsObj.finalOutFiles["CLOUD_MASK"] = paramsObj.cloudsImage

        # Calculate the proportion of the scene cover by cloud.
        paramsObj.propOfCloud = rsgislib.imagecalc.calcPropTrueExp('b1==1?1:b1==2?1:0', [rsgislib.imagecalc.BandDefn('b1', paramsObj.cloudsImage, 1)], paramsObj.validMaskImage)
        print("The scene is " + str(paramsObj.propOfCloud*100) + "% cloud.")
        paramsObj.calcdOutVals['ARCSI_CLOUD_COVER'] = paramsObj.propOfCloud

        if paramsObj.propOfCloud < 0.95: # Less than 95% cloud cover then process.
            print("Applying cloud masks to images...")
            paramsObj.processStageStr = paramsObj.processStageStr + "_mclds"
            outputRADImage = os.path.join(paramsObj.outFilePath, paramsObj.outBaseName + paramsObj.processStageStr + "_rad" + paramsObj.outFormatExt)
            rsgislib.imageutils.maskImage(paramsObj.radianceImage, paramsObj.cloudsImage, outputRADImage, paramsObj.outFormat, rsgisUtils.getRSGISLibDataTypeFromImg(paramsObj.radianceImage), 0, [1,2])
            paramsObj.radianceImage = outputRADImage
            paramsObj.sensorClass.setBandNames(paramsObj.radianceImage)
            paramsObj.finalOutFiles["RADIANCE"] = paramsObj.radianceImage
            if paramsObj.calcStatsPy:
                print("Calculating Statistics...")
                rsgislib.imageutils.popImageStats(paramsObj.radianceImage, True, 0.0, True)
            outputTOAImage = os.path.join(paramsObj.outFilePath, paramsObj.outBaseName + paramsObj.processStageStr + "_rad_toa" + paramsObj.outFormatExt)
            rsgislib.imageutils.maskImage(paramsObj.toaImage, paramsObj.cloudsImage, outputTOAImage, paramsObj.outFormat, rsgisUtils.getRSGISLibDataTypeFromImg(paramsObj.toaImage), 0, [1,2])
            paramsObj.toaImage = outputTOAImage
            paramsObj.sensorClass.setBandNames(paramsObj.toaImage)
            paramsObj.finalOutFiles["TOA"] = paramsObj.toaImage
            if paramsObj.calcStatsPy:
                print("Calculating Statistics...")
                rsgislib.imageutils.popImageStats(paramsObj.toaImage, True, 0.0, True)
        paramsObj.prodsCalculated["CLOUDS"] = True
        print("")

def performClearSkyMasking(paramsObj):
    if paramsObj.prodsToCalc["CLEARSKY"]:
        rsgisUtils = rsgislib.RSGISPyUtils()
        outName = paramsObj.outBaseName + "_clearsky" + paramsObj.outFormatExt
        paramsObj.clearskyImage = paramsObj.sensorClass.generateClearSkyMask(paramsObj.cloudsImage, paramsObj.validMaskImage, paramsObj.outFilePath, outName, paramsObj.outFormat, paramsObj.tmpPath, paramsObj.initClearSkyRegionDist, paramsObj.initClearSkyRegionMinSize, paramsObj.finalClearSkyRegionDist, paramsObj.clearSkyMorphSize)
        if paramsObj.calcStatsPy:
            print("Calculating Statistics...")
            rsgislib.rastergis.populateStats(paramsObj.clearskyImage, True, True)
        paramsObj.finalOutFiles["CLEARSKY_MASK"] = paramsObj.clearskyImage

        # Calculate the proportion of the scene which is clear sky.
        paramsObj.propOfClearSky = rsgislib.imagecalc.calcPropTrueExp('b1==1?1:0', [rsgislib.imagecalc.BandDefn('b1', paramsObj.clearskyImage, 1)], paramsObj.validMaskImage)
        print("The scene is " + str(paramsObj.propOfClearSky*100) + "% clear-sky.")
        paramsObj.calcdOutVals['ARCSI_CLEARSKY_COVER'] = paramsObj.propOfClearSky

        if paramsObj.propOfClearSky > 0.05: # Keep going if at least 5% of the scene is clear sky
            print("Applying clear-sky masks to images...")
            paramsObj.processStageStr = paramsObj.processStageStr + "_clearsky"
            outputRADImage = os.path.join(paramsObj.outFilePath, paramsObj.outBaseName + paramsObj.processStageStr + "_rad" + paramsObj.outFormatExt)
            rsgislib.imageutils.maskImage(paramsObj.radianceImage, paramsObj.clearskyImage, outputRADImage, paramsObj.outFormat, rsgisUtils.getRSGISLibDataTypeFromImg(paramsObj.radianceImage), 0, 0)
            paramsObj.radianceImage = outputRADImage
            paramsObj.sensorClass.setBandNames(paramsObj.radianceImage)
            paramsObj.finalOutFiles["RADIANCE"] = paramsObj.radianceImage
            if paramsObj.calcStatsPy:
                print("Calculating Statistics...")
                rsgislib.imageutils.popImageStats(paramsObj.radianceImage, True, 0.0, True)
            outputTOAImage = os.path.join(paramsObj.outFilePath, paramsObj.outBaseName + paramsObj.processStageStr + "_rad_toa" + paramsObj.outFormatExt)
            rsgislib.imageutils.maskImage(paramsObj.toaImage, paramsObj.clearskyImage, outputTOAImage, paramsObj.outFormat, rsgisUtils.getRSGISLibDataTypeFromImg(paramsObj.toaImage), 0, 0)
            paramsObj.toaImage = outputTOAImage
            paramsObj.sensorClass.setBandNames(paramsObj.toaImage)
            paramsObj.finalOutFiles["TOA"] = paramsObj.toaImage
            if paramsObj.calcStatsPy:
                print("Calculating Statistics...")
                rsgislib.imageutils.popImageStats(paramsObj.toaImage, True, 0.0, True)
        paramsObj.prodsCalculated["CLEARSKY"] = True
        print("")

def prepareDEM(paramsObj):
    if (paramsObj.demFile is not None) and (paramsObj.demFile != ""):
        rsgisUtils = rsgislib.RSGISPyUtils()
        outDEMNameTmp = os.path.join(paramsObj.outFilePath, (paramsObj.outBaseName + "_demtmp" + paramsObj.outFormatExt))
        rsgislib.imageutils.createCopyImage(paramsObj.radianceImage, outDEMNameTmp, 1, paramsObj.demNoDataVal, paramsObj.outFormat, rsgislib.TYPE_32FLOAT)

        inDEMDS = gdal.Open(paramsObj.demFile, gdal.GA_ReadOnly)
        outDEMDS = gdal.Open(outDEMNameTmp, gdal.GA_Update)

        print("Subset and reproject DEM...")
        gdal.ReprojectImage(inDEMDS, outDEMDS, None, None, gdal.GRA_CubicSpline)
        inDEMDS = None
        outDEMDS = None

        paramsObj.outDEMName = os.path.join(paramsObj.outFilePath, (paramsObj.outBaseName + "_dem" + paramsObj.outFormatExt))
        print("Output DEM: ", paramsObj.outDEMName)
        rsgislib.imageutils.maskImage(outDEMNameTmp, paramsObj.validMaskImage, paramsObj.outDEMName, paramsObj.outFormat, rsgislib.TYPE_32FLOAT, paramsObj.demNoDataVal, 0)

        paramsObj.outDEMNameMsk = os.path.join(paramsObj.outFilePath, (paramsObj.outBaseName + "_demmsk" + paramsObj.outFormatExt))
        if paramsObj.prodsToCalc["CLEARSKY"]:
            rsgislib.imageutils.maskImage(paramsObj.outDEMName, paramsObj.clearskyImage, paramsObj.outDEMNameMsk, paramsObj.outFormat, rsgislib.TYPE_32FLOAT, paramsObj.demNoDataVal, 0)
        elif paramsObj.prodsToCalc["CLOUDS"]:
            rsgislib.imageutils.maskImage(paramsObj.outDEMName, paramsObj.cloudsImage, paramsObj.outDEMNameMsk, paramsObj.outFormat, rsgislib.TYPE_32FLOAT, paramsObj.demNoDataVal, [1,2,3])
        else:
            paramsObj.outDEMNameMsk = paramsObj.outDEMName

        # Calculate DEM statistics and set no data value.
        if paramsObj.prodsToCalc["CLEARSKY"] or paramsObj.prodsToCalc["CLOUDS"]:
            rsgislib.imageutils.popImageStats(paramsObj.outDEMNameMsk, True, paramsObj.demNoDataVal, True)
        rsgislib.imageutils.popImageStats(paramsObj.outDEMName, True, paramsObj.demNoDataVal, True)

        # Remove tmp DEM file.
        rsgisUtils.deleteFileWithBasename(outDEMNameTmp)
        paramsObj.finalOutFiles["IMAGE_DEM"] = paramsObj.outDEMName

def calcTopoShadowMask(paramsObj):
    if paramsObj.prodsToCalc["TOPOSHADOW"]:
        rsgisUtils = rsgislib.RSGISPyUtils()
        outName = paramsObj.outBaseName + "_toposhad" + paramsObj.outFormatExt
        tmpDEMFile = paramsObj.outDEMNameMsk
        if paramsObj.fullImgOuts:
            tmpDEMFile = paramsObj.outDEMName
        paramsObj.topoShadowImage = paramsObj.sensorClass.generateTopoDirectShadowMask(tmpDEMFile, paramsObj.outFilePath, outName, paramsObj.outFormat, paramsObj.tmpPath)
        if paramsObj.calcStatsPy:
            print("Calculating Statistics...")
            rsgislib.rastergis.populateStats(paramsObj.topoShadowImage, True, True)
        paramsObj.finalOutFiles["TOPO_SHADOW_MASK"] = paramsObj.topoShadowImage

        paramsObj.processStageStr = paramsObj.processStageStr + "_topshad"
        if paramsObj.prodsToCalc["RAD"]:
            outputRADImage = os.path.join(paramsObj.outFilePath, paramsObj.outBaseName + paramsObj.processStageStr + "_rad"  + paramsObj.outFormatExt)
            rsgislib.imageutils.maskImage(paramsObj.radianceImage, paramsObj.topoShadowImage, outputRADImage, paramsObj.outFormat, rsgisUtils.getRSGISLibDataTypeFromImg(paramsObj.radianceImage), 0, 1)
            paramsObj.radianceImage = outputRADImage
            paramsObj.sensorClass.setBandNames(paramsObj.radianceImage)
            if paramsObj.calcStatsPy:
                print("Calculating Statistics...")
                rsgislib.imageutils.popImageStats(paramsObj.radianceImage, True, 0.0, True)
        if paramsObj.prodsToCalc["TOA"]:
            outputTOAImage = os.path.join(paramsObj.outFilePath, paramsObj.outBaseName + paramsObj.processStageStr + "_rad_toa" + paramsObj.outFormatExt)
            rsgislib.imageutils.maskImage(paramsObj.toaImage, paramsObj.topoShadowImage, outputTOAImage, paramsObj.outFormat, rsgisUtils.getRSGISLibDataTypeFromImg(paramsObj.toaImage), 0, 1)
            paramsObj.toaImage = outputTOAImage
            paramsObj.sensorClass.setBandNames(paramsObj.toaImage)
            if paramsObj.calcStatsPy:
                print("Calculating Statistics...")
                rsgislib.imageutils.popImageStats(paramsObj.toaImage, True, 0.0, True)

        paramsObj.prodsCalculated["TOPOSHADOW"] = True
        print("")

def performDOS(paramsObj):
    # Convert to an approximation of Surface Reflectance using a dark object subtraction
    if paramsObj.prodsToCalc["DOS"]:
        print("Convert to reflectance using dark object subtraction.")
        outName = paramsObj.outBaseName + paramsObj.processStageStr + "_rad_toa_dos" + paramsObj.outFormatExt
        outWholeName = paramsObj.outBaseName + paramsObj.processStageWholeImgStr + "_rad_toa_dos" + paramsObj.outFormatExt
        paramsObj.srefDOSImage, offVals = paramsObj.sensorClass.convertImageToReflectanceSimpleDarkSubtract(paramsObj.toaImage, paramsObj.outFilePath, outName, paramsObj.outFormat, paramsObj.dosOutRefl)
        paramsObj.calcdOutVals['ARCSI_DOS_OFFSETS'] = offVals
        if paramsObj.fullImgOuts:
            paramsObj.srefDOSWholeImage, offVals = paramsObj.sensorClass.convertImageToReflectanceSimpleDarkSubtract(paramsObj.toaImageWhole, paramsObj.outFilePath, outWholeName, paramsObj.outFormat, paramsObj.dosOutRefl, offVals)

        print("Setting Band Names...")
        paramsObj.sensorClass.setBandNames(paramsObj.srefDOSImage)
        paramsObj.finalOutFiles["SREF_DOS_IMG"] = paramsObj.srefDOSImage
        if paramsObj.fullImgOuts:
            paramsObj.sensorClass.setBandNames(paramsObj.srefDOSWholeImage)
            paramsObj.finalOutFiles["SREF_DOS_IMG_WHOLE"] = paramsObj.srefDOSWholeImage

        if paramsObj.calcStatsPy:
            print("Calculating Statistics...")
            rsgislib.imageutils.popImageStats(paramsObj.srefDOSImage, True, 0.0, True)
            if paramsObj.fullImgOuts:
                rsgislib.imageutils.popImageStats(paramsObj.srefDOSWholeImage, True, 0.0, True)
            print("")
        paramsObj.prodsCalculated["DOS"] = True


def estimateSceneAOT(paramsObj):
    # Use image to estimate AOT values
    if paramsObj.prodsToCalc["DOSAOTSGL"]:
        paramsObj.calcdOutVals['ARCSI_AOT_RANGE_MIN'] = paramsObj.minAOT
        paramsObj.calcdOutVals['ARCSI_AOT_RANGE_MAX'] = paramsObj.maxAOT
        paramsObj.aotVal = paramsObj.sensorClass.estimateSingleAOTFromDOS(paramsObj.radianceImage, paramsObj.toaImage, paramsObj.outDEMName, paramsObj.tmpPath, paramsObj.outBaseName, paramsObj.outFormat, paramsObj.aeroProfile, paramsObj.atmosProfile, paramsObj.grdRefl, paramsObj.minAOT, paramsObj.maxAOT, paramsObj.dosOutRefl)
        paramsObj.minAOT = paramsObj.aotVal - paramsObj.lowAOT
        if paramsObj.minAOT < 0.01:
            paramsObj.minAOT = 0.05
        paramsObj.maxAOT = paramsObj.aotVal + paramsObj.upAOT
        print("AOT Search Range = [" + str(paramsObj.minAOT) + ", " + str(paramsObj.maxAOT) + "]")
        paramsObj.calcdOutVals['ARCSI_AOT_VALUE'] = paramsObj.aotVal
        paramsObj.prodsCalculated["DOSAOTSGL"] = True

    if paramsObj.prodsToCalc["DDVAOT"]:
        outName = paramsObj.outBaseName + "_ddvaot" + paramsObj.outFormatExt
        paramsObj.aotFile = paramsObj.sensorClass.estimateImageToAODUsingDDV(paramsObj.radianceImage, paramsObj.toaImage, paramsObj.outDEMName, paramsObj.topoShadowImage, paramsObj.outFilePath, outName, paramsObj.outFormat, paramsObj.tmpPath, paramsObj.aeroProfile, paramsObj.atmosProfile, paramsObj.grdRefl, paramsObj.minAOT, paramsObj.maxAOT)
        dataset = gdal.Open(paramsObj.aotFile, gdal.GA_Update)
        dataset.GetRasterBand(1).SetDescription("AOT")
        dataset = None
        if paramsObj.calcStatsPy:
            print("Calculating Statistics...")
            rsgislib.imageutils.popImageStats(paramsObj.aotFile, True, 0.0, True)
        paramsObj.finalOutFiles["AOTIMG_DDV"] = paramsObj.aotFile
        paramsObj.calcdOutVals['ARCSI_AOT_RANGE_MIN'] = paramsObj.minAOT
        paramsObj.calcdOutVals['ARCSI_AOT_RANGE_MAX'] = paramsObj.maxAOT
        paramsObj.prodsCalculated["DDVAOT"] = True
        print("")

    if paramsObj.prodsToCalc["DOSAOT"]:
        outName = paramsObj.outBaseName + "_dosaot" + paramsObj.outFormatExt
        paramsObj.aotFile = paramsObj.sensorClass.estimateImageToAODUsingDOS(paramsObj.radianceImage, paramsObj.toaImage, paramsObj.outDEMName, paramsObj.topoShadowImage, paramsObj.outFilePath, outName, paramsObj.outFormat, paramsObj.tmpPath, paramsObj.aeroProfile, paramsObj.atmosProfile, paramsObj.grdRefl, paramsObj.minAOT, paramsObj.maxAOT, paramsObj.globalDOS, paramsObj.simpleDOS, paramsObj.dosOutRefl)
        dataset = gdal.Open(paramsObj.aotFile, gdal.GA_Update)
        dataset.GetRasterBand(1).SetDescription("AOT")
        dataset = None
        if paramsObj.calcStatsPy:
            print("Calculating Statistics...")
            rsgislib.imageutils.popImageStats(paramsObj.aotFile, True, 0.0, True)
        paramsObj.finalOutFiles["AOTIMG_DOS"] = paramsObj.aotFile
        paramsObj.calcdOutVals['ARCSI_AOT_RANGE_MIN'] = paramsObj.minAOT
        paramsObj.calcdOutVals['ARCSI_AOT_RANGE_MAX'] = paramsObj.maxAOT
        paramsObj.prodsCalculated["DOSAOT"] = True
        print("")

def calculateSREF(paramsObj):
    # Convert to Surface Reflectance using 6S Standard Models
    if paramsObj.prodsToCalc["SREF"]:
        arcsiUtils = ARCSIUtils()
        if (paramsObj.prodsToCalc["DDVAOT"] or paramsObj.prodsToCalc["DOSAOT"]) and (paramsObj.aotFile != ""):
            imgDS = gdal.Open(paramsObj.aotFile, gdal.GA_ReadOnly )
            imgBand = imgDS.GetRasterBand(1)
            (min,max,mean,stddev) = imgBand.ComputeStatistics(False)
            print("AOT Mean (Std Dev) = " + str(mean) + " (" + str(stddev) + ")")
            print("AOT [Min, Max] = [" + str(min) + "," + str(max) + "]")
            aotVal = mean
            if aotVal < 0.01:
                print("WARNING: Something has gone wrong as AOT value is 0 or below. Setting to 0.05")
                aotVal = 0.05
            imgDS = None

        if (paramsObj.aotVal is None) and (paramsObj.visVal is None) and (paramsObj.aotFile == ""):
            raise ARCSIException("Either the AOT or the visability need to specified.")
        elif (paramsObj.aotVal is None) and (paramsObj.aotFile == ""):
            print("Convert to vis to aot...")
            paramsObj.aotVal = arcsiUtils.convertVisabilityToAOT(paramsObj.visVal)

        if paramsObj.aotVal is not None:
            print("AOT Value: {}".format(paramsObj.aotVal))
            paramsObj.calcdOutVals['ARCSI_AOT_VALUE'] = paramsObj.aotVal

        if (paramsObj.demFile is None):
            paramsObj.processSREFStr = '_rad_sref'
            outName = paramsObj.outBaseName + paramsObj.processStageStr + paramsObj.processSREFStr + paramsObj.outFormatExt
            paramsObj.srefImage = paramsObj.sensorClass.convertImageToSurfaceReflSglParam(paramsObj.radianceImage, paramsObj.outFilePath, outName, paramsObj.outFormat, paramsObj.aeroProfile, paramsObj.atmosProfile, paramsObj.grdRefl, paramsObj.surfaceAltitude, paramsObj.aotVal, paramsObj.useBRDF, paramsObj.scaleFactor)
            if paramsObj.fullImgOuts:
                outName = paramsObj.outBaseName + paramsObj.processStageWholeImgStr + paramsObj.processSREFStr + paramsObj.outFormatExt
                paramsObj.sref6SWholeImage = paramsObj.sensorClass.convertImageToSurfaceReflSglParam(paramsObj.radianceImageWhole, paramsObj.outFilePath, outName, paramsObj.outFormat, paramsObj.aeroProfile, paramsObj.atmosProfile, paramsObj.grdRefl, paramsObj.surfaceAltitude, paramsObj.aotVal, paramsObj.useBRDF, paramsObj.scaleFactor)
            paramsObj.calcdOutVals['ARCSI_ELEVATION_VALUE'] = paramsObj.surfaceAltitude
        else:
            # Calc Min, Max Elevation for region intersecting with the image.
            statsElev = rsgislib.imagecalc.getImageStatsInEnv(paramsObj.outDEMName, 1, float(paramsObj.demNoDataVal), paramsObj.sensorClass.lonTL, paramsObj.sensorClass.lonBR, paramsObj.sensorClass.latBR, paramsObj.sensorClass.latTL)

            print("Minimum Elevation = ", statsElev[0])
            print("Maximum Elevation = ", statsElev[1])

            paramsObj.minElev = arcsiUtils.findMinimumElev(statsElev[0])
            paramsObj.maxElev = arcsiUtils.findMaximumElev(statsElev[1])

            paramsObj.calcdOutVals['ARCSI_LUT_ELEVATION_MIN'] = paramsObj.minElev
            paramsObj.calcdOutVals['ARCSI_LUT_ELEVATION_MAX'] = paramsObj.maxElev

            elevRange = (paramsObj.maxElev - paramsObj.minElev) / 100
            numElevSteps = math.ceil(elevRange) + 1
            print("Elevation Ranges from ", paramsObj.minElev, " to ", paramsObj.maxElev, " an LUT with ", numElevSteps, " will be created.")

            if (paramsObj.aotFile is None) or (paramsObj.aotFile == ""):
                print("Build an DEM LUT with AOT = " + str(paramsObj.aotVal) + "...")
                paramsObj.processSREFStr = '_rad_srefdem'
                outName = paramsObj.outBaseName + paramsObj.processStageStr + paramsObj.processSREFStr + paramsObj.outFormatExt
                paramsObj.srefImage, paramsObj.sixsLUTCoeffs = paramsObj.sensorClass.convertImageToSurfaceReflDEMElevLUT(paramsObj.radianceImage, paramsObj.outDEMName, paramsObj.outFilePath, outName, paramsObj.outFormat, paramsObj.aeroProfile, paramsObj.atmosProfile, paramsObj.grdRefl, paramsObj.aotVal, paramsObj.useBRDF, paramsObj.minElev, paramsObj.maxElev, paramsObj.scaleFactor)
                if paramsObj.fullImgOuts:
                    outName = paramsObj.outBaseName + paramsObj.processStageWholeImgStr + paramsObj.processSREFStr + paramsObj.outFormatExt
                    paramsObj.sref6SWholeImage, sixsLUTCoeffs = paramsObj.sensorClass.convertImageToSurfaceReflDEMElevLUT(paramsObj.radianceImageWhole, paramsObj.outDEMName, paramsObj.outFilePath, outName, paramsObj.outFormat, paramsObj.aeroProfile, paramsObj.atmosProfile, paramsObj.grdRefl, paramsObj.aotVal, paramsObj.useBRDF, paramsObj.minElev, paramsObj.maxElev, paramsObj.scaleFactor, paramsObj.sixsLUTCoeffs)
                #paramsObj.calcdOutVals['ARCSI_6S_COEFFICENTS'] = paramsObj.sixsLUTCoeffs
                paramsObj.aotLUT = False
            else:
                print("Build an AOT and DEM LUT...")
                statsAOT = rsgislib.imagecalc.getImageStatsInEnv(paramsObj.aotFile, 1, float(-9999), paramsObj.sensorClass.lonTL, paramsObj.sensorClass.lonBR, paramsObj.sensorClass.latBR, paramsObj.sensorClass.latTL)

                paramsObj.minAOT = arcsiUtils.findMinimumAOT(statsAOT[0])
                if paramsObj.minAOT < 0.01:
                    minAOT = 0.05
                paramsObj.maxAOT = arcsiUtils.findMaximumAOT(statsAOT[1])

                paramsObj.calcdOutVals['ARCSI_LUT_AOT_MIN'] = paramsObj.minAOT
                paramsObj.calcdOutVals['ARCSI_LUT_AOT_MAX'] = paramsObj.maxAOT

                aotRange = (paramsObj.maxAOT - paramsObj.minAOT) / 0.05
                numAOTSteps = math.ceil(aotRange) + 1
                print("AOT Ranges from ", paramsObj.minAOT, " to ", paramsObj.maxAOT, " an LUT with ", numAOTSteps, " will be created.")
                paramsObj.processSREFStr = '_rad_srefdemaot'
                outName = paramsObj.outBaseName + paramsObj.processStageStr + paramsObj.processSREFStr + paramsObj.outFormatExt
                paramsObj.srefImage, paramsObj.sixsLUTCoeffs = paramsObj.sensorClass.convertImageToSurfaceReflAOTDEMElevLUT(paramsObj.radianceImage, paramsObj.outDEMName, paramsObj.aotFile, paramsObj.outFilePath, outName, paramsObj.outFormat, paramsObj.aeroProfile, paramsObj.atmosProfile, paramsObj.grdRefl, paramsObj.useBRDF, paramsObj.minElev, paramsObj.maxElev, paramsObj.minAOT, paramsObj.maxAOT, paramsObj.scaleFactor)
                if paramsObj.fullImgOuts:
                    outName = paramsObj.outBaseName + paramsObj.processStageWholeImgStr + paramsObj.processSREFStr + paramsObj.outFormatExt
                    paramsObj.sref6SWholeImage, sixsLUTCoeffs = paramsObj.sensorClass.convertImageToSurfaceReflAOTDEMElevLUT(paramsObj.radianceImageWhole, paramsObj.outDEMName, paramsObj.aotFile, paramsObj.outFilePath, outName, paramsObj.outFormat, paramsObj.aeroProfile, paramsObj.atmosProfile, paramsObj.grdRefl, paramsObj.useBRDF, paramsObj.minElev, paramsObj.maxElev, paramsObj.minAOT, paramsObj.maxAOT, paramsObj.scaleFactor, paramsObj.sixsLUTCoeffs)
                #paramsObj.calcdOutVals['ARCSI_6S_COEFFICENTS'] = paramsObj.sixsLUTCoeffs
                paramsObj.aotLUT = True

        print("Setting Band Names...")
        paramsObj.sensorClass.setBandNames(paramsObj.srefImage)
        if paramsObj.fullImgOuts:
            paramsObj.sensorClass.setBandNames(paramsObj.sref6SWholeImage)
            paramsObj.finalOutFiles["SREF_DOS_IMG_WHOLE"] = paramsObj.sref6SWholeImage

        if paramsObj.calcStatsPy:
            print("Calculating Statistics...")
            rsgislib.imageutils.popImageStats(paramsObj.srefImage, True, 0.0, True)
            if paramsObj.fullImgOuts:
                rsgislib.imageutils.popImageStats(paramsObj.sref6SWholeImage, True, 0.0, True)
        paramsObj.finalOutFiles["SREF_6S_IMG"] = paramsObj.srefImage
        if paramsObj.fullImgOuts:
            paramsObj.finalOutFiles["SREF_6S_WHOLE_IMG"] = paramsObj.sref6SWholeImage
        paramsObj.prodsCalculated["SREF"] = True
        print("")

def calculateStandarisedSREF(paramsObj):
    if paramsObj.prodsToCalc["STDSREF"]:
        if paramsObj.sixsLUTCoeffs is None:
            raise ARCSIException("To calculate standardised reflectance \'STDSREF\' you need to have a LUT within 6S.\nTherefore, a DEM is required and optional an AOT surface should be calculated.")
        if not paramsObj.fullImgOuts:
            paramsObj.sref6SWholeImage = None
        paramsObj.processSREFStr = paramsObj.processSREFStr + '_stdsref'
        outName = paramsObj.outBaseName + paramsObj.processStageStr + paramsObj.processSREFStr + paramsObj.outFormatExt
        outNameWhole = paramsObj.outBaseName + paramsObj.processStageWholeImgStr + paramsObj.processSREFStr + paramsObj.outFormatExt
        paramsObj.stdSREFImg, paramsObj.stdSREFWholeImg = paramsObj.sensorClass.convertSREF2StdisedSREF(paramsObj.srefImage, paramsObj.sref6SWholeImage, paramsObj.outDEMName, paramsObj.topoShadowImage, paramsObj.outFilePath, outName, outNameWhole, paramsObj.outFormat, paramsObj.tmpPath, paramsObj.sixsLUTCoeffs, paramsObj.aotLUT, paramsObj.scaleFactor, brdfBeta=1, outIncidenceAngle=0, outExitanceAngle=0)

        print("Setting Band Names...")
        paramsObj.sensorClass.setBandNames(paramsObj.stdSREFImg)
        if paramsObj.fullImgOuts:
            paramsObj.sensorClass.setBandNames(paramsObj.stdSREFWholeImg)

        if paramsObj.calcStatsPy:
            print("Calculating Statistics...")
            rsgislib.imageutils.popImageStats(paramsObj.stdSREFImg, True, 0.0, True)
            if paramsObj.fullImgOuts:
                rsgislib.imageutils.popImageStats(paramsObj.stdSREFWholeImg, True, 0.0, True)
        paramsObj.finalOutFiles["STD_SREF_IMG"] = paramsObj.stdSREFImg
        if paramsObj.fullImgOuts:
            paramsObj.finalOutFiles["STD_SREF_WHOLE_IMG"] = paramsObj.stdSREFWholeImg
        paramsObj.prodsCalculated["STDSREF"] = True
        print("")

def exportMetaData(paramsObj):
    if paramsObj.prodsToCalc["METADATA"]:
        print("Exporting Meta-data file")
        outName = paramsObj.outBaseName + "_meta.json"
        paramsObj.finalOutFiles["METADATA"] = outName

        validMaskImagePath = ""
        if paramsObj.validMaskImage is not None:
            if paramsObj.reproject:
                validMaskImagePath = paramsObj.validMaskImageProj
            else:
                validMaskImagePath = paramsObj.validMaskImage

        paramsObj.sensorClass.generateMetaDataFile(paramsObj.outFilePath, outName, paramsObj.productsStr, validMaskImagePath, paramsObj.prodsToCalc["FOOTPRINT"], paramsObj.calcdOutVals, paramsObj.finalOutFiles)
        paramsObj.prodsCalculated["METADATA"] = True
        print("")

def runARCSI(inputHeader, inputImage, cloudMaskUsrImg, sensorStr, inWKTFile, outFormat, outFilePath, outBaseName, outWKTFile, outProj4File, projAbbv, xPxlResUsr, yPxlResUsr, productsStr, calcStatsPy, aeroProfileOption, atmosProfileOption, aeroProfileOptionImg, atmosProfileOptionImg,  grdReflOption, surfaceAltitude, atmosOZoneVal,atmosWaterVal, atmosOZoneWaterSpecified, aeroWaterVal, aeroDustVal, aeroOceanicVal, aeroSootVal, aeroComponentsSpecified, aotVal, visVal, tmpPath, minAOT, maxAOT, lowAOT, upAOT, demFile, demNoDataUsrVal, aotFile, globalDOS, dosOutRefl, simpleDOS, debugMode, scaleFactor, interpAlgor, interpAlgorResample, initClearSkyRegionDist, initClearSkyRegionMinSize, finalClearSkyRegionDist, clearSkyMorphSize, fullImgOuts, checkOutputs, classmlclouds, cloudtrainclouds, cloudtrainother, resample2LowResImg, fileEnding2Keep, cloud_methods):
    """
    A function contains the main flow of the software
    """
    try:
        # Initialise and parameters object.
        paramsObj = None
        paramsObj = prepParametersObj(inputHeader, inputImage, cloudMaskUsrImg, sensorStr, inWKTFile, outFormat, outFilePath, outBaseName, outWKTFile, outProj4File, projAbbv, xPxlResUsr, yPxlResUsr, productsStr, calcStatsPy, aeroProfileOption, atmosProfileOption, aeroProfileOptionImg, atmosProfileOptionImg,  grdReflOption, surfaceAltitude, atmosOZoneVal,atmosWaterVal, atmosOZoneWaterSpecified, aeroWaterVal, aeroDustVal, aeroOceanicVal, aeroSootVal, aeroComponentsSpecified, aotVal, visVal, tmpPath, minAOT, maxAOT, lowAOT, upAOT, demFile, demNoDataUsrVal, aotFile, globalDOS, dosOutRefl, simpleDOS, debugMode, scaleFactor, interpAlgor, interpAlgorResample, initClearSkyRegionDist, initClearSkyRegionMinSize, finalClearSkyRegionDist, clearSkyMorphSize, fullImgOuts, checkOutputs, classmlclouds, cloudtrainclouds, cloudtrainother, resample2LowResImg, fileEnding2Keep, cloud_methods)

        # Check Input image(s) is valid before proceeding.
        checkForValidInput(paramsObj)

        # Check if bands need resampling
        resampleBands(paramsObj)

        # Check if the image data needs mosaicking.
        mosaicInputImages(paramsObj)

        # Create valid image area mask and view angle images
        createValidMaskViewAngle(paramsObj)

        # Create Vector Footprint
        createFootprint(paramsObj)

        # Create Saturated image
        createSaturatedImage(paramsObj)

        # Convert imagery to radiance
        convertInputImageToRadiance(paramsObj)

        # Calculate Thermal Brightness
        calcThermalBrightness(paramsObj)

        # Calculate TOA Reflectance
        calcTOAReflectance(paramsObj)

        # Save the process stage string for using with whole image outputs.
        paramsObj.processStageWholeImgStr = paramsObj.processStageStr

        # Perform a cloud masking
        performCloudMasking(paramsObj)
        
        # Don't continue further if there is more than 95% cloud cover in the scene.
        if  (not paramsObj.prodsToCalc["CLOUDS"]) or (paramsObj.prodsToCalc["CLOUDS"] and paramsObj.propOfCloud < 0.95):
            # Perform clear sky masking
            performClearSkyMasking(paramsObj)
            
            # Don't continue further if there is less than 5% of the scene of clear sky
            if (not paramsObj.prodsToCalc["CLEARSKY"]) or (paramsObj.prodsToCalc["CLEARSKY"] and paramsObj.propOfClearSky > 0.05):
                # Prepare the DEM for later processing stages.
                prepareDEM(paramsObj)

                # Calculate Topographic shadow mask
                calcTopoShadowMask(paramsObj)

                # Perfrom Dark Object Subtraction (DOS)
                performDOS(paramsObj)

                # Estimate AOT for the scene
                estimateSceneAOT(paramsObj)

                # Calculate SREF
                calculateSREF(paramsObj)

                # Calculate Standarised SREF
                calculateStandarisedSREF(paramsObj)

            else:
                keys2Del = []
                for key in paramsObj.prodsToCalc.keys():
                    if paramsObj.prodsToCalc[key] is not paramsObj.prodsCalculated[key]:
                        if key in ['DDVAOT', 'DOSAOT', 'DOSAOTSGL', 'SREF', 'DOS', 'STDSREF', 'TOPOSHADOW']:
                            keys2Del.append(key)
                for key in keys2Del:
                    del paramsObj.prodsToCalc[key]
        else:
            keys2Del = []
            for key in paramsObj.prodsToCalc.keys():
                if paramsObj.prodsToCalc[key] is not paramsObj.prodsCalculated[key]:
                    if key in ['CLEARSKY', 'DDVAOT', 'DOSAOT', 'DOSAOTSGL', 'SREF', 'DOS', 'STDSREF', 'TOPOSHADOW']:
                        keys2Del.append(key)
            for key in keys2Del:
                del paramsObj.prodsToCalc[key]

        # Export metadata output file.
        exportMetaData(paramsObj)

        print('Clean up anything left over...')
        paramsObj.sensorClass.cleanFollowProcessing(paramsObj.outFilePath, paramsObj.fileEnding2Keep)
        
    except ARCSIException as e:
        print('Input Header: \'' + inputHeader + '\'', file=sys.stderr)
        if (paramsObj is not None) and (paramsObj.outBaseName is not None):
            print('Output Basename: \'' + paramsObj.outBaseName + '\'', file=sys.stderr)
        print("Error: {}".format(e), file=sys.stderr)
        if debugMode:
            raise
    except Exception as e:
        print('Input Header: \'' + inputHeader + '\'', file=sys.stderr)
        if (paramsObj is not None) and (paramsObj.outBaseName is not None):
            print('Output Basename: \'' + paramsObj.outBaseName + '\'', file=sys.stderr)
        print("Error: {}".format(e), file=sys.stderr)
        if debugMode:
            raise
    finally:
        if paramsObj is not None:
            failedProdsList = []
            # Check all requested products have been created
            for key in paramsObj.prodsToCalc.keys():
                if paramsObj.prodsToCalc[key] is not paramsObj.prodsCalculated[key]:
                    failedProdsList.append(key)
            if len(failedProdsList) > 0:
                print("Error: The following products were not generated:", file=sys.stderr)
                print(" ".join(failedProdsList), file=sys.stderr)
                print('Input Header: \'' + inputHeader + '\'', file=sys.stderr)
                if paramsObj.outBaseName is not None:
                    print('Output Basename: \'' + paramsObj.outBaseName + '\'', file=sys.stderr)

def _runARCSIPart1(paramsObj):
    try:
         # Check Input image(s) is valid before proceeding.
        checkForValidInput(paramsObj)

        # Check if bands need resampling
        resampleBands(paramsObj)

        # Check if the image data needs mosaicking.
        mosaicInputImages(paramsObj)

        # Create valid image area mask and view angle images
        createValidMaskViewAngle(paramsObj)

        # Create Vector Footprint
        createFootprint(paramsObj)

        # Create Saturated image
        createSaturatedImage(paramsObj)

        # Convert imagery to radiance
        convertInputImageToRadiance(paramsObj)

        # Calculate Thermal Brightness
        calcThermalBrightness(paramsObj)

        # Calculate TOA Reflectance
        calcTOAReflectance(paramsObj)

        # Save the process stage string for using with whole image outputs.
        paramsObj.processStageWholeImgStr = paramsObj.processStageStr

        # Perform a cloud masking
        performCloudMasking(paramsObj)
        
        # Don't continue further if there is more than 95% cloud cover in the scene.
        if  (not paramsObj.prodsToCalc["CLOUDS"]) or (paramsObj.prodsToCalc["CLOUDS"] and paramsObj.propOfCloud < 0.95):
            # Perform clear sky masking
            performClearSkyMasking(paramsObj)
            
            # Don't continue further if there is less than 5% of the scene of clear sky
            if (not paramsObj.prodsToCalc["CLEARSKY"]) or (paramsObj.prodsToCalc["CLEARSKY"] and paramsObj.propOfClearSky > 0.05):
                # Prepare the DEM for later processing stages.
                prepareDEM(paramsObj)

                # Calculate Topographic shadow mask
                calcTopoShadowMask(paramsObj)

                # Perfrom Dark Object Subtraction (DOS)
                performDOS(paramsObj)

                # Estimate AOT for the scene
                estimateSceneAOT(paramsObj)
            else:
                keys2Del = []
                for key in paramsObj.prodsToCalc.keys():
                    if paramsObj.prodsToCalc[key] is not paramsObj.prodsCalculated[key]:
                        if key in ['DDVAOT', 'DOSAOT', 'DOSAOTSGL', 'SREF', 'DOS', 'STDSREF', 'TOPOSHADOW']:
                            keys2Del.append(key)
                for key in keys2Del:
                    del paramsObj.prodsToCalc[key]
        else:
            keys2Del = []
            for key in paramsObj.prodsToCalc.keys():
                if paramsObj.prodsToCalc[key] is not paramsObj.prodsCalculated[key]:
                    if key in ['CLEARSKY', 'DDVAOT', 'DOSAOT', 'DOSAOTSGL', 'SREF', 'DOS', 'STDSREF', 'TOPOSHADOW']:
                        keys2Del.append(key)
            for key in keys2Del:
                del paramsObj.prodsToCalc[key]
        
    except ARCSIException as e:
        print("Error: {}".format(e), file=sys.stderr)
    except Exception as e:
        print("Error: {}".format(e), file=sys.stderr)

    return paramsObj

def _runARCSIPart2(paramsObj):
    try:
        # Don't continue further if there is more than 95% cloud cover in the scene.
        if  (not paramsObj.prodsToCalc["CLOUDS"]) or (paramsObj.prodsToCalc["CLOUDS"] and paramsObj.propOfCloud < 0.95):
            # Don't continue further if there is less than 5% of the scene of clear sky
            if (not paramsObj.prodsToCalc["CLEARSKY"]) or (paramsObj.prodsToCalc["CLEARSKY"] and paramsObj.propOfClearSky > 0.05):

                # Calculate SREF
                calculateSREF(paramsObj)

                # Calculate Standarised SREF
                calculateStandarisedSREF(paramsObj)

    except ARCSIException as e:
        print("Error: {}".format(e), file=sys.stderr)
    except Exception as e:
        print("Error: {}".format(e), file=sys.stderr)
    return paramsObj

def _runARCSIPart3(paramsObj):
    try:
        exportMetaData(paramsObj)
    except ARCSIException as e:
        print("Error: {}".format(e), file=sys.stderr)
    except Exception as e:
        print("Error: {}".format(e), file=sys.stderr)
    return paramsObj

def _runARCSIPart4(paramsObj):
    try:
        print('Clean up anything left over...')
        paramsObj.sensorClass.cleanFollowProcessing(paramsObj.outFilePath, paramsObj.fileEnding2Keep)
    except ARCSIException as e:
        print("Error: {}".format(e), file=sys.stderr)
    except Exception as e:
        print("Error: {}".format(e), file=sys.stderr)
    return paramsObj

def runARCSIMulti(inputHeaders, sensorStr, inWKTFile, outFormat, outFilePath, outBaseName, outWKTFile, outProj4File, projAbbv, xPxlResUsr, yPxlResUsr, productsStr, calcStatsPy, aeroProfileOption, atmosProfileOption, aeroProfileOptionImg, atmosProfileOptionImg,  grdReflOption, surfaceAltitude, atmosOZoneVal,atmosWaterVal, atmosOZoneWaterSpecified, aeroWaterVal, aeroDustVal, aeroOceanicVal, aeroSootVal, aeroComponentsSpecified, aotVal, visVal, tmpPath, minAOT, maxAOT, lowAOT, upAOT, demFile, demNoDataUsrVal, aotFile, globalDOS, dosOutRefl, simpleDOS, debugMode, scaleFactor, interpAlgor, interpAlgorResample, initClearSkyRegionDist, initClearSkyRegionMinSize, finalClearSkyRegionDist, clearSkyMorphSize, fullImgOuts, checkOutputs, classmlclouds, cloudtrainclouds, cloudtrainother, resample2LowResImg, ncores, fileEnding2Keep, cloud_methods):
    """
    A function contains the main flow of the software
    """
    try:
        rsgisUtils = rsgislib.RSGISPyUtils()
        inputHeadersLst = rsgisUtils.readTextFile2List(inputHeaders)
        paramsLst = []
        calc6SSREF = False
        exportMetaData = False
        calcAOT = False
        useAOTImage = False
        first = True
        for inputHeader in inputHeadersLst:
            print(inputHeader)
            # Initialise and parameters object.
            paramsObj = None
            paramsObj = prepParametersObj(inputHeader, None, None, sensorStr, inWKTFile, outFormat, outFilePath, outBaseName, outWKTFile, outProj4File, projAbbv, xPxlResUsr, yPxlResUsr, productsStr, calcStatsPy, aeroProfileOption, atmosProfileOption, aeroProfileOptionImg, atmosProfileOptionImg,  grdReflOption, surfaceAltitude, atmosOZoneVal,atmosWaterVal, atmosOZoneWaterSpecified, aeroWaterVal, aeroDustVal, aeroOceanicVal, aeroSootVal, aeroComponentsSpecified, aotVal, visVal, tmpPath, minAOT, maxAOT, lowAOT, upAOT, demFile, demNoDataUsrVal, aotFile, globalDOS, dosOutRefl, simpleDOS, debugMode, scaleFactor, interpAlgor, interpAlgorResample, initClearSkyRegionDist, initClearSkyRegionMinSize, finalClearSkyRegionDist, clearSkyMorphSize, fullImgOuts, checkOutputs, classmlclouds, cloudtrainclouds, cloudtrainother, resample2LowResImg, fileEnding2Keep, cloud_methods)
            paramsLst.append(paramsObj)
            if first:
                if paramsObj.prodsToCalc["DDVAOT"] or paramsObj.prodsToCalc["DOSAOT"] or paramsObj.prodsToCalc["DOSAOTSGL"]:
                    calcAOT = True
                    if paramsObj.prodsToCalc["DDVAOT"] or paramsObj.prodsToCalc["DOSAOT"]:
                        useAOTImage = True
                if paramsObj.prodsToCalc["SREF"]:
                    calc6SSREF = True
                if paramsObj.prodsToCalc["METADATA"]:
                    exportMetaData = True
                first = False

        plObj = Pool(ncores)
        paramsLst = plObj.map(_runARCSIPart1, paramsLst)
        
        if calcAOT:
            if useAOTImage:
                raise ARCSIException("Currently the --multi option does not support the merging of AOT images (i.e., from DDVAOT and DOSAOT) across multiple scenes.")
            else:
                # Replace the AOT value with the mean from all the scenes.
                aotSum = 0.0
                aotN = 0.0
                for paramsObj in paramsLst:
                    if paramsObj.aotVal is not None:
                        aotSum = aotSum + paramsObj.aotVal
                        aotN = aotN + 1
                if aotN > 0:
                    avgAOT = aotSum / aotN
                else:
                    avgAOT = 0.05
                for params in paramsLst:
                    paramsObj.aotVal = avgAOT
        
        if calc6SSREF:
            paramsLst = plObj.map(_runARCSIPart2, paramsLst)

        if exportMetaData:
            paramsLst = plObj.map(_runARCSIPart3, paramsLst)

        paramsLst = plObj.map(_runARCSIPart4, paramsLst)

    except ARCSIException as e:
        print("Error: {}".format(e), file=sys.stderr)
        if debugMode:
            raise
    except Exception as e:
        print("Error: {}".format(e), file=sys.stderr)
        if debugMode:
            raise

def print2ConsoleListSensors():
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

def print2ConsoleListProductDescription(product):
    """
    A function which lists the currently supported products
    and describes what that are and the parameters they require.
    """
    print("Hello World. this has not be written yet!")

def print2ConsoleListEnvVars():
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
    print("                       by a scale factor. So, for a scale factor")
    print("                       of 1000 a value of 20 is 2 % reflectance")
    print("ARCSI_USE_LOCALDOS     in place of the --localdos (variable ")
    print("                       values can be either `TRUE' or `FALSE') option")
    print("ARCSI_USE_SIMPLEDOS    in place of the --simpledos (variable ")
    print("                       values can be either `TRUE' or `FALSE') option")
    print("ARCSI_SCALE_FACTOR     in place of the --scalefac option")
    print("")
