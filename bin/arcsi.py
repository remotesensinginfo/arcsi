#! /usr/bin/env python

"""
Module that contains the ARSCI Main class.
"""

############################################################################
#  arcsi.py
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
# Purpose:  The 'main' interface for the ARCSI system controlling the
#           follow of the program and user interface.
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

# Import the print function (for Python 2)
from __future__ import print_function
# Import the system library
import sys
# Import the subprocess module
import subprocess
# Import the OS python module
import os
# Import the python Argument parser
import argparse
# Import the time module
import time
# Import the copy module
import copy
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

class ARCSI (object):
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
            
    
    def run(self, inputHeader, inputImage, sensorStr, inWKTFile, outFormat, outFilePath, outBaseName, 
            outWKTFile, projAbbv, productsStr, calcStatsPy, aeroProfileOption, atmosProfileOption, 
            aeroProfileOptionImg, atmosProfileOptionImg,  grdReflOption, surfaceAltitude, atmosOZoneVal,  
            atmosWaterVal, atmosOZoneWaterSpecified, aeroWaterVal, aeroDustVal, aeroOceanicVal,  
            aeroSootVal, aeroComponentsSpecified, aotVal, visVal, tmpPath, minAOT, maxAOT, lowAOT, upAOT, 
            demFile, aotFile, globalDOS, dosOutRefl, simpleDOS, debugMode, scaleFactor, interpAlgor):
        """
        A function contains the main flow of the software
        """
        
        startTime = time.time()
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
            if outBaseName is None:
                outBaseName = sensorClass.generateOutputBaseName()
                if not projAbbv is None:
                    outBaseNameProj = outBaseName + "_" + str(projAbbv)
            print("Image Base Name: " + outBaseName + "\n")
            
            # Step 4: Find the products which are to be generated.
            prodsToCalc["RAD"] = False
            prodsToCalc["TOA"] = False
            prodsToCalc["CLOUDS"] = False
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
            
            validMaskImage=None
            validMaskImageProj=""
            radianceImage=""
            saturateImage=""
            saturateImageProj=""
            thermalRadImage=""
            thermalBrightImage=""
            maskImage=""
            toaImage=""
            srefImage=""
            aotFile=""
            cloudsImage=""
            outDEMName=""
            topoShadowImage=""
            footprintShpFile=""
            metaDataFile=""
            
            # Check Input image(s) is valid before proceeding. 
            print('Checking Input Images are valid')
            sensorClass.checkInputImageValid()
            
            if sensorClass.imgNeedMosaicking():
            	print("Mosacking Input Image Tiles.")
            	sensorClass.mosaicImageTiles()
            
            # Get the valid image data maskImage
            outName = outBaseName + "_valid" + arcsiUtils.getFileExtension(outFormat)
            validMaskImage = sensorClass.generateValidImageDataMask(outFilePath, outName, "KEA")
            if not validMaskImage is None:
                rsgislib.rastergis.populateStats(validMaskImage, True, True)
            print("")
            
            if (not outWKTFile is None) and (not validMaskImage is None):
                outName = outBaseNameProj + "_valid" + arcsiUtils.getFileExtension(outFormat)
                validMaskImageProj = os.path.join(outFilePath, outName)
                validImgDS = gdal.Open(validMaskImage, gdal.GA_ReadOnly)
                if validImgDS is None:
                    raise ARCSIException('Could not open the Valid Image Mask ' + validMaskImage)
                geoTransform = validImgDS.GetGeoTransform()
                if geoTransform is None:
                    raise ARCSIException('Could read the geotransform from the Valid Image Mask ' + validMaskImage)
                xPxlRes = geoTransform[1]
                yPxlRes = geoTransform[5]
                validImgDS = None
                cmd = 'gdalwarp -t_srs ' + outWKTFile + ' -tr ' + str(xPxlRes) + ' ' + str(yPxlRes) + ' -ot Byte -wt Float32 ' \
                    + '-r near -tap -srcnodata 0 -dstnodata 0 -of ' + outFormat + ' -overwrite ' \
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
                print("")
                
            if prodsToCalc["FOOTPRINT"]:
                if validMaskImage is None:
                    raise ARCSIException("To generate a footprint a valid image mask is required - not supported by this sensor?")
                outFootprintLyrName = ''
                if not outWKTFile is None:
                    outFootprintLyrName = outBaseNameProj + "_footprint"
                    footprintShpFile = sensorClass.generateImageFootprint(validMaskImageProj, outFilePath, outFootprintLyrName)
                else:
                    outFootprintLyrName = outBaseName + "_footprint"
                    footprintShpFile = sensorClass.generateImageFootprint(validMaskImage, outFilePath, outFootprintLyrName)
                prodsCalculated["FOOTPRINT"] = True
                print("")
            
            if prodsToCalc["SATURATE"]:
                # Execute generation of the saturation image
                outName = outBaseName + "_sat" + arcsiUtils.getFileExtension(outFormat)
                saturateImage = sensorClass.generateImageSaturationMask(outFilePath, outName, outFormat)
                
                if not outWKTFile is None:
                    outNameProj = outBaseNameProj + "_sat" + arcsiUtils.getFileExtension(outFormat)
                    saturateImageProj = os.path.join(outFilePath, outNameProj)
                    cmd = 'gdalwarp -t_srs ' + outWKTFile + ' -tr ' + str(xPxlRes) + ' ' + str(yPxlRes) + ' -ot Byte ' \
                    + '-r near -tap -of ' + outFormat + ' -overwrite ' \
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
                    outImgName = outBaseName + "_rad_msk" + arcsiUtils.getFileExtension(outFormat)
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
                    outRadPathName = os.path.join(outFilePath, outBaseName + "_rad_vmsk" + arcsiUtils.getFileExtension(outFormat))
                    if sensorClass.maskInputImages() & (not maskImage is None):
                        outRadPathName = os.path.join(outFilePath, outBaseName + "_rad_msk_vmsk" + arcsiUtils.getFileExtension(outFormat))
                    rsgislib.imageutils.maskImage(radianceImage, validMaskImage, outRadPathName, outFormat, rsgisUtils.getRSGISLibDataTypeFromImg(radianceImage), 0.0, 0.0)
                    rsgisUtils.deleteFileWithBasename(radianceImage)
                    radianceImage = outRadPathName
                    if not thermalRadImage == None:
                        outThermPathName = os.path.join(outFilePath, outBaseName + "_therm_vmsk" + arcsiUtils.getFileExtension(outFormat))
                        rsgislib.imageutils.maskImage(thermalRadImage, validMaskImage, outThermPathName, outFormat, rsgisUtils.getRSGISLibDataTypeFromImg(thermalRadImage), 0.0, 0.0)
                        rsgisUtils.deleteFileWithBasename(thermalRadImage)
                        thermalRadImage = outThermPathName
                    if not outWKTFile is None:
                        rsgisUtils.deleteFileWithBasename(validMaskImage)
                        validMaskImage = validMaskImageProj
                
                if not outWKTFile is None:
                    if not radianceImage is None:
                        outName = outBaseNameProj + "_rad" + arcsiUtils.getFileExtension(outFormat)
                        if not validMaskImage is None:
                                outName = outBaseNameProj + "_rad_vmsk" + arcsiUtils.getFileExtension(outFormat)
                        if sensorClass.maskInputImages() & (not maskImage is None):
                            outName = outBaseNameProj + "_rad_msk" + arcsiUtils.getFileExtension(outFormat)
                            if not validMaskImage is None:
                                outName = outBaseNameProj + "_rad_msk_vmsk" + arcsiUtils.getFileExtension(outFormat)
                        
                        outRadImagePath = os.path.join(outFilePath, outName)
                        cmd = 'gdalwarp -t_srs ' + outWKTFile + ' -tr ' + str(xPxlRes) + ' ' + str(yPxlRes) + ' -ot Float32 -wt Float32 ' \
                        + '-r ' + interpAlgor + ' -tap -srcnodata 0 -dstnodata 0 -of ' + outFormat + ' -overwrite ' \
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
                        outName = outBaseNameProj + "_therm_rad" + arcsiUtils.getFileExtension(outFormat)
                        outThermRadImagePath = os.path.join(outFilePath, outName)
                        cmd = 'gdalwarp -t_srs ' + outWKTFile + ' -tr ' + str(xPxlRes) + ' ' + str(yPxlRes) + ' -ot Float32 -wt Float32 ' \
                        + '-r ' + interpAlgor + ' -srcnodata 0 -dstnodata 0 -of ' + outFormat + ' -overwrite ' \
                        + thermalRadImage + ' ' + outThermRadImagePath 
                        #print(cmd)
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
                prodsCalculated["RAD"] = True
                print("")
                        
            # Create and reproject the DEM file if not specified.
            if (not (demFile == None)) and (outDEMName == ""):
                # Interpolate image to output refl image resolution and convert projection to the same as the output image.
                outDEMName = os.path.join(outFilePath, (outBaseName + "_dem" + arcsiUtils.getFileExtension(outFormat)))
                print("Output DEM: ", outDEMName)
                rsgislib.imageutils.createCopyImage(radianceImage, outDEMName, 1, -32768.0, outFormat, rsgislib.TYPE_32FLOAT)  
                
                inDEMDS = gdal.Open(demFile, gdal.GA_ReadOnly)
                outDEMDS = gdal.Open(outDEMName, gdal.GA_Update)

                print("Subset and reproject DEM...")
                gdal.ReprojectImage(inDEMDS, outDEMDS, None, None, gdal.GRA_CubicSpline)
                
                inDEMDS = None
                outDEMDS = None
                if calcStatsPy:
                    print("Calculating Statistics...")
                    rsgislib.imageutils.popImageStats(outDEMName, True, -32768.0, True)  
                print("")
                            
            if prodsToCalc["TOPOSHADOW"]:
                # Execute generation of the saturation image
                outName = outBaseName + "_toposhad" + arcsiUtils.getFileExtension(outFormat)
                topoShadowImage = sensorClass.generateTopoDirectShadowMask(outDEMName, outFilePath, outName, outFormat, tmpPath)
                if calcStatsPy:
                    print("Calculating Statistics...")
                    rsgislib.rastergis.populateStats(topoShadowImage, True, True)
                prodsCalculated["TOPOSHADOW"] = True
                print("")
                        
            if prodsToCalc["THERMAL"]:
                # Execute calibrate thermal to brightness
                outName = outBaseName + "_thermal" + arcsiUtils.getFileExtension(outFormat)
                thermalBrightImage = sensorClass.convertThermalToBrightness(thermalRadImage, outFilePath, outName, outFormat, scaleFactor)
                if calcStatsPy:
                    print("Calculating Statistics...")
                    rsgislib.imageutils.popImageStats(thermalBrightImage, True, 0.0, True)
                prodsCalculated["THERMAL"] = True
                print("")
            
            # Step 6: Convert to TOA
            if prodsToCalc["TOA"]:
                # Execute conversion to top of atmosphere reflectance
                outName = outBaseName + "_rad_toa" + arcsiUtils.getFileExtension(outFormat)
                toaImage = sensorClass.convertImageToTOARefl(radianceImage, outFilePath, outName, outFormat, scaleFactor)
                print("Setting Band Names...")
                sensorClass.setBandNames(toaImage)
                if calcStatsPy:
                    print("Calculating Statistics...")
                    rsgislib.imageutils.popImageStats(toaImage, True, 0.0, True)
                prodsCalculated["TOA"] = True
                print("")
            
            # Step 7: Generate Cloud Masks
            if prodsToCalc["CLOUDS"]:
                # Execute conversion to top of atmosphere reflectance
                outName = outBaseName + "_clouds" + arcsiUtils.getFileExtension(outFormat)
                cloudsImage = sensorClass.generateCloudMask(toaImage, saturateImage, thermalBrightImage, validMaskImage, outFilePath, outName, outFormat, tmpPath, scaleFactor)
                print("Setting Band Names...")
                sensorClass.setBandNames(toaImage)
                if calcStatsPy:
                    print("Calculating Statistics...")
                    rsgislib.rastergis.populateStats(cloudsImage, True, True)
                print("Applying cloud masks to images...")
                outputRADImage = os.path.join(outFilePath, outBaseName + "_rad_mclds" + arcsiUtils.getFileExtension(outFormat))
                rsgislib.imageutils.maskImage(radianceImage, cloudsImage, outputRADImage, outFormat, rsgislib.TYPE_32FLOAT, 0, 2)
                radianceImage = outputRADImage
                sensorClass.setBandNames(radianceImage)
                if calcStatsPy:
                    print("Calculating Statistics...")
                    rsgislib.imageutils.popImageStats(radianceImage, True, 0.0, True)
                outputTOAImage = os.path.join(outFilePath, outBaseName + "_rad_toa_mclds" + arcsiUtils.getFileExtension(outFormat))
                rsgislib.imageutils.maskImage(toaImage, cloudsImage, outputTOAImage, outFormat, rsgislib.TYPE_16UINT, 0, 2)
                toaImage = outputTOAImage
                sensorClass.setBandNames(toaImage)
                if calcStatsPy:
                    print("Calculating Statistics...")
                    rsgislib.imageutils.popImageStats(toaImage, True, 0.0, True)
                prodsCalculated["CLOUDS"] = True
                print("")
            
            # Step 8: Convert to an approximation of Surface Reflectance using a dark object subtraction
            if prodsToCalc["DOS"]:
                print("Convert to reflectance using dark object subtraction.")
                outName = outBaseName + "_rad_toa_dos" + arcsiUtils.getFileExtension(outFormat)  
                if simpleDOS:
                    srefImage = sensorClass.convertImageToReflectanceSimpleDarkSubtract(toaImage, outFilePath, outName, outFormat, dosOutRefl)
                else:        
                    srefImage = sensorClass.convertImageToReflectanceDarkSubstract(toaImage, outFilePath, outName, outFormat, tmpPath, globalDOS, dosOutRefl)
            
                print("Setting Band Names...")
                sensorClass.setBandNames(srefImage)

                if calcStatsPy:
                    print("Calculating Statistics...")
                    rsgislib.imageutils.popImageStats(srefImage, True, 0.0, True)
                    print("")
                prodsCalculated["DOS"] = True
            
            # Step 9: Use image to estimate AOD values
            if prodsToCalc["DOSAOTSGL"]:
                aotVal = sensorClass.estimateSingleAOTFromDOS(radianceImage, toaImage, outDEMName, tmpPath, outBaseName, outFormat, aeroProfile, atmosProfile, grdRefl, minAOT, maxAOT, dosOutRefl)
                minAOT = aotVal - lowAOT
                if minAOT < 0.01:
                    minAOT = 0.05
                maxAOT = aotVal + upAOT
                print("AOT Search Range = [" + str(minAOT) + ", " + str(maxAOT) + "]")
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
                prodsCalculated["DOSAOT"] = True
                print("")
                        
            # Step 10: Convert to Surface Reflectance using 6S Standard Models
            if prodsToCalc["SREF"]:
                # Execute conversion to surface reflectance by applying 6S using a 'standard' modelled atmosphere.
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
                
                if (demFile == None):
                    outName = outBaseName + "_rad_sref" + arcsiUtils.getFileExtension(outFormat)
                    srefImage = sensorClass.convertImageToSurfaceReflSglParam(radianceImage, outFilePath, outName, outFormat, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF, scaleFactor)
                else:
                    # Calc Min, Max Elevation for region intersecting with the image.
                    statsElev = rsgislib.imagecalc.getImageStatsInEnv(demFile, 1, -32768.0, sensorClass.latTL, sensorClass.latBR, sensorClass.lonBR, sensorClass.lonTL)
                    
                    print("Minimum Elevation = ", statsElev[0])
                    print("Maximum Elevation = ", statsElev[1])
                    
                    minElev = self.findMinimumElev(statsElev[0])
                    maxElev = self.findMaximumElev(statsElev[1])
                    
                    elevRange = (maxElev - minElev) / 100
                    numElevSteps = math.ceil(elevRange) + 1
                    print("Elevation Ranges from ", minElev, " to ", maxElev, " an LUT with ", numElevSteps, " will be created.")
                        
                    if (aotFile == None) or (aotFile == ""):
                        print("Build an DEM LUT with AOT == " + str(aotVal) + "...")
                        outName = outBaseName + "_rad_srefdem" + arcsiUtils.getFileExtension(outFormat)
                        srefImage = sensorClass.convertImageToSurfaceReflDEMElevLUT(radianceImage, outDEMName, outFilePath, outName, outFormat, aeroProfile, atmosProfile, grdRefl, aotVal, useBRDF, minElev, maxElev, scaleFactor)
                    else:
                        print("Build an AOT and DEM LUT...")
                        statsAOT = rsgislib.imagecalc.getImageStatsInEnv(aotFile, 1, -9999, sensorClass.latTL, sensorClass.latBR, sensorClass.lonBR, sensorClass.lonTL)
                        
                        minAOT = self.findMinimumAOT(statsAOT[0])
                        if minAOT < 0.01:
                            minAOT = 0.05
                        maxAOT = self.findMaximumAOT(statsAOT[1])
                        
                        aotRange = (maxAOT - minAOT) / 0.05
                        numAOTSteps = math.ceil(aotRange) + 1
                        print("AOT Ranges from ", minAOT, " to ", maxAOT, " an LUT with ", numAOTSteps, " will be created.")
                        outName = outBaseName + "_rad_srefdemaot" + arcsiUtils.getFileExtension(outFormat)
                        srefImage = sensorClass.convertImageToSurfaceReflAOTDEMElevLUT(radianceImage, outDEMName, aotFile, outFilePath, outName, outFormat, aeroProfile, atmosProfile, grdRefl, useBRDF, minElev, maxElev, minAOT, maxAOT, scaleFactor)

                print("Setting Band Names...")
                sensorClass.setBandNames(srefImage)

                if calcStatsPy:
                    print("Calculating Statistics...")
                    rsgislib.imageutils.popImageStats(srefImage, True, 0.0, True)
                prodsCalculated["SREF"] = True
                print("")
                    
            if prodsToCalc["METADATA"]:
                print("Exporting Meta-data file")
                outName = outBaseName + "_meta.json"
                
                validMaskImagePath = ""
                if not validMaskImage is None:
                    if not outWKTFile is None:
                        validMaskImagePath = validMaskImageProj
                    else:
                        validMaskImagePath = validMaskImage
                sensorClass.generateMetaDataFile(outFilePath, outName, productsStr, validMaskImagePath, prodsToCalc["FOOTPRINT"])
                prodsCalculated["METADATA"] = True
                print("")
            
            print('Clean up anything left over...')
            sensorClass.cleanFollowProcessing()
                
        except ARCSIException as e:
            print("Error: {}".format(e), file=sys.stderr)
        except Exception as e:
            print("Error: {}".format(e), file=sys.stderr)
        finally:
            failedProdsList = []
            # Check all requested products have been created
            for key in prodsToCalc.keys():
                if prodsToCalc[key] is not prodsCalculated[key]:
                    failedProdsList.append(key)
            if len(failedProdsList) > 0:
                print("Error: The following products were not generated:",file=sys.stderr)
                print(" ".join(failedProdsList),file=sys.stderr)

            elapsedTime = (time.time() - startTime)
            hours = 0
            minutes = 0
            seconds = 0.0
            if elapsedTime > 60:
                minutes = round((elapsedTime/60),0)
                if minutes > 60:
                    hours = round((minutes/60),0)
                    minutes = round(minutes - (hours * 60),0)
                seconds = round(elapsedTime - (minutes * 60),2)
            else:
                seconds = round(elapsedTime,2)
                
            print("Processing took " + str(hours) + "h " + str(minutes) + "m " + str(seconds) + "s. Thank you for using ARCSI.\n\n")        
        
    def listSensors(self):
        """
        A function which lists the currently supported sensors
        and the names by which they should be specified to the 
        ARCSI command line argument.
        """
        print("Supported Sensors are:")
        print("\t-----------------------------------------------------------------------------------------------------")
        print("\tSensor        | Shorthand     | Functions")
        print("\t-----------------------------------------------------------------------------------------------------")
        print("\tLandsat 1 MSS | \'ls1\'       | RAD, TOA, DOSAOT, DOSAOTSGL, SREF, DOS, TOPOSHADOW, FOOTPRINT, METADATA")
        print("\tLandsat 2 MSS | \'ls2\'       | RAD, TOA, DOSAOT, DOSAOTSGL, SREF, DOS, TOPOSHADOW, FOOTPRINT, METADATA")
        print("\tLandsat 3 MSS | \'ls3\'       | RAD, TOA, DOSAOT, DOSAOTSGL, SREF, DOS, TOPOSHADOW, FOOTPRINT, METADATA")
        print("\tLandsat 4 MSS | \'ls4mss\'    | RAD, TOA, DOSAOT, DOSAOTSGL, SREF, DOS, TOPOSHADOW, FOOTPRINT, METADATA")
        print("\tLandsat 4 TM  | \'ls4tm\'     | RAD, TOA, DOSAOT, DDVAOT, DOSAOTSGL, SREF, DOS, THERMAL, TOPOSHADOW, FOOTPRINT, METADATA")
        print("\tLandsat 5 MSS | \'ls5mss\'    | RAD, TOA, DOSAOT, DOSAOTSGL, SREF, DOS, TOPOSHADOW, FOOTPRINT, METADATA")
        print("\tLandsat 5 TM  | \'ls5tm\'     | RAD, TOA, DOSAOT, DDVAOT, DOSAOTSGL, SREF, DOS, THERMAL, TOPOSHADOW, FOOTPRINT, METADATA")
        print("\tLandsat 7 ETM | \'ls7\'       | RAD, TOA, DOSAOT, DDVAOT, DOSAOTSGL, SREF, DOS, THERMAL, TOPOSHADOW, FOOTPRINT, METADATA")
        print("\tLandsat 8     | \'ls8\'       | RAD, TOA, DOSAOT, DDVAOT, DOSAOTSGL, SREF, DOS, THERMAL, TOPOSHADOW, FOOTPRINT, METADATA")
        print("\tRapideye      | \'rapideye\'  | RAD, TOA, DOSAOT, DOSAOTSGL, SREF, DOS, TOPOSHADOW, METADATA")
        print("\tWorldView2    | \'wv2\'       | RAD, TOA, DOSAOT, DOSAOTSGL, SREF, DOS, TOPOSHADOW, METADATA")
        print("\tSPOT5         | \'spot5\'     | RAD, TOA, DOSAOT, DOSAOTSGL, SREF, DOS, TOPOSHADOW, METADATA")
        print("\t-----------------------------------------------------------------------------------------------------")
        
    def listProductDescription(self, product):
        """
        A function which lists the currently supported products
        and describes what that are and the parameters they require.
        """
        print("Hello World. this has not be written yet!")
    
    def listEnvVars(self):
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


if __name__ == '__main__':
    """
    The command line user interface to ARCSI
    """
    
    print("ARCSI " + ARCSI_VERSION + " Copyright (C) " + ARCSI_COPYRIGHT_YEAR + " " + ARCSI_COPYRIGHT_NAMES)
    print("This program comes with ABSOLUTELY NO WARRANTY.")
    print("This is free software, and you are welcome to redistribute it")
    print("under certain conditions; See website (" + ARCSI_WEBSITE + ").")
    print("Bugs are to be reported to " + ARCSI_SUPPORT_EMAIL + ".\n")
    
    parser = argparse.ArgumentParser(prog='arcsi',
                                    description='''Software for the Atmospheric
                                                and Radiometric Correction of
                                                Satellite Imagery (ARCSI)''',
                                    epilog='''Using 6S and rsgislib this
                                            software (ARCSI) provides a unified
                                            set of scripts for taking \'RAW\'
                                            digital numbers and converting them
                                            into radiance, top of atmosphere
                                            reflectance, surface reflectance using
                                            modelled atmosphere via 6S or undertaking
                                            a dark object subtraction. New sensors
                                            are easy to add so get in touch if we
                                            don't currently support the sensor you 
                                            require.''')
    # Request the version number.
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version ' + ARCSI_VERSION)
    # Define the argument to define the debug mode for arcsi.
    parser.add_argument("--debug", action='store_true', default=False, 
                        help='''If define the debug mode will be activated, 
                        therefore intermediate files will not be deleted.''')
    # Define the argument for specifying the input images header file.
    parser.add_argument("-i", "--inputheader", type=str, 
                        help='''Specify the input image header file.''')
    # Define the argument for specifying the input image file - overriding the header file.
    parser.add_argument("--imagefile", type=str, 
                        help='''Specify the input image file, overriding the image header.''')
    # Define the argument for specifying the sensor.
    parser.add_argument("-s", "--sensor", choices=ARCSI_SENSORS_LIST,  
                        help='''Specify the sensor being processed.''')
    # Define the argument for requesting a list of the supported sensors.
    parser.add_argument("--sensorlist", action='store_true', default=False, 
                        help='''List the sensors which are supported 
                                and the require names.''')
    # Define the argument for requesting a list of the available environment variables.
    parser.add_argument("--envvars", action='store_true', default=False, 
                        help='''List the available environmental variables for ARCSI.''')
    # Define the argument for specifying the WKT projection file for the input file.
    parser.add_argument("--inwkt", type=str, 
                        help='''Specify the WKT projection of the input image with projection defined with WKT file.''')
    # Define the argument for specifying the WKT projection file for the outputs.                    
    parser.add_argument("--outwkt", type=str, 
                        help='''Transform the outputs to the projection defined with WKT file.''')
    # Define the argument for string added to the output files names indicating the projection.
    parser.add_argument("--projabbv", type=str, 
                        help='''Abbreviation or acronym for the project which will added to the file name.''')
    # Define the argument for specifying the image file format.
    parser.add_argument("-f", "--format", type=str, 
                        help='''Specify the image output format (GDAL name).''')
    # Define the argument for specifying the output image base file name if it is
    # not to be automatically generated.
    parser.add_argument("--outbasename", type=str,
                        help='''Specify the output file base name if it
                        is not to be generated by the system.''')
    # Define the argument for specifying the file path of the output images.
    parser.add_argument("-o", "--outpath", type=str,
                        help='''Specify the output file path.''')
                        
    # Define the argument for specifying the file path of the output images.
    parser.add_argument("-t", "--tmpath", type=str,
                        help='''Specify a tempory path for files to be written to temporarly during processing if required (DDVAOT, DOS and CLOUDS).''')
    
    # Define the argument which specifies the products which are to be generated.
    parser.add_argument("-p", "--prods", type=str, nargs='+', choices=ARCSI_PRODUCTS_LIST,
                        help='''Specify the output products which are to be
                        calculated, as a comma separated list.''')
    # Define the argument for requesting a list of products.
    parser.add_argument("--prodlist", action='store_true', default=False, 
                        help='''List the products which are supported and 
                        their input requirements.''')
    # Define the argument which specifies the standard aersol profile to use.
    parser.add_argument("--aeropro", type=str, choices=['NoAerosols', 'Continental', 
    'Maritime', 'Urban', 'Desert', 'BiomassBurning', 'Stratospheric'],
                        help='''Specify the 6S defined aersol profile to use. 
                        (NoAerosols, Continental, Maritime, Urban, Desert, BiomassBurning, Stratospheric)''')
    # Define the argument which specifies the standard atompheric profile to use.
    parser.add_argument("--atmospro", type=str, choices=['NoGaseousAbsorption', 'Tropical', 
    'MidlatitudeSummer', 'MidlatitudeWinter', 'SubarcticSummer', 'SubarcticWinter', 'USStandard1962'],
                        help='''Specify the 6S defined atmospheric profile to use. 
                        (NoGaseousAbsorption, Tropical, MidlatitudeSummer, MidlatitudeWinter, 
                        SubarcticSummer, SubarcticWinter, USStandard1962)''')  
    # Define the argument for specifying the file path for the image specifying the generic aerosol model.
    parser.add_argument("--aeroimg", type=str, help='''Specify the aerosol model image file path.''') 
    # Define the argument for specifying the file path for the image specifying the generic atmosphere model.
    parser.add_argument("--atmosimg", type=str, help='''Specify the atmosphere model image file path.''')                                              
    # Define the argument which specifies the amount of OZone in atmosphere
    parser.add_argument("--atmosozone", type=float, 
                        help='''Specify the total amount of ozone in a vertical path 
                        through the atmosphere (in cm-atm)''')
    # Define the argument which specifies the amount of water in the atmosphere
    parser.add_argument("--atmoswater", type=float,
                        help='''Specify the  total amount of water in a vertical path 
                        through the atmosphere (in g/cm^2)''')
    # Define the argument which specifies the proportion of water-like aerosols
    parser.add_argument("--aerowater", type=float,
                        help='''Specify the proportion of water-like aerosols
                        (water, dust, oceanic and soot proportions must add up 
                        to 1 although all do not been to be specified).''')
    # Define the argument which specifies the proportion of dust-like aerosols 
    parser.add_argument("--aerodust", type=float,
                        help='''Specify the proportion of dust-like aerosols
                        (water, dust, oceanic and soot proportions must add up 
                        to 1 although all do not been to be specified).''')
    # Define the argument which specifies the proportion of oceanic-like aerosols
    parser.add_argument("--aerooceanic", type=float,
                        help='''Specify the proportion of oceanic-like aerosols
                        (water, dust, oceanic and soot proportions must add up 
                        to 1 although all do not been to be specified).''')
    # Define the argument which specifies the proportion of soot-like aerosols
    parser.add_argument("--aerosoot", type=float,
                        help='''Specify the proportion of soot-like aerosols
                        (water, dust, oceanic and soot proportions must add up 
                        to 1 although all do not been to be specified).''')                        
    # Define the argument which specifies the ground reflectance.
    parser.add_argument("--grdrefl", type=str, default="GreenVegetation", choices=['GreenVegetation',  
    'ClearWater', 'Sand', 'LakeWater', 'BRDFHapke'], help='''Specify the ground reflectance used for the 
                                                             6S model. (GreenVegetation, ClearWater, Sand, 
                                                             LakeWater, BRDFHapke).''')
    # Define the argument for specifying the surface elevation for the scene.
    parser.add_argument("--surfacealtitude", type=float, default="0",
                        help='''Specify the altiude (in km) of the surface being sensed.''')
    # Define the argument for specifying the AOT value for the scene
    parser.add_argument("--aot", type=float, 
                        help='''Specifiy the AOT value for the scene. 
                                If the AOT is specified the visability is ignored.''')
    # Define the argument for specifying the visability value for the scene
    parser.add_argument("--vis", type=float, 
                        help='''Specifiy the visibility value for the scene. 
                                If the AOT is specified the visability is ignored.''')
    # Define the argument for specifying the AOT value for the scene
    parser.add_argument("--minaot", type=float, default=0.05,
                        help='''Specify the minimum AOT value within the search space 
                                used to identify AOT values for the scene''')
                        # Define the argument for specifying the AOT value for the scene
    parser.add_argument("--maxaot", type=float, default=0.5,
                        help='''Specify the maximum AOT value within the search space 
                                used to identify AOT values for the scene.''')
    # Define the argument for specifying the AOT value for the scene
    parser.add_argument("--lowaot", type=float, default=0.1,
                        help='''Specify the lower AOT amount to be removed from the AOT 
                                estimate for defining --minaot within search space. (Default 0.1)''')
                        # Define the argument for specifying the AOT value for the scene
    parser.add_argument("--upaot", type=float, default=0.4,
                        help='''Specify the upper AOT amount to be added to the AOT 
                                estimate for defining --maxaot within search space. (Default 0.4)''')
    # Define the argument for specifying the AOT image file for the scene
    parser.add_argument("--aotfile", type=str, 
                        help='''Specifiy an image file with AOT values for the
                                correction. An LUT for AOT and elevation will be generated.
                                Therefore, --dem needs to be provided alongside --aotfile.''')
    # Define the argument for specifying that statistics and pyramids should be built for 
    # all output images.
    parser.add_argument("--stats", action='store_true', default=False, 
                        help='''Specifies that the image statistics and
                        pyramids should be build for all output images.''')
    parser.add_argument("-d", "--dem", type=str,
                        help='''Specify a DEM which is to be used for building
                        an LUT and applying 6S coefficients with respect to elevation.''')
    parser.add_argument("--localdos", action='store_true', default=False, 
                        help='''Specifies that a local DOS should be applied
                        rather than a global DOS.''')
    parser.add_argument("--simpledos", action='store_true', default=False, 
                        help='''Specifies that a simple (basic) DOS should be applied
                        rather than the more complex variable global/local DOS methods.''')
    parser.add_argument("--dosout", type=float, default=20, 
                        help='''Specifies the reflectance value to which dark objects
                        are set to during the dark object subtraction. (Default is 20, 
                        which is equivalent to 2 percent reflectance.''')
    parser.add_argument("--scalefac", type=int, default=1000, 
                        help='''Specifies the scale factor for the reflectance 
                        products.''')
    parser.add_argument("--interp", type=str, default="cubic", 
                        choices=['near', 'bilinear', 'cubic', 'cubicspline', 'lanczos'], 
                        help='''Specifies interpolation algorithm when reprojecting the imagery
                                (Note. the options are those in gdalwarp).''')
    
                        
    # Call the parser to parse the arguments.
    args = parser.parse_args()
        
    arcsiObj = ARCSI()
    arcsiUtils = ARCSIUtils()
    
    if args.sensorlist:
        arcsiObj.listSensors()
    elif args.prodlist:
        arcsiObj.listProductDescription()
    elif args.envvars:
        arcsiObj.listEnvVars()
    else:
        # Check that the input header parameter has been specified.
        if args.inputheader == None:
            print("Error: No input header image file has been provided.\n")
            parser.print_help()
            sys.exit()
    
        # Check that the senor parameter has been specified.
        if args.sensor == None:
            print("Error: No sensor has been provided.\n")
            parser.print_help()
            sys.exit()
    
        # Check that the output image format has been specified.
        if args.format == None:
            # Print an error message if not and exit.
            envVar = arcsiUtils.getEnvironmentVariable("ARCSI_OUT_FORMAT")
            if envVar == None:
                print("Error: No output image format provided.\n")
                parser.print_help()
                sys.exit()
            else:
                print("Taking output format from environment variable.")
                args.format = envVar
        
        # Check that the output image format has been specified.
        if args.outpath == None:
            # Print an error message if not and exit.
            envVar = arcsiUtils.getEnvironmentVariable("ARCSI_OUTPUT_PATH")
            if envVar == None:
                print("Error: No output file path has been provided.\n")
                parser.print_help()
                sys.exit()
            else:
                print("Taking output file path from environment variable.")
                args.outpath = envVar
        
        if not os.path.exists(args.outpath):
            print("WARNING: Output directory does not exist so creating it...")
            os.makedirs(args.outpath)
        elif not os.path.isdir(args.outpath):
            print("ERROR: Output Path exists but is not a directory...\n")
            parser.print_help()
            sys.exit()
        
        if not args.outwkt is None:
            if not os.path.exists(args.outwkt):
                print("Error: The output WKT file does not exist.\n")
                sys.exit()
            elif args.projabbv == None:
                print("WARNING: It is recommended that a projection abbreviation or acronym is provided (--projabbv)...")
                 
        needAOD = False
        needAODMinMax = False
        needTmp = False
        needDEM = False
        for prod in args.prods:
            if prod == 'DDVAOT':
                needAODMinMax = True
                needTmp = True
                needDEM = True
            elif prod == 'SREF':
                needAOD = True
            elif prod == 'DOS':
                if not args.simpledos:
                    needTmp = True
            elif prod == 'CLOUDS':
                needTmp = True
            elif prod == 'DOSAOT':
                needAODMinMax = True
                needTmp = True
                needDEM = True
            elif prod == 'DOSAOTSGL':
                needAODMinMax = True
                needTmp = True
                needDEM = True
            elif prod == 'TOPOSHADOW':
                needTmp = True
                needDEM = True
            
        if needAODMinMax and (args.minaot == None) and (args.maxaot == None):
            envVarMinAOT = arcsiUtils.getEnvironmentVariable("ARCSI_MIN_AOT")
            if envVarMinAOT == None:
                print("Error: The min and max AOT values for the search should be specified.\n")
                parser.print_help()
                sys.exit()
            else:
                print("Taking min AOT from environment variable.")
                args.minaot = float(envVarMinAOT)
                
            envVarMaxAOT = arcsiUtils.getEnvironmentVariable("ARCSI_MAX_AOT")
            if envVarMaxAOT == None:
                print("Error: The min and max AOT values for the search should be specified.\n")
                parser.print_help()
                sys.exit()
            else:
                print("Taking max AOT from environment variable.")
                args.maxaot = float(envVarMaxAOT)
        
        if args.lowaot is None:     
            envVarLowAOT = arcsiUtils.getEnvironmentVariable("ARCSI_LOW_AOT")
            if not envVarLowAOT is None:
                args.lowaot = float(envVarLowAOT)
        
        if args.upaot is None:     
            envVarUpAOT = arcsiUtils.getEnvironmentVariable("ARCSI_UP_AOT")
            if not envVarUpAOT is None:
                args.upaot = float(envVarUpAOT)
            
        
        if needAOD and (args.aot == None) and (args.vis == None) and (not needAODMinMax) and (not args.aotfile):
            print("Error: Either the AOT or the Visability need to specified. Or --aotfile needs to be provided.\n")
            parser.print_help()
            sys.exit()
                    
        if needTmp and args.tmpath == None:
            envVar = arcsiUtils.getEnvironmentVariable("ARCSI_TMP_PATH")
            if envVar == None:
                print("Error: If the DDVAOT, DOS, DOSAOTSGL, CLOUDS or TOPOSHADOW product is set then a tempory path needs to be provided.\n")
                parser.print_help()
                sys.exit()
            else:
                print("Taking temp path from environment variable.")
                args.tmpath = envVar
        
        if needTmp: 
            if not os.path.exists(args.tmpath):
                print("WARNING: The temp path specified does not exist, it is being created.")
                os.makedirs(args.tmpath)
            if not os.path.isdir(args.tmpath):
                print("Error: The temp path specified is not a directory, please correct and run again.\n")
                parser.print_help()
                sys.exit()
            
        if args.dem == None:
            envVar = arcsiUtils.getEnvironmentVariable("ARCSI_DEM_PATH")
            if not envVar == None:
                args.dem = envVar
                print("Taking DEM path from environment variable.")
        
        if needDEM:
            if not os.path.exists(args.dem):
                print("Error: A file path to a DEM has either not been specified or does exist, please check it and run again.\n")
                parser.print_help()
                sys.exit()
        
        if args.aeroimg == None:
            envVar = arcsiUtils.getEnvironmentVariable("ARCSI_AEROIMG_PATH")
            if not envVar == None:
                args.aeroimg = envVar
                print("Taking aerosol profile image path from environment variable.")
            else:
                args.aeroimg = arcsilib.DEFAULT_ARCSI_AEROIMG_PATH

        if args.atmosimg == None:
            envVar = arcsiUtils.getEnvironmentVariable("ARCSI_ATMOSIMG_PATH")
            if not envVar == None:
                args.atmosimg = envVar
                print("Taking atmosphere profile image path from environment variable.")
            else:
                args.atmosimg = arcsilib.DEFAULT_ARCSI_ATMOSIMG_PATH


        if args.atmosimg == None:
            envVar = arcsiUtils.getEnvironmentVariable("ARCSI_ATMOSIMG_PATH")
            if not envVar == None:
                args.atmosimg = envVar
                print("Taking atmosphere profile image path from environment variable.")
                
        if args.outwkt == None:
            envVar = arcsiUtils.getEnvironmentVariable("ARCSI_OUTPUT_WKT")
            if not envVar == None:
                args.outwkt = envVar
                print("Taking output WKT from environment variable.")
        
        if args.projabbv == None:
            envVar = arcsiUtils.getEnvironmentVariable("ARCSI_PROJ_ABBV")
            if not envVar == None:
                args.projabbv = envVar
                print("Taking projection abbreviation from environment variable.")
        
        atmosOZoneWaterSpecified = False
        if (not args.atmosozone == None) and (args.atmoswater == None):
            print("Error: If the atmospheric ozone is defined then the atmospheric water needs to be specfied --atmoswater.\n")
            parser.print_help()
            sys.exit()
        elif (not args.atmoswater == None) and (args.atmosozone == None):
            print("Error: If the atmospheric water is defined then the atmospheric ozone needs to be specfied --atmosozone.\n")
            parser.print_help()
            sys.exit()
        elif (not args.atmoswater == None) and (not args.atmosozone == None):
            atmosOZoneWaterSpecified = True

        aeroComponentsSpecified = False
        if (not args.aerowater == None) or (not args.aerodust == None) or (not args.aerooceanic == None) or (not args.aerosoot == None):
            aeroComponentsSpecified = True
        
        if args.dosout == None:
            envVar = arcsiUtils.getEnvironmentVariable("ARCSI_OUTDOS_REFL")
            if not envVar == None:
                args.dosout = envVar
                print("Taking output DOS reflectance from environment variable.")
                
        envVar = arcsiUtils.getEnvironmentVariable("ARCSI_USE_LOCALDOS")
        if not envVar == None:
            if envVar == "TRUE":
                args.localdos = True
                print("Using local DOS method due to environment variable.")
            else:
                args.localdos = False
                print("Using global DOS method due to environment variable.")
        
        envVar = arcsiUtils.getEnvironmentVariable("ARCSI_USE_SIMPLEDOS")
        if not envVar == None:
            if envVar == "TRUE":
                args.simpledos = True
                print("Using simple DOS method due to environment variable.")
            else:
                args.simpledos = False
                print("Not using simple DOS method due to environment variable.")
        
        arcsiObj.run(args.inputheader, args.imagefile, args.sensor, args.inwkt, args.format, args.outpath, 
                     args.outbasename, args.outwkt, args.projabbv, args.prods, args.stats, args.aeropro,  
                     args.atmospro, args.aeroimg, args.atmosimg, args.grdrefl, args.surfacealtitude, 
                     args.atmosozone, args.atmoswater, atmosOZoneWaterSpecified, args.aerowater, 
                     args.aerodust, args.aerooceanic, args.aerosoot, aeroComponentsSpecified, 
                     args.aot, args.vis, args.tmpath, args.minaot, args.maxaot, args.lowaot, args.upaot, 
                     args.dem, args.aotfile, (not args.localdos), args.dosout, args.simpledos, args.debug, 
                     args.scalefac, args.interp)




