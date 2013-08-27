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

# Import the system library
import sys
# Import the python Argument parser
import argparse
# Import the time module
import time
# Import the ARCSI exception class
from arcsi.arcsiexception import ARCSIException
# Import the ARCSI utilities class
from arcsi.arcsiutils import ARCSIUtils
# Import the sensor classes
from arcsi.arcsisensor import ARCSIAbstractSensor
from arcsi.arcsisensorlandsat7 import ARCSILandsat7Sensor
from arcsi.arcsisensorlandsat5tm import ARCSILandsat5TMSensor
from arcsi.arcsisensorlandsat4tm import ARCSILandsat4TMSensor
from arcsi.arcsisensorlandsat5mss import ARCSILandsat5MSSSensor
from arcsi.arcsisensorlandsat4mss import ARCSILandsat4MSSSensor
from arcsi.arcsisensorlandsat3mss import ARCSILandsat3MSSSensor
from arcsi.arcsisensorlandsat2mss import ARCSILandsat2MSSSensor
from arcsi.arcsisensorlandsat1mss import ARCSILandsat1MSSSensor
from arcsi.arcsisensorlandsat8 import ARCSILandsat8Sensor
from arcsi.arcsisensorrapideye import ARCSIRapidEyeSensor
# Import the image utilities module from rsgislib
import rsgislib.imageutils
# Import the py6s module for running 6S from python.
import Py6S

class ARCSI (object):
    """
    The \'main\' class which executes the whole ARCSI package.
    """
    
    def sensorClassFactory(self, sensor):
        sensorClass = None
        if sensor == 'ls7':
            sensorClass = ARCSILandsat7Sensor()
        elif sensor == 'ls5tm':
            sensorClass = ARCSILandsat5TMSensor()
        elif sensor == 'ls4tm':
            sensorClass = ARCSILandsat4TMSensor()
        elif sensor == 'ls5mss':
            sensorClass = ARCSILandsat5MSSSensor()
        elif sensor == 'ls4mss':
            sensorClass = ARCSILandsat4MSSSensor()
        elif sensor == 'ls3':
            sensorClass = ARCSILandsat3MSSSensor()
        elif sensor == 'ls2':
            sensorClass = ARCSILandsat2MSSSensor()
        elif sensor == 'ls1':
            sensorClass = ARCSILandsat1MSSSensor()
        elif sensor == 'ls8':
            sensorClass = ARCSILandsat8Sensor()
        elif sensor == 'rapideye':
            sensorClass = ARCSIRapidEyeSensor()
        else:
            raise ARCSIException("Could not get a class representing the sensor specified from the factory.")
        return sensorClass
        
    def convertVisabilityToAOD(self, vis):
        return (3.9449/vis)+0.08498
    
    def run(self, inputHeader, sensorStr, inWKTFile, outFormat, outFilePath, outBaseName, productsStr, calcStatsPy, aeroProfileOption, atmosProfileOption, grdReflOption, surfaceAltitude, atmosOZoneVal, atmosWaterVal, aeroWaterVal, aeroDustVal, aeroOceanicVal, aeroSootVal, aotVal, visVal, tmpPath, minAOT, maxAOT):
        """
        A function contains the main flow of the software
        """
        print("ASCSI 0.1a Copyright (C) 2013  Peter Bunting")
        print("This program comes with ABSOLUTELY NO WARRANTY.")
        print("This is free software, and you are welcome to redistribute it")
        print("under certain conditions; See website (http://www.rsgislib.org/ascsi).")
        print("Bugs are to be reported to pfb@aber.ac.uk\n")
        
        startTime = time.time()
        arcsiUtils = ARCSIUtils()
        try:
            # Read WKT file if provided.
            wktStr = None
            if inWKTFile != None:
                wktStr = arcsiUtils.readTextFile(inWKTFile)
            # Step 1: Get the Sensor specific class from factory
            sensorClass = self.sensorClassFactory(sensorStr)
            # Step 2: Read header parameters
            sensorClass.extractHeaderParameters(inputHeader, wktStr)
            # Step 3: Get Output Image Base Name.
            if outBaseName == None:
                outBaseName = sensorClass.generateOutputBaseName()
            print("Image Base Name: " + outBaseName)
            # Step 4: Find the products which are to be generated.
            prodsToCalc = dict()
            prodsToCalc["RAD"] = False
            prodsToCalc["TOA"] = False
            prodsToCalc["DDVAOT"] = False
            prodsToCalc["SREFSTDMDL"] = False
            prodsToCalc["SREFSGLPARAM"] = False
            
            for prod in productsStr:
                if prod == 'RAD':
                    prodsToCalc["RAD"] = True
                elif prod == 'TOA':
                    prodsToCalc["RAD"] = True
                    prodsToCalc["TOA"] = True
                elif prod == 'DDVAOT':
                    prodsToCalc["RAD"] = True
                    prodsToCalc["TOA"] = True
                    prodsToCalc["DDVAOT"] = True
                elif prod == 'SREFSTDMDL':
                    prodsToCalc["RAD"] = True
                    prodsToCalc["SREFSTDMDL"] = True
                elif prod == 'SREFSGLPARAM':
                    prodsToCalc["RAD"] = True
                    prodsToCalc["SREFSGLPARAM"] = True
                
            
            radianceImage=""
            toaImage = ""
            srefImage = ""
            aodImage = ""
            # Step 5: Convert to Radiance
            if prodsToCalc["RAD"]:
                # Execute conversion to radiance
                outName = outBaseName + "_rad" + arcsiUtils.getFileExtension(outFormat)
                radianceImage = sensorClass.convertImageToRadiance(outFilePath, outName, outFormat)
                if calcStatsPy:
                    print("Calculating Statistics...")
                    #rsgislib.imageutils.popImageStats(radianceImage, True, 0.0, True)
            # Step 6: Convert to TOA
            if prodsToCalc["TOA"]:
                # Execute conversion to top of atmosphere reflectance
                outName = outBaseName + "_rad_toa" + arcsiUtils.getFileExtension(outFormat)
                toaImage = sensorClass.convertImageToTOARefl(radianceImage, outFilePath, outName, outFormat)
                if calcStatsPy:
                    print("Calculating Statistics...")
                    #rsgislib.imageutils.popImageStats(toaImage, True, 0.0, True)
            # Step 7: Use image to estimate AOD values 
            if prodsToCalc["DDVAOT"]:
                aeroProfile = None
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
                
                atmosProfile = None
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
                
                grdRefl = None
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
                else:
                    raise ARCSIException("The specified ground reflectance is unknown.")
                
                if (aotVal == None) and (visVal == None):
                    raise ARCSIException("Either the AOT or the visability need to specified.")
                elif (aotVal == None):
                    aotVal = self.convertVisabilityToAOD(visVal)
            
                outName = outBaseName + "_ddvaod" + arcsiUtils.getFileExtension(outFormat)
                aodImage = sensorClass.estimateImageToAOD(toaImage, outFilePath, outName, outFormat, tmpPath, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, minAOT, maxAOT)
                if calcStatsPy:
                    print("Calculating Statistics...")
                    #rsgislib.imageutils.popImageStats(aodImage, True, 0.0, True)
            # Step 8: Convert to Surface Reflectance using 6S Standard Models
            if prodsToCalc["SREFSTDMDL"]:
                # Execute conversion to surface reflectance by applying 6S using a 'standard' modelled atmosphere.
                outName = outBaseName + "_rad_srefstdmdl" + arcsiUtils.getFileExtension(outFormat)

                aeroProfile = None
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
                
                atmosProfile = None
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
                
                grdRefl = None
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
                else:
                    raise ARCSIException("The specified ground reflectance is unknown.")
                
                if (aotVal == None) and (visVal == None):
                    raise ARCSIException("Either the AOT or the visability need to specified.")
                elif (aotVal == None):
                    aotVal = self.convertVisabilityToAOD(visVal)
                
                print("AOT Value: "+ str(aotVal))
                
                srefImage = sensorClass.convertImageToSurfaceReflSglParam(radianceImage, outFilePath, outName, outFormat, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal)
                if calcStatsPy:
                    print("Calculating Statistics...")
                    rsgislib.imageutils.popImageStats(srefImage, True, 0.0, True)
            if prodsToCalc["SREFSGLPARAM"]:
                # Execute conversion to surface reflectance by applying 6S using a 'standard' modelled atmosphere.
                outName = outBaseName + "_rad_srefsglparam" + arcsiUtils.getFileExtension(outFormat)
                
                if calcStatsPy:
                    print("Calculating Statistics...")
                    #rsgislib.imageutils.popImageStats(srefImage, True, 0.0, True)
                
        except ARCSIException as e:
            print("Error: " + str(e))
        except Exception as e:
            print("Error: " + str(e))
        finally:
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
                
            print("\nProcessing took " + str(hours) + "h " + str(minutes) + "m " + str(seconds) + "s. Thank you for using ARCSI.")        
        
    def listSensors(self):
        """
        A function which lists the currently supported sensors
        and the names by which they should be specified to the 
        ARCSI command line argument.
        """
        print("Supported Sensors are:")
        print("\t-------------------------------------------------------")
        print("\tSenor         | Shorthand   | Functions")
        print("\t-------------------------------------------------------")
        print("\tLandsat 1 MSS | \'ls1\'       | RAD, TOA, SREFSTDMDL")
        print("\tLandsat 2 MSS | \'ls2\'       | RAD, TOA, SREFSTDMDL")
        print("\tLandsat 3 MSS | \'ls3\'       | RAD, TOA, SREFSTDMDL")
        print("\tLandsat 4 MSS | \'ls4mss\'    | RAD, TOA, SREFSTDMDL")
        print("\tLandsat 4 TM  | \'ls5tm\'     | RAD, TOA, SREFSTDMDL")
        print("\tLandsat 5 MSS | \'ls5mss\'    | RAD, TOA, SREFSTDMDL")
        print("\tLandsat 5 TM  | \'ls5tm\'     | RAD, TOA, SREFSTDMDL")
        print("\tLandsat 7 ETM | \'ls7\'       | RAD, TOA, SREFSTDMDL")
        print("\tLandsat 8     | \'ls8\'       | RAD, TOA, DDVAOT, SREFSTDMDL")
        print("\tRapideye      | \'rapideye\'  | RAD, TOA")
        #print("\tSPOT 5        | \'spot5\'     | RAD, TOA")
        #print("\tASTER         | \'aster\'     | RAD, TOA")
        #print("\tIRS P6        | \'irsp6\'     | RAD, TOA")
        print("\t-------------------------------------------------------")
        
    def listProductDescription(self):
        """
        A function which lists the currently supported products
        and describes what that are and the parameters they require.
        """
        print("Hello World.")


if __name__ == '__main__':
    """
    The command line user interface to ARCSI
    """
    parser = argparse.ArgumentParser(prog='arcsi',
                                    description='''Software for the Atmospheric
                                                and Radiometric Correction of
                                                Satellite Imagery (ARCSI)''',
                                    epilog='''Using 6S and rsgislib this
                                            software (ARCSI) provides a unified
                                            set of scripts for taking \'RAW\'
                                            digital numbers and converting them
                                            into radiance, top of atmosphere
                                            reflectance or surface reflectance,
                                            including options for BRDF
                                            correction.''')
    # Request the version number.
    parser.add_argument('-v', '--version', action='version', version='%(prog)s Version 0.1a')
    # Define the argument for specifying the input images header file.
    parser.add_argument("-i", "--inputheader", type=str, 
                        help='''Specify the input image header file.''')
    # Define the argument for specifying the sensor.
    parser.add_argument("-s", "--sensor", choices=['ls1', 'ls2', 'ls3', 'ls4mss', 'ls4tm',
                                                   'ls5mss', 'ls5tm', 'ls7', 
                                                   'ls8', 'spot5', 'aster', 'irsp6',
                                                   'rapideye'],  
                        help='''Specify the sensor being processed.''')
    # Define the argument for requesting a list of the supported sensors.
    parser.add_argument("--sensorlist", action='store_true', default=False, 
                        help='''List the sensors which are supported 
                                and the require names.''')
    # Define the argument for specifying the WKT projection file 
    # for the intput file.
    parser.add_argument("--inwkt", type=str, 
                        help='''Specify the WKT projection of the input image''')
    # Define the argument for specifying the image file format.
    parser.add_argument("-f", "--format", type=str, default="KEA",
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
                        help='''Specify a tempory path for files to be written to temporarly during processing if required (--ddvaod).''')
    
    # Define the argument which specifies the products which are to be generated.
    parser.add_argument("-p", "--prods", type=str, nargs='+', choices=['RAD', 'TOA', 'DDVAOT', 'SREFSTDMDL', 'SREFSGLPARAM'],
                        help='''Specify the output products which are to be
                        calculated, as a comma separated list. (RAD, TOA, DDVAOT, SREFSTDMDL, SREFSGLPARAM)''')
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
    'ClearWater', 'Sand', 'LakeWater'], help='''Specify the ground reflectance used for the 
                                                6S model. (GreenVegetation, ClearWater, Sand, LakeWater).''')
    # Define the argument for specifying the surface elevation for the scene.
    parser.add_argument("--surfacealtitude", type=float, default="0",
                        help='''Specify the altiude (in km) of the surface being sensed.''')
    # Define the argument for specifying the AOT value for the scene
    parser.add_argument("--aot", type=float, 
                        help='''Specifiy the AOT or visability value for the scene. If the AOT is specified the visability is ignored.''')
    # Define the argument for specifying the visability value for the scene
    parser.add_argument("--vis", type=float, 
                        help='''Specifiy the AOT or visability value for the scene. If the AOT is specified the visability is ignored.''')
    # Define the argument for specifying the AOT value for the scene
    parser.add_argument("--minaot", type=float, default=0.05,
                        help='''Specifiy the AOT or visability value for the scene. If the AOT is specified the visability is ignored.''')
                        # Define the argument for specifying the AOT value for the scene
    parser.add_argument("--maxaot", type=float, default=0.5,
                        help='''Specifiy the AOT or visability value for the scene. If the AOT is specified the visability is ignored.''')
    # Define the argument for specifying that statistics and pyramids should be built for 
    # all output images.
    parser.add_argument("--stats", action='store_true', default=False, 
                        help='''Specifies that the image statistics and
                        pyramids should be build for all output images.''')
    # Call the parser to parse the arguments.
    args = parser.parse_args()
    
    arcsiObj = ARCSI()
    
    if args.sensorlist:
        arcsiObj.listSensors()
    if args.prodlist:
        arcsiObj.listProductDescription()
    else:
        # Check that the input header parameter has been specified.
        if args.inputheader == None:
            print("Error: No input header image file has been provided.")
            sys.exit()
    
        # Check that the senor parameter has been specified.
        if args.sensor == None:
            print("Error: No sensor has been provided.")
            sys.exit()
    
        # Check that the output image format has been specified.
        if args.format == None:
            # Print an error message if not and exit.
            print("Error: No output image format provided.")
            sys.exit()
        
        # Check that the output image format has been specified.
        if args.outpath == None:
            # Print an error message if not and exit.
            print("Error: No output file path has been provided.")
            sys.exit()
        
        needAOD = False
        needAODMinMax = False
        needTmp = False
        for prod in args.prods:
            if prod == 'DDVAOT':
                needAODMinMax = True
                needTmp = True
            elif prod == 'SREFSTDMDL':
                needAOD = True
            elif prod == 'SREFSGLPARAM':
                needAOD = True
        if needAOD and (args.aot == None) and (args.vis == None):
            print("Error: Either the AOT or the Visability need to specified.")
            sys.exit()
            
        if needAODMinMax and (args.minaot == None) and (args.maxaot == None):
            print("Error: The min and max AOT values for the search should be specified.")
            sys.exit()
                
        if needTmp and args.tmpath == None:
                print("Error: If the DDVAOT product is set then a tempory path needs to be provided.")
                sys.exit()
            
        arcsiObj.run(args.inputheader, args.sensor, args.inwkt, args.format, args.outpath, args.outbasename, args.prods, args.stats, args.aeropro, args.atmospro, args.grdrefl, args.surfacealtitude, args.atmosozone, args.atmoswater, args.aerowater, args.aerodust, args.aerooceanic, args.aerosoot, args.aot, args.vis, args.tmpath, args.minaot, args.maxaot)




