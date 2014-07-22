#! /usr/bin/env python

"""
Module that contains the ARSCI command to build lists of commands for arcsi.py.
"""

############################################################################
#  arcsibuildcmdslist.py
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
# Purpose:  A script to build lists of commands for the arcsi.py script
#           and export them as a shell script which can either be executed
#           on its own or used as import into the GNU parallel command. 
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 29/01/2014
# Version: 1.0
#
# History:
# Version 1.0 - Created.
#
############################################################################

# Import the python os.path module
import os.path
# Import the python sys module
import sys
# Import the python glob module
import glob
# Import the python Argument parser
import argparse

class ARCSIBuildCommands (object):
        
    def buildCmds(self, inputDIR, outputFile, headerEnding, noFolders, sensor, inwkt, 
                     format, outpath, prods, stats, aeropro, atmospro, 
                     aeroimg, atmosimg, grdrefl, surfacealtitude, 
                     atmosozone, atmoswater, aerowater, 
                     aerodust, aerooceanic, aerosoot, aot, vis, tmpath, 
                     minaot, maxaot, dem, localdos, dosout):
        
        inputDIR = os.path.abspath(inputDIR)
        outputFile = os.path.abspath(outputFile)
        if not aeroimg == None:
            aeroimg = os.path.abspath(aeroimg)
        if not atmosimg == None:
            atmosimg = os.path.abspath(atmosimg)
        if not dem == None:
            dem = os.path.abspath(dem)
        if not tmpath == None:
            tmpath = os.path.abspath(tmpath)
        if not outpath == None:
            outpath = os.path.abspath(outpath)
        if not inwkt == None:
            inwkt = os.path.abspath(inwkt)
            
        headersFilesList = list()
        searchPath = ""
        if noFolders:
            searchPath = os.path.join(inputDIR, "*" + headerEnding)
        else:
            searchPath = os.path.join(inputDIR, "*", "*" + headerEnding)
        print("Search Path: ", searchPath)
        headersFilesList = glob.glob(searchPath)
                
        prodsStr = ""
        first = True
        for prod in prods:
            if first:
                prodsStr = prod
                first = False
            else:
                prodsStr = prodsStr + " " + prod

        outFile = open(outputFile, 'w+')
        for hFile in headersFilesList:
            print("Processing :", hFile)
            cmd = "arcsi.py -s " + sensor + " -p " + prodsStr + " -i " + hFile
            if not outpath == None:
                cmd = cmd + " --outpath " + outpath
            if stats:
               cmd = cmd + " --stats"
            if not format == None:
               cmd = cmd + " --format " + format
            else:
                cmd = cmd + " --format KEA"
            if not tmpath == None:
               cmd = cmd + " --tmpath " + tmpath
            if not dem == None:
               cmd = cmd + " --dem " + dem
            if not aeroimg == None:
               cmd = cmd + " --aeroimg " + aeroimg
            elif not aeropro == None:
               cmd = cmd + " --aeropro " + aeropro
            if not atmosimg == None:
               cmd = cmd + " --atmosimg " + atmosimg
            elif not atmospro == None:
               cmd = cmd + " --atmospro " + atmospro
            if not surfacealtitude == None:
               cmd = cmd + " --surfacealtitude " + str(surfacealtitude)
            if (not minaot == None) and (not maxaot == None):
                cmd = cmd + " --minaot " + str(minaot) + " --maxaot " + str(maxaot)
            elif not aot == None:
               cmd = cmd + " --aot " + str(aot)
            elif not vis == None:
               cmd = cmd + " --vis " + str(vis)
            if not inwkt == None:
               cmd = cmd + " --inwkt " + inwkt
            if not grdrefl == None:
                cmd = cmd + " --grdrefl " + grdrefl
            if not atmosozone == None:
                cmd = cmd + " --atmosozone " + str(atmosozone)
            if not atmoswater == None:
                cmd = cmd + " --atmoswater " + str(atmoswater)
            if not aerowater == None:
                cmd = cmd + " --aerowater " + str(aerowater)
            if not aerodust == None:
                cmd = cmd + " --aerodust " + str(aerodust)
            if not aerooceanic == None:
                cmd = cmd + " --aerooceanic " + str(aerooceanic)
            if not aerosoot == None:
                cmd = cmd + " --aerosoot " + str(aerosoot)
            if not dosout == None:
                cmd = cmd + " --dosout " + str(dosout)
            if localdos:
                cmd = cmd + " --localdos "
            
            print(cmd)
            outFile.write(cmd + "\n")
        outFile.flush()
        outFile.close()

if __name__ == '__main__':
    """
    The command line user interface to ARCSI Data Extraction Tool.
    """
    parser = argparse.ArgumentParser(prog='arcsiextractdata',
                                    description='''ARCSI command extract data
                                                   from tar or tar.gz archives.''',
                                    epilog='''A tools to extract data
                                              from tar or tar.gz archives into
                                              individual directories per image''')
    # Request the version number.
    parser.add_argument('-v', '--version', action='version', version='%(prog)s Version 0.1a')
    
    parser.add_argument("-i", "--input", type=str, required=True, 
                        help='''Input directory containing the data to be processed''')

    parser.add_argument("-o", "--output", type=str, required=True, 
                        help='''Output text file (shell script) with the list of commands.''')
                        
    parser.add_argument("-e", "--header", type=str, required=True, 
                        help='''The extension / unquie file ending for the input header files.''')
    
    parser.add_argument("--nofolders", action='store_true', default=False, 
                        help='''Specifies whether that the the commands should only look for files
                                within the specified input directory and not go down the directory tree.''')
                                                           
    # Define the argument for specifying the sensor.
    parser.add_argument("-s", "--sensor", choices=['ls1', 'ls2', 'ls3', 'ls4mss', 'ls4tm',
                                                   'ls5mss', 'ls5tm', 'ls7', 
                                                   'ls8', 'rapideye', 'wv2'],  
                        help='''Specify the sensor being processed.''')
                        
    # Define the argument for specifying the WKT projection file 
    # for the intput file.
    parser.add_argument("--inwkt", type=str, 
                        help='''Specify the WKT projection of the input image''')
                        
    # Define the argument for specifying the image file format.
    parser.add_argument("-f", "--format", type=str, 
                        help='''Specify the image output format (GDAL name).''')
                        
    parser.add_argument("--outpath", type=str, 
                        help='''Specifiy the output directory for the products from ARCSI.''')
                        
    # Define the argument for specifying the file path of the output images.
    parser.add_argument("-t", "--tmpath", type=str,
                        help='''Specify a tempory path for files to be written to temporarly during processing if required (DDVAOT, DOSUB and CLOUDS).''')
    
    # Define the argument which specifies the products which are to be generated.
    parser.add_argument("-p", "--prods", type=str, nargs='+', choices=['RAD', 'SATURATE', 'TOA', 'CLOUDS', 'DDVAOT', 'DOSAOT', 'SREFSTDMDL', 'DOSUB', 'THERMAL'],
                        help='''Specify the output products which are to be
                        calculated, as a comma separated list. (RAD, SATURATE, TOA, CLOUDS, DDVAOT, DOSAOT, SREFSTDMDL, DOSUB, THERMAL)''')
                        
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
                                                6S model. (GreenVegetation, ClearWater, Sand, LakeWater, BRDFHapke).''')
                                                
    # Define the argument for specifying the surface elevation for the scene.
    parser.add_argument("--surfacealtitude", type=float, 
                            help='''Specify the altiude (in km) of the surface being sensed.''')
                        
    # Define the argument for specifying the AOT value for the scene
    parser.add_argument("--aot", type=float, 
                        help='''Specifiy the AOT or visability value for the scene. 
                                If the AOT is specified the visability is ignored.''')
                                
    # Define the argument for specifying the visability value for the scene
    parser.add_argument("--vis", type=float, 
                        help='''Specifiy the AOT or visability value for the scene. 
                                If the AOT is specified the visability is ignored.''')
                                
    # Define the argument for specifying the AOT value for the scene
    parser.add_argument("--minaot", type=float,
                        help='''Specifiy the AOT or visibility value for the scene. 
                                If the AOT is specified the visibility is ignored.''')
                                
    # Define the argument for specifying the AOT value for the scene
    parser.add_argument("--maxaot", type=float,
                        help='''Specifiy the AOT or visability value for the scene. 
                                If the AOT is specified the visability is ignored.''')
    
    
    # Define the argument for specifying that statistics and pyramids should be built for 
    # all output images.
    parser.add_argument("--stats", action='store_true', default=False, 
                        help='''Specifies that the image statistics and
                        pyramids should be build for all output images.''')
                        
    parser.add_argument("--dem", type=str,
                        help='''Specify a DEM which is to be used for building
                        an LUT and applying 6S coefficients with respect to elevation.''')
    
    parser.add_argument("--localdos", action='store_true', default=False, 
                        help='''Specifies that a local DOS should be applied
                        rather than a global DOS.''')
    
    parser.add_argument("--dosout", type=float, 
                        help='''Specifies the reflectance value to which dark objects
                        are set to during the dark object subtraction. (Default is 20, 
                        which is equivalent to 2 % reflectance.''')
                        
    # Call the parser to parse the arguments.
    args = parser.parse_args()
    
    if args.input == None:
        print("An input directory was not specified.")
        parser.print_help()
        sys.exit()
        
    if args.output == None:
        print("An output file name and path was not specified.")
        parser.print_help()
        sys.exit()
        
    if args.header == None:
        print("The unique file name ending of the header files was not specified.")
        parser.print_help()
        sys.exit()
    
    if args.sensor == None:
        print("The sensor needs to be specified.")
        parser.print_help()
        sys.exit()
        
    if args.prods == None:
        print("The list of produced to be generated must be specified.")
        parser.print_help()
        sys.exit()
    
    arcsiObj = ARCSIBuildCommands()
    
    arcsiObj.buildCmds(args.input, args.output, args.header, args.nofolders, args.sensor, args.inwkt, 
                     args.format, args.outpath, args.prods, args.stats, args.aeropro, args.atmospro, 
                     args.aeroimg, args.atmosimg, args.grdrefl, args.surfacealtitude, 
                     args.atmosozone, args.atmoswater, args.aerowater, 
                     args.aerodust, args.aerooceanic, args.aerosoot, args.aot, args.vis, args.tmpath, 
                     args.minaot, args.maxaot, args.dem, args.localdos, args.dosout)
    
    