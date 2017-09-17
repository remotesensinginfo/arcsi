#! /usr/bin/env python

"""
Module that contains the ARSCI Main class.
"""

############################################################################
#  arcsimpi.py
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
# Date: 16/09/2017
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
# Import the os.path python module
import os.path
# Import the python Argument parser
import argparse
# Import ARCSI library
import arcsilib
# Import ARCSI execution class
from arcsilib.arcsirun import *
# Import the ARCSI utilities class
from arcsilib.arcsiutils import ARCSIUtils
# Import the ARCSI enum
from arcsilib.arcsiutils import ARCSIEnum
# Import the ARCSI exception class
from arcsilib.arcsiexception import ARCSIException
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
# Import the list of gdal file formats arcsi supports
from arcsilib import ARCSI_GDALFORMATS_LIST
# Import rsgislib library
import rsgislib
# Import the MPI library
from mpi4py import MPI

# Define MPI message tags
mpiTags = ARCSIEnum('READY', 'DONE', 'EXIT', 'START')

arcsiStages = ARCSIEnum('ARCSIPART1', 'ARCSIPART2', 'ARCSIPART3', 'ARCSIPART4')

# Initializations and preliminaries
mpiComm = MPI.COMM_WORLD      # get MPI communicator object
mpiSize = mpiComm.size           # total number of processes
mpiRank = mpiComm.rank           # rank of this process
mpiStatus = MPI.Status()      # get MPI status object


if __name__ == '__main__':
    """
    The command line user interface to ARCSI
    """
    print("ARCSI " + ARCSI_VERSION + " Copyright (C) " + ARCSI_COPYRIGHT_YEAR + " " + ARCSI_COPYRIGHT_NAMES)
    print("This program comes with ABSOLUTELY NO WARRANTY.")
    print("This is free software, and you are welcome to redistribute it")
    print("under certain conditions; See website (" + ARCSI_WEBSITE + ").")
    print("Bugs are to be reported to " + ARCSI_SUPPORT_EMAIL + ".\n")

    if len(sys.argv) == 1:
        print('arcsimpi.py [options]')
        print('help : arcsimpi.py --help')
        print('  or : arcsimpi.py -h')
        print('')
    else:
        parser = argparse.ArgumentParser(prog='arcsimpi.py',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                        description='''Software for the Atmospheric and Radiometric Correction of Satellite Imagery (ARCSI)''',
                                        epilog='''Using 6S and rsgislib this software (ARCSI) provides a unified set of scripts for taking \'RAW\' digital numbers and converting them into radiance, top of atmosphere reflectance, surface reflectance using modelled atmosphere via 6S or undertaking a dark object subtraction. New sensors are easy to add so get in touch if we don't currently support the sensor you require.''')

        # Request the version number.
        parser.add_argument('-v', '--version', action='version', version='%(prog)s version ' + ARCSI_VERSION)
        # Define the argument to define the debug mode for arcsi.
        parser.add_argument("--debug", action='store_true', default=False,
                            help='''If defined the debug mode will be activated,
                            therefore intermediate files will not be deleted.''')
        parser.add_argument("-i", "--inputheaders", type=str,
                            help='''Specify the text file listing the image header files for the images to be processed.''')
        # Define the argument for specifying the sensor.
        parser.add_argument("-s", "--sensor", choices=ARCSI_SENSORS_LIST,
                            help='''Specify the sensor being processed.''')
        # Define the argument for requesting a list of the supported sensors.
        parser.add_argument("--sensorlist", action='store_true', default=False,
                            help='''List the sensors which are supported
                                    and the require names.''')
        # Define the argument for specifying the file path of the output images.
        parser.add_argument("-o", "--outpath", type=str,
                            help='''Specify the output file path.''')
        # Define the argument for specifying the file path of the output images.
        parser.add_argument("-t", "--tmpath", type=str,
                            help='''Specify a tempory path for files to be written to temporarly 
                                    during processing if required (DDVAOT, DOS and CLOUDS).''')
        # Define the argument which specifies the products which are to be generated.
        parser.add_argument("-p", "--prods", type=str, nargs='+', choices=ARCSI_PRODUCTS_LIST,
                            help='''Specify the output products which are to be
                            calculated, as a comma separated list.''')
        # Define the argument for requesting a list of products.
        parser.add_argument("--prodlist", action='store_true', default=False,
                            help='''List the products which are supported and
                            their input requirements.''')
        # Define the argument for specifying that statistics and pyramids should be built for
        # all output images.
        parser.add_argument("--stats", action='store_true', default=False,
                            help='''Specifies that the image statistics and
                            pyramids should be build for all output images.''')
        parser.add_argument("-d", "--dem", type=str,
                            help='''Specify a DEM which is to be used for building
                            an LUT and applying 6S coefficients with respect to elevation.''')
        parser.add_argument("--demnodata", type=float,
                            help='''Specify a no data value for the input DEM image file.''')
        # Define the argument for specifying the output image base file name if it is
        # not to be automatically generated.
        parser.add_argument("--outbasename", type=str,
                            help='''Specify the output file base name if it
                            is not to be generated by the system.''')
        # Define the argument for specifying the image file format.
        parser.add_argument("-f", "--format", type=str, choices=ARCSI_GDALFORMATS_LIST, default='KEA',
                            help='''Specify the image output format (Note. Current just the KEA file format is supported, 
                            use gdal_translate to convert to other formats (e.g., GeoTIFF) following completion.).''')
        # Define the argument stating that alongsided the masked products should none masked products.
        parser.add_argument("--fullimgouts", action='store_true', default=False,
                            help='''If set then alongside the masked outputs (e.g., clouds) then SREF (DOS and/or modelled) 
                                    versions of the full images (i.e., without mask applied) will also be outputted.''')
        # Define the argument for specifying the WKT projection file for the input file.
        parser.add_argument("--inwkt", type=str,
                            help='''Specify the WKT projection of the input image with projection defined with WKT file.''')
        # Define the argument for specifying the WKT projection file for the outputs.
        parser.add_argument("--outwkt", type=str,
                            help='''Transform the outputs to the projection defined with WKT file.''')
        # Define the argument for specifying the proj4 projection string for the outputs.
        parser.add_argument("--outproj4", type=str,
                            help='''Transform the outputs to the projection defined using a proj4 string and provided within a text file.''')
        # Define the argument for string added to the output files names indicating the projection.
        parser.add_argument("--projabbv", type=str,
                            help='''Abbreviation or acronym for the project which will added to the file name.''')
        # Define the argument for the output x pixel resolution (if image re-projected).
        parser.add_argument("--ximgres", type=float,
                            help='''Float for the output image pixel x resolution (if re-projected). 
                                    Optional, if not provided the input image resolution is used.''')
        # Define the argument for the output y pixel resolution (if image re-projected).
        parser.add_argument("--yimgres", type=float,
                            help='''Float for the output image pixel y resolution (if re-projected). 
                                    Optional, if not provided the input image resolution is used.''')
        # Define the argument which specifies the standard aersol profile to use.
        parser.add_argument("--aeropro", type=str, choices=['NoAerosols', 'Continental',
                            'Maritime', 'Urban', 'Desert', 'BiomassBurning', 'Stratospheric'],
                            help='''Specify the 6S defined aersol profile to use.
                            (NoAerosols, Continental, Maritime, Urban, Desert, BiomassBurning, Stratospheric)''')
        # Define the argument which specifies the standard atompheric profile to use.
        parser.add_argument("--atmospro", type=str, choices=['NoGaseousAbsorption', 'Tropical',
                            'MidlatitudeSummer', 'MidlatitudeWinter', 'SubarcticSummer', 'SubarcticWinter', 
                            'USStandard1962'],
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
                                    estimate for defining --minaot within search space.''')
                            # Define the argument for specifying the AOT value for the scene
        parser.add_argument("--upaot", type=float, default=0.4,
                            help='''Specify the upper AOT amount to be added to the AOT
                                    estimate for defining --maxaot within search space.''')
        # Define the argument for specifying the AOT image file for the scene
        parser.add_argument("--aotfile", type=str,
                            help='''Specifiy an image file with AOT values for the
                                    correction. An LUT for AOT and elevation will be generated.
                                    Therefore, --dem needs to be provided alongside --aotfile.''')
        parser.add_argument("--localdos", action='store_true', default=False,
                            help='''Specifies that a local DOS should be applied
                            rather than a global DOS.''')
        parser.add_argument("--simpledos", action='store_true', default=False,
                            help='''Specifies that a simple (basic) DOS should be applied
                            rather than the more complex variable global/local DOS methods.''')
        parser.add_argument("--dosout", type=float, default=20,
                            help='''Specifies the reflectance value to which dark objects
                            are set to during the dark object subtraction. The default is 20,
                            which is equivalent to 2 percent reflectance.''')
        parser.add_argument("--scalefac", type=int, default=1000,
                            help='''Specifies the scale factor for the reflectance
                            products.''')
        
        parser.add_argument("--interp", type=str, default="cubic",
                            choices=['near', 'bilinear', 'cubic', 'cubicspline', 'lanczos', 'average', 'mode', 'max', 'min', 'med'],
                            help='''Specifies interpolation algorithm when reprojecting the imagery
                                    (Note. the options are those in gdalwarp).''')
        
        parser.add_argument("--interpresamp", type=str, default="near",
                            choices=['near', 'bilinear', 'cubic', 'cubicspline', 'lanczos', 'average'],
                            help='''Specifies interpolation algorithm when resampling image bands to a new resolution (e.g., Sentinel-2)
                                    (Note. the options are those in gdalwarp).''')

        parser.add_argument("--cs_initdist", type=int, default=3000,
                                         help='''When clear-sky regions are being defined this parameter
                                                   is the initial distance (m) from cloud objects to generate the initial
                                                   clear-sky regions.''')
        parser.add_argument("--cs_initminsize", type=int, default=3000,
                                         help='''When clear-sky regions are being defined this parameter
                                                   is the minimum size (in pixels) of the initial objects.''')
        parser.add_argument("--cs_finaldist", type=int, default=1000,
                                         help='''When clear-sky regions are being defined this parameter
                                                   is final distance (m) from the cloud objects defining clear sky
                                                   regions.''')
        parser.add_argument("--cs_morphop", type=int, default=21,
                                         help='''When clear-sky regions are being defined this parameter
                                                   is the size of the morphological opening operator used
                                                   to finalise the result.''')
        parser.add_argument("--checkouts", action='store_true', default=False,
                            help='''Specifies that the output path should be checked for files with the same base name.
                            If a file with the same base name is found then processing will not proceed - i.e., files will
                            not be overwritten.''')
        # Define the argument for requesting a list of the available environment variables.
        parser.add_argument("--envvars", action='store_true', default=False,
                            help='''List the available environmental variables for ARCSI.''')

        parser.add_argument("--classmlclouds", action='store_true', default=False,
                            help='''Specifies that the generic machine learning based clouds classification process
                            should be used. Note. --cloudtrainclouds and --cloudtrainother need to be specified if this 
                            option is used. ''')
        # Define the argument for specifying a HDF5 file with training samples for classifying clouds using the generic cloud 
        # masking process using the extra random forest algorithm.
        parser.add_argument("--cloudtrainclouds", type=str,
                            help='''Specify a hdf5 file with the training for classifying clouds.''')
        # Define the argument for specifying a HDF5 file with training samples for classifying non-clouds using the generic cloud 
        # masking process using the extra random forest algorithm.
        parser.add_argument("--cloudtrainother", type=str,
                            help='''Specify a hdf5 file with the training for classifying non-clouds''')
        parser.add_argument("--resample2lowres", action='store_true', default=False, 
                            help='''If image data is provided at multiple image spatial resolutions then this
                                    switch specifies that the higher resolution images should be resampled to the 
                                    same resolution as the lower resolution images (Default: lower resolution are 
                                    resampled to the higher resolution). Example, using this switch will mean Sentinel-2
                                    imagery outputted at 20m rather than 10m resolution.''')
        # Provide a list of file ends which are to be kept, all others will be deleted from the output directory.
        parser.add_argument("-k", "--keepfileends", type=str, nargs='+', default=None,
                            help='''Provide a list of file endings which are to be kept following the completion of the processing.''')

        # Call the parser to parse the arguments.
        args = parser.parse_args()

        arcsiUtils = ARCSIUtils()

        if args.sensorlist:
            arcsilib.arcsirun.print2ConsoleListSensors()
        elif args.prodlist:
            arcsilib.arcsirun.print2ConsoleListProductDescription()
        elif args.envvars:
            arcsilib.arcsirun.print2ConsoleListEnvVars()
        else:
            # Check that the input header parameter has been specified.
            if args.inputheaders == None:
                print("Error: No list of input header files has been provided.\n")
                sys.exit()

            # Check that the senor parameter has been specified.
            if args.sensor == None:
                print("Error: No sensor has been provided.\n")
                sys.exit()

            # Check that the output image path has been specified.
            if args.outpath == None:
                # Print an error message if not and exit.
                envVar = arcsiUtils.getEnvironmentVariable("ARCSI_OUTPUT_PATH")
                if envVar == None:
                    print("Error: No output file path has been provided.\n")
                    sys.exit()
                else:
                    print("Taking output file path from environment variable.")
                    args.outpath = envVar

            if not os.path.exists(args.outpath):
                print("WARNING: Output directory does not exist so creating it...")
                os.makedirs(args.outpath)
            elif not os.path.isdir(args.outpath):
                print("ERROR: Output Path exists but is not a directory...\n")
                sys.exit()

            if not args.outwkt is None:
                if not os.path.exists(args.outwkt):
                    print("Error: The output WKT file does not exist.\n")
                    sys.exit()
                elif args.projabbv == None:
                    print("WARNING: It is recommended that a projection abbreviation or acronym is provided (--projabbv)...")

            if not args.outproj4 is None:
                if not os.path.exists(args.outproj4):
                    print("Error: The output proj4 file does not exist.\n")
                    sys.exit()
                elif args.projabbv == None:
                    print("WARNING: It is recommended that a projection abbreviation or acronym is provided (--projabbv)...")

            if args.classmlclouds:
                if (args.cloudtrainclouds == None) or (args.cloudtrainother == None):
                    print("Error: If --classmlclouds if specified then --cloudtrainclouds and --cloudtrainother are also required.\n")
                    sys.exit()

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
                    sys.exit()
                else:
                    print("Taking min AOT from environment variable.")
                    args.minaot = float(envVarMinAOT)

                envVarMaxAOT = arcsiUtils.getEnvironmentVariable("ARCSI_MAX_AOT")
                if envVarMaxAOT == None:
                    print("Error: The min and max AOT values for the search should be specified.\n")
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
                sys.exit()

            if needTmp and args.tmpath == None:
                envVar = arcsiUtils.getEnvironmentVariable("ARCSI_TMP_PATH")
                if envVar == None:
                    print("Error: If the DDVAOT, DOS, DOSAOTSGL, CLOUDS or TOPOSHADOW product is set then a tempory path needs to be provided.\n")
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
                    sys.exit()

            if args.dem == None:
                envVar = arcsiUtils.getEnvironmentVariable("ARCSI_DEM_PATH")
                if not envVar == None:
                    args.dem = envVar
                    print("Taking DEM path from environment variable.")

            if needDEM:
                if (args.dem == None) or (not os.path.exists(args.dem)):
                    print("Error: A file path to a DEM has either not been specified or does exist, please check it and run again.\n")
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
                    print("Taking output WKT file from environment variable.")

            if args.outproj4 == None:
                envVar = arcsiUtils.getEnvironmentVariable("ARCSI_OUTPUT_PROJ4")
                if not envVar == None:
                    args.outproj4 = envVar
                    print("Taking output proj4 file from environment variable.")


            if args.projabbv == None:
                envVar = arcsiUtils.getEnvironmentVariable("ARCSI_PROJ_ABBV")
                if not envVar == None:
                    args.projabbv = envVar
                    print("Taking projection abbreviation from environment variable.")

            atmosOZoneWaterSpecified = False
            if (not args.atmosozone == None) and (args.atmoswater == None):
                print("Error: If the atmospheric ozone is defined then the atmospheric water needs to be specfied --atmoswater.\n")
                sys.exit()
            elif (not args.atmoswater == None) and (args.atmosozone == None):
                print("Error: If the atmospheric water is defined then the atmospheric ozone needs to be specfied --atmosozone.\n")
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

            envVar = arcsiUtils.getEnvironmentVariable("ARCSI_SCALE_FACTOR")
            if not envVar == None:
                args.scalefac = int(envVar)

            runTimer = rsgislib.RSGISTime()
            runTimer.start(True)
            
            try:
                rsgisUtils = rsgislib.RSGISPyUtils()
                inputHeadersLst = rsgisUtils.readTextFile2List(args.inputheaders)
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
                    paramsObj = prepParametersObj(inputHeader, None, None, args.sensor, args.inwkt, args.format, args.outpath, args.outbasename, args.outwkt, args.outproj4, args.projabbv, args.ximgres, args.yimgres, args.prods, args.stats, args.aeropro, args.atmospro, args.aeroimg, args.atmosimg, args.grdrefl, args.surfacealtitude, args.atmosozone, args.atmoswater, atmosOZoneWaterSpecified, args.aerowater, args.aerodust, args.aerooceanic, args.aerosoot, aeroComponentsSpecified, args.aot, args.vis, args.tmpath, args.minaot, args.maxaot, args.lowaot, args.upaot, args.dem, args.demnodata, args.aotfile, (not args.localdos), args.dosout, args.simpledos, args.debug, args.scalefac, args.interp, args.interpresamp, args.cs_initdist, args.cs_initminsize, args.cs_finaldist, args.cs_morphop, args.fullimgouts, args.checkouts, args.classmlclouds, args.cloudtrainclouds, args.cloudtrainother, args.resample2lowres, args.keepfileends)
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

                paramsLstTmp = []
                nTasks = len(paramsLst)
                taskIdx = 0
                nWorkers = mpiSize - 1
                completedTasks = 0
                while completedTasks < nTasks:
                    rtnParamsObj = mpiComm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=mpiStatus)
                    source = mpiStatus.Get_source()
                    tag = mpiStatus.Get_tag()
                    if tag == mpiTags.READY:
                        # Worker is ready, so send it a task
                        if taskIdx < nTasks:
                            mpiComm.send([arcsiStages.ARCSIPART1, paramsLst[taskIdx]], dest=source, tag=mpiTags.START)
                            print("Sending task %d to worker %d" % (taskIdx, source))
                            taskIdx += 1
                        #else: Do nothing
                    elif tag == mpiTags.DONE:
                        print("Got data from worker %d" % source)
                        paramsLstTmp.append(rtnParamsObj)
                        ++completedTasks
                    elif tag == tags.EXIT:
                        raise ARCSIException("MPI worker was closed - worker was still needed so there is a bug here somewhere... Please report to mailing list.")
                paramsLst = paramsLstTmp
                
                if calcAOT:
                    if useAOTImage:
                        raise ARCSIException("Currently the --multi option does not support the merging of AOT images (i.e., from DDVAOT and DOSAOT) across multiple scenes.")
                    else:
                        # Replace the AOT value with the mean from all the scenes.
                        aotSum = 0.0
                        aotN = 0.0
                        for paramsObj in paramsLst:
                            aotSum = aotSum + paramsObj.aotVal
                            aotN = aotN + 1
                        avgAOT = aotSum / aotN
                        for params in paramsLst:
                            paramsObj.aotVal = avgAOT
                
                """
                if calc6SSREF:
                    paramsLst = plObj.map(_runARCSIPart2, paramsLst)

                if exportMetaData:
                    paramsLst = plObj.map(_runARCSIPart3, paramsLst)

                paramsLst = plObj.map(_runARCSIPart4, paramsLst)
                """

                paramsLstTmp = list()
                nTasks = len(paramsLst)
                taskIdx = 0
                nWorkers = mpiSize - 1
                closedWorkers = 0
                while closedWorkers < nWorkers:
                    rtnParamsObj = mpiComm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=mpiStatus)
                    source = mpiStatus.Get_source()
                    tag = mpiStatus.Get_tag()
                    if tag == mpiTags.READY:
                        # Worker is ready, so send it a task
                        if taskIdx < nTasks:
                            mpiComm.send([arcsiStages.ARCSIPART4, paramsLst[taskIdx]], dest=source, tag=mpiTags.START)
                            print("Sending task %d to worker %d" % (taskIdx, source))
                            taskIdx += 1
                        else:
                            mpiComm.send(None, dest=source, tag=mpiTags.EXIT)
                    elif tag == mpiTags.DONE:
                        print("Got data from worker %d" % source)
                        paramsLstTmp.append(rtnParamsObj)
                    elif tag == tags.EXIT:
                        print("Worker %d exited." % source)
                        closedWorkers += 1

            except ARCSIException as e:
                print("Error: {}".format(e), file=sys.stderr)
                if args.debug:
                    raise
            except Exception as e:
                print("Error: {}".format(e), file=sys.stderr)
                if args.debug:
                    raise

            runTimer.end(True, "ARCSI took ", " to process. Thank you for using ARCSI.")
            print("\n\n")
else:
    # Worker processes execute code below
    while True:
        mpiComm.send(None, dest=0, tag=mpiTags.READY)
        tskData = mpiComm.recv(source=0, tag=MPI.ANY_TAG, status=mpiStatus)
        tag = mpiStatus.Get_tag()
        paramsObj = None

        if tag == mpiTags.START:
            # Do work!
            if tskData[0] == arcsiStages.ARCSIPART1:
                paramsObj = _runARCSIPart1(tskData[1])
            elif tskData[0] == arcsiStages.ARCSIPART2:
                paramsObj = _runARCSIPart2(tskData[1])
            elif tskData[0] == arcsiStages.ARCSIPART2:
                paramsObj = _runARCSIPart3(tskData[1])
            elif tskData[0] == arcsiStages.ARCSIPART4:
                paramsObj = _runARCSIPart4(tskData[1])
            else:
                raise ARCSIException("Don't recognise processing stage")

            mpiComm.send(paramsObj, dest=0, tag=mpiTags.DONE)
        elif tag == mpiTags.EXIT:
            break

    mpiComm.send(None, dest=0, tag=tags.EXIT)
