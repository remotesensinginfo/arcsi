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

import os
import sys
import fnmatch
import argparse
from arcsilib import ARCSI_VERSION
from arcsilib import ARCSI_SENSORS_LIST
from arcsilib import ARCSI_GDALFORMATS_LIST
from arcsilib import ARCSI_PRODUCTS_LIST
from arcsilib import ARCSI_CLOUD_METHODS_LIST
from arcsilib.arcsiexception import ARCSIException


class ARCSIBuildCommands(object):
    def getListOfFiles(self, searchDIR, searchStr, depth):
        inDIRCount = searchDIR.count(os.path.sep)
        outFiles = []
        for dirName, subdirList, fileList in os.walk(searchDIR):
            for fname in fileList:
                fname = str(fname)
                if fnmatch.fnmatch(fname, searchStr):
                    outFiles.append(os.path.abspath(os.path.join(dirName, fname)))
            dirLevel = dirName.count(os.path.sep) - inDIRCount
            if dirLevel >= depth:
                subdirList[:] = []
        return outFiles

    def buildCmds(
        self,
        inputPath,
        inputIsDIR,
        outputFile,
        headerSearchStr,
        searchDepth,
        sensor,
        inwkt,
        format,
        outpath,
        prods,
        stats,
        aeropro,
        atmospro,
        aeroimg,
        atmosimg,
        grdrefl,
        surfacealtitude,
        atmosozone,
        atmoswater,
        aerowater,
        aerodust,
        aerooceanic,
        aerosoot,
        aot,
        vis,
        tmpath,
        minaot,
        maxaot,
        dem,
        dem_no_data_val,
        localdos,
        dosout,
        simpledos,
        scalefac,
        outwkt,
        outproj4,
        projabbv,
        interp,
        interpresamp,
        checkouts,
        fullimgouts,
        cloudmethods,
        classmlclouds,
        cloudtrainclouds,
        cloudtrainother,
        resample2lowres,
        keepfileends,
        multi,
        ncores,
    ):

        inputPath = os.path.abspath(inputPath)
        outputFile = os.path.abspath(outputFile)

        headersFilesList = []
        if inputIsDIR:
            if headerSearchStr.count("*") == 0:
                raise ARCSIException(
                    "The search string you have provided does not have any '*' - which is needed for searching."
                )
            headersFilesList = self.getListOfFiles(
                inputPath, headerSearchStr, searchDepth
            )
        else:

            headersFilesList = rsgislib.tools.utils.read_text_file_no_new_lines2List(
                inputPath
            )

        prodsStr = ""
        first = True
        for prod in prods:
            if first:
                prodsStr = prod
                first = False
            else:
                prodsStr = prodsStr + " " + prod

        outFile = open(outputFile, "w+")
        for hFile in headersFilesList:
            print("Processing :", hFile)

            sensorOUT = sensor
            if sensor == "LANDSAT":
                basefilename = os.path.basename(hFile)
                filePrefix3 = basefilename[:3]
                filePrefix4 = basefilename[:4]

                if filePrefix3 == "LM1" or filePrefix4 == "LM01":
                    sensorOUT = "ls1"
                elif filePrefix3 == "LM2" or filePrefix4 == "LM02":
                    sensorOUT = "ls2"
                elif filePrefix3 == "LM3" or filePrefix4 == "LM03":
                    sensorOUT = "ls3"
                elif filePrefix3 == "LM4" or filePrefix4 == "LM04":
                    sensorOUT = "ls4mss"
                elif filePrefix3 == "LM5" or filePrefix4 == "LM05":
                    sensorOUT = "ls5mss"
                elif (
                    filePrefix3 == "LT4"
                    or filePrefix4 == "LS04"
                    or filePrefix4 == "LE04"
                    or filePrefix4 == "LT04"
                ):
                    sensorOUT = "ls4tm"
                elif (
                    filePrefix3 == "LT5"
                    or filePrefix4 == "LS05"
                    or filePrefix4 == "LE05"
                    or filePrefix4 == "LT05"
                ):
                    sensorOUT = "ls5tm"
                elif (
                    filePrefix3 == "LE7"
                    or filePrefix4 == "LS07"
                    or filePrefix4 == "LE07"
                    or filePrefix4 == "LT07"
                ):
                    sensorOUT = "ls7"
                elif (
                    filePrefix3 == "LC8"
                    or filePrefix4 == "LS08"
                    or filePrefix4 == "LC08"
                ):
                    sensorOUT = "ls8"
                else:
                    raise ARCSIException(
                        'Sensor was not recognised for file: "' + hFile + '"'
                    )

            cmd = "arcsi.py -s " + sensorOUT + " -p " + prodsStr + ' -i "' + hFile + '"'
            if outpath is not None:
                cmd = cmd + ' --outpath "' + os.path.abspath(outpath) + '"'
            if stats:
                cmd = cmd + " --stats"
            if not format == None:
                cmd = cmd + " --format " + format
            else:
                cmd = cmd + " --format KEA"
            if tmpath is not None:
                cmd = cmd + ' --tmpath "' + os.path.abspath(tmpath) + '"'
            if dem is not None:
                cmd = cmd + ' --dem "' + os.path.abspath(dem) + '"'
            if dem_no_data_val is not None:
                cmd = cmd + " --demnodata " + str(dem_no_data_val)
            if aeroimg is not None:
                cmd = cmd + ' --aeroimg "' + os.path.abspath(aeroimg) + '"'
            elif aeropro is not None:
                cmd = cmd + " --aeropro " + aeropro
            if atmosimg is not None:
                cmd = cmd + ' --atmosimg "' + os.path.abspath(atmosimg) + '"'
            elif atmospro is not None:
                cmd = cmd + " --atmospro " + atmospro
            if surfacealtitude is not None:
                cmd = cmd + " --surfacealtitude " + str(surfacealtitude)
            if (minaot is not None) and (maxaot is not None):
                cmd = cmd + " --minaot " + str(minaot) + " --maxaot " + str(maxaot)
            elif aot is not None:
                cmd = cmd + " --aot " + str(aot)
            elif vis is not None:
                cmd = cmd + " --vis " + str(vis)
            if inwkt is not None:
                cmd = cmd + ' --inwkt "' + os.path.abspath(inwkt) + '"'
            if grdrefl is not None:
                cmd = cmd + " --grdrefl " + grdrefl
            if atmosozone is not None:
                cmd = cmd + " --atmosozone " + str(atmosozone)
            if atmoswater is not None:
                cmd = cmd + " --atmoswater " + str(atmoswater)
            if aerowater is not None:
                cmd = cmd + " --aerowater " + str(aerowater)
            if aerodust is not None:
                cmd = cmd + " --aerodust " + str(aerodust)
            if aerooceanic is not None:
                cmd = cmd + " --aerooceanic " + str(aerooceanic)
            if aerosoot is not None:
                cmd = cmd + " --aerosoot " + str(aerosoot)
            if dosout is not None:
                cmd = cmd + " --dosout " + str(dosout)
            if scalefac is not None:
                cmd = cmd + " --scalefac " + str(scalefac)
            if outwkt is not None:
                cmd = cmd + ' --outwkt "' + os.path.abspath(outwkt) + '"'
            if outproj4 is not None:
                cmd = cmd + ' --outproj4 "' + os.path.abspath(outproj4) + '"'
            if projabbv is not None:
                cmd = cmd + " --projabbv " + str(projabbv)
            if interp is not None:
                cmd = cmd + " --interp " + str(interp)
            if interpresamp is not None:
                cmd = cmd + " --interpresamp " + str(interpresamp)
            if localdos:
                cmd = cmd + " --localdos "
            if simpledos:
                cmd = cmd + " --simpledos "
            if checkouts:
                cmd = cmd + " --checkouts "
            if fullimgouts:
                cmd = cmd + " --fullimgouts "
            if cloudmethods is not None:
                cmd = cmd + " --cloudmethods " + str(cloudmethods)
            if classmlclouds:
                cmd = cmd + " --classmlclouds "
                if not cloudtrainclouds == None:
                    cmd = cmd + ' --cloudtrainclouds "' + str(cloudtrainclouds) + '"'
                if not cloudtrainother == None:
                    cmd = cmd + ' --cloudtrainother "' + str(cloudtrainother) + '"'
            if resample2lowres:
                cmd = cmd + " --resample2lowres "
            if keepfileends is not None:
                cmd = cmd + " --keepfileends "
                for fileEnd in keepfileends:
                    cmd = cmd + ' "' + fileEnd + '"'
            if multi:
                cmd = cmd + " --multi --ncores " + str(ncores)

            print(cmd)
            outFile.write(cmd + "\n")
        outFile.flush()
        outFile.close()


if __name__ == "__main__":
    """
    The command line user interface to ARCSI Data Extraction Tool.
    """
    parser = argparse.ArgumentParser(
        prog="arcsibuildcmdslist.py",
        description="""ARCSI command to build arcsi.py commands
                                                for a set of input images using the same options.""",
        epilog="""A tools to build arcsi.py commands
                                           for a set of input images using
                                           the same options""",
    )
    # Request the version number.
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s version " + ARCSI_VERSION
    )

    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="""Input directory containing the data to be processed or text file listing paths to header files (don't set --header if text file)""",
    )

    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="""Output text file (shell script) with the list of commands.""",
    )

    parser.add_argument(
        "-e",
        "--header",
        type=str,
        required=False,
        help="""A \'UNIX\' search string for identifying the image headers. Note, multiple \'*\' can be used for the search string. If no \'*\' is provided then """,
    )

    parser.add_argument(
        "-d",
        "--depth",
        type=int,
        required=False,
        default=1000,
        help="""The depth within the directory tree from the input path which should be searched for image header files.""",
    )

    # Define the argument for specifying the sensor.
    sensorList = ARCSI_SENSORS_LIST
    sensorList.append("LANDSAT")
    parser.add_argument(
        "-s",
        "--sensor",
        required=True,
        choices=ARCSI_SENSORS_LIST,
        help="""Specify the sensor being processed. Note, if the 
                        option 'LANDSAT' is selected then the sensor will 
                        automatically be identified from the file name.""",
    )

    # Define the argument for specifying the WKT projection file
    # for the intput file.
    parser.add_argument(
        "--inwkt", type=str, help="""Specify the WKT projection of the input image"""
    )

    # Define the argument for specifying the image file format.
    parser.add_argument(
        "-f",
        "--format",
        type=str,
        choices=ARCSI_GDALFORMATS_LIST,
        default="KEA",
        help="""Specify the image output format (Note. Current just the KEA file format is supported, 
                        use gdal_translate to convert to other formats (e.g., GeoTIFF) following completion.).""",
    )

    # Define the argument stating that alongsided the masked products should none masked products.
    parser.add_argument(
        "--fullimgouts",
        action="store_true",
        default=False,
        help="""If set then alongside the masked outputs (e.g., clouds) then SREF (DOS and/or modelled) 
                            versions of the full images (i.e., without mask applied) will also be outputted.""",
    )

    parser.add_argument(
        "--outpath",
        type=str,
        help="""Specifiy the output directory for the products from ARCSI.""",
    )

    # Define the argument for specifying the file path of the output images.
    parser.add_argument(
        "-t",
        "--tmpath",
        type=str,
        help="""Specify a tempory path for files to be written to temporarly during processing if required (DDVAOT, DOSUB and CLOUDS).""",
    )

    # Define the argument which specifies the products which are to be generated.
    parser.add_argument(
        "-p",
        "--prods",
        type=str,
        required=True,
        nargs="+",
        choices=ARCSI_PRODUCTS_LIST,
        help="""Specify the output products which are to be
                        calculated, as a comma separated list.""",
    )

    # Define the argument which specifies the standard aersol profile to use.
    parser.add_argument(
        "--aeropro",
        type=str,
        choices=[
            "NoAerosols",
            "Continental",
            "Maritime",
            "Urban",
            "Desert",
            "BiomassBurning",
            "Stratospheric",
        ],
        help="""Specify the 6S defined aersol profile to use.
                        (NoAerosols, Continental, Maritime, Urban, Desert, BiomassBurning, Stratospheric)""",
    )

    # Define the argument which specifies the standard atompheric profile to use.
    parser.add_argument(
        "--atmospro",
        type=str,
        choices=[
            "NoGaseousAbsorption",
            "Tropical",
            "MidlatitudeSummer",
            "MidlatitudeWinter",
            "SubarcticSummer",
            "SubarcticWinter",
            "USStandard1962",
        ],
        help="""Specify the 6S defined atmospheric profile to use.
                        (NoGaseousAbsorption, Tropical, MidlatitudeSummer, MidlatitudeWinter,
                        SubarcticSummer, SubarcticWinter, USStandard1962)""",
    )

    # Define the argument for specifying the file path for the image specifying the generic aerosol model.
    parser.add_argument(
        "--aeroimg", type=str, help="""Specify the aerosol model image file path."""
    )

    # Define the argument for specifying the file path for the image specifying the generic atmosphere model.
    parser.add_argument(
        "--atmosimg", type=str, help="""Specify the atmosphere model image file path."""
    )

    # Define the argument which specifies the amount of OZone in atmosphere
    parser.add_argument(
        "--atmosozone",
        type=float,
        help="""Specify the total amount of ozone in a vertical path
                        through the atmosphere (in cm-atm)""",
    )

    # Define the argument which specifies the amount of water in the atmosphere
    parser.add_argument(
        "--atmoswater",
        type=float,
        help="""Specify the  total amount of water in a vertical path
                        through the atmosphere (in g/cm^2)""",
    )

    # Define the argument which specifies the proportion of water-like aerosols
    parser.add_argument(
        "--aerowater",
        type=float,
        help="""Specify the proportion of water-like aerosols
                        (water, dust, oceanic and soot proportions must add up
                        to 1 although all do not been to be specified).""",
    )

    # Define the argument which specifies the proportion of dust-like aerosols
    parser.add_argument(
        "--aerodust",
        type=float,
        help="""Specify the proportion of dust-like aerosols
                        (water, dust, oceanic and soot proportions must add up
                        to 1 although all do not been to be specified).""",
    )

    # Define the argument which specifies the proportion of oceanic-like aerosols
    parser.add_argument(
        "--aerooceanic",
        type=float,
        help="""Specify the proportion of oceanic-like aerosols
                        (water, dust, oceanic and soot proportions must add up
                        to 1 although all do not been to be specified).""",
    )

    # Define the argument which specifies the proportion of soot-like aerosols
    parser.add_argument(
        "--aerosoot",
        type=float,
        help="""Specify the proportion of soot-like aerosols
                        (water, dust, oceanic and soot proportions must add up
                        to 1 although all do not been to be specified).""",
    )

    # Define the argument which specifies the ground reflectance.
    parser.add_argument(
        "--grdrefl",
        type=str,
        choices=["GreenVegetation", "ClearWater", "Sand", "LakeWater", "BRDFHapke"],
        help="""Specify the ground reflectance used for the
                                                6S model. (GreenVegetation, ClearWater, Sand, LakeWater, BRDFHapke).""",
    )

    # Define the argument for specifying the surface elevation for the scene.
    parser.add_argument(
        "--surfacealtitude",
        type=float,
        help="""Specify the altiude (in km) of the surface being sensed.""",
    )

    # Define the argument for specifying the AOT value for the scene
    parser.add_argument(
        "--aot",
        type=float,
        help="""Specifiy the AOT or visability value for the scene.
                                If the AOT is specified the visability is ignored.""",
    )

    # Define the argument for specifying the visability value for the scene
    parser.add_argument(
        "--vis",
        type=float,
        help="""Specifiy the AOT or visability value for the scene.
                                If the AOT is specified the visability is ignored.""",
    )

    # Define the argument for specifying the AOT value for the scene
    parser.add_argument(
        "--minaot",
        type=float,
        help="""Specifiy the AOT or visibility value for the scene.
                                If the AOT is specified the visibility is ignored.""",
    )

    # Define the argument for specifying the AOT value for the scene
    parser.add_argument(
        "--maxaot",
        type=float,
        help="""Specifiy the AOT or visability value for the scene.
                                If the AOT is specified the visability is ignored.""",
    )

    # Define the argument for specifying that statistics and pyramids should be built for
    # all output images.
    parser.add_argument(
        "--stats",
        action="store_true",
        default=False,
        help="""Specifies that the image statistics and
                        pyramids should be build for all output images.""",
    )

    parser.add_argument(
        "--dem",
        type=str,
        help="""Specify a DEM which is to be used for building
                        an LUT and applying 6S coefficients with respect to elevation.""",
    )

    parser.add_argument(
        "--demnodata",
        type=float,
        help="""Specify a no data value for the input DEM image file.""",
    )

    parser.add_argument(
        "--localdos",
        action="store_true",
        default=False,
        help="""Specifies that a local DOS should be applied
                        rather than a global DOS.""",
    )

    parser.add_argument(
        "--simpledos",
        action="store_true",
        default=False,
        help="""Specifies that a simple (basic) DOS should be applied
                        rather than the more complex variable global/local DOS methods.""",
    )

    parser.add_argument(
        "--dosout",
        type=float,
        help="""Specifies the reflectance value to which dark objects
                        are set to during the dark object subtraction. (Default is 20,
                        which is equivalent to 2 % reflectance.""",
    )

    parser.add_argument(
        "--scalefac",
        type=int,
        help="""Specifies the scale factor for the reflectance
                        products.""",
    )

    parser.add_argument(
        "--outwkt",
        type=str,
        help="""Transform the outputs to the projection defined with WKT file.""",
    )

    parser.add_argument(
        "--projabbv",
        type=str,
        help="""Abbreviation or acronym for the project which will added to the file name.""",
    )

    # Define the argument for specifying the proj4 projection string for the outputs.
    parser.add_argument(
        "--outproj4",
        type=str,
        help="""Transform the outputs to the projection defined using a proj4 string and provided within a text file.""",
    )

    parser.add_argument(
        "--interp",
        type=str,
        default="cubic",
        choices=[
            "near",
            "bilinear",
            "cubic",
            "cubicspline",
            "lanczos",
            "average",
            "mode",
            "max",
            "min",
            "med",
        ],
        help="""Specifies interpolation algorithm when reprojecting the imagery
                                (Note. the options are those in gdalwarp).""",
    )

    parser.add_argument(
        "--interpresamp",
        type=str,
        default="near",
        choices=["near", "bilinear", "cubic", "cubicspline", "lanczos", "average"],
        help="""Specifies interpolation algorithm when resampling image bands to a new resolution (e.g., Sentinel-2)
                                (Note. the options are those in gdalwarp).""",
    )

    parser.add_argument(
        "--checkouts",
        action="store_true",
        default=False,
        help="""Specifies that the output path should be checked for files with the same base name.
                    If a file with the same base name is found then processing will not proceed - i.e., files will
                    not be overwritten.""",
    )

    parser.add_argument(
        "--cloudmethods",
        type=str,
        default=None,
        choices=ARCSI_CLOUD_METHODS_LIST,
        help="""Specify the method(s) of cloud masking. Current Sentinel-2 and Landsat have options).
                            Sentinel-2: FMASK, FMASK_DISP or S2CLOUDLESS. Landsat: FMASK or LSMSK, FMASK is current
                            the default for both.""",
    )

    parser.add_argument(
        "--classmlclouds",
        action="store_true",
        default=False,
        help="""Specifies that the generic machine learning based clouds classification process
                        should be used. Note. --cloudtrainclouds and --cloudtrainother need to be specified if this 
                        option is used. """,
    )

    parser.add_argument(
        "--cloudtrainclouds",
        type=str,
        help="""Specify a hdf5 file with the training for classifying clouds.""",
    )

    parser.add_argument(
        "--cloudtrainother",
        type=str,
        help="""Specify a hdf5 file with the training for classifying non-clouds""",
    )

    parser.add_argument(
        "--resample2lowres",
        action="store_true",
        default=False,
        help="""If image data is provided at multiple image spatial resolutions then this
                                switch specifies that the higher resolution images should be resampled to the 
                                same resolution as the lower resolution images (Default: lower resolution are 
                                resampled to the higher resolution). Example, using this switch will mean Sentinel-2
                                imagery outputted at 20m rather than 10m resolution.""",
    )

    parser.add_argument(
        "-k",
        "--keepfileends",
        type=str,
        nargs="+",
        default=None,
        help="""Provide a list of file endings which are to be kept following the completion of the processing.""",
    )

    parser.add_argument(
        "--ncores",
        type=int,
        default=1,
        help="""Number of cores available for processing when using the --multi option.
                                If a value of -1 is provided then all available cores will be used.""",
    )

    parser.add_argument(
        "--multi",
        action="store_true",
        default=False,
        help="""If defined the multiple mode will be activated,
                        therefore --inputheader needs to be a text file listing
                        a series of input image header files.""",
    )

    # Call the parser to parse the arguments.
    args = parser.parse_args()

    inputIsDIR = True
    if args.header == None:
        print(
            "Header ending not specified therefore expecting input to be a text file listing header files."
        )
        inputIsDIR = False
        if not os.path.isfile(args.input):
            print("Error: Input is expected to be a file but is not, please check.")
            sys.exit()

    arcsiObj = ARCSIBuildCommands()

    arcsiObj.buildCmds(
        args.input,
        inputIsDIR,
        args.output,
        args.header,
        args.depth,
        args.sensor,
        args.inwkt,
        args.format,
        args.outpath,
        args.prods,
        args.stats,
        args.aeropro,
        args.atmospro,
        args.aeroimg,
        args.atmosimg,
        args.grdrefl,
        args.surfacealtitude,
        args.atmosozone,
        args.atmoswater,
        args.aerowater,
        args.aerodust,
        args.aerooceanic,
        args.aerosoot,
        args.aot,
        args.vis,
        args.tmpath,
        args.minaot,
        args.maxaot,
        args.dem,
        args.demnodata,
        args.localdos,
        args.dosout,
        args.simpledos,
        args.scalefac,
        args.outwkt,
        args.outproj4,
        args.projabbv,
        args.interp,
        args.interpresamp,
        args.checkouts,
        args.fullimgouts,
        args.cloudmethods,
        args.classmlclouds,
        args.cloudtrainclouds,
        args.cloudtrainother,
        args.resample2lowres,
        args.keepfileends,
        args.multi,
        args.ncores,
    )
