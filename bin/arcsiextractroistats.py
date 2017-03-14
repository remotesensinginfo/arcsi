#!/usr/bin/env python


############################################################################
#  arcsiextractroistats.py
#
#  Copyright 2014 ARCSI.
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
# Purpose:  A class extracts data for an ROI for a set of ARCSI outputted
#           image files given an input shapefile.
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 05/11/2014
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
# Import the glob module.
import glob
# Import python system library
import sys
#Import the python file paths module
import os.path
# Import the datetime module
import datetime
# Import the python Argument parser
import argparse
# Import the ARCSI exception class
from arcsilib.arcsiexception import ARCSIException
# Import the arcsi version number
from arcsilib import ARCSI_VERSION
# Import the RSGISLib Zonal Stats module
from rsgislib import zonalstats


class ARCSIExtractROIStats (object):
    """
    A class extracts data for an ROI for a set of ARCSI outputted image files
    given an input shapefile.
    """

    def parseDateFromFileName(self, imageFilePath):
        try:
            imageFile = os.path.basename(imageFilePath)
            tokens = imageFile.split("_")
            dateStr = tokens[1]
            dateObj = datetime.datetime.strptime(dateStr, "%Y%m%d")
            return dateObj
        except ARCSIException as e:
            raise e
        except Exception as e:
            raise e

    def importZonalStatsData(self, statsFile):
        try:
            statsData = []
            stats = open(statsFile, 'r')
            row = 0
            for statsRow in stats:
                statsRow = statsRow.strip()
                if row > 0:
                    data = statsRow.split(',')
                    statsData.append(data)
                row+=1
            stats.close()
            return statsData
        except ARCSIException as e:
            raise e
        except Exception as e:
            raise e

    def importZonalStatsFields(self, statsFile):
        try:
            statsData = []
            stats = open(statsFile, 'r')
            row = 0
            for statsRow in stats:
                statsRow = statsRow.strip()
                if row == 0:
                    statsData = statsRow.split(',')
                    break
            stats.close()
            return statsData
        except ARCSIException as e:
            raise e
        except Exception as e:
            raise e

    def exportStats2TextFile(self, outputFilePath, headerFields, statsData):
        try:
            outputFile = open(outputFilePath, 'w')
            rowStr = ''
            first = True
            for field in headerFields:
                if first:
                    rowStr = str(field)
                    first = False
                else:
                    rowStr = rowStr + ", " + str(field)
            rowStr = rowStr + "\n"
            outputFile.write(rowStr)

            for data in statsData:
                first = True
                for field in data:
                    if first:
                        rowStr = field.strftime("%Y%m%d")
                        first = False
                    else:
                        rowStr = rowStr + ", " + str(field)
                rowStr = rowStr + "\n"
                outputFile.write(rowStr)

        except ARCSIException as e:
            raise e
        except Exception as e:
            raise e

    def extractSingleFileStats(self, tmpOutFile, roiFile, imageFile, pixelIntersectMethod):
        try:
            zonalattributes = zonalstats.ZonalAttributes(minThreshold=0, maxThreshold=10000, calcCount=False, calcMin=True, calcMax=True, calcMean=True, calcStdDev=True, calcMode=False, calcSum=False)
            zonalstats.pixelStats2TXT(imageFile, roiFile, tmpOutFile, zonalattributes, True, True, zonalstats.METHOD_POLYCONTAINSPIXELCENTER, False)
        except ARCSIException as e:
            raise e
        except Exception as e:
            raise e

    def sortByDate(self, dateVals, statsData):
        try:
            dateVals.sort()
            statsDataSorted = []
            for dateVal in dateVals:
                for statVals in statsData:
                    if statVals[0].strftime("%Y%m%d") == dateVal.strftime("%Y%m%d"):
                        statsDataSorted.append(statVals)

            return statsDataSorted
        except ARCSIException as e:
            raise e
        except Exception as e:
            raise e

    def extractImageFileStats(self, inputImagesLoc, outputFile, roiFile):
        try:
            inputImages = glob.glob(inputImagesLoc)
            if len(inputImages) == 0:
                raise ARCSIException("No input images were found.")
            statsData2Export = []
            statsData2ExportFields = []
            dateVals = []
            numFields = 0
            first = True
            for imageFile in inputImages:
                print("Processing " + imageFile)
                if not os.path.exists(imageFile):
                    raise ARCSIException("Image file does not exist.")
                self.extractSingleFileStats(outputFile, roiFile, imageFile, zonalstats.METHOD_POLYCONTAINSPIXELCENTER)
                imageDate = self.parseDateFromFileName(imageFile)
                dateVals.append(imageDate)
                if first:
                    statsData2ExportFields = self.importZonalStatsFields(outputFile)
                    numFields = len(statsData2ExportFields)
                    statsData2ExportFields.insert(0, "Date")
                    first = False
                statsData = self.importZonalStatsData(outputFile)
                for stats in statsData:
                    if len(stats) != numFields:
                        raise ARCSIException("The number of fields is incorrect.")
                    stats.insert(0, imageDate)
                    statsData2Export.append(stats)
                print("Completed image " + imageDate.strftime("%Y%m%d"))
            print("Sorting data into date order")
            statsData2Export = self.sortByDate(dateVals, statsData2Export)
            print("Exporting Data")
            self.exportStats2TextFile(outputFile, statsData2ExportFields, statsData2Export)
            print("Completed Processing.")
        except ARCSIException as e:
            raise e
        except Exception as e:
            raise e


if __name__ == '__main__':
    """
    The command line user interface to ARCSI Extract ROI Stats tool.
    """
    parser = argparse.ArgumentParser(prog='arcsiextractroistats.py',
                                    description='''ARCSI command for extracting
                                                   image pixels stats.''',
                                    epilog='''A tool to extracts data for an ROI for
                                              a set of ARCSI outputted image files
                                              given an input shapefile.''')
    # Request the version number.
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version ' + ARCSI_VERSION)

    parser.add_argument("-i", "--input", type=str,
                        help='''Specify the directory containing the input image files
                                and files selected from the input directory, glob
                                supported wild characters .''')

    parser.add_argument("-r", "--roi", type=str,
                        help='''String to specify the region of interest (ROI). Must be
                        a shapefile of the image projection as the input images.''')

    parser.add_argument("-o", "--output", type=str,
                        help='''An output text file with the zonal stats results.''')




    # Call the parser to parse the arguments.
    args = parser.parse_args()

    if args.input == None:
        print("Input file information was not specified.")
        parser.print_help()
        sys.exit()

    if args.output == None:
        print("An output file was not specified.")
        parser.print_help()
        sys.exit()

    if args.roi == None:
        print("An ROI shapefile was not specified.")
        parser.print_help()
        sys.exit()

    arcsiObj = ARCSIExtractROIStats()
    arcsiObj.extractImageFileStats(args.input, args.output, args.roi)











