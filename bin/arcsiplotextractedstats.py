#!/usr/bin/env python


############################################################################
#  arcsiplotextractedstats.py
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
# Purpose:  A class to plot the statistics extracted using the
#           arcsiextractroistats.py command.
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

# Import the future functionality (for Python 2)
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals
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
# Import matplotlib module
import matplotlib.pyplot as plt
import matplotlib.dates as matdate
# Import numpy
import numpy

class ARCSIPlotExtractedStats (object):
    """
    A class to plot the statistics extracted using the
    arcsiextractroistats.py command.
    """

    def parseInputFile(self, inputStatsFile, featureID):
        try:
            statsData = []
            fieldNames = []
            inStats = open(inputStatsFile, 'r')
            row = 0
            for statsRow in inStats:
                statsRow = statsRow.strip()
                if row == 0:
                    fieldNamesTmp = statsRow.split(',')
                    for field in fieldNamesTmp:
                        fieldNames.append(field.strip())
                else:
                    rowStrData = statsRow.split(',')
                    fieldIdx = 0
                    rowData = []
                    for data in rowStrData:
                        data = data.strip()
                        if fieldIdx == 0:
                            rowData.append(datetime.datetime.strptime(data, "%Y%m%d"))
                        elif fieldIdx == 1:
                            rowData.append(int(data))
                        elif fieldIdx == (len(rowStrData)-1):
                            rowData.append(int(data))
                        else:
                            rowData.append(float(data))
                        fieldIdx += 1
                    if rowData[1] == featureID:
                        statsData.append(rowData)
                row+=1

            inStats.close()
            return statsData, fieldNames
        except ARCSIException as e:
            raise e
        except Exception as e:
            raise e

    def generatePlotStdDevRegion(self, datesVals, minVals, maxVals, avgVals, stdVals, outputPlotFile, fieldBase, plotTitle, nStdDev):
        try:
            overallMean = numpy.mean(avgVals)
            overallStdDev = numpy.std(avgVals)
            overallMeanStdDev = numpy.mean(stdVals)

            print("Mean = ", overallMean)
            print("Std Dev of Mean = ", overallStdDev)
            print("Mean Std Dev = ", overallMeanStdDev)

            low1StdDev = overallMean - (overallStdDev * nStdDev)
            upp1StdDev = overallMean + (overallStdDev * nStdDev)

            #print("Lower Bound = ", low1StdDev)
            #print("Upper Bound = ", upp1StdDev)

            lowerErrVals = []
            upperErrVals = []

            lowerVals = []
            upperVals = []

            for i in range(len(datesVals)):
                lowVal = avgVals[i] - stdVals[i]
                if lowVal < minVals[i]:
                    lowVal = minVals[i]
                lowerErrVals.append(avgVals[i]-lowVal)
                lowerVals.append(low1StdDev)

                upVal = avgVals[i] + stdVals[i]
                if upVal > maxVals[i]:
                    upVal = maxVals[i]
                upperErrVals.append(upVal-avgVals[i])
                upperVals.append(upp1StdDev)

            asymmetric_error = [lowerErrVals, upperErrVals]
            dateFlt = matdate.date2num(datesVals)

            fig = plt.figure(figsize=(10, 5), dpi=80)
            ax1 = fig.add_subplot(111)

            ax1.plot_date(datesVals, avgVals, 'k-', label='Mean', zorder=10)
            ax1.errorbar(dateFlt, avgVals, yerr=asymmetric_error, ecolor='k', label='1 Sample Std. Dev.')
            ax1.fill_between(datesVals, lowerVals, upperVals, alpha=0.2, linewidth=1.0, facecolor=[0.70,0.70,0.70], edgecolor=[0.70,0.70,0.70], label='1 Mean Std. Dev.', zorder=-1)

            ax1Range = ax1.axis('tight')

            plt.grid(color='k', linestyle='--', linewidth=0.5)
            plt.title(plotTitle)
            plt.xlabel("Date")
            ax1.set_ylabel(fieldBase)

            plt.savefig(outputPlotFile, format='PDF')

        except ARCSIException as e:
            raise e
        except Exception as e:
            raise e


    def generatePlotSimple(self, datesVals, minVals, maxVals, avgVals, stdVals, outputPlotFile, fieldBase, plotTitle):
        try:
            lowerErrVals = []
            upperErrVals = []

            for i in range(len(datesVals)):
                lowVal = avgVals[i] - stdVals[i]
                if lowVal < minVals[i]:
                    lowVal = minVals[i]
                lowerErrVals.append(avgVals[i]-lowVal)

                upVal = avgVals[i] + stdVals[i]
                if upVal > maxVals[i]:
                    upVal = maxVals[i]
                upperErrVals.append(upVal-avgVals[i])

            asymmetric_error = [lowerErrVals, upperErrVals]

            dateFlt = matdate.date2num(datesVals)

            fig = plt.figure(figsize=(10, 5), dpi=80)
            ax1 = fig.add_subplot(111)

            ax1.plot_date(datesVals, avgVals, 'k-', label='Mean', zorder=10)
            ax1.errorbar(dateFlt, avgVals, yerr=asymmetric_error, ecolor='k')

            ax1Range = ax1.axis('tight')

            plt.grid(color='k', linestyle='--', linewidth=0.5)
            plt.title(plotTitle)
            plt.xlabel("Date")
            ax1.set_ylabel(fieldBase)

            plt.savefig(outputPlotFile, format='PDF')

        except ARCSIException as e:
            raise e
        except Exception as e:
            raise e



    def processAndPlotData(self, inputStatsFile, outputPlotFile, fieldBase, plotTitle, simplePlot, featureID):
        try:
            print("Processing Plot \'" + plotTitle + "\'")
            fieldBase = fieldBase.strip()
            statsData, fieldNames = self.parseInputFile(inputStatsFile, featureID)
            #print(statsData)

            minField = fieldBase + "Min"
            maxField = fieldBase + "Max"
            avgField = fieldBase + "Avg"
            stdField = fieldBase + "Std"
            dateField = "Date"

            minFieldIdx = -1
            maxFieldIdx = -1
            avgFieldIdx = -1
            stdFieldIdx = -1
            dateFieldIdx = -1

            idx = 0
            for field in fieldNames:
                if field == minField:
                    minFieldIdx = idx
                elif field  == maxField:
                    maxFieldIdx = idx
                elif field == avgField:
                    avgFieldIdx = idx
                elif field == stdField:
                    stdFieldIdx = idx
                elif field == dateField:
                    dateFieldIdx = idx
                idx += 1

            #print("Date Idx = ", dateFieldIdx)
            #print("Min Idx = ", minFieldIdx)
            #print("Max Idx = ", maxFieldIdx)
            #print("Avg Idx = ", avgFieldIdx)
            #print("Std Idx = ", stdFieldIdx)

            if dateFieldIdx == -1:
                raise ARCSIException("Date Field was not found.")
            if minFieldIdx == -1:
                raise ARCSIException("Min Field was not found.")
            if maxFieldIdx == -1:
                raise ARCSIException("Max Field was not found.")
            if avgFieldIdx == -1:
                raise ARCSIException("Avg Field was not found.")
            if stdFieldIdx == -1:
                raise ARCSIException("Std Field was not found.")

            datesVals = []
            minVals = []
            maxVals = []
            avgVals = []
            stdVals = []

            for dataRow in statsData:
                datesVals.append(dataRow[dateFieldIdx])
                minVals.append(dataRow[minFieldIdx])
                maxVals.append(dataRow[maxFieldIdx])
                avgVals.append(dataRow[avgFieldIdx])
                stdVals.append(dataRow[stdFieldIdx])

            #print(datesVals)
            #print(minVals)
            #print(maxVals)
            #print(avgVals)
            #print(stdVals)

            if simplePlot:
                self.generatePlotSimple(datesVals, minVals, maxVals, avgVals, stdVals, outputPlotFile, fieldBase, plotTitle)
            else:
                self.generatePlotStdDevRegion(datesVals, minVals, maxVals, avgVals, stdVals, outputPlotFile, fieldBase, plotTitle, 1.0)
            print("Completed")
        except ARCSIException as e:
            raise e
        except Exception as e:
            raise e


if __name__ == '__main__':
    """
    The command line user interface to ARCSI Plot Extracted Data.
    """
    parser = argparse.ArgumentParser(prog='arcsiplotextractedstats.py',
                                    description='''ARCSI command to plot
                                                   data.''',
                                    epilog='''A tool to plot data extracted
                                              using arcsi.''')
    # Request the version number.
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version ' + ARCSI_VERSION)

    parser.add_argument("-i", "--input", type=str,
                        help='''Specify the input statistics file.''')

    parser.add_argument("-f", "--field", type=str,
                        help='''String to specify the field to be plotted.
                                Note, this should be the base name, i.e.,
                                Red not RedAvg or RedMax''')

    parser.add_argument("-e", "--feature", type=int, default=0,
                        help='''Integer specifying the feature from the
                        shapefile defining the ROI to be plotted. If just
                        1 feature then value is zero (default).''')

    parser.add_argument("-t", "--title", type=str,
                        help='''String for the title to be added to the plot.''')

    parser.add_argument("-o", "--output", type=str,
                        help='''Output PDF file for the plot.''')

    parser.add_argument("--simple", action='store_true', default=False,
                        help='''Just a simple line plot with error bars''')

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

    if args.field == None:
        print("An input field needs to be specified.")
        parser.print_help()
        sys.exit()

    if args.title == None:
        print("A title for the plot is required.")
        parser.print_help()
        sys.exit()

    arcsiObj = ARCSIPlotExtractedStats()
    arcsiObj.processAndPlotData(args.input, args.output, args.field, args.title, args.simple, args.feature)











