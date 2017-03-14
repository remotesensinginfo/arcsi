#! /usr/bin/env python

"""
Module that contains the ARSCI command to extract data from archives.
"""

############################################################################
#  arcsiextractdata.py
#
#  Copyright 2015 ARCSI.
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
# Purpose:  A script to unarchive data from tar / tar.gz files into
#           a directory structure with a directory for contents of
#           each archive.
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 14/02/2015
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
# Import the arcsi version number
from arcsilib import ARCSI_VERSION


class ARCSISensitivityAnalysis (object):




    def run(self):
        print("Hello World...")





if __name__ == '__main__':
    """
    The command line user interface to ARCSI Sensitivity Analysis Tool.
    """
    parser = argparse.ArgumentParser(prog='arcsisensitivity.py',
                                    description='''ARCSI command to run a
                                                    sensitivity analysis on 6S.''',
                                    epilog='''A tools to run a sensitivity
                                              analysis on 6S for a range of
                                              parameter inputs.''')
    # Request the version number.
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version ' + ARCSI_VERSION)

    # Call the parser to parse the arguments.
    args = parser.parse_args()


    arcsiObj = ARCSISensitivityAnalysis()

    arcsiObj.run()
