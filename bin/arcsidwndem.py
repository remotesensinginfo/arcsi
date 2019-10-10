#! /usr/bin/env python
############################################################################
#  arcsidwndem.py
#
#  Copyright 2019 ARCSI.
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
# Date: 10/10/2019
# Version: 1.0
#
# History:
# Version 1.0 - Created.
#
############################################################################

import rsgislib
import os.path
import elevation
import argparse
# Import the arcsi version number
from arcsilib import ARCSI_VERSION
# Import the list of sensors arcsi supports
from arcsilib import ARCSI_SENSORS_LIST
# Import the ARCSI exception class
from arcsilib.arcsiexception import ARCSIException
# Import the ARCSI sensor factory class
from arcsilib.arcsiutils import ARCSISensorFactory

if __name__ == '__main__':
    """
    The command line user interface to ARCSI DEM download tool.
    """
    parser = argparse.ArgumentParser(prog='arcsidwndem.py',
                                    description='''Download DEM for area of interest defined by input image.''')
    # Request the version number.
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version ' + ARCSI_VERSION)
    # Define the argument for specifying the sensor.
    parser.add_argument("-s", "--sensor", required=True, choices=ARCSI_SENSORS_LIST, help='''Specify the sensor of the input image.''')
    # Define the argument for specifying the input directory to be processed.
    parser.add_argument("-i", "--inputheader", type=str, required=True, help='''Specify the input image header file.''')
    # Define the argument for specifying the output file.
    parser.add_argument("-o", "--output", type=str, required=True, help='''The output DEM file - outputs as a GeoTIFF.''')
    parser.add_argument("-b", "--buffer", type=float, default=0.5, help='''Specify the buffer around the image for which the DEM will be downloaded.''')
    parser.add_argument("-l", "--limit", type=int, default=10, help="A limit on the number of tiles which can be downloaded.")
    # Call the parser to parse the arguments.
    args = parser.parse_args()
    
    output_img = os.path.abspath(args.output)
    
    sensorFact = ARCSISensorFactory()
    sensor_cls_obj = sensorFact.getSensorClassFromName(args.sensor, False, None)
    sensor_cls_obj.extractHeaderParameters(args.inputheader, None)
    image_bbox_latlon = sensor_cls_obj.getBBOXLatLon()
    
    bounds_ext = (image_bbox_latlon[0]-args.buffer, image_bbox_latlon[2]-args.buffer, image_bbox_latlon[1]+args.buffer, image_bbox_latlon[3]+args.buffer)
    print(bounds_ext)
    try:
        elevation.clip(bounds=bounds_ext, output=output_img, max_download_tiles=args.limit)
    except Exception as e:
        print("An error has occurred when downloading and processing the DEM data. Try re-running as data is cached.")
        raise e
    elevation.clean()







