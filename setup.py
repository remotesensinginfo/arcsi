#!/usr/bin/env python
"""
Setup script for ARCSI. Use like this for Unix:

$ python setup.py install

"""
#  This file is part of 'ARCSI' - Atmospheric and Radiometric Correction of Satellite Imagery
#  Copyright (C) 2013  Pete Bunting
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with ARCSI.  If not, see <http://www.gnu.org/licenses/>.
#
#
# Purpose:  Installation of the ARCSI software
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 27/08/2013
# Version: 1.0
#
# History:
# Version 1.0 - Created.

import os
from setuptools import setup
import glob

import arcsilib

setup(name='arcsi',
    version=arcsilib.ARCSI_VERSION,
    description='Atmospheric and Radiometric Correction of Satellite Imagery',
    author='Pete Bunting and Dan Clewley',
    author_email='pfb@aber.ac.uk, daniel.clewley@gmail.com',
    scripts=glob.glob("bin/*.py"),
    packages=['arcsilib', 'arcsilib/s2cloudless'],
    #package_dir={'arcsilib': 'arcsilib', 'arcsilib/s2cloudless': 'arcsilib/s2cloudless'},
    data_files=[(os.path.join('share','arcsi'),
                [os.path.join('data','WorldAerosolParams.kea'),
                 os.path.join('data','WorldAtmosphereParams.kea'),
                 os.path.join('data','pixel_s2_cloud_detector_lightGBM_v0.1.txt')])],
    license='LICENSE.txt',
    url='https://github.com/remotesensinginfo/arcsi',
    classifiers=['Intended Audience :: Developers',
                 'Intended Audience :: Remote Sensing Scientists',
                 'Intended Audience :: Atmospheric Scientists',
                 'Operating System :: OS Independent',
                 'Programming Language :: Python :: 3.7',
                 'Programming Language :: Python :: 3.8',
                 'Programming Language :: Python :: 3.9',
                 'Programming Language :: Python :: 3.10',
                 'Programming Language :: Python :: 3.11'])
