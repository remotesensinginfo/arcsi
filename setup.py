#!/usr/bin/env python
"""
Setup script for ARCSI. Use like this for Unix:

$ python setup.py install

"""
#  This file is part of 'ARCSI' - Atmospheric and Radiometic Correction of Satellite Imagery
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


from distutils.core import setup

setup(name='ACRSI', 
    version='1.0', 
    description='Atmospheric and Radiometic Correction of Satellite Imagery',
    author='Pete Bunting',
    author_email='pfb@aber.ac.uk',
    scripts=['bin/arcsi.py', 'bin/arcsisolarirradiance.py', 'bin/arcsispecresponsefuncs.py', 'bin/arcsiextractdata.py', 'bin/arcsibuildcmdslist.py', 'bin/arcsisortlandsat.py'],
    packages=['arcsilib'],
    license='LICENSE.txt',
    url='https://bitbucket.org/petebunting/arcsi',
    classifiers=['Intended Audience :: Developers',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.2',
          'Programming Language :: Python :: 3.3'])
