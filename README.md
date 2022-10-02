# README #

The Atmospheric and Radiometric Correction of Satellite Imagery (ARCSI) software provides a command line tool for the atmospheric correction of Earth Observation imagery. The aim of ARCSI is to provide as automatic as possible method of retrieving the atmospheric correction parameters and using them to parameterise 6S.

To install ARCSI tot eh default location following command is to be used from within the ARCSI source directory:

``
pip install .
``


if you want to install into another location using the --prefix option.

``
pip install . --prefix=/to/install/path
``


To run ARCSI you need up to date versions of:

* RSGISLib
* KEALib
* GDAL
* RIOS 
* Py6S
* 6S

## Need support? ##

If you need support using ARCSI or think you've found a bug please email us on rsgislib-support@googlegroups.com.

## Contribution guidelines ##

If you would like to contribute to the library with bug fixes or new functionality (e.g. new sensors) please fork the repository and send a pull request. Or email me Dr Pete Bunting (pfb@aber.ac.uk). 

