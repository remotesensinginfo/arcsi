Atmospheric and Radiometric Correction of Satellite Imagery (ARCSI)
====================================================================

The Atmospheric and Radiometric Correction of Satellite Imagery (ARCSI) software provides a command line tool for the generation of Analysis Ready Data (ARD) optical data including atmospheric correction, cloud masking, topographic correction etc. of Earth Observation optical imagery (Blue-SWIR). The aim of ARCSI is to provide as automatic as possible method of generating Analysis Ready Data (ARD).


ARCSI is available to download from https://bitbucket.org/petebunting/arcsi/

Installation
-------------

To install ARCSI, we'd recommend you use our the builds available thought the conda-forge conda channel. Once you have installed miniconda (https://conda.io/miniconda.html; it must be a 64bit version) then you can install using the following command:

::

    conda install -c conda-forge arcsi 

However, we would recommend you create a new conda environment, which you do with the following command:

::

    conda create -n arcsienv python=3.5
    source activate arcsienv
    conda install -c conda-forge arcsi


Supported Sensors
-----------------

Currently we support:

* Sentinel-2
* Landsat 4, 5, 7, 8 (TM, ETM+, OLI)
* Landsat MSS 1, 2, 3, 4, 5
* Rapideye
* WorldView 2 and 3
* Pleiades
* SPOT 5, 6 and 7


Getting Started
---------------

To get started the basic command is:

::

    arcsi.py -s ls5tm -p CLOUDS DOSAOTSGL STDSREF SATURATE TOPOSHADOW FOOTPRINT METADATA \
    -o ./Outputs/ --stats --format KEA --tmpath ./tmp --dem ./UKSRTM_90m.kea \
    --k  clouds.kea meta.json sat.kea toposhad.kea valid.kea stdsref.kea \
    -i LT05_L1TP_203024_19950815_20180217_01_T1/LT05_L1TP_203024_19950815_20180217_01_T1_MTL.txt

This is processing a Landsat-5 scene from the file downloaded from the USGS to standardised surface reflectance (i.e., topographically corrected surface reflectance) which is masked for clouds, cloud shadows and topographic shadows. The aerosols (AOT) is also automatically derived.

The following command is the same but for Sentinel-2:

::
 
    arcsi.py -s sen2 --stats --format KEA \
    -p CLOUDS DOSAOTSGL STDSREF SATURATE TOPOSHADOW FOOTPRINT METADATA SHARP \
    -o ./Outputs  --dem ./UKSRTM_90m.kea --tmpath ./tmp  \
    --k  clouds.kea meta.json sat.kea toposhad.kea valid.kea stdsref.kea \
    -i S2A_MSIL1C_20170617T113321_N0205_R080_T30UVD_20170617T113319.SAFE/MTD_MSIL1C.xml



Information
---------------

.. toctree::
   :maxdepth: 2

   about
   tutorials
   cmdtools


Library Documentation
-----------------------

.. toctree::
   :maxdepth: 2
 
   arcsi_run
   arcsi_sensor


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

