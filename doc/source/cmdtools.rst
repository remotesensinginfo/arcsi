


Command Line Tools
===================

ARCSI is command line driven. While it is written in python and could be run directly from python using the arcsilib.arcsirun module this is not documented and arcsi.py also does some checking and setting up of default values so it would take a little effort and reference to the arcsi.py source code to setup your own code to call arcsilib.arcsirun. This is something we plan to improve at a later date.


Running ARCSI
--------------

These commands are for running ARCSI to generate analysis ready data (ARD).


arcsi.py
~~~~~~~~~

This is the main command that may well be the only command you use within ARCSI. This is the command which you use to process your data to retrieve ARD.

::

    arcsi.py -s ls5tm -p CLOUDS DOSAOTSGL STDSREF SATURATE TOPOSHADOW FOOTPRINT METADATA \
    -o ./Outputs/ --stats --format KEA --tmpath ./tmp --dem ./UKSRTM_90m.kea \
    --k  clouds.kea meta.json sat.kea toposhad.kea valid.kea stdsref.kea \
    -i LT05_L1TP_203024_19950815_20180217_01_T1/LT05_L1TP_203024_19950815_20180217_01_T1_MTL.txt


arcsimpi.py
~~~~~~~~~~~	

This command supports the identical functions to arcsi.py but can be used with MPI on computational clusters to process a group of input images using multiple processing cores (see arcsibuildmultifilelists.py).


Downloading Data
-----------------

These commands allow you to automate the processing of downloading your EO data, currently Landsat and Sentinel-2 using the Google Cloud.


arcsisetuplandsatdb.py
~~~~~~~~~~~~~~~~~~~~~~
This command creates a local copy of the Google database of landsat acquasitions as a SQLite database. 

::

    arcsisetuplandsatdb.py -f landsatdb_20180216.db


arcsisetupsen2db.py
~~~~~~~~~~~~~~~~~~~~

This command creates a local copy of the Google database of Sentinel-2 acquasitions as a SQLite database. 

::

    arcsisetupsen2db.py -f sen2db_20180216.db


arcsigenlandsatdownlst.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This command queries the database setup by arcsisetuplandsatdb.py generating a shell script for downloading the scenes meeting the query parameters.


::

    arcsigenlandsatdownlst.py -f landsatdb_20180216.db -p 204 -r 24 --collection T1  \
    --cloudcover 70 --startdate 2000-01-01 --outpath ./Downloads/ --multi \
    --lstcmds -o DwnldCmds.sh



arcsigensen2downlst.py
~~~~~~~~~~~~~~~~~~~~~~~

This command queries the database setup by arcsisetupsen2db.py generating a shell script for downloading the scenes meeting the query parameters.


::

    arcsigensen2downlst.py --source GOOG -f sen2db_20180216.db -t 30UVD \
    --cloudcover 70 --startdate 2017-01-01 --enddate 2017-12-31  \
    --outpath ./Downloads/ --multi --lstcmds -o DwnldCmds.sh




arcsidwnldgoog.py
~~~~~~~~~~~~~~~~~~

If you don't use the --lstcmds you can use this command to download the list of URLs and it will check whether the download already exists, skipping those that do. 



Batch Processing
----------------

A useful feature is being able to generate the arcsi.py commands for a large number of input files making it significantly easier to process a large number of input files. In addition to the commands listed here, arcsiextractdata.py may also be useful.


arcsibuildcmdslist.py
~~~~~~~~~~~~~~~~~~~~~~

This command allows you to build the arcsi.py commands for a large number of input files. This is really useful and can save a lot of time and error! You need to specify a search string and wildcard for finding the input image header files.

::

    arcsibuildcmdslist.py -s ls5tm -f KEA --stats -p CLOUDS DOSAOTSGL STDSREF \
    --outpath ./Outputs --dem ../UKSRTM_90m.kea \
    --keepfileends stdsref.kea clouds.kea \
    --tmpath ./tmp -i ./Inputs -e "*MTL.txt" -o LSARCSICmds.sh

arcsibuildmultifilelists.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~

A unique option within ARCSI is to process images files are a path/scene. For example, Sentinel-2 provides individual granules and landsat is cut into rows. ARCSI allows these images to be processed as a single job and tries to ensure that there isn't a boundary step between the individual input scenes/granules due to changes in parameterisation. To undertake this analysis you still use arcsi.py or arcsimpi.py (and you can use the mulitple core options) but requires a text file listing the header files for the scenes to be passed. This command (arcsibuildmultifilelists.py) provides the functionality to automatically build those text files from a directory of scenes you have downloaded (e.g., from Google).

arcsibuildfilenameslu.py
~~~~~~~~~~~~~~~~~~~~~~~~~~

ARCSI creates a standard unique name for each sensor using the meta-data for that image. However, if you have batch processed a large number of scenes then knowledge of which input dataset results in a particular output file might be lost. This command can build a look up table with this information.


Sentinel-2
-----------

These commands provide specific functionality for the Sentinel-2 sensor.

arcsichecksen2ver.py
~~~~~~~~~~~~~~~~~~~~~~

This command can be used to check that multiple versions of the same granule have been downloaded. ESA do reprocess some data on ocassion and Google keep all versions so you can end up with more than one copy of the same image on your system. This can cause a lot of problems using processing where file are over written. Therefore, this command was created to identify those images and select the most recent.


arcsisplitsen2granules.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first releases of Sentinel-2 data were provided by ESA as whole scenes (i.e., multiple granules). However, due to the size of the download ESA soon split the granules into indivdual downloads. ARCSI process granules and therefore these older files need to be split into individual granule files for processing. This is the functionality provided by this command. Please note that if you download your data from Google this will already have been split into granules as this is how Google have stored the data.


Landsat
--------

These commands provide specific functionality for the Landsat sensors.

arcsisortlandsat.py
~~~~~~~~~~~~~~~~~~~~

This command can sort a set of Landsat archives into a directory structure based on the sensor (i.e., Landsat 1, Landsat 2, ... Landsat 5 MSS, Landsat 5 TM, ... Landsat 8).


Error Checking
---------------

This functions can be useful for double checking that everything is processed through correctly particularly when you have a very large dataset that you cannot manually check.


arcsicheckfilespresent.py
~~~~~~~~~~~~~~~~~~~~~~~~~~

This command using the input data to check whether all the expected outputs have been produced.


arcsifindnotprocessed.py
~~~~~~~~~~~~~~~~~~~~~~~~~~
This comamnd will do a quick check as to whether an output file with the basename of input data has been created. 


arcsiremoveduplicates.py
~~~~~~~~~~~~~~~~~~~~~~~~~
This command aims to find duplicate input files and remove them before you undertake any processing. This is similar to arcsichecksen2ver.py.


Other Utilities
-----------------

arcsiextractdata.py
~~~~~~~~~~~~~~~~~~~~

Unless you have downloaded your data from Google you will probably have a set of archives (e.g., tar.gz, zip) for your images. It is useful to extract these into their own directories. This command provides functionality for this and can extract all archives within a directory or just a single input file.

::

    arcsiextractdata.py -i ./InputDIR -o ./OutputDIR




arcsiextractroistats.py
~~~~~~~~~~~~~~~~~~~~~~~~~~

This command can be used to extract, via zonal statistics, values for a region defined by a shapefile from a directory input images. Can you useful for QA checking your ARCSI outputs over expected invarient features.


Development Utilties
---------------------

These tools are not expected to be useful for the average user but are very useful for some creating a new sensor to be added to the ARCSI source code.


arcsispecresponsefuncs.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This command can resample a set of spectral response functions to a new sample interval. 6S requires spectral response functions to be sampled at 2.5 nm.


arcsisolarirradiance.py
~~~~~~~~~~~~~~~~~~~~~~~~

This command can use the spectral response functions to calculate the solar irradiance for the input band. This value is required for converting at sensor radiance to at sensor reflectance (also called top of atmosphere reflectance; TOA).


arcsicreatepy6scall.py
~~~~~~~~~~~~~~~~~~~~~~~

This command generates the py6S syntax for defining the sensor response function rather than using, if available, the in-built sensor response functions.





* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

