#!/usr/bin/env python
"""
Script to generate .rst file for scripts documentation. Runs scripts to get output and
inserts into text.

Inspired by answer suggesting modifying Makefile here:

http://stackoverflow.com/questions/7250659/python-code-to-generate-part-of-sphinx-documentation-is-it-possible

"""

import subprocess
import os

def get_command_out(command):
   """ Get output from command """

   out = subprocess.Popen(command,stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE,stderr=subprocess.PIPE)

   (stdout, stderr) = out.communicate()

   stdout_lines = stdout.decode().split('\n')
   tabbed_lines = ''
   for line in stdout_lines:
       tabbed_lines += '    {}\n'.format(line)

   return tabbed_lines

outfile = os.path.join(os.path.split(__file__)[0],'scripts.rst')

arcsi_out = get_command_out(['arcsi.py','-h'])
arcsibuildcmdslist_out = get_command_out(['arcsibuildcmdslist.py','-h'])
arcsiextractdata_out = get_command_out(['arcsiextractdata.py','-h'])
arcsiextractroistats_out = get_command_out(['arcsiextractroistats.py','-h'])
arcsiplotextractedstats_out = get_command_out(['arcsiplotextractedstats.py','-h'])
arcsisolarirradiance_out = get_command_out(['arcsisolarirradiance.py','-h'])
arcsisortlandsat_out = get_command_out(['arcsisortlandsat.py','-h'])
arcsispecresponsefuncs_out = get_command_out(['arcsispecresponsefuncs.py','-h'])

scripts_text = '''

Scripts
========

arcsi.py
-------------------

Main ARCSI script

.. code-block:: text 

{}


arcsibuildcmdslist.py
-------------------

.. code-block:: text 

{}

arcsiextractdata.py
-------------------

.. code-block:: text 

{}

arcsiextractroistats.py
----------------------

.. code-block:: text 

{}

arcsiplotextractedstats.py
------------------------

.. code-block:: text 

{}


arcsisolarirradiance.py
----------------------

.. code-block:: text 

{}

arcsisortlandsat.py
----------------------

.. code-block:: text 

{}


arcsispecresponsefuncs.py
------------------------

.. code-block:: text 

{}



'''.format(arcsi_out, 
        arcsibuildcmdslist_out, 
        arcsiextractdata_out, 
        arcsiextractroistats_out,
        arcsiplotextractedstats_out,
        arcsisolarirradiance_out, 
        arcsisortlandsat_out,
        arcsispecresponsefuncs_out)

with open(outfile,'w') as f:
    f.write(scripts_text)
