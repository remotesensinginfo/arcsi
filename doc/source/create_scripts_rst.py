#!/usr/bin/env python
"""
Script to generate .rst file for scripts documentation. Runs scripts to get output and
inserts into text.

Inspired by answer suggesting modifying Makefile here:

http://stackoverflow.com/questions/7250659/python-code-to-generate-part-of-sphinx-documentation-is-it-possible

"""

import subprocess
import os
import glob

def get_script_help(script):
    """ Get output from command """

    out = subprocess.check_output([script, "-h"])

    stdout_lines = out.decode().split('\n')
    tabbed_lines = ''
    for line in stdout_lines:
        tabbed_lines += '    {}\n'.format(line)
    return tabbed_lines


if __name__ == "__main__":

    outfile = os.path.join(os.path.split(__file__)[0],'scripts.rst')
    scripts_list = glob.glob("../../bin/*py")


    scripts_text = '''

Scripts
========

'''
    for script in scripts_list:
        try:
            script_out = get_script_help(script)
        except Exception as err:
            print("Couldn't run {}\n{}".format(script, err))

        scripts_text = '''{}

{}
------------------------------

.. code-block:: text

{}


'''.format(scripts_text, os.path.basename(script), script_out)

    with open(outfile,'w') as f:
        f.write(scripts_text)
