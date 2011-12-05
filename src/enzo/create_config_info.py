"""
This script gets enzo version information and all compile options and 
generates a source file to be included in enzo that will print out this 
information when the code is run.  All information is written to the 
Enzo_Build file, which can be used as a loadable compiler settings 
file.

Author: Britton Smith <brittonsmith@gmail.com>
Date: December, 2011
"""

import os
import time

def get_hg_info():
    try:
        from mercurial import hg, ui, commands 
        u = ui.ui() 
        u.pushbuffer() 
        repo = hg.repository(u, os.path.expanduser('../..')) 
        commands.identify(u, repo) 
        my_info = u.popbuffer().strip().split()
        return my_info[:2]
    except ImportError:
        print "WARNING: could not get version information."
        return ('unknown', 'unknown')

def get_options(filename, my_options=None, get_list_order=False):
    if my_options is None: my_options = {}
    if get_list_order: option_list = []
    if not os.path.exists(filename): return my_options
    for line in open(filename).readlines():
        if line.startswith('#'): continue
        if '=' in line:
            key, val = line.split('=', 2)
            my_options[key.strip()] = val.strip()
            if get_list_order: option_list.append(key.strip())
    if get_list_order:
        return (my_options, option_list)
    else:
        return my_options

if __name__ == "__main__":
    my_options, option_list = get_options('Make.config.settings', get_list_order=True)
    my_options = get_options('Make.config.override', my_options=my_options)
    machine_info = get_options('Make.config.machine')
    changeset, branch = get_hg_info()
    output = open('auto_show_compile_options.C', 'w')
    output.write('#include <stdio.h>\n')
    output.write('void auto_show_compile_options(void) {\n')
    output.write('   FILE *opf;\n')
    output.write('   opf = fopen("Enzo_Build", "w");\n')
    output.write('   fprintf(opf, "### Enzo build information:\\n");\n')
    output.write('   fprintf(opf, "### Compiled: %s\\n");\n' % time.ctime())
    output.write('   fprintf(opf, "### Machine name: %s\\n");\n' % machine_info['CONFIG_MACHINE'])
    output.write('   fprintf(opf, "### Changeset: %s\\n");\n' % changeset)
    output.write('   fprintf(opf, "### Branch: %s\\n");\n' % branch)
    output.write('   fprintf(opf, "###\\n");\n')
    output.write('   fprintf(opf, "### Use this as a compile settings file by renaming\\n");\n')
    output.write('   fprintf(opf, "### it Make.settings.<keyword> and moving it to\\n");\n')
    output.write('   fprintf(opf, "### either src/enzo or .enzo in your home directory.\\n");\n')
    output.write('   fprintf(opf, "### Then, type \\\"make load-config-<keyword>\\\" to load\\n");\n')
    output.write('   fprintf(opf, "### compiler options.\\n");\n')
    output.write('   fprintf(opf, "\\n");\n')
    for key in option_list:
        output.write('   fprintf(opf, "%s = %s\\n");\n' % (key, my_options[key]))
    output.write('   fclose(opf);\n')
    output.write('}\n')
    output.close()
