"""
This script gets enzo version information and all compile options and 
generates a source file (auto_show_compile_options.C) to be included 
in enzo that will print out this information when the code is run.  
All information is written to the Enzo_Build file, which can be used 
as a loadable compiler settings file.  If uncommitted changes exist 
when Enzo is built, the output of 'hg diff' is written to 
Enzo_Build_Diff.

Author: Britton Smith <brittonsmith@gmail.com>
Date: December, 2011
"""

import os
import time

def get_hg_info():
    try:
        #from mercurial import hg, ui, commands 
        #from mercurial.error import RepoError
        import hglib
    except ImportError:
        print("WARNING: could not get version information.  Please install mercurial and/or hglib.")
        return ('unknown', 'unknown', None)
    
    try:
        client = hglib.open('../..')
        changeset = str(client.tip().node)
        branch = str(client.tip().branch)
        my_diff = client.diff().decode('utf-8')
        
        return (changeset, branch, my_diff)
    
    except RepoError:
        print("WARNING: could not get version information.")
        return ('unknown', 'unknown', None)

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
    # get default compile options
    my_options, option_list = get_options('Make.config.settings', get_list_order=True)
    # override defaults with custom settings
    my_options = get_options('Make.config.override', my_options=my_options)
    # get machine name
    machine_info = get_options('Make.config.machine')
    # get mercurial version information and any uncommitted changes
    changeset, branch, diff = get_hg_info()
    output_file = 'Enzo_Build'
    diff_file = '%s_Diff' % output_file
    output = open('auto_show_compile_options.C', 'w')
    output.write('#include <stdio.h>\n')
    output.write('void auto_show_compile_options(void) {\n')
    output.write('   FILE *opf;\n')
    output.write('   opf = fopen("%s", "w");\n' % output_file)
    output.write('   fprintf(opf, "### Enzo build information:\\n");\n')
    output.write('   fprintf(opf, "### Compiled: %s\\n");\n' % time.ctime())
    output.write('   fprintf(opf, "### Machine name: %s\\n");\n' % machine_info['CONFIG_MACHINE'])
    output.write('   fprintf(opf, "### Branch: %s\\n");\n' % branch)
    output.write('   fprintf(opf, "### Changeset: %s\\n");\n' % changeset)
    if diff:
        output.write('   fprintf(opf, "###\\n");\n')
        output.write('   fprintf(opf, "### Enzo was built with uncommitted changes.\\n");\n')
        output.write('   fprintf(opf, "### Check %s to see the output of hg diff.\\n");\n' % diff_file)
    output.write('   fprintf(opf, "###\\n");\n')
    output.write('   fprintf(opf, "### Use this as a compile settings file by renaming\\n");\n')
    output.write('   fprintf(opf, "### it Make.settings.<keyword> and moving it to\\n");\n')
    output.write('   fprintf(opf, "### either src/enzo or .enzo in your home directory.\\n");\n')
    output.write('   fprintf(opf, "### Then, type \\\"make load-config-<keyword>\\\" to load\\n");\n')
    output.write('   fprintf(opf, "### compiler options.\\n");\n')
    output.write('   fprintf(opf, "###\\n");\n')
    output.write('   fprintf(opf, "### The current compile options are shown below.\\n");\n')
    output.write('   fprintf(opf, "### For more information on the meaning of the\\n");\n')
    output.write('   fprintf(opf, "### options, consult src/enzo/Make.config.settings.\\n");\n')
    output.write('   fprintf(opf, "\\n");\n')
    for key in option_list:
        output.write('   fprintf(opf, "%s = %s\\n");\n' % (key, my_options[key]))
    output.write('   fclose(opf);\n')
    if diff:
        output.write('\n')
        output.write('   opf = fopen("%s", "w");\n' % diff_file)
        for line in diff.split('\n'):
            line = line.replace('\\', '\\\\')
            line = line.replace('"', '\\"')
            line = line.replace('\r', '\\n')
            line = line.replace('%', '%%')
            output.write('   fprintf(opf, "')
            output.write(r"%s" % line)
            output.write('\\n");\n')
        output.write('   fclose(opf);\n')
    output.write('}\n')
    output.close()
