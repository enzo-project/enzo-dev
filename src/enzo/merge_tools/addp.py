#!/usr/bin/env python

"""
addp.py
add a global parameter to the 4 enzo files that need it:
global_data.h, SetDefaultGlobalValues.C, ReadParameterFile.C, WriteParameterFile.C

addp.py type name initial_value
"""

import shutil
import sys
def insert(filename,output_line):
    """looks for the line
    //MHDCT variables
    inserts *line* after"""
    backup_name = filename + '.old'
    shutil.copy(filename,backup_name)
    file_in = open(backup_name,'r')
    file_out = open(filename,'w')
    FoundIt = False
    print_string = filename + " "
    for line in file_in:
        file_out.write(line)
        try:
            line.index('//MHDCT variables')
            FoundIt = True
            file_out.write(output_line+"\n");
            print_string += output_line
        except:
            pass
    file_in.close()
    file_out.close()
    if not FoundIt:
        print_string += "//MHDCT variables not found"
    print print_string
if len(sys.argv) != 4 and len(sys.argv) != 5:
    print """addp.py
    add a global parameter to the 4 enzo files that need it:
    global_data.h, SetDefaultGlobalValues.C, ReadParameterFile.C, WriteParameterFile.C

    addp.py type name initial_value comment
    The comment (if present) goes in global_data.h"""

else:    
    type = sys.argv[1]
    name = sys.argv[2]
    value_in = sys.argv[3]
    comment = None
    if len(sys.argv) == 5:
        comment = sys.argv[4]

    Ctypes = {'int':'int','float':'float'}
    OutTypes = {'int':'ISYM','float':'GSYM'}
    InTypes = {'int':'ISYM','float':'GSYM'}
    if type not in ['float','int']:
        print "type (%s) not recognized."%type
        sys.exit(1)


    global_string = "EXTERN %s %s;"%(Ctypes[type],name)
    if comment:
        global_string += "// %s"%comment

    Set_Default_string = "  %s = %s;"%(name,value_in)
    input_control = '%"'+InTypes[type]
    ReadParameter_string = '    ret += sscanf(line,"%s = %s,&%s);'%(name,input_control,name)
    output_control = '%"'+OutTypes[type]+'"\\n"'
    output_string = name + " "*(30 - len(name))
    WriteParameter_string = '  fprintf(fptr,"%s=%s,%s);'%(output_string, output_control,name)
    insert('global_data.h',global_string)
    insert('SetDefaultGlobalValues.C',Set_Default_string)
    insert('ReadParameterFile.C',ReadParameter_string)
    insert('WriteParameterFile.C',WriteParameter_string)
#end
