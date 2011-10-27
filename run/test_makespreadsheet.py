#!/usr/bin/env python
#
#  Looks through the run directory for "*.enzotest" and makes a csv spreadsheet of all answer test scripts and their properties.
#
#  Author: David Collins (dcollins4096@gmail.com), 2011-06-14 11:19 AM.  It's a bright sunny day here in Los Alamos.
#

import fnmatch
import os


#Hunt for enzotest files.
dbg = 0
matches = []
for root, dirnames, filenames in os.walk('.'):
  for filename in fnmatch.filter(filenames, '*.enzotest'):
      matches.append(os.path.join(root, filename))

#Generate dictionary, make list of attributes
tests={}
attribute_list=['name','nprocs','max_time_minutes','dimensionality','runtime','answer_testing_script','hydro','gravity','cooling','chemistry','cosmology','author','mhd','radiation','AMR']
for file in matches:
    if dbg > 0:
        print file
    lines = open(file,'r').readlines()
    tests[file]={}
    for line in lines:
        if line.strip():
            key, value = line.split("=")
            if value[-1]=='\n':
                value = value[:-1]
            if value.count("#") > 0:
                value = value[0:value.index("#")]
            key=key.strip()
            value=value.strip()
            if key not in attribute_list:  
                attribute_list.append(key)
            tests[file][key]=value

#make csv
csv = open('test_spreadsheet.csv','w')
dir( csv)

#head row
csv.write("filename"+",%s"*len(attribute_list)%tuple(attribute_list)+"\n")

#lines
for file in matches:
    csv.write(file)
    for attr in attribute_list:
        csv.write(", %s"%(tests[file].get(attr,'')))
    csv.write("\n")
csv.close()


#end
