#!/bin/csh -f

# This script backs up enzo data hierarchies.  Execute it in the directory
# where enzo has dumped all of its datasets.  Procedurally, this script does
# the following:
#
# 1.  Determine if a given dataset exists by ensuring that there's a parameter
#     file and a non-zero-sized hierarchy file
#
# 2.  Creating a directory called $filename.dir where $filename is the same
#     as the restart parameter file
#
# 3.  Moving all enzo parameter, hierarchy, boundary and grid files into that
#     directory.
#
# 4.  Tarring up that directory into a file called $filename.tar -- when 
#     untarred, the grid files will untar into $filename.dir.  
#
# 5.  Optionally remove the directory $filename.dir -- NOTE:  This should
#     only be done if you REALLY trust the filesystem.
#
# Usage:   ./example_backup_script.bat <file_root> [<erase_directory>]
#
#    <file_root> is the first part of the parameter name.  For example, if
#      you want to back up all hierarchies that are named hier_0000, hier_0001,
#      etc. then <file_root> is hier_
#
#    <erase_directory> a boolean flag to tell the script whether or not you 
#      want to erase the directory that the grid files are placed into prior
#      to tarring.  The default is 0, which means that the script will not 
#      erase the directories.  A value of 1 means to erase them.
#
#   Only one argument MUST be used - <file_root>.  <erase_directory> is assumed
#   to be 0 if not explicitly stated.
#
# Example usage: ./example_backup_script.bat RedshiftOutput 1
#
#  This backups up all hierarchies named RedshiftOutput0000, 
#  RedshiftOutpu0001, etc. and erases the directories.
#
# NOTES:
#   This script is only good for up to 10,000 total grid files.
#   If you have more than 10,000 grid files you'll have to modify
#   this code by uncommenting a set of lines in the large foreach loop
#   (marked below), which will allow you to backup up to 100,000 grid files.
#
# Written by:  Brian O'Shea, February 2004
#
#

# are there only one or two arguments?  If not, error message
if( ($#argv != 2) && ($#argv != 1) ) then

  echo 'error:  script must have one or two arguments.'
  echo 'example_backup_script.bat <file_root> [<erase_directory>]'
  echo 'you have ' $#argv 'arguments.  Exiting...'
  echo ' '

  exit

endif

# if there are two arguments, first one is file root and
# second one is the erase_directory flag.  Read them in.
if( $#argv == 2 ) then
  set fileroot = ${1}
  set erase_directory = ${2}
  echo 'fileroot, erase_directory = '$fileroot', '$erase_directory

    # make sure erase_directory flag is a sane value
    if( ($erase_directory != 0) && ($erase_directory != 1) ) then
      echo 'erase_directory can only be 0 or 1!  Exiting...'
      echo ' '
      exit
    endif

endif

# if only one argument, set erase_directory to 0
if( $#argv == 1 ) then
  set fileroot = ${1}
  set erase_directory = 0   # off by default
  echo 'fileroot, erase_directory = '$fileroot', '$erase_directory
endif

# loop from 0 to 10000
set startfile = 0
set endfile = 10000
    

set i = $startfile
    
while ($i < $endfile)

  # get filename based on hierarchy number    
  if ( $i < 10 ) set filename = $fileroot'000'$i
  if ( $i >= 10 && $i < 100) set filename = $fileroot'00'$i
  if ( $i >= 100 && $i < 1000) set filename = $fileroot'0'$i
  if ( $i >= 1000 && $i < 10000) set filename = $fileroot$i
    
  set hierfile = $filename.hierarchy

  # if the parameter file ($filename) exists and there is a 
  # non-zero length hierarchy file ($hierfile) then this particular
  # enzo dump is done writing out, so it's safe to back it up    
  if ( -e $filename && !(-z $hierfile) ) then

    echo 'moving hierarchy ' $filename 'into '$filename'.dir'

    # create directory and move all of the hierarchy, boundary, 
    # parameter files into it
    mkdir $filename.dir

    mv $filename $filename.dir/.
    mv $filename.boundary $filename.dir/.
    mv $filename.boundary.hdf $filename.dir/.
    mv $filename.hierarchy $filename.dir/. 
	
    # loop over all possible grid filenames for this hierarchy.
    # If that grid exists, put it into the subdirectory.

    foreach m (0 1 2 3 4 5 6 7 8 9)
      foreach j (0 1 2 3 4 5 6 7 8 9)
        foreach k (0 1 2 3 4 5 6 7 8 9)
	  foreach l (0 1 2 3 4 5 6 7 8 9)
	  
	    # takes care of files 0-9999
	    if(-e $filename.grid$m$j$k$l) then
	      mv $filename.grid$m$j$k$l $filename.dir/.  >& /dev/null
	    endif

            ## If you want to back up files 10,000 - 99,999,
            ## uncomment the five lines after this line:

	    #foreach n (1 2 3 4 5 6 7 8 9)
	    #  if(-e $filename.grid$n$m$j$k$l) then
            #      mv $filename.grid$n$m$j$k$l $filename.dir/.  >& /dev/null
	    #  endif
	    #end
 
	  end
        end
      end
    end

    echo 'done moving all files into '$filename'.dir'
    echo 'tarring '$filename'.dir into '$filename'.tar'

    # tar directory
    tar cf $filename.tar $filename.dir

    # remove directory if user wants it.
    if($erase_directory == 1) then
      echo 'removing directory '$filename'.dir'
      rm -rf $filename.dir
    endif

    echo 'done backing up hierarchy '$filename 'into file '$filename'.tar'
    
  endif 
    
  @ i++
end
