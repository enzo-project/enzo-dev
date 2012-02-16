#!/usr/bin/env python

# chronos.py
# date: 10.13.11
# author: Cameron Hummels
# description: polls enzo performance information from XXX. 

# Use the cool OptionParser (only works with Python 2.6, broken with python 2.7
from optparse import OptionParser

usage = "usage: %prog [options] <XXX file>"
parser = OptionParser(usage)
(opts, args) = parser.parse_args()
if len(args) != 1:
        parser.error("incorrect number of arguments")

# If we're just writing the figure (and not displaying it), then don't use
# x-windows backend to render the image (or else it will break if it doesn't
# have $display set).
#if opts.write_fig:
import matplotlib as mpl
#matplotlib.use('Agg')
import pylab as pl
import numpy as np
from matplotlib.patches import Polygon
from matplotlib import cm

filename = args[0]

### Create empty lists
cycle_list = []
data_dict = {}

### Open the file, and pull out the data.  
### Comments are ignored.  Each block of multiline text is treated individually.
### The first line is expected to be the cycle_number, which is stored into an
### array.  Each subsequent line uses the first string as a dictionary key, 
### with the remaining columns the data for that key.  

input = open(filename, "r")
for line in input:
    if not line.startswith("#") and line.strip():
        linelist = line.split()
        if linelist[0] != 'Cycle_Number':
            exit("Expected Cycle_Number at beginning of output block")
        cycle = linelist[1]
        cycle_list.append(cycle)
        cycle = [int(cycle)]
        # Start to parse individual lines of data block until empty line
        for line in input:
            if line.strip() == '':
                break
            if not line.startswith("#"):
                linelist = line.split()
                line_key = linelist[0]
                ### Convert line_key's underscores to spaces
                line_key = " ".join(line_key.split('_'))
                line_value = np.array(cycle + linelist[1:],dtype='float64')
                ### If existing key, then append this, otherwise make new entry in dict
                if data_dict.__contains__(line_key):
                    data_dict[line_key] = np.column_stack((data_dict[line_key],line_value))
                else:
                    data_dict[line_key] = line_value

input.close()

times = np.array(cycle_list,dtype='int32')

### data_dict is a dictionary with keys representing different timing processes 
### (e.g. level information, rebuild_hierarchy, timed chunks of codes, etc.) 

### the payload of the dictionary for each key is an (N+1)xM array, where N = the
### number of columns of data associated with each timing process (the first is 
### the cycle number), and M represents the number of cycles.

pl.figure(1)

# can remove in final, but makes more readable.  for all levels.
cycle = data_dict['Total'][0,:]
mean = data_dict['Total'][1,:]
stddev = data_dict['Total'][2,:]
min = data_dict['Total'][3,:]
max = data_dict['Total'][4,:]

# plot cycle vs mean walltime / proc
# with min/max plotted over
pl.plot(cycle, mean)
minmax = pl.fill_between(cycle,min,max,facecolor='g')
pl.xlim([np.min(cycle),np.max(cycle)])
pl.ylim([0,1.1*np.max(max)])
minmax.set_alpha(0.5)
#sigmax = pl.fill_between(cycle,sigmin,sigmax,facecolor='b')
#sigmax.set_alpha(0.5)
pl.xlabel("Cycle Number")
pl.ylabel("Average Walltime Spent (seconds)")
pl.suptitle("Aggragate Walltime Spent with min/max per Processor") 
#pl.savefig("time_total.png")
pl.show()
pl.clf()

# plot total cells/sec/processor processed over time
pl.plot(cycle, data_dict['Total Cells/Sec/Processor:'][1,:])
pl.xlabel("Cycle Number")
pl.ylabel("Total Cells/Sec/Processor")
pl.suptitle("Cells Processed")
#pl.savefig("time_total.png")
pl.show()
pl.clf()

### Now plot up time taken per level, stacked to get cumulative time
cumulate_grids = np.zeros((2,len(cycle_list)))

### Figure out how many Levels of refinement have been recorded
i = 0
legend_list = []
while data_dict.has_key('Level %i' % i):
    cumulate_grids[1,:] += data_dict['Level %i' % i]
    i += 1
#XXX
    # plot timestep vs num_cells on diff levels
    pl.fill_between(times,cumulate_grids[0,:],cumulate_grids[1,:],color=cm.jet(0))
    pl.plot(times,cumulate_grids[1,:],color=cm.jet(0))
    legend_list.append('Level %i' % 0)
    cumulate_grids[0,:] = cumulate_grids[1,:]

#    pl.plot(times, data[:,i,0])
    #pl.fill_between(ax,times,data[:,i-1,0],data[:,i,0])
    cumulate_grids[1,:] += data[:,i,0]
    pl.fill_between(times,cumulate_grids[0,:],cumulate_grids[1,:],color=cm.jet(1.*(i-1)/len(levels)))
    pl.plot(times,cumulate_grids[1,:],color=cm.jet(1.*(i-1)/len(levels)))
    cumulate_grids[0,:] = cumulate_grids[1,:]
    legend_list.append('Level %i' % (i-1))
pl.xlabel("Cycle Number")
pl.ylabel("Number of Grids")
leg = pl.legend(legend_list)
pl.suptitle("Number of Grids vs Cycle") 

#pl.show()
pl.savefig("grids.png")
pl.clf()



# plot total timestep vs mean cells / sec / proc
# with min/max/stddev plotted over
cumulate_grids = np.zeros((2,len(cycle_list)))

cumulate_grids[1,:] += data[:,1,1]
legend_list = []
pl.fill_between(times,cumulate_grids[0,:],cumulate_grids[1,:],color=cm.jet(0))
pl.plot(times,cumulate_grids[1,:],color=cm.jet(0))
legend_list.append('Level %i' % 0)
cumulate_grids[0,:] = cumulate_grids[1,:]

for i in levels:
#    pl.plot(times, data[:,i,1])
    cumulate_grids[1,:] += data[:,i,1]
    pl.fill_between(times,cumulate_grids[0,:],cumulate_grids[1,:],color=cm.jet(1.*(i-1)/len(levels)))
    pl.plot(times,cumulate_grids[1,:],color=cm.jet(1.*(i-1)/len(levels)))
    cumulate_grids[0,:] = cumulate_grids[1,:]
    legend_list.append('Level %i' % (i-1))

#    minmax = pl.fill_between(times,min,max,facecolor='g')
    pl.xlim([np.min(times),np.max(times)])
#    pl.ylim([0,np.max(max)])
#    minmax.set_alpha(0.5)
#    sigmax = pl.fill_between(times,sigmin,sigmax,facecolor='b')
#    sigmax.set_alpha(0.5)

pl.xlabel("Cycle Number")
pl.ylabel("Walltime Spent (seconds)")
leg = pl.legend(legend_list)
pl.suptitle("Walltime Spent on Levels") 
#pl.show()
pl.savefig("time_levels.png")
pl.clf()

pl.close()
