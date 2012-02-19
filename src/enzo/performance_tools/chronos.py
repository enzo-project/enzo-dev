#!/usr/bin/env python

# chronos.py
# date: 10.13.11
# author: Cameron Hummels
# description: plots performance information from chronos.out

from optparse import OptionParser
usage = "usage: %prog <.out file>"
parser = OptionParser(usage)
(opts, args) = parser.parse_args()
if len(args) != 1:
        parser.error("incorrect number of arguments")

#if opts.write_fig:
import matplotlib as mpl
#matplotlib.use('Agg')
import pylab as pl
import numpy as np
from matplotlib import cm

filename = args[0]

def is_iterable(obj):
    """
    Checks to see if an object is iterable (i.e. if it is a list-like object)
    """
    return isinstance(obj, basestring) or getattr(obj, '__iter__', False)

def preserve_extrema(extrema, xdata, ydata):
    """
    Keep track of the universal extrema over multiple x-y datasets
    """
    minx = np.min([extrema[0], np.min(xdata)])
    maxx = np.max([extrema[1], np.max(xdata)])
    miny = np.min([extrema[2], np.min(ydata)])
    maxy = np.max([extrema[3], np.max(ydata)])
    return [minx,maxx,miny,maxy]

def plot_quantity(data, field_label, y_field_output, y_field_axis_label="",
                  x_field_output=0, x_field_axis_label="Cycle Number",
                  display=True, filename=""):
    """
    Produce a plot for the given label/output from the input file.

    Example:
    plot_quantity(data, "RebuildHierarchy", 1, "Mean Time (sec)")
                  
    """
    extrema = np.zeros(4)
    legend_list = []
    if is_iterable(field_label) and is_iterable(y_field_output):
        assert len(field_label) == len(y_field_output)
        for i in range(len(field_label)):
            xdata = data[field_label[i]][x_field_output]
            ydata = data[field_label[i]][y_field_output[i]]
            pl.plot(xdata,ydata)
            extrema = preserve_extrema(extrema,xdata,ydata)
            legend_list.append(field_label[i])
    else:
        xdata = data[field_label][x_field_output]
        ydata = data[field_label][y_field_output]
        pl.plot(xdata,ydata)
        extrema = preserve_extrema(extrema,xdata,ydata)
    pl.xlim(extrema[0:2])
    pl.ylim([extrema[2],1.1*extrema[3]])
    if len(legend_list) > 0:
        pl.legend(legend_list)
    #zerodata = np.zeros(len(ydata))
    #fillin = pl.fill_between(cycles,zerodata,ydata,facecolor='g')
    #fillin.set_alpha(0.5)
    pl.xlabel(x_field_axis_label)
    pl.ylabel(y_field_axis_label)
    if filename:
        pl.savefig(filename)
    if display:
        pl.show()
    pl.clf()
    
def plot_stack(data, field_label, y_field_output, y_field_axis_label="",
               x_field_output=0, x_field_axis_label="Cycle Number",
               display=True, filename="", repeated_field=""):
    """
    Produce a plot for the given label/outputs where each quantity is stacked
    on top of the previous quantity.

    Example:
    plot_stack(data, ["Level 0", "Level 1", "Level 2"], 1, "Mean Time (sec)")
    """
    extrema = np.zeros(4)
    legend_list = []

    ### If a repeated_field, figure out how many repeated fields there are.
    ### including any that were defined in the original field_label arg.
    if repeated_field:
        i = 0
        while data.has_key(repeated_field + ' %i' % i):
            field_label.append(repeated_field + ' %i' % i)
            i += 1
    num_fields = len(field_label)

    xdata = data[field_label[0]][x_field_output]
    cumulate_fields = np.zeros((2,len(xdata)))

    for i in range(0,num_fields):
        cumulate_fields[1] += data[field_label[i]][y_field_output]
        color = cm.jet(1.*i/num_fields)
        pl.fill_between(xdata,cumulate_fields[0],cumulate_fields[1],color=color)
        pl.plot(xdata,cumulate_fields[1],color=color)
        legend_list.append(field_label[i])
        extrema = preserve_extrema(extrema,xdata,cumulate_fields[1])
        # Move level down and repeat for new level
        cumulate_fields[0] = cumulate_fields[1]

    pl.xlim(extrema[0:2])
    pl.ylim([extrema[2],1.1*extrema[3]])
    if len(legend_list) > 0:
        pl.legend(legend_list)
    pl.xlabel(x_field_axis_label)
    pl.ylabel(y_field_axis_label)
    if filename:
        pl.savefig(filename)
    if display:
        pl.show()
    pl.clf()
    
### Create empty data structure which will store all input
data = {}

########
### XXX
### Need to play with more efficient manner of loading in data:
### Perhaps reading it all into memory first, then creating the dicts
### and the arrays before filling them.  Also, look into recarrays for
### giving specific handles on data entries.

### XXX 
### Need to add overplotting min/max over means.
### Better docstrings and better examples.

### Will all user-defined datasets give mean, sigma, min, max by default?
### Anything else?
########

### Open the file, and pull out the data.  
### Comments are ignored.  Each block of multiline text is treated individually.
### The first line of each text block is expected to be the cycle_number, 
### which is stored into an
### array.  Each subsequent line uses the first string as a dictionary key, 
### with the remaining columns the data for that key.  
### A blank line signifies the end of a text block

new_block = True
input = open(filename, "r")

for line in input:

# parse individual lines of data block until empty line, which signifies new
# block of text
    if line.strip() == '':
        new_block = True
        continue

    if not line.startswith("#") and line.strip():
        linelist = line.split()

        # if we're starting a new block of text
        if new_block:
            if linelist[0] != 'Cycle_Number':
                exit("Expected Cycle_Number at beginning of output block")
            else:
                cycle = int(linelist[1])
                new_block = False
                continue

        # if we've made it this far, cycle is defined and we're
        # populating our dictionary with data for cycle #x.
        line_key = linelist[0]

        ### Convert line_key's underscores to spaces
        line_key = " ".join(line_key.split('_'))

        line_value = np.array([cycle] + linelist[1:],dtype='float64') 
        line_value = np.nan_to_num(line_value)  # error checking

        ### If line_key is not already in our dictionary, then we create
        ### a new array as the dictionary payload

        if not data.__contains__(line_key):
            data[line_key] = line_value

        ### Otherwise, add a new row to the existing dataset
        else:
            data[line_key] = np.column_stack((data[line_key],
                                                   line_value))

input.close()

# Time to plot something with our data.  Data is now stored as a single
# dictionary, and each key represents a different "function" (ie the first
# string for a line in the chronos output file).  The payload for
# each dictionary key is an NxM array where N is the number of cycles
# and M is the number of fields outputted for that "function".

pl.figure(1)

# can remove in final, but makes more readable.  for all levels (ie Total).
cycles = data['Total'][0]
mean = data['Total'][1]
stddev = data['Total'][2]
min = data['Total'][3]
max = data['Total'][4]

#####
# testing 
#####
plot_quantity(data, ['Total', 'Level 0', 'Level 1', 'Level 2'], [1,1,1,1], "Total Time")
plot_quantity(data, 'Total', 1, "Total Time")
plot_stack(data, [], 1, "Mean Time (sec)", repeated_field="Level")
plot_stack(data, ['Total'], 1, "Mean Time (sec)", repeated_field="Level")
plot_stack(data, ['RebuildHierarchy','SolveHydroEquations'], 1, "Mean Time (sec)" )

######
# plot cycle vs mean walltime / proc
# with min/max plotted over in translucent green
######
pl.plot(cycles, mean)
minmax = pl.fill_between(cycles,min,max,facecolor='g')
pl.xlim([np.min(cycles),np.max(cycles)])
pl.ylim([0,1.1*np.max(max)])
minmax.set_alpha(0.5)
pl.xlabel("Cycle Number")
pl.ylabel("Average Walltime Spent (seconds)")
pl.suptitle("Aggregate Walltime Spent with min/max per Processor") 
#pl.savefig("time_total.png")
pl.show()
pl.clf()

######
# plot total cells/sec/processor processed over time
######
xdata = data['Total Cells/Sec/Processor'][0]
ydata = data['Total Cells/Sec/Processor'][1]
pl.plot(xdata, ydata)
pl.xlim([np.min(xdata),np.max(xdata)])
pl.ylim([np.min(ydata),np.max(ydata)*1.1])
pl.xlabel("Cycle Number")
pl.ylabel("Total Cells/Sec/Processor")
pl.suptitle("Cells Processed")
#pl.savefig("cell_total.png")
pl.show()
pl.clf()

######
# Now plot up time taken per level, stacked to get cumulative time
######

### Figure out how many Levels of refinement have been recorded
### And cycle through them, plotting a fill_between from the 0th to the 1st
### level, then the 1st to the 2nd, etc.
i = 0
legend_list = []
while data.has_key('Level %i' % i):
    i += 1
num_levels = i

# this temporarily stores the data from the ith to the i+1th level to be plotted
# plots from the cumulate_grids[0] to the cumulate_grids[1]
# starts at 0.
cumulate_grids = np.zeros((2,len(cycles)))

for i in range(0,num_levels):
    cumulate_grids[1] += data['Level %i' % i][6]
    color = cm.jet(1.*i/num_levels)
    pl.fill_between(cycles,cumulate_grids[0],cumulate_grids[1],color=color)
    pl.plot(cycles,cumulate_grids[1],color=color)
    legend_list.append('Level %i' % i)
    # Move level down and repeat for new level
    cumulate_grids[0] = cumulate_grids[1]

pl.xlabel("Cycle Number")
pl.ylabel("Grids")
pl.xlim([np.min(cycles),np.max(cycles)])
pl.ylim([0.,np.max(cumulate_grids[1])*1.1])
pl.legend(legend_list)
pl.suptitle("Grids vs Cycle") 

#pl.savefig("grids.png")
pl.show()
pl.clf()


######
# Now plot total timestep vs mean cells / sec / proc
######

### Figure out how many Levels of refinement have been recorded
### And cycle through them, plotting a fill_between from the 0th to the 1st
### level, then the 1st to the 2nd, etc.
i = 0
legend_list = []
while data.has_key('Level %i' % i):
    i += 1
num_levels = i

# this temporarily stores the data from the ith to the i+1th level to be plotted
# plots from the cumulate_grids[0] to the cumulate_grids[1]
# starts at 0.
cumulate_grids = np.zeros((2,len(cycles)))

for i in range(0,num_levels):
    cumulate_grids[1] += data['Level %i' % i][1]
    color = cm.jet(1.*i/num_levels)
    pl.fill_between(cycles,cumulate_grids[0],cumulate_grids[1],color=color)
    pl.plot(cycles,cumulate_grids[1],color=color)
    legend_list.append('Level %i' % i)
    # Move level down and repeat for new level
    cumulate_grids[0] = cumulate_grids[1]

pl.xlabel("Cycle Number")
pl.ylabel("Mean Time Spent")
pl.xlim([np.min(cycles),np.max(cycles)])
pl.ylim([0.,np.max(cumulate_grids[1])*1.1])
pl.legend(legend_list)
pl.suptitle("Mean Time Spent vs Cycle") 

#pl.savefig("grids.png")
pl.show()
pl.clf()
