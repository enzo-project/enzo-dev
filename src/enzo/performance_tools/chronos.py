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

def is_listlike(obj):
    """
    Checks to see if an object is listlike (but not a string)
    """
    import types
    return not (isinstance(obj, basestring) or type(obj) is types.IntType)

def preserve_extrema(extrema, xdata, ydata):
    """
    Keep track of the universal extrema over multiple x-y datasets

    Parameters
    ----------
    extrema : 5-element list
        This keeps track of the current min/max for x and y datasets.  It
        has a 5th element to keep track whether the extrema have yet been
        set. [xmin, xmax, ymin, ymax, already_set_bit]
        For already_set_bit, 0=no, 1=yes.
    xdata,ydata : array-like
        The dataset you wish to set add to your minima/maxima
    """
    ### If setting extrema for the first time, just use xdata/ydata min/max
    if extrema[4] == 0:
        minx = np.min(xdata)
        maxx = np.max(xdata)
        miny = np.min(ydata[np.nonzero(ydata)])
        maxy = np.max(ydata[np.nonzero(ydata)])
    ### Otherwise, preserve the existing extrema
    else:
        minx = np.min([extrema[0], np.min(xdata)])
        maxx = np.max([extrema[1], np.max(xdata)])
        miny = np.min([extrema[2], np.min(ydata[np.nonzero(ydata)])])
        maxy = np.max([extrema[3], np.max(ydata[np.nonzero(ydata)])])
    return [minx,maxx,miny,maxy,1]

def plot_quantity(data, field_label, y_field_output, y_field_axis_label="",
                  x_field_output=0, x_field_axis_label="Cycle Number",
                  display=True, filename="", repeated_field="", 
                  log_y_axis="Auto"):
    """
    Produce a plot for the given quantity(s) from the input file.

    Parameters
    ----------
    data : dictionary of recarrays
        The dictionary containing all of the data for use in plotting.
    field_label : string or array_like of strings
        The label of the field you wish to plot.  If you wish to plot
        multiple fields, enumerate them in an array or tuple. Ex: "Level 0"
    y_field_output : int or array_like of ints
        The index of the field you wish to plot. Ex: 1 for mean time.
    y_field_axis_label : string, optional
        The y axis label on the resulting plot. Default = ""
    x_field_output : int, optional
        The index of the x data you wish to plot. Default = 0 (Cycles)
    x_field_axis_label : string, optional
        The x axis label on the resulting plot. Default = "Cycle Number"
    display : bool, optional
        Do you wish to display the data to the screen?
    filename : string, optional
        The filename where I will store your plotted data.
    repeated_field : string, optional
        If you have a regularly named set of fields you wish to plot 
        against each other (e.g. "Level 0", "Level 1", "Level 2"), then
        include the string here and they will all be included automatically
        and in order (e.g. "Level").
    log_y_axis : string, optional
        This controls whether the plot will use logarithmic units for the
        y axis.  Valid settings are "Auto", "On", and "Off".  When "Auto" is
        used, the code automatically recognizes when you have a maximum 
        y value more than 3 orders of magnitude greater than your minimum y
        value (for non-zero values) at which point it plots the y axis in 
        log units.

    See Also
    --------
    plot_stack

    Examples
    --------
    To produce a simple plot of the mean time taken over the course of the
    simulation and display it to the screen:

    >>> plot_quantity(data, "RebuildHierarchy", 1, "Mean Time (sec)")

    To produce a plot comparing the RebuildHiearchy and SolveHydroEquations
    maximum time taken over the course of the simulation and save it 
    to a file: "test.png"

    >>> plot_quantity(data, ["RebuildHierarchy", "SolveHydroEquations"],
    4, "Maximum Time (sec)", display=False, filename="test.png")

    To produce a plot comparing the maximum time from RebuildHiearchy and 
    the minimum time from SolveHydroEquations taken over the course of the 
    simulation and save it to a file: "test.png"

    >>> plot_quantity(data, ["RebuildHierarchy", "SolveHydroEquations"],
    [4,3], "Time (sec)", display=False, filename="test.png")

    To produce a plot comparing the mean time taken by all of the different
    levels over the course of the simulation and save it to a file: "test.png"

    >>> plot_quantity(data, [], 1, "Mean Time (sec)", display=False, 
                      filename="test.png", repeated_field="Level")
    """
    extrema = np.zeros(5)
    legend_list = []

    ### Convert plots of single quantities to list format for homogenous 
    ### processing.
    if not is_listlike(field_label):
        field_label = [field_label]

    ### If there is a repeated_field, figure out how many fields 
    ### there are including any that were defined in the original 
    ### field_label argument.
    if repeated_field:
        i = 0
        while data.has_key(repeated_field + ' %i' % i):
            field_label.append(repeated_field + ' %i' % i)
            i += 1
    num_fields = len(field_label)

    ### If y_field_output is a single index, then replace it with a list of
    ### identical indices
    if not is_listlike(y_field_output):
        y_field_output = num_fields*[y_field_output]

    ### Loop through the y datasets to figure out the extrema
    for i in range(len(field_label)):
        xdata = data[field_label[i]][x_field_output]
        ydata = data[field_label[i]][y_field_output[i]]
        extrema = preserve_extrema(extrema,xdata,ydata)
    if log_y_axis=="Auto":
        if extrema[3]/extrema[2] > 1e3:
            log_y_axis="On"
        else:
            log_y_axis="Off"

    ### Now for the actual plotting
    for i in range(len(field_label)):
        xdata = data[field_label[i]][x_field_output]
        ydata = data[field_label[i]][y_field_output[i]]
        if log_y_axis=="On":
            pl.semilogy(xdata,ydata)
        else:
            pl.plot(xdata,ydata)
        legend_list.append(field_label[i])

    pl.xlim(extrema[0:2])
    if log_y_axis=="On":
        y_log_range = 1.2*np.log10(extrema[3]/extrema[2])
        pl.ylim([extrema[2],extrema[2]*10**y_log_range])
    else:
        pl.ylim([0.,1.2*extrema[3]])
    pl.legend(legend_list,2)
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
               display=True, filename="", repeated_field="",
               log_y_axis="Auto"):
    """
    Produce a plot for the given label/outputs where each quantity is stacked
    on top of the previous quantity.

    Parameters
    ----------
    data : dictionary of recarrays
        The dictionary containing all of the data for use in plotting.
    field_label : array_like of strings
        The labels of the fields you wish to plot.  If you wish to plot
        a single field, then use plot_quantity instead.
    y_field_output : int or array_like of ints
        The index of the field you wish to plot. Ex: 1 for mean time.
        If you have a single value for many field_labels, it is assumed
        that such a value will index all of them.  If you have an array_like
        structure of ints, each one will index the corresponding key in 
        field_lable.  
    y_field_axis_label : string, optional
        The y axis label on the resulting plot. Default = ""
    x_field_output : int, optional
        The index of the x data you wish to plot. Default = 0 (Cycles)
    x_field_axis_label : string, optional
        The x axis label on the resulting plot. Default = "Cycle Number"
    display : bool, optional
        Do you wish to display the data to the screen?
    filename : string, optional
        The filename where I will store your plotted data.
    repeated_field : string, optional
        If you have a regularly named set of fields you wish to plot 
        against each other (e.g. "Level 0", "Level 1", "Level 2"), then
        include the string here and they will all be included automatically
        and in order (e.g. "Level").
    log_y_axis : string, optional
        This controls whether the plot will use logarithmic units for the
        y axis.  Valid settings are "Auto", "On", and "Off".  When "Auto" is
        used, the code automatically recognizes when you have a maximum 
        y value more than 3 orders of magnitude greater than your minimum y
        value (for non-zero values) at which point it plots the y axis in 
        log units.

    Example:
    plot_stack(data, ["Level 0", "Level 1", "Level 2"], 1, "Mean Time (sec)")
    """
    extrema = np.zeros(5)
    legend_list = []
    
    ### If a repeated_field, figure out how many repeated fields there are.
    ### including any that were defined in the original field_label arg.
    if repeated_field:
        i = 0
        while data.has_key(repeated_field + ' %i' % i):
            field_label.append(repeated_field + ' %i' % i)
            i += 1
    num_fields = len(field_label)

    ### If y_field_output is a single index, then replace it with a list of
    ### identical indices
    if not is_listlike(y_field_output):
        y_field_output = num_fields*[y_field_output]

    ### Since we're stacking (and making plots from a bottom bound to an
    ### upper bound for each y value, we need to cumulate it as we stack.
    ### ydata_cum is the same size as xdata, but it has two indices, one
    ### for the bottom bound of each stacked quantity and one for the top bound.
    xdata = data[field_label[0]][x_field_output]
    ydata_cum = np.zeros((2,len(xdata)))

    ### Loop through the y datasets to figure out the extrema
    for i in range(len(field_label)):
        xdata = data[field_label[i]][x_field_output]
        ydata_cum[1] += data[field_label[i]][y_field_output[i]]
        extrema = preserve_extrema(extrema,xdata,ydata_cum[1])
    if log_y_axis=="Auto":
        if extrema[3]/extrema[2] > 1e3:
            log_y_axis="On"
        else:
            log_y_axis="Off"

    ### Reset the cumulative ydata array to have the lowest value found
    ### or zero (for linear plots)
    if log_y_axis == "On":
        ydata_cum = np.zeros((2,len(xdata)))+extrema[2]
    else:
        ydata_cum = np.zeros((2,len(xdata)))

    for i in range(0,num_fields):
        ydata_cum[1] += data[field_label[i]][y_field_output[i]]
        color = cm.jet(1.*i/num_fields)
        pl.fill_between(xdata,ydata_cum[0],ydata_cum[1],color=color)
        if log_y_axis=="On":
            pl.semilogy(xdata,ydata_cum[1],color=color)
        else:
            pl.plot(xdata,ydata_cum[1],color=color)
        legend_list.append(field_label[i])
        # Move our top bound to the bottom bound for our next iteration
        ydata_cum[0] = ydata_cum[1]

    pl.xlim(extrema[0:2])
    if log_y_axis=="On":
        y_log_range = 1.2*np.log10(extrema[3]/extrema[2])
        pl.ylim([extrema[2],extrema[2]*10**y_log_range])
    else:
        pl.ylim([0.,1.2*extrema[3]])
    if len(legend_list) > 0:
        pl.legend(legend_list,2)
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

#####
# testing 
#####
plot_quantity(data, ['Total', 'Level 0', 'Level 1', 'Level 2'], [1,1,1,1], "Total Time", log_y_axis="On")
plot_quantity(data, 'Total', 1, "Total Time")
plot_stack(data, [], 1, "Mean Time (sec)", repeated_field="Level",log_y_axis="Off")
plot_stack(data, ['RebuildHierarchy','SolveHydroEquations'], 1, "Mean Time (sec)",log_y_axis="On" )
