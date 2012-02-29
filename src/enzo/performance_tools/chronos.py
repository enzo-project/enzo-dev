#!/usr/bin/env python
### Chronos.py
### Date: 10.13.11
### Author: Cameron Hummels and Sam Skillman
### Description: 

### Chronos is a module for plotting the performance information which 
### comes out of Enzo.  It consists of a few universal helper functions
### and a new class: chronos.  You can use it one of two ways.  You can
### import this library in python, create a chronos object, and then
### create a few plots using the plot_quantity and plot_stack functions,
### like this:

### $ python
### >>> import chronos as c
### >>> object = c.chronos('chronos.out')
### >>> object.plot_quantity('Total','mean time','Mean Time (sec)')
### >>> object.plot_stack([],'mean time','Mean Time (sec)', \ 
###                       repeated_field='Level')

### Or you can just call this python file from the command line after
### editing the source file to have it print out whatever plots you want
### it to generate (some defaults are included) like this:

### $ python chronos.py chronos.out

import matplotlib as mpl
import pylab as pl
import numpy as np
from matplotlib import cm

def is_listlike(obj):
    """
    Checks to see if an object is listlike (but not a string)
    """
    return not isinstance(obj, basestring)
    
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

def smooth(x, window_len=11, window='hanning'):
    """
    Smooth a 1D array using a window with requested size.

    This function taken from http://www.scipy.org/Cookbook/SignalSmooth
    
    CBH modified it to cutoff the (window_len-1)/2 values at each end of the
    final 1D array (the reflected part's remnants), so that the final
    array size is the same as the original array size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the beginning and end part of the output signal.
    
    Parameters
    ----------
    x : 1D array_like
        The input signal 
    window_len: int, optional
        The dimension of the smoothing window; should be an odd integer
    window: string
        the type of smoothing window to use: 'flat', 'hanning', 'hamming', 
        'bartlett', 'blackman'. flat window will produce a moving 
        average smoothing.

    Returns
    -------
    out : 1D array
        The smoothed version of input x

    See Also
    --------
    np.hanning, np.hamming, np.bartlett, np.blackman, np.convolve
    scipy.signal.lfilter
 
    Examples
    --------
    >>> t = np.linspace(-4,4,100)
    >>> x = np.sin(t)
    >>> x_noisy = x + pl.randn(len(t))*0.1
    >>> x_smooth = smooth(x_noisy)
    """
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."

    if window_len<3:
        return x

    if window_len%2 == 0:
        raise ValueError, "smooth requires an odd-valued window_len." 

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    s = np.r_[x[window_len-1:0:-1], x, x[-1:-window_len:-1]]
    if window == 'flat': # moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.'+window+'(window_len)')

    y = np.convolve(w/w.sum(), s, mode='valid')
    clip_len = (window_len-1)/2
    return y[clip_len:len(y)-clip_len]


class chronos:
    """
    This simple class stores the performance information collected by Enzo.
    Enzo puts this information in an ASCII file typically called "chronos.out".
    That file is called by the constructor for this class and converted
    into the internal structure called "data".  
    
    "data" is a dictionary, where each key is one of the field labels 
    from "chronos.out" (i.e. the first word from a line, e.g. "Level 0").  
    The payload for each value is a python record array of N dimension, 
    where N is the number of cycles for that key.

    For each cycle, there are M entries. In the case of the "Level X" and 
    "Total" keys, these M entries are:

    "cycle"
    "mean time", 
    "std_dev time", 
    "min time", 
    "max time", 
    "cells", 
    "grids", 
    "cells/processor/sec"

    In all other cases, these M entries are:
    "cycle"
    "mean time", 
    "std_dev time", 
    "min time", 
    "max time"

    Since this is a record array, you can use these entries as the indices
    when indexing the array (e.g. data['Total']['cycle']).

    """
    def __init__(self, filename):
        self.filename = filename
        self.data = self.build_struct(filename)
        self.fields = self.data.keys()
 
    def build_struct(self, filename):
        """
        Build the internal data structure which holds all of the data from
        your external textfile outputted by enzo.  It allocates the memory 
        all at once instead of building the structure line by line.

        Parameters
        ----------
        filename : string
            The name of the file used as input

        Returns
        -------
        out : dictionary
            A dictionary of recarrays.  Each "label" (e.g. "Total") within the 
            input text file is represented as a key in the dictionary.  The 
            payload for each key is an N-dim record array where N is the 
            number of cycles over which that label appeared. See the
            docstrings for the chronos class for more information.
        """

        ### Open the file, and pull out the data.  
        ### Comments are ignored.  Each block of multiline text is treated 
        ### individually.  The first line of each text block is expected to 
        ### be the cycle_number, which is stored into an array.  Each 
        ### subsequent line uses the first string as a dictionary key, 
        ### with the remaining columns the data for that key.  
        ### A blank line signifies the end of a text block
        key_tally = {}
        key_counter = {}
        data = {}

        input = open(filename, "r")
        input_list = input.readlines()
        input.close()

        ### First run through the file contents counting up all instances of 
        ### each key (first word of each line) so as to allocate sufficiently
        ### big arrays in the structure to house the data.
        for line in input_list:
            if not line.startswith("#") and line.strip():
                line_list = line.split()
                line_key = line_list[0]
                line_key = " ".join(line_key.split('_'))
                if not key_tally.__contains__(line_key):
                    key_tally[line_key] = 0
                key_tally[line_key] += 1
        del key_tally['Cycle Number'] # Don't want it as a key

        ### Now build the dictionary with the appropriate memory and 
        ### retraverse the input file to fill it with the input information
        for key in key_tally.keys():
            key_counter[key] = 0
            if key == "Total" or key.startswith('Level'):
                records = [('cycle', 'float'), ('mean time', 'float'),
                           ('std_dev time', 'float'), ('min time', 'float'),
                           ('max time', 'float'), ('cells', 'float'),
                           ('grids', 'float'), ('cells/processor/sec', 'float')]
            else:
                records = [('cycle', 'float'), ('mean time', 'float'),
                           ('std_dev time', 'float'), ('min time', 'float'),
                           ('max time', 'float')]

            data[key] = np.zeros(key_tally[key], dtype=records)

        for line in input_list:
            if not line.startswith("#") and line.strip():
                line_list = line.split()
   
                if line_list[0] == 'Cycle_Number':
                    cycle = int(line_list[1])
                    continue
    
                ### If we've made it this far, cycle is defined and we're
                ### populating our dictionary with data for cycle #x.
                line_key = line_list[0]
    
                ### Convert line_key's underscores to spaces
                line_key = " ".join(line_key.split('_'))
        
                line_value = np.array([cycle] + line_list[1:],dtype='float64') 
                line_value = np.nan_to_num(line_value)  # error checking
    
                data[line_key][key_counter[line_key]] = line_value
                key_counter[line_key] += 1
        return data

    def plot_quantity(self, field_label, y_field_index, 
                      y_field_axis_label="", x_field_index='cycle', 
                      x_field_axis_label="Cycle Number",
                      filename="chronos.png", repeated_field="", 
                      log_y_axis="Auto", smooth_len=0, bounds="Off",
                      fractional=False):
        """
        Produce a plot for the given quantity(s) from the chronos data.
    
        Parameters
        ----------
        field_label : string or array_like of strings
            The label of the field you wish to plot.  If you wish to plot
            multiple fields, enumerate them in an array or tuple. Ex: "Level 0"
        y_field_index : string
            The index of the field you wish to plot on the y axis. 
            Ex: "mean time", "std_dev time", "min time", "max time",
            "cells", "grids", "cells/processor/sec".
            If you have a single value for many field_labels, it is assumed
            that such a value will index all of them.  If you have an array_like
            structure of ints, each one will index the corresponding key in 
            field_lable.  
        y_field_axis_label : string, optional
            The y axis label on the resulting plot. Default = the y_field_index
            of the recarray.
        x_field_index : string, optional
            The index of the field you wish to plot on the x axis.
            Default = "cycle"
        x_field_axis_label : string, optional
            The x axis label on the resulting plot. Default = "Cycle Number"
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
        smooth_len : int, optional
            This value controls the amount by which smoothing occurs over
            N consecutive cycles of data.  Default = 0 (i.e. None). 
            Must be an odd number (recommended 5-11)
        bounds : string, optional
            This controls whether to overplot additional bounding data over
            the existing plotted quantities.  Valid values of this variable
            are "minmax", "sigma" and "Off".  "minmax" overplots the minima and
            maxima bounds, whereas "sigma" plots the mean +/- 1 sigma bounds.
        fractional : bool, optional
            When set to true, the plotted values shown in fractions of the 
            equivalent field in "Total".
    
        See Also
        --------
        plot_stack
    
        Examples
        --------
        To produce a simple plot of the mean time taken over the course of the
        simulation and save it to chronos.png:
    
        >>> plot_quantity("RebuildHierarchy", "mean time", "Mean Time (sec)")
    
        To produce a plot comparing the RebuildHiearchy and SolveHydroEquations
        maximum time taken over the course of the simulation and save it 
        to file "test.png":
    
        >>> plot_quantity(["RebuildHierarchy", "SolveHydroEquations"],
        "max time", "Maximum Time (sec)", filename="test.png")
    
        To produce a plot comparing the maximum time from RebuildHiearchy and 
        the minimum time from SolveHydroEquations taken over the course of the 
        simulation and save it to file "test.png":
    
        >>> plot_quantity(["RebuildHierarchy", "SolveHydroEquations"],
        ["max time", "min time"], "Time (sec)", filename="test.png")
    
        To produce a plot comparing the mean time taken by all of the different
        levels over the course of the simulation and save it to file "test.png": 
        >>> plot_quantity([], "mean time", "Mean Time (sec)", 
        filename="test.png", repeated_field="Level")
        """
        data = self.data
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
    
        ### If y_field_index is a single index, then replace it with a list of
        ### identical indices
        if not is_listlike(y_field_index):
            y_field_index = num_fields*[y_field_index]
    
        ### Create a normalization vector to use on each vector
        ### before plotting.  In non-fractional case, this vector is 1.
        if fractional:
            norm = data['Total']
        else:
            records = [('cycle', 'float'), ('mean time', 'float'),
                       ('std_dev time', 'float'), ('min time', 'float'),
                       ('max time', 'float'), ('cells', 'float'),
                       ('grids', 'float'), ('cells/processor/sec', 'float')]
            norm = np.ones(data['Total'].shape, dtype=records)    

        ### Loop through the y datasets to figure out the extrema
        for i in range(len(field_label)):
            xdata = data[field_label[i]][x_field_index]
            ydata = data[field_label[i]][y_field_index[i]] / \
                    norm[y_field_index[i]]
            if smooth_len:
                ydata = smooth(ydata,smooth_len)
            extrema = preserve_extrema(extrema,xdata,ydata)
        if log_y_axis=="Auto":
            if extrema[3]/extrema[2] > 1e3:
                log_y_axis="On"
            else:
                log_y_axis="Off"
    
        ### Now for the actual plotting
        for i in range(len(field_label)):
            color = cm.jet(1.*i/num_fields)
            xdata = data[field_label[i]][x_field_index]
            ydata = data[field_label[i]][y_field_index[i]] / \
                    norm[y_field_index[i]]
            if smooth_len:
                ydata = smooth(ydata,smooth_len)
            if log_y_axis=="On":
                pl.semilogy(xdata,ydata,color=color)
            else:
                pl.plot(xdata,ydata,color=color)
            if not bounds == "Off":
                zerodata = np.zeros(len(ydata))
                if bounds == "minmax":
                    min_bound = data[field_label[i]]["min time"] / \
                    norm["min time"]
                    max_bound = data[field_label[i]]["max time"] / \
                    norm["max time"]
                else:
                    min_bound = ydata - data[field_label[i]]["std_dev time"] / \
                    norm["std_dev time"]
                    max_bound = ydata + data[field_label[i]]["std_dev time"] / \
                    norm["std_dev time"]
                if smooth_len:
                    min_bound = smooth(min_bound, smooth_len)
                    max_bound = smooth(max_bound, smooth_len)
                fillin = pl.fill_between(xdata,min_bound,max_bound,
                                         facecolor='0.5')
                fillin.set_alpha(0.5)
            legend_list.append(field_label[i])
    
        pl.xlim(extrema[0:2])
        if log_y_axis=="On":
            y_log_range = 1.2*np.log10(extrema[3]/extrema[2])
            ### To assure there is a labeled tick mark on the y-axis
            if y_log_range < 1.:
                y_log_range = 1.
            pl.ylim([extrema[2],extrema[2]*10**y_log_range])
        else:
            pl.ylim([0.,1.2*extrema[3]])

        ### Set up the legend with smaller text to allow for more entries
        legend = pl.legend(legend_list,loc='best')
        ltext  = legend.get_texts()
        pl.setp(ltext, fontsize='xx-small')

        pl.xlabel(x_field_axis_label)
        if not y_field_axis_label:
            y_field_axis_label=y_field_index
        if fractional:
            pl.ylabel(y_field_axis_label + ", as Fraction of Total")
        else:
            pl.ylabel(y_field_axis_label)
        pl.savefig(filename)
        pl.clf()
        
    def plot_stack(self, field_label, y_field_index, 
                   y_field_axis_label="", x_field_index='cycle', 
                   x_field_axis_label="Cycle Number", filename="chronos.png", 
                   repeated_field="", log_y_axis="Auto", smooth_len=0, 
                   fractional=False):
        """
        Produce a plot for the given label/indices where each quantity is 
        stacked on top of the previous quantity.
    
        Parameters
        ----------
        field_label : array_like of strings
            The labels of the fields you wish to plot.  If you wish to plot
            a single field, then use plot_quantity instead.
        y_field_index : string or array_like of strings
            The index of the field you wish to plot on the y axis. 
            Ex: "mean time", "std_dev time", "min time", "max time",
            "cells", "grids", "cells/processor/sec".
            If you have a single value for many field_labels, it is assumed
            that such a value will index all of them.  If you have an array_like
            structure of ints, each one will index the corresponding key in 
            field_lable.  
        y_field_axis_label : string, optional
            The y axis label on the resulting plot. Default = the y_field_index
            of the recarray.
        x_field_index : string, optional
            The index of the field you wish to plot on the x axis.
            Default = "cycle"
        x_field_axis_label : string, optional
            The x axis label on the resulting plot. Default = "Cycle Number"
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
        smooth_len : int, optional
            This value controls the amount by which smoothing occurs over
            N consecutive cycles of data.  Default = 0 (i.e. None)
            Must be an odd number (recommended 5-11)
        fractional : bool, optional
            When set to true, the plotted values shown in fractions of the 
            equivalent field in "Total".
    
        See Also
        --------
        plot_quantity
    
        Examples
        --------
        >>> plot_stack(["Level 0", "Level 1", "Level 2"], "mean time", 
                       "Mean Time (sec)")
        """
        data = self.data
        extrema = np.zeros(5)
        legend_list = []

        ### Convert plots of single quantities to list format for homogenous 
        ### processing.
        if not is_listlike(field_label):
            field_label = [field_label]
        
        ### If a repeated_field, figure out how many repeated fields there are.
        ### including any that were defined in the original field_label arg.
        if repeated_field:
            i = 0
            while data.has_key(repeated_field + ' %i' % i):
                field_label.append(repeated_field + ' %i' % i)
                i += 1
        num_fields = len(field_label)
    
        ### If y_field_index is a single index, then replace it with a list of
        ### identical indices
        if not is_listlike(y_field_index):
            y_field_index = num_fields*[y_field_index]
    
        ### Since we're stacking (and making plots from a bottom bound to an
        ### upper bound for each y value, we need to cumulate it as we stack.
        ### ydata_cum is the same size as xdata, but it has two indices, one
        ### for the bottom bound of each stacked quantity and one for the top 
        ### bound.  It is cumulatively added as we loop through the plotted
        ### quantities.
        xdata = data[field_label[0]][x_field_index]
        ydata_cum = np.zeros((2,len(xdata)))

        ### Create a normalization vector to use on each vector
        ### before plotting.  In non-fractional case, this vector is 1.
        if fractional:
            norm = data['Total']
        else:
            records = [('cycle', 'float'), ('mean time', 'float'),
                       ('std_dev time', 'float'), ('min time', 'float'),
                       ('max time', 'float'), ('cells', 'float'),
                       ('grids', 'float'), ('cells/processor/sec', 'float')]
            norm = np.ones(data['Total'].shape, dtype=records)    
    
        ### Loop through the y datasets to figure out the extrema
        for i in range(len(field_label)):
            xdata = data[field_label[i]][x_field_index]
            ydata = data[field_label[i]][y_field_index[i]] / \
                    norm[y_field_index[i]]
            if smooth_len:
                ydata = smooth(ydata, smooth_len)
            ydata_cum[1] += ydata
            extrema = preserve_extrema(extrema,xdata,ydata_cum[1])
        if log_y_axis=="Auto":
            if extrema[3]/extrema[2] > 1e3:
                log_y_axis="On"
            else:
                log_y_axis="Off"
    
        ### Reset the cumulative ydata array to have the lowest value found
        ### or zero (for linear plots)
        ydata_cum = np.zeros((2,len(xdata)))
        if log_y_axis == "On":
            ydata_cum[0] += extrema[2]
    
        for i in range(0,num_fields):
            ydata = data[field_label[i]][y_field_index[i]] / \
                    norm[y_field_index[i]]
            if smooth_len:
                ydata = smooth(ydata, smooth_len)
            ydata_cum[1] += ydata
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
            ### To assure there is a labeled tick mark on the y-axis
            if y_log_range < 1.:
                y_log_range = 1.
            pl.ylim([extrema[2],extrema[2]*10**y_log_range])
        else:
            pl.ylim([0.,1.2*extrema[3]])
        if len(legend_list) > 0:
            legend = pl.legend(legend_list,loc='best')
            ltext  = legend.get_texts()
            pl.setp(ltext, fontsize='xx-small')

        pl.xlabel(x_field_axis_label)
        if not y_field_axis_label:
            y_field_axis_label=y_field_index
        if fractional:
            pl.ylabel(y_field_axis_label + ", as Fraction of Total")
        else:
            pl.ylabel(y_field_axis_label)
        pl.savefig(filename)
        pl.clf()
        
### If chronos.py is invoked from the command line, these are its default
### behaviors:
### -- Build a chronos object from the provided filename
### -- Make some handy plots and write them out

if __name__ == "__main__":
    from optparse import OptionParser
    usage = "usage: %prog <.out file>"
    parser = OptionParser(usage)
    parser.add_option("-s","--smooth",dest="nsmooth",type='int',
                      default=0,help="Set number of cycles over which to smooth (odd)")
    (opts, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")
    filename = args[0]

    ### Build a chronos object from the data and generate some default plots
    c = chronos(filename)
    c.plot_quantity('Total', 'mean time', "Mean Time (sec)", 
                    repeated_field="Level", filename='c1.png',
                    smooth_len=opts.nsmooth, bounds='minmax')
    c.plot_quantity('Total', 'mean time', "Mean Time (sec)", 
                    repeated_field="Level", filename='c2.png',
                    smooth_len=opts.nsmooth, bounds='minmax', fractional=True)
    c.plot_stack([], 'mean time', "Mean Time (sec)", repeated_field="Level", 
                 filename='c3.png', smooth_len=opts.nsmooth)
    c.plot_stack(['RebuildHierarchy','SolveHydroEquations'], 'mean time', 
                 "Mean Time", filename='c4.png', smooth_len=opts.nsmooth,
                 fractional=True)
    c.plot_stack([], 'cells' ,'Number of Cells', repeated_field="Level",
                 filename='c5.png', smooth_len=opts.nsmooth)
    c.plot_quantity(['Total'], 'cells/processor/sec', 'Cells/sec/processor', 
                    repeated_field='Level', filename='c6.png', 
                    smooth_len=opts.nsmooth)
