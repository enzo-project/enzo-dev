#!/usr/bin/env python
### chronos.py
### date: 10.13.11
### author: Cameron Hummels
### description: plots performance information from chronos.out

import matplotlib as mpl
import pylab as pl
import numpy as np
from matplotlib import cm

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
    in the begining and end part of the output signal.
    
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

    def __init__(self, filename):
        self.filename = filename
        self.data = self.build_struct(filename)
        self.fields = self.data.keys()
 
    def build_struct(self, filename):
        """
        Build the internal data structure which holds all of the data from
        your external textfile outputted by chronos.

        Parameters
        ----------
        filename : string
            The name of the file used as input

        Returns
        -------
        out : dictionary
            A dictionary of arrays.  Each "label" (e.g. "Total") within the 
            input text file is represented as a key in the dictionary.  The 
            payload for each key is an NxM array where N is the number of 
            cycles over which that label appeared, and M is the number of 
            fields for each label.  The first field for each label is the 
            cycle number.
        """

        ### Open the file, and pull out the data.  
        ### Comments are ignored.  Each block of multiline text is treated 
        ### individually.  The first line of each text block is expected to 
        ### be the cycle_number, which is stored into an array.  Each 
        ### subsequent line uses the first string as a dictionary key, 
        ### with the remaining columns the data for that key.  
        ### A blank line signifies the end of a text block

        ### Create empty data structure which will store all input
        data = {}

        new_block = True
        input = open(filename, "r")

        ### Parse individual lines of data block until empty line, which 
        ### signifies new block of text
        for line in input:
            if line.strip() == '':
                new_block = True
                continue
    
            if not line.startswith("#") and line.strip():
                linelist = line.split()
        
                ### If we're starting a new block of text
                if new_block:
                    if linelist[0] != 'Cycle_Number':
                        exit("Expected Cycle_Number at start of output block")
                    else:
                        cycle = int(linelist[1])
                        new_block = False
                        continue
        
                ### If we've made it this far, cycle is defined and we're
                ### populating our dictionary with data for cycle #x.
                line_key = linelist[0]
        
                ### Convert line_key's underscores to spaces
                line_key = " ".join(line_key.split('_'))
        
                line_value = np.array([cycle] + linelist[1:],dtype='float64') 
                line_value = np.nan_to_num(line_value)  # error checking
        
                ### If line_key is not already in our dictionary, then we 
                ### create a new array as the dictionary payload
        
                if not data.__contains__(line_key):
                    data[line_key] = line_value
        
                ### Otherwise, add a new row to the existing dataset
                else:
                    data[line_key] = np.column_stack((data[line_key],
                                                      line_value))
        input.close()
        return data

    def plot_quantity(self, field_label, y_field_output, y_field_axis_label="",
                      x_field_output=0, x_field_axis_label="Cycle Number",
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
        y_field_output : int or array_like of ints
            The index of the field you wish to plot. Ex: 1 for mean time.
        y_field_axis_label : string, optional
            The y axis label on the resulting plot. Default = ""
        x_field_output : int, optional
            The index of the x data you wish to plot. Default = 0 (Cycles)
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
    
        >>> plot_quantity("RebuildHierarchy", 1, "Mean Time (sec)")
    
        To produce a plot comparing the RebuildHiearchy and SolveHydroEquations
        maximum time taken over the course of the simulation and save it 
        to file "test.png":
    
        >>> plot_quantity(["RebuildHierarchy", "SolveHydroEquations"],
        4, "Maximum Time (sec)", filename="test.png")
    
        To produce a plot comparing the maximum time from RebuildHiearchy and 
        the minimum time from SolveHydroEquations taken over the course of the 
        simulation and save it to file "test.png":
    
        >>> plot_quantity(["RebuildHierarchy", "SolveHydroEquations"],
        [4,3], "Time (sec)", filename="test.png")
    
        To produce a plot comparing the mean time taken by all of the different
        levels over the course of the simulation and save it to file "test.png": 
    
        >>> plot_quantity([], 1, "Mean Time (sec)", filename="test.png", 
        repeated_field="Level")
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
    
        ### If y_field_output is a single index, then replace it with a list of
        ### identical indices
        if not is_listlike(y_field_output):
            y_field_output = num_fields*[y_field_output]
    
        ### Create a normalization vector to use on each vector
        ### before plotting.  In non-fractional case, this vector is 1.
        if fractional:
            norm = 1./data['Total']
        else:
            norm = np.ones(data['Total'].shape)    

        ### Loop through the y datasets to figure out the extrema
        for i in range(len(field_label)):
            xdata = data[field_label[i]][x_field_output]
            ydata = data[field_label[i]][y_field_output[i]] * \
                    norm[y_field_output[i]]
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

            xdata = data[field_label[i]][x_field_output]
            ydata = data[field_label[i]][y_field_output[i]] * \
                    norm[y_field_output[i]]
            if smooth_len:
                ydata = smooth(ydata,smooth_len)
            if log_y_axis=="On":
                pl.semilogy(xdata,ydata)
            else:
                pl.plot(xdata,ydata)
            if not bounds == "Off":
                zerodata = np.zeros(len(ydata))
                if bounds == "minmax":
                    min_bound = data[field_label[i]][3] * \
                    norm[3]
                    max_bound = data[field_label[i]][4] * \
                    norm[4]
                else:
                    min_bound = ydata - data[field_label[i]][2] * \
                    norm[2]
                    max_bound = ydata + data[field_label[i]][2] * \
                    norm[2]
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
            pl.ylim([extrema[2],extrema[2]*10**y_log_range])
        else:
            pl.ylim([0.,1.2*extrema[3]])

        ### Set up the legend with smaller text to allow for more entries
        legend = pl.legend(legend_list,2)
        ltext  = legend.get_texts()
        pl.setp(ltext, fontsize='small')

        pl.xlabel(x_field_axis_label)
        if fractional:
            pl.ylabel(y_field_axis_label + ", as Fraction of Total")
        else:
            pl.ylabel(y_field_axis_label)
        pl.savefig(filename)
        pl.clf()
        
    def plot_stack(self, field_label, y_field_output, y_field_axis_label="",
                   x_field_output=0, x_field_axis_label="Cycle Number",
                   filename="chronos.png", repeated_field="", 
                   log_y_axis="Auto", smooth_len=0, fractional=False):
        """
        Produce a plot for the given label/outputs where each quantity is 
        stacked on top of the previous quantity.
    
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
        >>> plot_stack(["Level 0", "Level 1", "Level 2"], 1, "Mean Time (sec)")
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
    
        ### If y_field_output is a single index, then replace it with a list of
        ### identical indices
        if not is_listlike(y_field_output):
            y_field_output = num_fields*[y_field_output]
    
        ### Since we're stacking (and making plots from a bottom bound to an
        ### upper bound for each y value, we need to cumulate it as we stack.
        ### ydata_cum is the same size as xdata, but it has two indices, one
        ### for the bottom bound of each stacked quantity and one for the top 
        ### bound.  It is cumulatively added as we loop through the plotted
        ### quantities.
        xdata = data[field_label[0]][x_field_output]
        ydata_cum = np.zeros((2,len(xdata)))

        ### Create a normalization vector to use on each vector
        ### before plotting.  In non-fractional case, this vector is 1.
        if fractional:
            norm = 1./data['Total']
        else:
            norm = np.ones(data['Total'].shape)    
    
        ### Loop through the y datasets to figure out the extrema
        for i in range(len(field_label)):
            xdata = data[field_label[i]][x_field_output]
            ydata = data[field_label[i]][y_field_output[i]] * \
                    norm[y_field_output[i]]
            if smooth_len:
                ydata_cum[1] += smooth(ydata, smooth_len)
            else:
                ydata_cum[1] += ydata
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
            ydata = data[field_label[i]][y_field_output[i]] * \
                    norm[y_field_output[i]]
            if smooth_len:
                ydata_cum[1] += smooth(ydata, smooth_len)
            else:
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
            pl.ylim([extrema[2],extrema[2]*10**y_log_range])
        else:
            pl.ylim([0.,1.2*extrema[3]])
        if len(legend_list) > 0:
            legend = pl.legend(legend_list,2)
            ltext  = legend.get_texts()
            pl.setp(ltext, fontsize='small')

        pl.xlabel(x_field_axis_label)
        if fractional:
            pl.ylabel(y_field_axis_label + ", as Fraction of Total")
        else:
            pl.ylabel(y_field_axis_label)
        pl.savefig(filename)
        pl.clf()
        
########
### XXX -- Things left to do:

### Perhaps reading the data all into memory first, then creating the dicts
### and the arrays before filling them  
    
### Look into recarrays for giving specific handles on data entries

### Better docstrings and better examples
########

### If chronos.py is invoked from the command line, these are its default
### behaviors:
### -- Build a chronos object from the provided filename
### -- Make some handy plots and write them out

if __name__ == "__main__":
    from optparse import OptionParser
    usage = "usage: %prog <.out file>"
    parser = OptionParser(usage)
    (opts, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")
    filename = args[0]

    ### Build a chronos object from the data and generate some plots
    c = chronos(filename)
    c.plot_quantity(['Total'], 1, "Total Time", repeated_field="Level", 
                    log_y_axis="On", filename='c1.png')
    c.plot_quantity(['Total'], 1, "Total Time", repeated_field="Level",
                    log_y_axis="On", filename='c1f.png', 
                    fractional=True)
    c.plot_quantity(['Total'], 1, "Total Time", repeated_field="Level",
                    log_y_axis="Auto", filename='c1s.png',
                    smooth_len=11, bounds="minmax", fractional=True)
    c.plot_quantity('Total', 1, "Total Time", filename='c2.png', 
                    bounds="minmax", fractional=True)
    c.plot_quantity('Total', 1, "Total Time", filename='c2s.png',
                    smooth_len=15, bounds="minmax", fractional=True)
    c.plot_stack([], 1, "Mean Time (sec)", repeated_field="Level", 
                 log_y_axis="Off", filename='c3.png')
    c.plot_stack([], 1, "Mean Time (sec)", repeated_field="Level", 
                 log_y_axis="Off", filename='c3f.png', fractional=True)
    c.plot_stack([], 1, "Mean Time (sec)", repeated_field="Level", 
                 log_y_axis="Off", filename='c3s.png', smooth_len=11)
    c.plot_stack(['RebuildHierarchy','SolveHydroEquations'], 1, 
                 "Mean Time (sec)",log_y_axis="On", filename='c4.png')
    c.plot_stack(['Total','RebuildHierarchy','SolveHydroEquations'], 1, 
                 "Mean Time (sec)",log_y_axis="On", filename='c4s.png',
                 smooth_len=19)
    c.plot_stack(['Total','RebuildHierarchy','SolveHydroEquations'], 1, 
                 "Mean Time (sec)",log_y_axis="On", filename='c4sf.png',
                 smooth_len=19, fractional=True)
    c.plot_stack(['Total','RebuildHierarchy','SolveHydroEquations'],1, "Mean Time",
                 filename='c4sfa.png',smooth_len=19,fractional=True, log_y_axis="On",
                 repeated_field="Level")
