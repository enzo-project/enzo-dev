#!/usr/bin/env python

# chronos.py
# date: 10.13.11
# author: Cameron Hummels
# description: polls enzo performance information from XXX. 

# Use the cool OptionParser (only works with Python 2.6, broken with python 2.7
from optparse import OptionParser

usage = "usage: %prog [options] <XXX file>"
parser = OptionParser(usage)
#parser.add_option("-t", "--title", dest="title", action="store", 
#                  default=None, type="string", help="set title in plot")
#parser.add_option("-x", "--x-title", dest="x_title", action="store", 
#                  default=None, type="string", help="set x-axis title in plot")
#parser.add_option("-y", "--y-title", dest="y_title", action="store", 
#                  default=None, type="string", help="set y-axis title in plot")
#parser.add_option("-l", "--log-type", dest="log_type", action="store", 
#                  help="set plot-type 0=lin, 1=log-lin, 2=lin-log, 3=log-log",
#                  default=0, type="int")
#parser.add_option("-w", "--write", dest="write_fig", action="store", 
#                  default=None, type="string", 
#                  help="write fig to an output .png file instead of to screen")
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

#def fill_between(ax, x, y1, y2, **kwargs):
#    # add x,y2 in reverse order for proper polygon filling
#    verts = zip(x,y1) + [(x[i], y2[i]) for i in range(len(x)-1,-1,-1)]
#    poly = Polygon(verts, **kwargs)
#    ax.add_patch(poly)
#    ax.autoscale_view()
#    return poly

filename = args[0]

### Create empty lists
#x = []
#y = []
time_list = []
data_list = []

### Open the filename, and for each non-comment line pull the x and y 
### column values out and append them to our lists: x and y
input = open(filename, "r")

#import pdb; pdb.set_trace()
for line in input:
    data_sublist = []
    if not line.startswith("#") and line.strip():
        linelist = line.split()
        if linelist[0] == 'Total':
            break
        time_list.append(linelist[1])
        data_sublist.append(linelist[2:7])
        ### XXX Poll data for individual processors in verbose case
        for line in input:
            if line.strip() == '':
                break
            if not line.startswith("#"):
                linelist = line.split()
                data_sublist.append(linelist[2:7])
        data_list.append(data_sublist)            

### XXX Poll data for last paragraph "Total" 
#for line in input:
#    if not line.startswith("#"):
#        linelist = line.split()
#        data_sublist.append(linelist[2:7])
#    data_list.append(data_sublist)            


input.close()

times = np.array(time_list,dtype='int32')
data = np.array(data_list,dtype='float32')

### data is an array of dimension timesteps x levels+1 x 5 

### the first dimension is over timesteps
### the second dimension is over levels (but the first entry is the sum over all levels)

### the 5 values in the 3rd dimension are: 
### (num_cells, mean_cells/sec/proc, stddev_cells/sec/proc, min_cells/sec/proc, max_cells/sec/proc)

pl.figure(1)

# plot total timestep vs mean cells / sec / proc
# with min/max/stddev plotted over
#pl.plot(times, data[:,0,1])
#    pylab.xlabel(opts.x_title)

# plot time on each level vs timestep
#for i in range(data.shape[1]):  
#    pl.plot(times, data[:,i,0])

# plot time on each processsor vs timestep
#for i in range(data.shape[2]):  
#    pl.plot(times, data[:,0,i])

#for level in range(data.shape[2])
#    pl.plot(timesteps, data[,:,level], 'b-')

# Add extras to the plot and save or display it?
#if not opts.title:
#    opts.title = "%s -- Column %i vs Column %i\n" % (filename,x_column,y_column)
#pylab.suptitle(opts.title) 
# perhaps I need to check if x_title is defined before doing the following line
#if opts.x_title:
#    pylab.xlabel(opts.x_title)
#if opts.y_title:
#    pylab.ylabel(opts.y_title)
#
#if opts.write_fig:
#    pylab.savefig(opts.write_fig)
#else:

input.close()


### data is an array of dimension timesteps x levels+1 x 5 

### the first dimension is over timesteps
### the second dimension is over levels (but the first entry is the sum over all levels)

### the 5 values in the 3rd dimension are: 
### (num_cells, mean_cells/sec/proc, stddev_cells/sec/proc, min_cells/sec/proc, max_cells/sec/proc)

#pl.figure(1)

cumulate_grids = np.zeros((2,len(time_list)))

cumulate_grids[1,:] += data[:,1,0]
#import pdb;pdb.set_trace()
# plot timestep vs num_cells on diff levels
legend_list = []
levels = range(2,data.shape[1])
pl.fill_between(times,cumulate_grids[0,:],cumulate_grids[1,:],color=cm.jet(0))
pl.plot(times,cumulate_grids[1,:],color=cm.jet(0))
legend_list.append('Level %i' % 0)
cumulate_grids[0,:] = cumulate_grids[1,:]

for i in levels:
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

# can remove in final, but makes more readable.  for all levels.
mean = data[:,0,1]
stddev = data[:,0,2]
min = data[:,0,3]
max = data[:,0,4]
sigmin = mean - stddev
sigmax = mean + stddev


# plot total timestep vs mean cells / sec / proc
# with min/max/stddev plotted over
cumulate_grids = np.zeros((2,len(time_list)))

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

# plot total timestep vs mean cells / sec / proc
# with min/max/stddev plotted over
pl.plot(times, mean)
minmax = pl.fill_between(times,min,max,facecolor='g')
pl.xlim([np.min(times),np.max(times)])
pl.ylim([0,1.1*np.max(max)])
minmax.set_alpha(0.5)
#sigmax = pl.fill_between(times,sigmin,sigmax,facecolor='b')
#sigmax.set_alpha(0.5)
pl.xlabel("Cycle Number")
pl.ylabel("Walltime Spent (seconds)")
pl.suptitle("Aggragate Walltime Spent with min/max for Processors") 
pl.savefig("time_total.png")
pl.clf()
#pl.show()

# plot times on each level vs timestep
#for i in range(data.shape[1]):  
#    pl.plot(timesteps, data[:,i,0])

# plot times on each processsor vs timestep
#for i in range(data.shape[2]):  
#    pl.plot(timesteps, data[:,0,i])

#for level in range(data.shape[2])
#    pl.plot(timesteps, data[,:,level], 'b-')

# Add extras to the plot and save or display it?
#if not opts.title:
#    opts.title = "%s -- Column %i vs Column %i\n" % (filename,x_column,y_column)
#pylab.suptitle(opts.title) 
# perhaps I need to check if x_title is defined before doing the following line
#if opts.x_title:
#    pylab.xlabel(opts.x_title)
#if opts.y_title:
#    pylab.ylabel(opts.y_title)
#
#if opts.write_fig:
#    pylab.savefig(opts.write_fig)
#else:
pl.close()
