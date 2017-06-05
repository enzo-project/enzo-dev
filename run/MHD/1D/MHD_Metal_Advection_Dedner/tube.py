import yt
import numpy as np
import pdb
import pylab
import matplotlib.pyplot as plt
import matplotlib
import h5py
import copy
nar = np.array

        
def tube(ds_list,fields=None,times=None, points=[(0.0,0.0,0.0),(1.0,0.0,0.0)],width=2,
         filename='FILE.png',fig_args={},plot_args={},legend=False,
         renorm=False,
         offsets=None,
         axis=0,coord=(0.505,0.505),
         debug=0,delta=False,
         labels=None, xlim = None): 
    """Plots 1d profiles for each ds in *ds_list*, for each field in *fields* and time in *times*.
    One plot per field, with oober and time overplotted in each field window.
    Uses the AMRRay from yt, with start and end points given by *points*.
    If *fields* is not given, the minimum set is taken from oober.fields for each oober in *ooberlist*.
    If *times* is not given, each oober uses the values in oober.frames
    *renorm* removes the mean of each field.
    *offsets*, if not None, should be a list of multiplicitive offsets for each oober
    *axis*,*coord* gives the arguments to ortho_ray
    *delta* = True takes the difference of oober[0] with oober[1:]
    """
    #global plt
    #ds_list = ensure_list(ds_list)
    fig = plt.figure()
    ray_set = {}

    # Set up times and fields to loop over

    if fields == None:
        fields = ['density']
    if legend == 1: 
        extra_for_legend = 1
    else:
        extra_for_legend = 0
    n_rows = int(1.*(len(fields)+extra_for_legend)/width+0.5)
    # Make list of rays
    #pdb.set_trace()
    if 0:
        """this is a dirty hack."""
        if labels is None:
            labs = range(len(ds_list))
        else:
            labs = labels
        for ds, lab in zip(ds_list, labs):
            ray_set[lab]= ds.ortho_ray(axis,coord) 
        return ray_set
    else:
        for ds in ds_list:
            ray_set[ds]= ds.ortho_ray(axis,coord) 

    #make the figure
    fig = plt.figure(**fig_args)
    first_ax = None
    Handels=[]
    
    units={'density':'code_mass/code_length**3','pressure':'code_mass/(code_length*code_time**2)'}
    for x in 'xyz':
        units['B'+x] = 'code_magnetic'
        units[x+'-velocity']='code_length/code_time'

    for i,field in enumerate(fields):

        fiducial_field = []
        ax = fig.add_subplot(n_rows,width,i+1)
        if debug > 0:
            print "n_rows %d width %d i %d field %s"%(n_rows,width,i+1, field)
        counter = -1
        for n_ds,ds in enumerate(ds_list):
            counter += 1
            this_ray = ray_set[ds]
            this_x = copy.copy(this_ray['xyz'[axis]].v[:])
            sort_x = np.argsort(this_x)
            this_x = this_x[sort_x]
            this_y = copy.copy(this_ray[field])
            if units.has_key(field):
                this_y = this_y.in_units(units[field])
            this_y =this_y.v[sort_x]
            if delta:
                if n_ds==0:
                    fiducial_field = this_y
                    continue
                else:
                    this_y -= fiducial_field
                    nonzero = fiducial_field != 0
                    this_y[nonzero] /= fiducial_field[nonzero]


            if hasattr(renorm,'has_key'):
                if renorm.has_key(field):
                    this_y -= renorm[field]
            if offsets is not None:
                this_y *= offsets[n_ds]


            if first_ax is None:
                plot_args['label'] =  '%s %0.2f'%(labels[n_ds],ds['InitialTime'])
            L = ax.plot(this_x,this_y,**plot_args)
            if this_y.min() > 0:
                ax.set_yscale('log')

        if first_ax is None:
            first_ax = ax
            Handels, LabelsToUse = first_ax.get_legend_handles_labels()
        if xlim is not None:
            ax.set_xlim(xlim)
        ax.set_ylabel(field)

    if legend == 1:
        ax = fig.add_subplot(n_rows,width,i+2)
        ax.legend( Handels, LabelsToUse,loc = 0)
    elif legend == 2:
        first_ax.legend(Handels,LabelsToUse)
    fig.savefig(filename)
    plt.close(fig)
    print filename
