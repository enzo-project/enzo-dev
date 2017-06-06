# numpy-based error-checking script for Turner & Stone equilibration tests
# Daniel R. Reynolds, reynolds@smu.edu

# imports
import numpy as np

# set the total number of snapshots
nt = 100

# set the solution tolerance
tol = 0.01

# load reference solution from disk
t_ref = np.load('times.npy')
e_ref = np.load('e_sol.npy')

# define some helpful functions
def get_params(file):
    """Returns t, dUnit, tUnit, lUnit from a given parameter file"""
    import shlex
    f = open(file)
    for line in f:
        text = shlex.split(line)
        if ("InitialTime" in text):
            tval = float(text[len(text)-1])
        elif ("DensityUnits" in text):
            dUnit = float(text[len(text)-1])
        elif ("TimeUnits" in text):
            tUnit = float(text[len(text)-1])
        elif ("LengthUnits" in text):
            lUnit = float(text[len(text)-1])
    return [tval, dUnit, tUnit, lUnit]

def load_vals(tdump):
    """Returns Eg, etot, tval from a given data dump"""
    import h5py
    sdump = repr(tdump).zfill(4)
    pfile = 'DD' + sdump + '/data' + sdump
    hfile = pfile + '.cpu0000'
    tval, dUnit, tUnit, lUnit = get_params(pfile)
    f = h5py.File(hfile,'r')
    Eg3D = f.get('/Grid00000001/Grey_Radiation_Energy')
    et3D = f.get('/Grid00000001/Total_Energy')
    nx, ny, nz = Eg3D.shape
    Egval = sum(sum(sum(Eg3D)))/ny/nz/nx*dUnit*lUnit*lUnit/tUnit/tUnit
    etval = sum(sum(sum(et3D)))/ny/nz/nx*lUnit*lUnit/tUnit/tUnit
    result = [Egval, etval, tval*tUnit]
    return result


# initialize the data arrays
etot, times = [], []

# loop over snapshots, loading values and times
for it in range(0,nt+1):
    Egval, etotval, timeval = load_vals(it)
    etot.append(etotval)
    times.append(timeval)

# compare solutions at same times (interpolate if needed)
ntimes = t_ref.size
e_err = []
for it in range(0,nt+1):
    e_val = etot[it]
    t_val = times[it]
    for i in range(0,ntimes-1):
        if ((t_val >= t_ref[i]) and (t_val < t_ref[i+1])):
            t0 = t_ref[i]
            t1 = t_ref[i+1]
            e0 = e_ref[i]
            e1 = e_ref[i+1]
    e_sol = e0 + (e1-e0)*(t_val-t0)/(t1-t0)
    e_err.append((e_val - e_sol)/e_sol)

# compute error norm
err_norm = (np.sum(np.multiply(e_err,e_err))/nt)**(0.5)
if (err_norm < tol):
    print('Error of ',err_norm,' is below tolerance ',tol)
    print('PASS')
else:
    print('Error of ',err_norm,' is above tolerance ',tol)
    print('FAIL')
