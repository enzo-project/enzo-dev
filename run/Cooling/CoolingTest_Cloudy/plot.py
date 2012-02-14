from matplotlib import pyplot
from yt.mods import *
import numpy as na

def _H_NumberDensity(field, data):
    "Get total Hydrogen nuclei number density."
    if data.pf['MultiSpecies'] == 0:
        return (0.76 * data['Density'])

    fieldData = na.zeros(data['HI_Density'].shape, 
                         dtype=data['HI_Density'].dtype)

    if data.pf['MultiSpecies'] > 0:
        fieldData += data["HI_Density"]
        fieldData += data["HII_Density"]
    if data.pf["MultiSpecies"] > 1:
        fieldData += data["HM_Density"]
        fieldData += data["H2I_Density"]
        fieldData += data["H2II_Density"]
    if data.pf["MultiSpecies"] > 2:
        fieldData += data["HDI_Density"] / 3.0
    return fieldData
def _ConvertHNumberDensity(data):
    return (1 / 1.67e-24)
add_field("H_NumberDensity", units=r"\rm{cm}^{-3}",
          function=_H_NumberDensity,
          convert_function=_ConvertHNumberDensity)

def plot_cooling_rate(input_file, coordinates, axes, labels=None):
    "Plot cooling rate vs. T for various densities and metallicities."

    pf = load(input_file)
    grid = pf.h.grids[0]

    cooling = grid['Gas_Energy'] * grid['Density'] / grid['Cooling_Time'] / \
        grid['H_NumberDensity']**2

    for q, coord in enumerate(coordinates):
        if labels is None:
            my_coord = list(coord)
            my_coord.append(0)
            my_coord = tuple(my_coord)
            label = "log(n$_{\mathrm{H}}$/cm$^{-3}$) = %.1f, log(Z/Z$_{\odot}$) = %.1f" % \
                (na.log10(grid['H_NumberDensity'][my_coord]),
                 na.log10(grid['Metallicity'][my_coord]))
        else:
            label = labels[q]
        axes.loglog(grid['Temperature'][coord], cooling[coord], label=label)

def plot_cooling_solutions(axes):
    """
    Plot some known cooling rates:
    1. CIE atomic H/He (Black 1981).
    2. Z = 0.5, 1 Z_sun (Sarazin & White 1987).
    """

    black1981 = file("primordial_cie.dat")
    t_hhe = []
    c_hhe = []
    for line in black1981:
        if not line.startswith('#') and len(line) > 1:
            online = line.split()
            t_hhe.append(float(online[0]))
            c_hhe.append(float(online[1]))
    t_hhe = na.power(10, t_hhe)
    c_hhe = na.power(10, c_hhe)

    sz1987 = file("cool_rates.in")
    t_sz = []
    c1_sz = []
    c2_sz = []
    for line in sz1987:
        if not line.startswith('#') and len(line) > 1:
            online = line.split()
            t_sz.append(float(online[0]))
            c1_sz.append(float(online[1]))
            c2_sz.append(float(online[2]))
    t_sz = na.power(10, t_sz)
    c1_sz = na.power(10, c1_sz)
    c2_sz = na.power(10, c2_sz)

    #axes.loglog(t_sz, c2_sz, label='Z = 0.5 Z$_{\odot}$ (Sarazin & White 1987)')
    axes.loglog(t_sz, c1_sz, label='Z = Z$_{\odot}$ (Sarazin & White 1987)')
    axes.loglog(t_hhe, c_hhe, label='H/He (Black 1981)')

pyplot.clf()
axes = pyplot.axes()
axes.set_xlabel('T [K]')
axes.set_ylabel('$\Lambda/n_{H}^{2}$ [erg s$^{-1}$ cm$^{3}$]')
plot_cooling_rate('DD0001/DD0001', [(1, 4)], axes, 
                  labels=['Cloudy, Z = Z$_{\odot}$'])
plot_cooling_solutions(axes)
axes.set_xlim(10, 1e8)
axes.legend(prop=dict(size=10), loc='best')
pyplot.savefig('cooling_rates.png')
