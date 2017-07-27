from matplotlib import pyplot

from yt.mods import *

# Create image slices and plot the evolution of the total angular momentum of the system. 

es = simulation("RotatingCylinder.enzo", "Enzo")
es.get_time_series(initial_time=0.0, final_time=26)

AngMom = []
time = []

for pf in es:
    # calculate total angular momentum of simulation volume
    data = pf.h.all_data()
    my_ang = np.array(data.quantities["TotalQuantity"](["AngularMomentumX",
                                                        "AngularMomentumY",
                                                        "AngularMomentumZ"]))
    AngMom.append(np.sqrt((my_ang**2).sum()))

    # get current time
    time.append(pf.current_time * pf.time_units['Myr'])

    # create density slices
    plot = SlicePlot(pf, 'x', 'Density')
    plot.annotate_text((-0.1, -0.05), "t = %f" % pf.current_time)
    plot.save()

# plot the percentage change in angular momentum
AngMom = np.array(AngMom)
AngMomInitial = AngMom[0]
AngMomPercentageChange = 100. * (AngMom - AngMom[0]) / AngMom[0]
AngMomPercentageChangeEachMyr = 100. * (AngMom[1:] / AngMom[:-1] - 1) / 2.
AngMomPercentageChangeEachMyr = np.concatenate([[0.0], AngMomPercentageChangeEachMyr])
    
pyplot.xlabel("Time [Myr]")
pyplot.ylabel("Percentage Increase in Total Angular Momentum")
pyplot.plot(time, AngMomPercentageChange, 'b', label='Net change')
pyplot.plot(time, AngMomPercentageChangeEachMyr, 'r', label='Change per Myr')
pyplot.legend()
pyplot.savefig('AngMom.png')

print("Net L change, Maximum L change: ", max(AngMomPercentageChange), \
  max(AngMomPercentageChangeEachMyr))


