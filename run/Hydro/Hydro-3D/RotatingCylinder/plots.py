from yt.mods import *
import yt.extensions.EnzoSimulation as ES

# Create image slices and plot the evolution of the total angular momentum of the system. 

es = ES.EnzoSimulation("RotatingCylinder.enzo", initial_time=0.0, final_time=26)

AngMom = []
time = []
i = 0

for output in es.allOutputs:
    # load up a dataset
    pf = load(output['filename'])

    # calculate total angular momentum of simulation volume
    data = pf.h.all_data()
    AngMom.append(data.quantities["TotalQuantity"]("AngularMomentum")[0])    

    # calculate the time of output
    time.append(2.0*i)
    i = i + 1

    # create density slices
    pc = PlotCollection(pf, center=[0.5,0.5,0.5])
    pc.add_slice("Density", 0)
    pc.save("t%s" % str(time[i-1]))


# plot the percentage change in angular momentum
AngMomInitial = AngMom[0]
AngMomPercentageChange = []
AngMomPercentageChangeEachMyr = []

for i, item in enumerate(AngMom):
    AngMomPercentageChange.append(100.0*(item - AngMomInitial)/AngMomInitial)
    AngMomPercentageChangeEachMyr.append(100.0*(item - AngMom[i-1])/AngMom[i-1]/2.0)

pylab.xlabel("Time [Myr]")
pylab.ylabel("Percentage Increase in Total Angular Momentum")
pylab.plot(time, AngMomPercentageChange, 'b', label='Net change')
pylab.plot(time, AngMomPercentageChangeEachMyr, 'r', label='Change per Myr')
pylab.legend()
pylab.savefig('AngMom.png')
pylab.show()

print "Net L change, Maximum L change: ", max(AngMomPercentageChange), max(AngMomPercentageChangeEachMyr)


