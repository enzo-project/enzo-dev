from yt.config import ytcfg
ytcfg["yt","loglevel"] = '50'
ytcfg["yt","suppressStreamLogging"] = 'True'

from yt.mods import *
from yt.extensions.enzo_test import YTStaticOutputTest

class TestGravityTest(YTStaticOutputTest):
    name = "GravityTest_plot"

    def run(self):
        Data = na.loadtxt("TestGravityCheckResults.out")
        radius = Data[:,0]
        ForceRadialComputed = Data[:,2]
        ForceRadialTrue = Data[:,3]

        Error = (ForceRadialComputed-ForceRadialTrue)/ForceRadialTrue
        indices = na.where((radius > 1.0) & (radius < 8.0))

        rmsError = na.std(Error[indices])
        print "rms error = "+str(rmsError)
        self.result = rmsError

    def compare(self, old_result):
        self.compare_value_delta(self.result, old_result, 0.01)

    def plot(self):
        # There's not much to plot, so we just return an empty list.
        return []

if __name__ == "__main__":
    run()
