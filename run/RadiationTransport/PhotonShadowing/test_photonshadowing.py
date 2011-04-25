from yt.config import ytcfg
ytcfg["yt","loglevel"] = '50'
ytcfg["yt","suppressStreamLogging"] = 'True'

from yt.mods import *
from yt.utilities.answer_testing.api import YTStaticOutputTest, run_main, create_test
import matplotlib.pyplot as plt

class TestPhotonShadowing1(YTStaticOutputTest):

    def run(self):
        # self.pf already exists
        sl = self.pf.h.slice(2,0.5)
        frb = FixedResolutionBuffer(sl, (0,1,0,1), (200,200))
        self.result = frb["HII_Density"]

    def compare(self, old_result):
        current_buffer = self.result
        old_buffer = old_result

        # We want our arrays to agree to some delta
        self.compare_array_delta(current_buffer, old_buffer, 5e-3)

    def plot(self):
        plt.clf()
        plt.imshow(self.result, interpolation='nearest',
                   origin='lower')
        fn = '%s_%s.png' % (self.pf, self.field)
        plt.savefig(fn)
        return [fn]

class TestPhotonShadowing2(YTStaticOutputTest):
    name = "photon_shadow_temp_plot"

    def run(self):
        # self.pf already exists
        sl = self.pf.h.slice(2,0.5)
        frb = FixedResolutionBuffer(sl, (0,1,0,1), (200,200))
        self.result = frb["Temperature"]

    def compare(self, old_result):
        current_buffer = self.result
        old_buffer = old_result

        # We want our arrays to agree to some delta
        self.compare_array_delta(current_buffer, old_buffer, 5e-3)

    def plot(self):
        return []

for f in ['HI_Density', 'HII_Density', 'Temperature', 'HI_kph']:
    create_test(TestPhotonShadowing1, 'photon_shadow_test_%s' % f, field=f)
