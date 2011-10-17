from yt.config import ytcfg
ytcfg["yt","loglevel"] = '50'
ytcfg["yt","suppressStreamLogging"] = 'True'

from yt.mods import *
from yt.utilities.answer_testing.api import \
    YTStaticOutputTest, create_test, run_main
import matplotlib.pyplot as plt

class TestAMRPhotonTest(YTStaticOutputTest):

    def run(self):
        # self.pf already exists
        sl = self.pf.h.slice(2,1.0/64)
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

for f in ['Density', 'HI_Fraction', 'HII_Fraction', 'HI_kph']:
    create_test(TestAMRPhotonTest, 'amr_photon_test_%s' % f, field=f)
