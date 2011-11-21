from yt.mods import *
from yt.utilities.answer_testing.api import YTStaticOutputTest, create_test
import pylab
class TestGardinerImage(YTStaticOutputTest):
    field = None

    def run(self):
        # self.pf already exists
        sl = self.pf.h.slice(2, 0.5)
        frb = FixedResolutionBuffer(sl, (0.0, 1.0, 0.0, 1.0), (200,200))
        self.result = frb[self.field]

    def compare(self, old_result):
        current_buffer = self.result
        old_buffer = old_result

        # We want our arrays to agree to some delta
        self.compare_array_delta(current_buffer, old_buffer, 3e-5)

    def plot(self):
        pylab.clf()
        pylab.imshow(self.result,
            interpolation='nearest', origin='lower')
        fn = "%s_%s_projection.png" % (self.pf, self.field)
        pylab.savefig(fn)
        return [fn]
for field in ['Density','Pressure','Bx','By']:
    create_test(TestGardinerImage,'test_gardiner_%s'%field,field=field)
