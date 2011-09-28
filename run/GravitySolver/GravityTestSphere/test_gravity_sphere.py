from yt.config import ytcfg
ytcfg["yt","loglevel"] = '50'
ytcfg["yt","suppressStreamLogging"] = 'True'

from yt.mods import *
from yt.utilities.answer_testing.api import YTStaticOutputTest

class TestGravitySphere(YTStaticOutputTest):
    name = "TestGravitySphere_plot"

    def run(self):
        # self.pf already exists
        ray = self.pf.h.ray([0.5,0.5,0.5], [1.0,1.0,1.0])
        self.result = ray["x-velocity"]

    def compare(self, old_result):
        current_buffer = self.result
        old_buffer = old_result

        # We want our arrays to agree to some delta
        self.compare_array_delta(current_buffer, old_buffer, 5e-3)

    def plot(self):
        # There's not much to plot, so we just return an empty list.
        return []

if __name__ == "__main__":
    run()
