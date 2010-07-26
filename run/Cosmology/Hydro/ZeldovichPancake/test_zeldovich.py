from yt.config import ytcfg
ytcfg["yt","loglevel"] = '50'
ytcfg["yt","suppressStreamLogging"] = 'True'

from yt.mods import *
from yt.extensions.enzo_test import YTStaticOutputTest, run

class TestZeldovich(YTStaticOutputTest):
    name = "zeldovich_plot"

    def run(self):
        # self.pf already exists
        ray = self.pf.h.ray([0.0,0.5,0.5], [1.0,0.5,0.5])
        self.result = ray["Density"]

    def compare(self, old_result):
        current_buffer = self.result
        old_buffer = old_result

        # We want our arrays to agree to some delta
        self.compare_array_delta(current_buffer, old_buffer, 5e-3)

    def plot(self):
        # There's not much to plot, so we just return an empty list.
        return []

class TestZeldovichMax(YTStaticOutputTest):
    name = "zeldovich_max"

    def run(self):
        # self.pf already exists
        self.result = self.pf.h.find_max("Density")[0]

    def compare(self, old_result):
        self.compare_value_delta(self.result, old_result, 5e-3)

    def plot(self):
        # There's not much to plot, so we just return an empty list.
        return []

if __name__ == "__main__":
    run()
