from yt.mods import *
from yt.utilities.answer_testing.api import YTStaticOutputTest

class TestShockImage(YTStaticOutputTest):
    name = "shock_image"

    def run(self):
        # self.pf already exists
        sl = self.pf.h.slice(2, 0.5)
        frb = FixedResolutionBuffer(sl, (0.0, 1.0, 0.0, 1.0), (150, 150))
        self.result = frb["Density"]

    def compare(self, old_result):
        current_buffer = self.result
        old_buffer = old_result

        # We want our arrays to agree to some delta
        self.compare_array_delta(current_buffer, old_buffer, 5e-3)

    def plot(self):
        # There's not much to plot, so we just return an empty list.
        return []

class TestMaxValue(YTStaticOutputTest):
    name = "shock_max"

    def run(self):
        # self.pf already exists
        self.result = self.pf.h.find_max("Density")[0]

    def compare(self, old_result):
        self.compare_value_delta(self.result, old_result, 5e-3)

    def plot(self):
        # There's not much to plot, so we just return an empty list.
        return []
