from yt.mods import *
from yt.utilities.answer_testing.api import YTStaticOutputTest

class TestCollapseProjection(YTStaticOutputTest):
    name = "collapse_projection"

    def run(self):
        # self.pf already exists
        proj = self.pf.h.proj(0, 'Density', center=[0.5, 0.5, 0.5])
        frb = FixedResolutionBuffer(proj, (0.0, 1.0, 0.0, 1.0), (150, 150))
        self.result = frb["Density"]

    def compare(self, old_result):
        current_buffer = self.result
        old_buffer = old_result

        # We want our arrays to agree to some delta
        self.compare_array_delta(current_buffer, old_buffer, 5e-3)

    def plot(self):
        # There's not much to plot, so we just return an empty list.
        return []

class TestCollapseMaxValue(YTStaticOutputTest):
    name = "collapse_max"

    def run(self):
        # self.pf already exists
        self.result = self.pf.h.find_max("Density")[0]

    def compare(self, old_result):
        self.compare_value_delta(self.result, old_result, 5e-3)

    def plot(self):
        # There's not much to plot, so we just return an empty list.
        return []
