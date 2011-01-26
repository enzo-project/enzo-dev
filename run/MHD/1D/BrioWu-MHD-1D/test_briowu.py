from yt.mods import *
from yt.extensions.enzo_test import YTStaticOutputTest

class TestBrioWuDensity(YTStaticOutputTest):
    name = "briowu_density"

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

class TestBrioWuBy(YTStaticOutputTest):
    name = "briowu_By"

    def run(self):
        # self.pf already exists
        ray = self.pf.h.ray([0.0,0.5,0.5], [1.0,0.5,0.5])
        self.result = ray["By"]

    def compare(self, old_result):
        current_buffer = self.result
        old_buffer = old_result

        # We want our arrays to agree to some delta
        self.compare_array_delta(current_buffer, old_buffer, 5e-3)

    def plot(self):
        # There's not much to plot, so we just return an empty list.
        return []

