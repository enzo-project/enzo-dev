from yt.mods import *
from yt.utilities.answer_testing.api import YTStaticOutputTest

class TestCoolingValues(YTStaticOutputTest):
    name = "cooling_values"

    def run(self):
        # self.pf already exists
        self.result = self.pf.h.grids[0]['Cooling_Time']

    def compare(self, old_result):
        current_buffer = self.result
        old_buffer = old_result

        # We want our arrays to agree to some delta
        self.compare_array_delta(current_buffer, old_buffer, 1e-6)

    def plot(self):
        # There's not much to plot, so we just return an empty list.
        return []
