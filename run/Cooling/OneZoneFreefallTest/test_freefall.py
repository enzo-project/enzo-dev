from yt.mods import *
from yt.utilities.answer_testing.api import YTStaticOutputTest

class TestFreefallTemperatures(YTStaticOutputTest):
    name = "freefall_temperatures"

    def run(self):
        # self.pf already exists
        self.result = self.pf.h.all_data()['Temperature']

    def compare(self, old_result):
        current_buffer = self.result
        old_buffer = old_result

        # We want our arrays to agree to some delta
        self.compare_array_delta(current_buffer, old_buffer, 1e-6)

    def plot(self):
        # There's not much to plot, so we just return an empty list.
        return []

class TestFreefallDustTemperatures(YTStaticOutputTest):
    name = "freefall_dust_temperatures"

    def run(self):
        # self.pf already exists
        self.result = self.pf.h.all_data()['Dust_Temperature']

    def compare(self, old_result):
        current_buffer = self.result
        old_buffer = old_result

        # We want our arrays to agree to some delta
        self.compare_array_delta(current_buffer, old_buffer, 1e-6)

    def plot(self):
        # There's not much to plot, so we just return an empty list.
        return []
