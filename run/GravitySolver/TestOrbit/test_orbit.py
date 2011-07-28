from yt.mods import *
from yt.utilities.answer_testing.api import YTStaticOutputTest

class TestOrbitPosition(YTStaticOutputTest):
    name = "test_orbit_position"

    def run(self):
        # self.pf already exists
        dd = self.pf.h.all_data()
        particle_position = na.array([dd['particle_position_x'][1], dd['particle_position_y'][1]])
        print particle_position
        self.result = particle_position

    def compare(self, old_result):
        current_buffer = self.result
        old_buffer = old_result

        # We want our arrays to agree to some delta
        self.compare_array_delta(current_buffer, old_buffer, 5e-3)

    def plot(self):
        # There's not much to plot, so we just return an empty list.
        return []

