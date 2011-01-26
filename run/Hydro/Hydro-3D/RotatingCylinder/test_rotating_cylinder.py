from yt.mods import *
from yt.extensions.enzo_test import YTStaticOutputTest

class TestCyclinderImage(YTStaticOutputTest):
    name = "cylinder_image"

    def run(self):
        # self.pf already exists
        sl = self.pf.h.proj(2, 'Density', center=[0.5,0.5,0.5])
        frb = FixedResolutionBuffer(sl, (0.0, 1.0, 0.0, 1.0), (64, 64))
        self.result = frb["Density"]

    def compare(self, old_result):
        current_buffer = self.result
        old_buffer = old_result

        # We want our arrays to agree to some delta
        self.compare_array_delta(current_buffer, old_buffer, 5e-3)

    def plot(self):
        # There's not much to plot, so we just return an empty list.
        return []

class TestLVariation(YTStaticOutputTest):
    name = "deltaL"

    def run(self):
        # self.pf already exists
        data = self.pf.h.all_data()
        AngMom = []
        AngMom.append(data.quantities["TotalQuantity"]("AngularMomentum")[0])
        AngMomInitial = AngMom[0]
        AngMomPercentageChange = []
        for i, item in enumerate(AngMom):
            AngMomPercentageChange.append(100.0*(item - AngMomInitial)/AngMomInitial)
        self.result = max(AngMomPercentageChange)

    def compare(self, old_result):
        self.compare_value_delta(self.result, old_result, 5e-3)

    def plot(self):
        # There's not much to plot, so we just return an empty list.
        return []

if __name__ == "__main__":
    run()
