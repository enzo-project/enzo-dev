import glob

from yt.mods import *
from yt.testing import *
from yt.utilities.answer_testing.framework import \
    AnswerTestingTest, \
    sim_dir_load
from yt.frontends.enzo.answer_testing_support import \
    requires_outputlog

_pf_name = os.path.basename(os.path.dirname(__file__)) + ".enzo"
_dir_name = os.path.dirname(__file__)


class TestAPMGravityTest(AnswerTestingTest):
    _type_name = "APMGravityTest"
    _attrs = ()

    def __init__(self):
        self.ds = None

    def run(self):

        rmsError_list = []
        for f in glob.glob('TestGravityCheckResults*'):
            Data = np.loadtxt(f)

            radius = Data[:, 4]
            ForceRadialComputed = Data[:, 6]
            ForceRadialTrue = Data[:, 7]

            Error = (ForceRadialComputed-ForceRadialTrue)/ForceRadialTrue
            indices = np.where(radius > 0.4)
            rmsError_list.append(np.std(Error[indices]))

        return sum(sigma**2 for sigma in rmsError_list)**0.5


@requires_outputlog(_dir_name, _pf_name)
def test_gravity_test():
    yield TestAPMGravityTest()
