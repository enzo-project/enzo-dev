from yt.mods import *
from yt.testing import *
from yt.utilities.answer_testing.framework import \
     AnswerTestingTest, \
     sim_dir_load
from yt.frontends.enzo.answer_testing_support import \
     requires_outputlog

tolerance = ytcfg.getint("yt", "answer_testing_tolerance")
     
_pf_name = os.path.basename(os.path.dirname(__file__)) + ".enzo"
_dir_name = os.path.dirname(__file__)

class TestGravityTest(AnswerTestingTest):
    _type_name = "GravityTest"
    _attrs = ()

    def __init__(self):
        self.pf = None
    
    def run(self):
        Data = np.loadtxt("TestGravityCheckResults.out")
        radius = Data[:,0]
        ForceRadialComputed = Data[:,2]
        ForceRadialTrue = Data[:,3]

        Error = (ForceRadialComputed-ForceRadialTrue)/ForceRadialTrue
        indices = np.where((radius > 1.0) & (radius < 8.0))

        rmsError = np.std(Error[indices])
        return rmsError

    def compare(self, new_result, old_result):
        assert_allclose(new_result, old_result, rtol=10.**-tolerance, atol=0)

@requires_outputlog(_dir_name, _pf_name)
def test_gravity_test():
    yield TestGravityTest()
