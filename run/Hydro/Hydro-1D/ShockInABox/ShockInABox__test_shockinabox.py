import os
import yt
from yt.funcs import *
from yt.testing import *
from yt.frontends.enzo.answer_testing_support import \
    requires_outputlog

# Verifies that OutputLog exists
@requires_outputlog(os.path.dirname(__file__), "ShockInABox.enzo")
def test_shockinabox():
    pf = yt.load('DD0010/data0010')
    mach = pf.find_max(('enzo', 'Mach'))[0]
    myname = 'ShockInABox_Mach'
    yield assert_allclose, mach, 3.0, 1.0e-2
