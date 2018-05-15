import yt
import numpy as np
from yt.utilities.answer_testing.framework import assert_equal
import os,glob

test_data_dir   = os.environ.get("COSMO_TEST_DATA_DIR", None)
compare_answers = int(os.environ.get("COSMO_TEST_COMPARE",0))

def test_output_number():
    ds = yt.load('DD0000/DD0000')

    DDnum = len(glob.glob('DD????/DD????'))
    RDnum = len(glob.glob('RD????/RD????'))

    output_data = {'number_of_files' : np.array([DDnum,RDnum])}

    filename = "outputnum_results.h5"
    yt.save_as_dataset(ds, filename, output_data)

    if compare_answers:
        compare_filename = os.path.join(test_data_dir, filename)
        ds_comp = yt.load(compare_filename)

        assert_equal(output_data['num_files'], ds_comp.data['num_files'])

    return

