import numpy as np

from snputils.snp.io.read import SNPReader


def test_auto_reader(data_path, snpobj_pgen):
    reader = SNPReader(data_path + "/pgen/subset.pgen")
    snpobj = reader.read(phased=True)

    assert np.array_equal(snpobj.calldata_gt, snpobj_pgen.calldata_gt)
