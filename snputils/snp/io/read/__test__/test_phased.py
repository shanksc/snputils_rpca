import numpy as np


def test_vcf_summed_strands(snpobj_vcf, snpobj_vcf_summed_strands):
    assert np.array_equal(snpobj_vcf.calldata_gt.sum(axis=2, dtype=np.int8), snpobj_vcf_summed_strands.calldata_gt)


def test_bed_summed_strands(snpobj_bed, snpojb_bed_summed_strands):
    assert np.array_equal(snpobj_bed.calldata_gt.sum(axis=2, dtype=np.int8), snpojb_bed_summed_strands.calldata_gt)


def test_pgen_summed_strands(snpobj_pgen, snpobj_pgen_summed_strands):
    assert np.array_equal(snpobj_pgen.calldata_gt.sum(axis=2, dtype=np.int8), snpobj_pgen_summed_strands.calldata_gt)
