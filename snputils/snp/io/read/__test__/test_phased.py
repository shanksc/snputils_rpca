import numpy as np


def test_vcf_unphased(snpobj_vcf, snpobj_vcf_unphased):
    assert np.array_equal(snpobj_vcf.calldata_gt.sum(axis=2, dtype=np.int8), snpobj_vcf_unphased.calldata_gt)


def test_bed_unphased(snpobj_bed, snpojb_bed_unphased):
    assert np.array_equal(snpobj_bed.calldata_gt.sum(axis=2, dtype=np.int8), snpojb_bed_unphased.calldata_gt)


def test_pgen_unphased(snpobj_pgen, snpobj_pgen_unphased):
    assert np.array_equal(snpobj_pgen.calldata_gt.sum(axis=2, dtype=np.int8), snpobj_pgen_unphased.calldata_gt)
