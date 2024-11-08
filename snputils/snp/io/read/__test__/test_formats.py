import numpy as np

# # VCF - BED


def test_vcf_bed_samples(snpobj_vcf, snpobj_bed):
    assert snpobj_vcf.samples is not None
    assert snpobj_bed.samples is not None
    assert np.array_equal(snpobj_vcf.samples, snpobj_bed.samples)


def test_vcf_bed_variants_ref(snpobj_vcf, snpobj_bed):
    assert snpobj_vcf.variants_ref is not None
    assert snpobj_bed.variants_ref is not None
    assert np.array_equal(snpobj_vcf.variants_ref, snpobj_bed.variants_ref)


# TODO
# def test_vcf_bed_variants_alt(snpobj_vcf, snpobj_bed):
#     assert snpobj_vcf.variants_alt is not None
#     assert snpobj_bed.variants_alt is not None
#     assert np.array_equal(snpobj_vcf.variants_alt, snpobj_bed.variants_alt)


def test_vcf_bed_variants_chrom(snpobj_vcf, snpobj_bed):
    assert snpobj_vcf.variants_chrom is not None
    assert snpobj_bed.variants_chrom is not None
    assert np.array_equal(snpobj_vcf.variants_chrom, snpobj_bed.variants_chrom)


def test_vcf_bed_variants_id(snpobj_vcf, snpobj_bed):
    assert snpobj_vcf.variants_id is not None
    assert snpobj_bed.variants_id is not None
    assert np.array_equal(snpobj_vcf.variants_id, snpobj_bed.variants_id)


def test_vcf_bed_variants_pos(snpobj_vcf, snpobj_bed):
    assert snpobj_vcf.variants_pos is not None
    assert snpobj_bed.variants_pos is not None
    assert np.array_equal(snpobj_vcf.variants_pos, snpobj_bed.variants_pos)


# # BED - PGEN


def test_bed_pgen_samples(snpobj_bed, snpobj_pgen):
    assert snpobj_bed.samples is not None
    assert snpobj_pgen.samples is not None
    assert np.array_equal(snpobj_bed.samples, snpobj_pgen.samples)


def test_bed_pgen_variants_ref(snpobj_bed, snpobj_pgen):
    assert snpobj_bed.variants_ref is not None
    assert snpobj_pgen.variants_ref is not None
    assert np.array_equal(snpobj_bed.variants_ref, snpobj_pgen.variants_ref)


# TODO
# def test_bed_pgen_variants_alt(snpobj_bed, snpobj_pgen):
#     assert snpobj_bed.variants_alt is not None
#     assert snpobj_pgen.variants_alt is not None
#     assert np.array_equal(snpobj_bed.variants_alt, snpobj_pgen.variants_alt)


def test_bed_pgen_variants_chrom(snpobj_bed, snpobj_pgen):
    assert snpobj_bed.variants_chrom is not None
    assert snpobj_pgen.variants_chrom is not None
    assert np.array_equal(snpobj_bed.variants_chrom, snpobj_pgen.variants_chrom)


def test_bed_pgen_variants_id(snpobj_bed, snpobj_pgen):
    assert snpobj_bed.variants_id is not None
    assert snpobj_pgen.variants_id is not None
    assert np.array_equal(snpobj_bed.variants_id, snpobj_pgen.variants_id)


def test_bed_pgen_variants_pos(snpobj_bed, snpobj_pgen):
    assert snpobj_bed.variants_pos is not None
    assert snpobj_pgen.variants_pos is not None
    assert np.array_equal(snpobj_bed.variants_pos, snpobj_pgen.variants_pos)


# VCF - PGEN


def test_vcf_pgen_calldata_gt(snpobj_vcf, snpobj_pgen):
    assert np.array_equal(snpobj_vcf.calldata_gt, snpobj_pgen.calldata_gt)


def test_vcf_pgen_samples(snpobj_vcf, snpobj_pgen):
    assert snpobj_vcf.samples is not None
    assert snpobj_pgen.samples is not None
    assert np.array_equal(snpobj_vcf.samples, snpobj_pgen.samples)


def test_vcf_pgen_variants_ref(snpobj_vcf, snpobj_pgen):
    assert snpobj_vcf.variants_ref is not None
    assert snpobj_pgen.variants_ref is not None
    assert np.array_equal(snpobj_vcf.variants_ref, snpobj_pgen.variants_ref)


# TODO
# def test_vcf_pgen_variants_alt(snpobj_vcf, snpobj_pgen):
#     assert snpobj_vcf.variants_alt is not None
#     assert snpobj_pgen.variants_alt is not None
#     assert np.array_equal(snpobj_vcf.variants_alt, snpobj_pgen.variants_alt)


def test_vcf_pgen_variants_chrom(snpobj_vcf, snpobj_pgen):
    assert snpobj_vcf.variants_chrom is not None
    assert snpobj_pgen.variants_chrom is not None
    assert np.array_equal(snpobj_vcf.variants_chrom, snpobj_pgen.variants_chrom)


def test_vcf_pgen_variants_id(snpobj_vcf, snpobj_pgen):
    assert snpobj_vcf.variants_id is not None
    assert snpobj_pgen.variants_id is not None
    assert np.array_equal(snpobj_vcf.variants_id, snpobj_pgen.variants_id)


def test_vcf_pgen_variants_pos(snpobj_vcf, snpobj_pgen):
    assert snpobj_vcf.variants_pos is not None
    assert snpobj_pgen.variants_pos is not None
    assert np.array_equal(snpobj_vcf.variants_pos, snpobj_pgen.variants_pos)
