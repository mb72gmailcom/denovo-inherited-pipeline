from pathlib import Path

import pytest

from inherited.af import is_rare, load_af_json
from inherited.analyze import get_position
from inherited.classify import (
    classify_father_son,
    classify_mother_son,
    classify_trio,
)
from inherited.families import load_family_relations, normalize_sex
from inherited.genotype import get_good_site, is_good, is_hom_ref, is_mendelian_diploid
from inherited.xchrom import chrom_mode_for, is_x_chrom, is_y_chrom, male_x_bucket, x_region


FIXTURES = Path(__file__).parent / "fixtures"


def test_get_position_parses_without_sample_fields():
    assert get_position("chr22\t12345\tid\tA\tG\t.\t.\t.\tGT\t0/1\n") == 12345


def test_get_position_rejects_malformed_record():
    with pytest.raises(ValueError, match="Malformed VCF"):
        get_position("chr22\t12345")


def test_load_af_json_scalar_and_object(tmp_path):
    table = load_af_json(FIXTURES / "tiny_af.json")
    assert table["var_rare"] == 0.001
    assert is_rare(table, "var_common") is False
    assert is_rare(table, "var_rare") is True


def test_load_af_json_object_value(tmp_path):
    path = tmp_path / "af.json"
    path.write_text('{"k1": {"AF": 0.2, "AF_EUR": 0.01}, "k2": {"AF": 0.005}}')
    table = load_af_json(path)
    assert table["k1"] == 0.01
    assert table["k2"] == 0.005


def test_load_family_relations():
    rel = load_family_relations(FIXTURES / "families.tsv")
    assert rel.trio_cl["child1"] == ("ma1", "fa1")
    assert rel.trios_ids == [["child1", "fa1", "ma1"]]
    assert rel.family_size["fam1"] == 3
    assert "child1" in rel.female_children
    assert not rel.male_children


def test_load_family_relations_splits_sex_and_skips_unknown():
    rel = load_family_relations(FIXTURES / "families_x.tsv")
    assert rel.female_children == {"girl1"}
    assert rel.male_children == {"boy1"}
    assert "unknown1" not in rel.trio_cl


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("1", "male"),
        ("Male", "male"),
        ("male", "male"),
        ("2", "female"),
        ("Female", "female"),
        ("female", "female"),
        (".", None),
        ("", None),
        ("other", None),
    ],
)
def test_normalize_sex(value, expected):
    assert normalize_sex(value) is expected


def test_is_good_rejects_low_dp():
    assert not is_good("0/1", "5", "2,3", "0,0,0,0", "30", 1)


def test_is_good_rejects_missing_ad():
    assert not is_good("0/1", "30", ".", "0,0,0,0", "30", 1)
    assert not is_good("0/1", "30", "15,.", "0,0,0,0", "30", 1)


def test_is_good_allows_homozygous_reference():
    assert is_good("0/0", "30", "30,0", "0,0,0,0", "30", 1)


def test_is_good_multiallelic_cleans_missing_ad_as_zero():
    assert not is_good("0/1", "30", "10,4,.,6", "0,0,0,0", "30", 2, clean_missing_ad_as_zero=True)
    assert not is_good("0/1", "30", "10,4,.,6", "0,0,0,0", "30", 1, clean_missing_ad_as_zero=True)
    assert is_good("0/3", "30", "10,4,.,6", "0,0,0,0", "30", 3, clean_missing_ad_as_zero=True)


def test_is_good_multiallelic_strict_without_clean():
    assert not is_good("0/3", "30", "10,4,.,6", "0,0,0,0", "30", 3)


def test_get_good_site_handles_missing_ad():
    sample = "0/1:30:.:0,0,0,0:30:0,30,30:."
    ac, gt, gq = get_good_site(sample, 1)
    assert ac == -1


def test_get_good_site_homozygous_reference_returns_zero():
    sample = "0/0:30:30,0:0,0,0,0:30:0,30,30:."
    assert get_good_site(sample, 1) == (0, "0/0", "30")


def test_get_good_site_skip_qc_if_no_alt_bypasses_depth():
    # Low DP would fail normal QC, but non-carrier children skip QC entirely.
    sample = "0/0:1:1,0:0,0,0,0:1:0,1,1:."
    assert get_good_site(sample, 1) == (-1, ".", "0")
    assert get_good_site(sample, 1, skip_qc_if_no_alt=True) == (0, "0/0", "1")


def test_get_good_site_skip_qc_if_no_alt_still_filters_carriers():
    sample = "0/1:5:2,3:0,0,0,0:30:0,30,30:."
    assert get_good_site(sample, 1, skip_qc_if_no_alt=True) == (-1, ".", "0")


def test_get_good_site_haploid_thresholds():
    low_ab = "1/1:10:5,5:0,0,0,0:30:0,30,30:."
    assert get_good_site(low_ab, 1, haploid=True) == (-1, ".", "0")
    high_ab = "1/1:10:1,9:0,0,0,0:30:0,30,30:."
    assert get_good_site(high_ab, 1, haploid=True) == (2, "1/1", "30")


def test_get_good_site_multiallelic_clean_ad():
    sample = "0/3:30:10,4,.,6:0,0,0,0:30:0,30,30:."
    assert get_good_site(sample, 2, clean_ad=True) == (-1, ".", "0")
    assert get_good_site(sample, 3, clean_ad=True) == (1, "0/3", "30")
    assert get_good_site(sample, 3, clean_ad=False) == (-1, ".", "0")


def test_get_good_site_counts_alt():
    sample = "0/1:30:15,15:0,0,0,0:30:0,30,30:."
    ac, gt, gq = get_good_site(sample, 1)
    assert ac == 1
    assert gt == "0/1"
    assert gq == "30"


def test_get_good_site_returns_negative_one_when_not_good():
    sample = "0/1:5:2,3:0,0,0,0:30:0,30,30:."
    ac, gt, gq = get_good_site(sample, 1)
    assert ac == -1
    assert gt == "."
    assert gq == "0"


@pytest.mark.parametrize(
    ("ac", "mac", "fac", "m_gt", "f_gt", "c_gt", "alt_index", "expected"),
    [
        (1, 2, 2, "1/1", "1/1", "0/1", 1, "mendelian_bad"),
        (2, 2, 2, "1/1", "1/1", "1/1", 1, "inherited"),
        (2, 1, 0, "0/1", "0/0", "1/1", 1, "mendelian_bad"),
        (1, 1, 0, "0/1", "0/0", "0/1", 1, "inherited"),
        (2, 0, 1, "0/0", "0/1", "1/1", 1, "mendelian_bad"),
        (1, 0, 1, "0/0", "0/1", "0/1", 1, "inherited"),
        (1, 0, 0, "0/0", "0/0", "0/1", 1, "denovo"),
        (2, 0, 0, "0/0", "0/0", "1/1", 1, "denovo"),
        (1, -1, 0, "0/1", "0/0", "0/1", 1, None),
        (1, 0, -1, "0/0", "0/1", "0/1", 1, None),
        (0, 0, 0, "0/0", "0/0", "0/0", 1, None),
        # Multiallelic: allele-2 counts look inherited, but child has unexplained allele 1.
        (1, 1, 1, "0/2", "0/2", "1/2", 2, "mendelian_bad"),
        (1, 1, 1, "0/2", "0/2", "0/2", 2, "inherited"),
        (1, 0, 0, "0/2", "0/2", "0/3", 3, None),
        (1, 1, 1, "0/2", "0/1", "1/2", 2, "inherited"),
    ],
)
def test_classify_trio(ac, mac, fac, m_gt, f_gt, c_gt, alt_index, expected):
    assert classify_trio(ac, mac, fac, m_gt, f_gt, c_gt, alt_index) == expected


@pytest.mark.parametrize(
    ("ac", "mac", "expected"),
    [
        (1, 1, "inherited"),
        (1, 0, "denovo"),
        (1, -1, None),
        (0, 0, None),
    ],
)
def test_classify_mother_son(ac, mac, expected):
    assert classify_mother_son(ac, mac) == expected


@pytest.mark.parametrize(
    ("ac", "fac", "expected"),
    [
        (1, 1, "inherited"),
        (1, 0, "denovo"),
        (1, -1, None),
        (0, 0, None),
    ],
)
def test_classify_father_son(ac, fac, expected):
    assert classify_father_son(ac, fac) == expected


@pytest.mark.parametrize(
    ("gt", "expected"),
    [
        ("0/0", True),
        ("0|0", True),
        ("0", True),
        ("0/1", False),
        ("0/2", False),
        ("1/1", False),
        (".", False),
        ("0/.", False),
    ],
)
def test_is_hom_ref(gt, expected):
    assert is_hom_ref(gt) is expected


@pytest.mark.parametrize(
    ("m_gt", "f_gt", "c_gt", "expected"),
    [
        ("0/0", "0/1", "0/1", True),
        ("0/1", "0/0", "0/1", True),
        ("0/1", "0/1", "0/1", True),
        ("0/1", "0/1", "1/1", True),
        ("1/1", "1/1", "1/1", True),
        ("1/1", "1/1", "0/1", False),
        ("0/0", "0/0", "0/1", False),
        ("0/2", "0/2", "0/2", True),
        ("0/2", "0/2", "1/2", False),
        ("0/2", "0/1", "1/2", True),
        ("0/1", ".", "0/1", False),
        ("0", "0/1", "0/1", False),
    ],
)
def test_is_mendelian_diploid(m_gt, f_gt, c_gt, expected):
    assert is_mendelian_diploid(m_gt, f_gt, c_gt) is expected


def test_x_region_and_buckets():
    assert is_x_chrom("X")
    assert is_x_chrom("chrX")
    assert not is_x_chrom("22")
    assert is_y_chrom("Y")
    assert is_y_chrom("chrY")
    assert chrom_mode_for("chrX") == "x"
    assert chrom_mode_for("chrY") == "y"
    assert chrom_mode_for("22") == "autosomal"
    assert x_region(10001) == "par1"
    assert x_region(2781479) == "par1"
    assert x_region(2781480) == "nonPar"
    assert x_region(155701383) == "par2"
    assert male_x_bucket(15000) == "males_par1"
    assert male_x_bucket(5_000_000) == "males_nonPar"
    assert male_x_bucket(155800000) == "males_par2"
