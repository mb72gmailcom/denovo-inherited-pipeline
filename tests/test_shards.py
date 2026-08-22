from pathlib import Path

import pytest

from inherited.shards import discover_vcf_shards


def _touch(path: Path) -> Path:
    path.write_text("")
    return path


def test_discover_vcf_shards_accepts_pattern_that_includes_chrom(tmp_path):
    _touch(tmp_path / "SPARK.WGS.2026_08.gatk.chr21_37500001_40000000.vcf.gz")
    _touch(tmp_path / "SPARK.WGS.2026_08.gatk.chr21_40000001_42500000.vcf.gz")

    for pattern in (
        "SPARK.WGS.2026_08.gatk",
        "SPARK.WGS.2026_08.gatk.chr21",
        "SPARK.WGS.2026_08.gatk.chr21_",
    ):
        shards = discover_vcf_shards(tmp_path, pattern)
        assert [shard.start for shard in shards] == [37500001, 40000001]
        assert shards[0].chrom == "chr21"


def test_discover_vcf_shards_sorts_and_parses_coordinates(tmp_path):
    _touch(tmp_path / "SPARK.WGS.2026_08.gatk.chr2_95000001_97500000.vcf.gz")
    _touch(tmp_path / "SPARK.WGS.2026_08.gatk.chr2_92500001_95000000.vcf.gz")
    _touch(tmp_path / "other.chr2_1_1000.vcf.gz")

    shards = discover_vcf_shards(tmp_path, "SPARK.WGS.2026_08.gatk")
    assert [shard.start for shard in shards] == [92500001, 95000001]
    assert shards[0].end == 95000000
    assert shards[0].chrom == "chr2"
    assert shards[1].path.name.endswith("97500000.vcf.gz")


def test_discover_vcf_shards_rejects_mixed_chromosomes(tmp_path):
    _touch(tmp_path / "callset.chr2_1_1000.vcf.gz")
    _touch(tmp_path / "callset.chr21_1_1000.vcf.gz")

    with pytest.raises(ValueError, match="more than one contig"):
        discover_vcf_shards(tmp_path, "callset")


def test_discover_vcf_shards_rejects_overlaps(tmp_path):
    _touch(tmp_path / "callset.chr2_1_1000.vcf")
    _touch(tmp_path / "callset.chr2_1000_2000.vcf")

    with pytest.raises(ValueError, match="Overlapping"):
        discover_vcf_shards(tmp_path, "callset")


def test_discover_vcf_shards_requires_matches(tmp_path):
    with pytest.raises(FileNotFoundError, match="No VCF shards"):
        discover_vcf_shards(tmp_path, "callset")
