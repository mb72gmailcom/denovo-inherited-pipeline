from pathlib import Path

import pytest

from inherited.shards import chrom_filename_tokens, discover_vcf_shards


def _touch(path: Path) -> Path:
    path.write_text("")
    return path


def test_chrom_filename_tokens_accepts_chr_prefix():
    assert chrom_filename_tokens("chr2") == ("chr2", "2")
    assert chrom_filename_tokens("2") == ("2", "chr2")
    assert chrom_filename_tokens("X") == ("X", "chrX")


def test_discover_vcf_shards_sorts_and_parses_coordinates(tmp_path):
    _touch(tmp_path / "SPARK.WGS.2026_08.gatk.chr2_95000001_97500000.vcf.gz")
    _touch(tmp_path / "SPARK.WGS.2026_08.gatk.chr2_92500001_95000000.vcf.gz")
    _touch(tmp_path / "SPARK.WGS.2026_08.gatk.chr21_1_1000.vcf.gz")
    _touch(tmp_path / "other.chr2_1_1000.vcf.gz")

    shards = discover_vcf_shards(tmp_path, "SPARK.WGS.2026_08.gatk", "chr2")
    assert [shard.start for shard in shards] == [92500001, 95000001]
    assert shards[0].end == 95000000
    assert shards[1].chrom == "chr2"
    assert shards[1].path.name.endswith("97500000.vcf.gz")


def test_discover_vcf_shards_does_not_match_chr21_when_asking_chr2(tmp_path):
    _touch(tmp_path / "callset.chr2_1_1000.vcf.gz")
    _touch(tmp_path / "callset.chr21_1_1000.vcf.gz")
    _touch(tmp_path / "callset.chr22_1_1000.vcf.gz")

    shards = discover_vcf_shards(tmp_path, "callset", "2")
    assert len(shards) == 1
    assert shards[0].chrom == "chr2"


def test_discover_vcf_shards_rejects_overlaps(tmp_path):
    _touch(tmp_path / "callset.chr2_1_1000.vcf")
    _touch(tmp_path / "callset.chr2_1000_2000.vcf")

    with pytest.raises(ValueError, match="Overlapping"):
        discover_vcf_shards(tmp_path, "callset", "chr2")


def test_discover_vcf_shards_requires_matches(tmp_path):
    with pytest.raises(FileNotFoundError, match="No VCF shards"):
        discover_vcf_shards(tmp_path, "callset", "chr2")
