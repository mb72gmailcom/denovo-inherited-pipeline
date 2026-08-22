import pytest

from inherited.cli import build_parser
from inherited.constants import (
    DEFAULT_AB,
    DEFAULT_AB_HOM,
    DEFAULT_DP,
    DEFAULT_GQ,
    DEFAULT_HAPLO_AB,
    DEFAULT_HAPLO_DP,
)

_ANALYZE_MIN = [
    "analyze",
    "--vcf",
    "a.vcf",
    "--af-json",
    "af.json",
    "--family-file",
    "fam.tsv",
    "-o",
    "out",
]


def test_analyze_parser_qc_defaults():
    args = build_parser().parse_args(_ANALYZE_MIN)
    assert args.gq_threshold == DEFAULT_GQ
    assert args.dp_threshold == DEFAULT_DP
    assert args.dp_haploid_threshold == DEFAULT_HAPLO_DP
    assert args.ab_threshold == DEFAULT_AB
    assert args.ab_hom_threshold == DEFAULT_AB_HOM
    assert args.ab_haploid_threshold == DEFAULT_HAPLO_AB
    assert args.vcf_dir is None
    assert args.vcf_pattern is None
    assert args.chrom is None


def test_analyze_parser_qc_overrides():
    args = build_parser().parse_args(
        [
            *_ANALYZE_MIN,
            "--gq-threshold",
            "30",
            "--dp-threshold",
            "15",
            "--dp-haploid-threshold",
            "8",
            "--ab-threshold",
            "0.3",
            "--ab-hom-threshold",
            "0.95",
            "--ab-haploid-threshold",
            "0.8",
        ]
    )
    assert args.gq_threshold == 30
    assert args.dp_threshold == 15
    assert args.dp_haploid_threshold == 8
    assert args.ab_threshold == 0.3
    assert args.ab_hom_threshold == 0.95
    assert args.ab_haploid_threshold == 0.8


def test_analyze_parser_vcf_dir_mode():
    args = build_parser().parse_args(
        [
            "analyze",
            "--vcf-dir",
            "shards",
            "--vcf-pattern",
            "SPARK.WGS.2026_08.gatk",
            "--chr",
            "chr2",
            "--af-json",
            "af.json",
            "--family-file",
            "fam.tsv",
            "-o",
            "out",
        ]
    )
    assert args.vcf is None
    assert args.vcf_dir.as_posix() == "shards"
    assert args.vcf_pattern == "SPARK.WGS.2026_08.gatk"
    assert args.chrom == "chr2"


def test_analyze_parser_rejects_vcf_and_vcf_dir_together():
    with pytest.raises(SystemExit):
        build_parser().parse_args(
            [
                *_ANALYZE_MIN,
                "--vcf-dir",
                "shards",
            ]
        )
