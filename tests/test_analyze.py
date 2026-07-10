import json
from pathlib import Path

import pytest

from inherited.analyze import analyze_vcf
from inherited.checkpoint import Checkpoint, CumulativeStats, load_checkpoint, save_checkpoint
from inherited.output import glob_result_tsvs, read_result_tsv, serialize_payload

FIXTURES = Path(__file__).parent / "fixtures"


def test_analyze_vcf_writes_short_format_single_file(tmp_path):
    stats = analyze_vcf(
        vcf_path=FIXTURES / "tiny.vcf",
        af_json_path=FIXTURES / "tiny_af.json",
        family_file=FIXTURES / "families.tsv",
        output_dir=tmp_path / "out",
        multiallelic=True,
        block_size=1,
        segment_size=0,
    )

    inherited_records = read_result_tsv(tmp_path / "out" / "inherited.tsv", short_format=True)
    bad_records = read_result_tsv(tmp_path / "out" / "mendelian_bad.tsv", short_format=True)

    assert len(inherited_records) == 1
    chrom, pos, ref, alt, patient_ids = inherited_records[0]
    assert chrom == "22" and pos == "3000" and patient_ids == ["child1"]
    assert len(bad_records) == 1
    assert stats.inherited_variants == 1
    assert json.loads((tmp_path / "out" / "inherited_per_variant.json").read_text())["var_inh"] == 1


def test_analyze_vcf_writes_segmented_output(tmp_path):
    stats = analyze_vcf(
        vcf_path=FIXTURES / "tiny.vcf",
        af_json_path=FIXTURES / "tiny_af.json",
        family_file=FIXTURES / "families.tsv",
        output_dir=tmp_path / "out",
        segment_size=1,
        block_size=1,
    )

    inherited_files = glob_result_tsvs(tmp_path / "out", "inherited")
    bad_files = glob_result_tsvs(tmp_path / "out", "mendelian_bad")
    assert len(inherited_files) >= 1
    assert len(bad_files) >= 1
    assert inherited_files[0].name == "inherited_00000.tsv"
    assert (tmp_path / "out" / "checkpoint.json").is_file()
    assert (tmp_path / "out" / "stats_cumulative.json").is_file()
    assert stats.inherited_variants == 1


def test_analyze_vcf_skips_variants_in_repeat_intervals(tmp_path):
    stats = analyze_vcf(
        vcf_path=FIXTURES / "tiny.vcf",
        af_json_path=FIXTURES / "tiny_af.json",
        family_file=FIXTURES / "families.tsv",
        output_dir=tmp_path / "out",
        segment_size=0,
        repeats_path=FIXTURES / "tiny_repeats.bed",
    )

    inherited_records = read_result_tsv(tmp_path / "out" / "inherited.tsv", short_format=True)
    bad_records = read_result_tsv(tmp_path / "out" / "mendelian_bad.tsv", short_format=True)

    assert len(inherited_records) == 0
    assert len(bad_records) == 1
    assert bad_records[0][1] == "4000"
    assert stats.inherited_variants == 0
    assert stats.mendelian_bad_variants == 1


def test_resume_continues_from_checkpoint(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    save_checkpoint(
        out,
        Checkpoint(
            chrom="22",
            last_pos=3000,
            segment_index=0,
            cumulative=CumulativeStats(
                variants_seen=3,
                alleles_tested=1,
                inherited_entries=1,
                inherited_variants=1,
                inherited_per_person={"child1": 1},
            ),
            completed=False,
        ),
    )

    stats = analyze_vcf(
        vcf_path=FIXTURES / "tiny.vcf",
        af_json_path=FIXTURES / "tiny_af.json",
        family_file=FIXTURES / "families.tsv",
        output_dir=out,
        segment_size=1000,
        block_size=1,
        resume=True,
    )

    assert stats.inherited_variants == 1
    assert stats.mendelian_bad_variants == 1
    assert stats.variants_seen == 4
    assert len(read_result_tsv(glob_result_tsvs(out, "mendelian_bad")[0])) == 1


def test_resume_rejects_completed_checkpoint(tmp_path):
    out = tmp_path / "out"
    analyze_vcf(
        vcf_path=FIXTURES / "tiny.vcf",
        af_json_path=FIXTURES / "tiny_af.json",
        family_file=FIXTURES / "families.tsv",
        output_dir=out,
        segment_size=1000,
    )
    assert load_checkpoint(out).completed is True

    with pytest.raises(ValueError, match="completed"):
        analyze_vcf(
            vcf_path=FIXTURES / "tiny.vcf",
            af_json_path=FIXTURES / "tiny_af.json",
            family_file=FIXTURES / "families.tsv",
            output_dir=out,
            segment_size=1000,
            resume=True,
        )


def test_serialize_payload_short_and_full():
    hits = {"child1": ("0/0", "0/1", "0/1", "30")}
    assert serialize_payload(hits, short_format=True) == "child1"
    assert serialize_payload(hits, short_format=False) == "child1=0/0|0/1|0/1|30"
