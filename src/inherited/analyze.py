from __future__ import annotations

import gzip
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from inherited import __version__
from inherited.af import is_rare, load_af_json
from inherited.checkpoint import load_checkpoint
from inherited.constants import (
    DEFAULT_AB,
    DEFAULT_AF_THRESHOLD,
    DEFAULT_BLOCK_SIZE,
    DEFAULT_DP,
    DEFAULT_GQ,
    DEFAULT_HAPLO_AB,
    DEFAULT_HAPLO_DP,
    DEFAULT_MEMORY_BLOCK,
    DEFAULT_SEGMENT_SIZE,
)
from inherited.debug import log_memory_if_due
from inherited.families import build_trio_indices, load_family_relations
from inherited.genotype import get_good_site
from inherited.output import ResultWriter
from inherited.repeats import RepeatIntervalFilter


def get_nfields(line: str, n: int) -> list[str]:
    return line.rstrip().split("\t")[:n]


@dataclass
class AnalysisStats:
    variants_seen: int = 0
    alleles_tested: int = 0
    inherited_entries: int = 0
    mendelian_bad_entries: int = 0
    inherited_variants: int = 0
    mendelian_bad_variants: int = 0


def analyze_vcf(
    vcf_path: Path,
    af_json_path: Path,
    family_file: Path,
    output_dir: Path,
    *,
    multiallelic: bool = True,
    af_threshold: float = DEFAULT_AF_THRESHOLD,
    debug: bool = False,
    memory_block: int = DEFAULT_MEMORY_BLOCK,
    block_size: int = DEFAULT_BLOCK_SIZE,
    segment_size: int = DEFAULT_SEGMENT_SIZE,
    short_format: bool = True,
    resume: bool = False,
    repeats_path: Path | None = None,
) -> AnalysisStats:
    """Scan a VCF, classify trios, and stream results to segmented TSV files."""
    if resume and segment_size <= 0:
        raise ValueError("--resume requires --segment-size > 0")

    af_table = load_af_json(af_json_path)
    relations = load_family_relations(family_file)

    checkpoint = load_checkpoint(output_dir) if resume else None
    if resume and checkpoint is None:
        raise FileNotFoundError(f"No checkpoint found in {output_dir}")
    if resume and checkpoint.completed:
        raise ValueError(f"Checkpoint in {output_dir} is already marked completed")

    if checkpoint is not None:
        writer = ResultWriter.from_checkpoint(
            output_dir,
            checkpoint,
            block_size=block_size,
            segment_size=segment_size,
            short_format=short_format,
        )
        resume_last_pos = checkpoint.last_pos
    else:
        writer = ResultWriter(
            output_dir,
            block_size=block_size,
            segment_size=segment_size,
            short_format=short_format,
        )
        resume_last_pos = -1

    repeat_filter: RepeatIntervalFilter | None = None
    if repeats_path is not None:
        repeat_filter = RepeatIntervalFilter(repeats_path)
        if resume_last_pos >= 0:
            repeat_filter.advance_past(resume_last_pos)

    try:
        opener = gzip.open if str(vcf_path).endswith(".gz") else open
        with opener(vcf_path, "rt", encoding="utf-8") as handle:
            trios_ind: list[tuple[int, int, int]] = []
            sample_header: list[str] = []

            for line in handle:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    sample_header = line.strip().split("\t")[9:]
                    _, trios_ind = build_trio_indices(sample_header, relations.trio_cl)
                    continue

                fields = line.rstrip().split("\t")
                pos = int(fields[1])
                if pos <= resume_last_pos:
                    continue
                if repeat_filter is not None and repeat_filter.in_repeat(pos):
                    continue

                if multiallelic:
                    _process_multiallelic_line(
                        line,
                        af_table,
                        af_threshold,
                        trios_ind,
                        sample_header,
                        writer,
                    )
                else:
                    _process_biallelic_line(
                        line,
                        af_table,
                        af_threshold,
                        trios_ind,
                        sample_header,
                        writer,
                    )

                log_memory_if_due(
                    writer.cumulative.variants_seen,
                    debug=debug,
                    memory_block=memory_block,
                )
    finally:
        if repeat_filter is not None:
            repeat_filter.close()
        writer.close()

    writer.finalize(completed=True)
    return _stats_from_writer(writer)


def _stats_from_writer(writer: ResultWriter) -> AnalysisStats:
    return AnalysisStats(
        variants_seen=writer.cumulative.variants_seen,
        alleles_tested=writer.cumulative.alleles_tested,
        inherited_entries=writer.inherited_entries,
        mendelian_bad_entries=writer.mendelian_bad_entries,
        inherited_variants=writer.inherited_variants,
        mendelian_bad_variants=writer.mendelian_bad_variants,
    )


def _process_multiallelic_line(
    line: str,
    af_table: dict[str, float],
    af_threshold: float,
    trios_ind: list[tuple[int, int, int]],
    sample_header: list[str],
    writer: ResultWriter,
) -> None:
    chrom, pos, keys, ref, alts = get_nfields(line, 5)
    if len(ref) > 1:
        return

    skeys = keys.split(";")
    salts = alts.split(",")
    sample_fields = line.rstrip().split("\t")[9:]
    writer.cumulative.variants_seen += 1

    for alt_index, key in enumerate(skeys, start=1):
        if alt_index > len(salts):
            break
        if not is_rare(af_table, key, af_threshold):
            continue

        alt = salts[alt_index - 1]
        if len(alt) > 1:
            continue

        writer.cumulative.alleles_tested += 1
        _process_trios_for_allele(
            chrom,
            pos,
            ref,
            alt,
            key,
            alt_index,
            sample_fields,
            sample_header,
            trios_ind,
            writer,
            clean_ad=True,
        )


def _process_biallelic_line(
    line: str,
    af_table: dict[str, float],
    af_threshold: float,
    trios_ind: list[tuple[int, int, int]],
    sample_header: list[str],
    writer: ResultWriter,
) -> None:
    chrom, pos, key, ref, alt = get_nfields(line, 5)
    if len(ref) > 1 and len(alt) > 1:
        return
    if not is_rare(af_table, key, af_threshold):
        return

    writer.cumulative.variants_seen += 1
    writer.cumulative.alleles_tested += 1
    sample_fields = line.rstrip().split("\t")[9:]
    _process_trios_for_allele(
        chrom,
        pos,
        ref,
        alt,
        key,
        1,
        sample_fields,
        sample_header,
        trios_ind,
        writer,
    )


def _process_trios_for_allele(
    chrom: str,
    pos: str,
    ref: str,
    alt: str,
    variant_key: str,
    alt_index: int,
    sample_fields: list[str],
    sample_header: list[str],
    trios_ind: list[tuple[int, int, int]],
    writer: ResultWriter,
    *,
    clean_ad: bool = False,
) -> None:
    parents_cache: dict[int, list[object]] = {}
    inherited_hits: dict[str, tuple[str, str, str, str]] = {}
    bad_hits: dict[str, tuple[str, str, str, str]] = {}

    for child_idx, mother_idx, father_idx in trios_ind:
        child_sample = sample_fields[child_idx]
        ac, child_gt, child_gq = get_good_site(child_sample, alt_index, clean_ad=clean_ad)
        if ac <= 0:
            continue

        if mother_idx in parents_cache:
            mac, mother_gt, mother_gq = parents_cache[mother_idx]
        else:
            mac, mother_gt, mother_gq = get_good_site(
                sample_fields[mother_idx], alt_index, clean_ad=clean_ad
            )
            parents_cache[mother_idx] = [mac, mother_gt, mother_gq]

        if father_idx in parents_cache:
            fac, father_gt, father_gq = parents_cache[father_idx]
        else:
            fac, father_gt, father_gq = get_good_site(
                sample_fields[father_idx], alt_index, clean_ad=clean_ad
            )
            parents_cache[father_idx] = [fac, father_gt, father_gq]

        pid = sample_header[child_idx]
        bucket: str | None = None
        record: tuple[str, str, str, str]

        if mac > 0 and fac > 0:
            if mac == 2 and fac == 2 and ac < 2:
                bucket = "mendelian_bad"
            else:
                bucket = "inherited"
            record = (mother_gt, father_gt, child_gt, child_gq)
        elif mac > 0:
            bucket = "mendelian_bad" if ac == 2 else "inherited"
            record = (mother_gt, "0/0", child_gt, child_gq)
        elif fac > 0:
            bucket = "mendelian_bad" if ac == 2 else "inherited"
            record = ("0/0", father_gt, child_gt, child_gq)
        else:
            continue

        if bucket == "inherited":
            inherited_hits[pid] = record
        else:
            bad_hits[pid] = record

    if inherited_hits:
        writer.write_inherited(chrom, pos, ref, alt, variant_key, inherited_hits)
    if bad_hits:
        writer.write_mendelian_bad(chrom, pos, ref, alt, bad_hits)


def save_run_params(
    output_dir: Path,
    *,
    vcf_path: Path,
    af_json_path: Path,
    family_file: Path,
    multiallelic: bool,
    af_threshold: float,
    debug: bool = False,
    memory_block: int = DEFAULT_MEMORY_BLOCK,
    block_size: int = DEFAULT_BLOCK_SIZE,
    segment_size: int = DEFAULT_SEGMENT_SIZE,
    short_format: bool = True,
    resume: bool = False,
    repeats_path: Path | None = None,
) -> Path:
    """Write the parameters for this run into the chromosome output directory."""
    output_dir.mkdir(parents=True, exist_ok=True)
    params_path = output_dir / "params.json"

    payload: dict[str, Any] = {
        "package_version": __version__,
        "run_at": datetime.now(timezone.utc).isoformat(),
        "vcf": str(vcf_path.resolve()),
        "af_json": str(af_json_path.resolve()),
        "family_file": str(family_file.resolve()),
        "output_dir": str(output_dir.resolve()),
        "multiallelic": multiallelic,
        "af_threshold": af_threshold,
        "debug": debug,
        "memory_block": memory_block,
        "block_size": block_size,
        "segment_size": segment_size,
        "short_format": short_format,
        "resume": resume,
        "remove_repeats": str(repeats_path.resolve()) if repeats_path is not None else None,
        "quality_filters": {
            "gq": DEFAULT_GQ,
            "dp": DEFAULT_DP,
            "ab": DEFAULT_AB,
            "haplo_dp": DEFAULT_HAPLO_DP,
            "haplo_ab": DEFAULT_HAPLO_AB,
        },
    }
    with params_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
    return params_path
