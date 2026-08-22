from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path


_SHARD_NAME = re.compile(
    r"^(?P<stem>.+)\.(?P<chrom>[^._]+)_(?P<start>\d+)_(?P<end>\d+)\.vcf(?:\.gz)?$"
)


@dataclass(frozen=True)
class VcfShard:
    path: Path
    chrom: str
    start: int
    end: int


def discover_vcf_shards(vcf_dir: Path, pattern: str) -> list[VcfShard]:
    """Find shards named ``{stem}.{chr}_{start}_{end}.vcf.gz`` in one directory.

    ``pattern`` may be the callset stem (``SPARK.WGS.2026_08.gatk``) or that
    stem plus the contig (``...gatk.chr21`` or ``...gatk.chr21_``). The
    directory is expected to hold one chromosome. Intervals are inclusive
    ``[start, end]``.
    """
    if not pattern:
        raise ValueError("--vcf-pattern must be a non-empty filename prefix")
    if not vcf_dir.is_dir():
        raise FileNotFoundError(f"VCF directory not found: {vcf_dir}")

    prefix = pattern.rstrip("_")
    shards: list[VcfShard] = []
    for path in vcf_dir.iterdir():
        if not path.is_file():
            continue
        match = _SHARD_NAME.fullmatch(path.name)
        if match is None:
            continue
        stem = match.group("stem")
        chrom = match.group("chrom")
        if prefix not in (stem, f"{stem}.{chrom}"):
            continue
        start = int(match.group("start"))
        end = int(match.group("end"))
        if start > end:
            raise ValueError(
                f"Invalid shard coordinates in {path.name}: start {start} > end {end}"
            )
        shards.append(VcfShard(path=path, chrom=chrom, start=start, end=end))

    if not shards:
        expected = f"{prefix}.{{chr}}_{{start}}_{{end}}.vcf.gz"
        raise FileNotFoundError(f"No VCF shards matching {expected} in {vcf_dir}")

    shards.sort(key=lambda shard: (shard.start, shard.end, shard.path.name))
    chroms = {shard.chrom for shard in shards}
    if len(chroms) > 1:
        raise ValueError(
            f"Matched shards for more than one contig in {vcf_dir}: {sorted(chroms)}"
        )
    for prev, current in zip(shards, shards[1:]):
        if current.start <= prev.end:
            raise ValueError(
                f"Overlapping VCF shards: {prev.path.name} [{prev.start}, {prev.end}] "
                f"and {current.path.name} [{current.start}, {current.end}]"
            )
    return shards
