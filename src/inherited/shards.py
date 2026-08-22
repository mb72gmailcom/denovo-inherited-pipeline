from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class VcfShard:
    path: Path
    chrom: str
    start: int
    end: int


def discover_vcf_shards(vcf_dir: Path, pattern: str) -> list[VcfShard]:
    """Find ``{pattern}.{chr}_{start}_{end}.vcf.gz`` shards in one directory.

    The directory is expected to hold one chromosome. The contig token is
    taken from the filenames. Intervals are treated as inclusive ``[start, end]``.
    """
    if not pattern:
        raise ValueError("--vcf-pattern must be a non-empty filename prefix")
    if not vcf_dir.is_dir():
        raise FileNotFoundError(f"VCF directory not found: {vcf_dir}")

    rx = re.compile(
        "^"
        + re.escape(pattern)
        + r"\.(.+)_(\d+)_(\d+)\.vcf(?:\.gz)?$"
    )
    shards: list[VcfShard] = []
    for path in vcf_dir.iterdir():
        if not path.is_file():
            continue
        match = rx.fullmatch(path.name)
        if match is None:
            continue
        chrom = match.group(1)
        start = int(match.group(2))
        end = int(match.group(3))
        if start > end:
            raise ValueError(
                f"Invalid shard coordinates in {path.name}: start {start} > end {end}"
            )
        shards.append(VcfShard(path=path, chrom=chrom, start=start, end=end))

    if not shards:
        expected = f"{pattern}.{{chr}}_{{start}}_{{end}}.vcf.gz"
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
