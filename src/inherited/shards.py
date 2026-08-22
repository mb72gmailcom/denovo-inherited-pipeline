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


def chrom_filename_tokens(chrom: str) -> tuple[str, ...]:
    """Return contig tokens that may appear in a split pVCF filename."""
    raw = chrom.strip()
    if not raw:
        raise ValueError("--chr must be a non-empty contig name")
    if raw.lower().startswith("chr"):
        bare = raw[3:]
        tokens = (raw, bare, f"chr{bare}")
    else:
        tokens = (raw, f"chr{raw}")
    seen: list[str] = []
    for token in tokens:
        if token and token not in seen:
            seen.append(token)
    return tuple(seen)


def discover_vcf_shards(vcf_dir: Path, pattern: str, chrom: str) -> list[VcfShard]:
    """Find ``{pattern}.{chr}_{start}_{end}.vcf.gz`` shards for one contig.

    Intervals are treated as inclusive ``[start, end]``. Matching is exact on
    the contig token so ``chr2`` does not pick up ``chr21`` / ``chr22``.
    """
    if not pattern:
        raise ValueError("--vcf-pattern must be a non-empty filename prefix")
    if not vcf_dir.is_dir():
        raise FileNotFoundError(f"VCF directory not found: {vcf_dir}")

    shards: list[VcfShard] = []
    seen: set[Path] = set()
    for token in chrom_filename_tokens(chrom):
        rx = re.compile(
            "^"
            + re.escape(pattern)
            + r"\."
            + re.escape(token)
            + r"_(\d+)_(\d+)\.vcf(?:\.gz)?$"
        )
        for path in vcf_dir.iterdir():
            if not path.is_file() or path in seen:
                continue
            match = rx.fullmatch(path.name)
            if match is None:
                continue
            start = int(match.group(1))
            end = int(match.group(2))
            if start > end:
                raise ValueError(
                    f"Invalid shard coordinates in {path.name}: start {start} > end {end}"
                )
            shards.append(VcfShard(path=path, chrom=token, start=start, end=end))
            seen.add(path)

    if not shards:
        expected = f"{pattern}.{chrom}_{{start}}_{{end}}.vcf.gz"
        raise FileNotFoundError(f"No VCF shards matching {expected} in {vcf_dir}")

    shards.sort(key=lambda shard: (shard.start, shard.end, shard.path.name))
    chroms = {shard.chrom for shard in shards}
    if len(chroms) > 1:
        raise ValueError(
            f"Matched shards for more than one contig token in {vcf_dir}: {sorted(chroms)}"
        )
    for prev, current in zip(shards, shards[1:]):
        if current.start <= prev.end:
            raise ValueError(
                f"Overlapping VCF shards: {prev.path.name} [{prev.start}, {prev.end}] "
                f"and {current.path.name} [{current.start}, {current.end}]"
            )
    return shards
