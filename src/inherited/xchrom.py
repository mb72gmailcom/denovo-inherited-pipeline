from __future__ import annotations

from inherited.constants import PAR1_END, PAR1_START, PAR2_END, PAR2_START

# hg38 PARs on chrX; intervals are closed on both ends: [start, end].
PAR1 = (PAR1_START, PAR1_END)
PAR2 = (PAR2_START, PAR2_END)

X_BUCKET_FEMALES = "females"
X_BUCKET_MALE_PAR1 = "male_par1"
X_BUCKET_MALE_PAR2 = "male_par2"
X_BUCKET_MALE_NONPAR = "male_nonPar"

X_BUCKETS = (
    X_BUCKET_FEMALES,
    X_BUCKET_MALE_PAR1,
    X_BUCKET_MALE_NONPAR,
    X_BUCKET_MALE_PAR2,
)


def is_x_chrom(chrom: str) -> bool:
    return chrom in ("X", "chrX")


def x_region(pos: int) -> str:
    """Return ``par1``, ``par2``, or ``nonPar`` for a chrX position."""
    if PAR1[0] <= pos <= PAR1[1]:
        return "par1"
    if PAR2[0] <= pos <= PAR2[1]:
        return "par2"
    return "nonPar"


def male_x_bucket(pos: int) -> str:
    region = x_region(pos)
    if region == "par1":
        return X_BUCKET_MALE_PAR1
    if region == "par2":
        return X_BUCKET_MALE_PAR2
    return X_BUCKET_MALE_NONPAR
