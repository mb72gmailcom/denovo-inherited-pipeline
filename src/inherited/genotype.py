from __future__ import annotations

from inherited.constants import (
    DEFAULT_AB,
    DEFAULT_DP,
    DEFAULT_GQ,
    DEFAULT_HAPLO_AB,
    DEFAULT_HAPLO_DP,
)


def _parse_int_field(value: str) -> int | None:
    value = value.strip()
    if value in ("", "."):
        return None
    try:
        return int(value)
    except ValueError:
        return None


def _parse_ad_field(ad: str, clean_missing_as_zero: bool = False) -> list[int] | None:
    if not ad or ad.strip() == ".":
        return None

    values: list[int] = []
    for part in ad.split(","):
        part = part.strip()
        if part in ("", "."):
            if clean_missing_as_zero:
                values.append(0)
                continue
            return None
        try:
            values.append(int(part))
        except ValueError:
            return None
    return values or None


def _gt_alleles(gt: str) -> list[str]:
    return gt.replace("|", "/").split("/")


def is_good(
    gt: str,
    dp: str,
    ad: str,
    sb: str,
    gq: str,
    alt_index: int = 1,
    *,
    clean_missing_ad_as_zero: bool = False,
    dp_min: int = DEFAULT_DP,
    ab_min: float = DEFAULT_AB,
    gq_min: int = DEFAULT_GQ,
    alleles: list[str] | None = None,
) -> bool:
    if "." in gt:
        return False

    dp_value = _parse_int_field(dp)
    if dp_value is None or dp_value < dp_min:
        return False

    gq_value = _parse_int_field(gq)
    if gq_value is None or gq_value < gq_min:
        return False

    ads = _parse_ad_field(ad, clean_missing_as_zero=clean_missing_ad_as_zero)
    if ads is None or len(ads) <= alt_index:
        return False
    if sum(ads) == 0:
        return False

    # Require sufficient alt support only when the genotype includes that alt.
    # Homozygous-reference calls are allowed through with DP/GQ alone.
    parsed = alleles if alleles is not None else _gt_alleles(gt)
    has_alt = any(
        allele.isdigit() and int(allele) == alt_index for allele in parsed
    )
    if has_alt and ads[alt_index] / sum(ads) < ab_min:
        return False
    return True


def get_good_site(
    sample_field: str,
    alt_index: int = 1,
    *,
    clean_ad: bool = False,
    haploid: bool = False,
    skip_qc_if_no_alt: bool = False,
) -> tuple[int, str, str]:
    """Return allele count, GT, and GQ for one alternate allele.

    Returns ``-1`` when the site fails quality filters or FORMAT is incomplete.
    Returns ``0`` for a good-quality homozygous-reference call.

    When ``clean_ad`` is True, missing AD values (``.``) are treated as zero depth
    for that allele. Used on the multiallelic path so one missing ALT depth does
    not reject all alleles at the site.

    When ``haploid`` is True, use haploid depth/AB thresholds while still counting
    alternate alleles from the GT field.

    When ``skip_qc_if_no_alt`` is True and the genotype carries no copies of the
    queried alt, return ``(0, gt, gq)`` without DP/AD/GQ checks. Used for child
    lookups where non-carriers are discarded immediately.
    """
    parts = sample_field.split(":")
    if len(parts) < 6:
        return -1, ".", "0"

    gt, dp, ad, sb, gq = parts[0], parts[1], parts[2], parts[3], parts[4]
    if "." in gt:
        return -1, ".", "0"

    alleles = _gt_alleles(gt)
    ac = sum(allele.isdigit() and int(allele) == alt_index for allele in alleles)
    if skip_qc_if_no_alt and ac == 0:
        return 0, gt, gq

    dp_min = DEFAULT_HAPLO_DP if haploid else DEFAULT_DP
    ab_min = DEFAULT_HAPLO_AB if haploid else DEFAULT_AB
    if is_good(
        gt,
        dp,
        ad,
        sb,
        gq,
        alt_index,
        clean_missing_ad_as_zero=clean_ad,
        dp_min=dp_min,
        ab_min=ab_min,
        alleles=alleles,
    ):
        return ac, gt, gq
    return -1, ".", "0"
