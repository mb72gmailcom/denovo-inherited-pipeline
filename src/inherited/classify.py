from __future__ import annotations

from inherited.genotype import (
    alleles_only_ref_and_alt,
    is_hom_ref,
    is_mendelian_diploid,
)


def classify_trio(
    ac: int,
    mac: int,
    fac: int,
    m_gt: str,
    f_gt: str,
    c_gt: str,
    alt_index: int = 1,
) -> str | None:
    """Classify a diploid trio for one alternate allele.

    Returns:
        ``inherited``, ``mendelian_bad``, ``denovo``, or ``None`` to skip.

    Both parents must pass QC (``mac >= 0`` and ``fac >= 0``). Denovo requires
    homozygous-reference parental genotypes and a child genotype that uses only
    ref and the queried alt (so ``0/0 × 0/0 → 1/2`` is not denovo). When
    genotypes only use ref and the queried alt, cheap allele-count dosage rules
    are enough; full Mendelian checking runs only when another allele appears.
    """
    if ac <= 0:
        return None
    if mac < 0 or fac < 0:
        return None

    if mac == 0 and fac == 0:
        if not (is_hom_ref(m_gt) and is_hom_ref(f_gt)):
            return None
        if alleles_only_ref_and_alt(c_gt, str(alt_index)):
            return "denovo"
        return "mendelian_bad"

    # Dosage patterns that can never be Mendelian for this alt.
    if mac == 2 and fac == 2 and ac == 1:
        return "mendelian_bad"
    if ac == 2 and mac > 0 and fac == 0:
        return "mendelian_bad"
    if ac == 2 and fac > 0 and mac == 0:
        return "mendelian_bad"

    # Common path: site is effectively biallelic for this alt → inherited.
    alt_token = str(alt_index)
    if (
        alleles_only_ref_and_alt(m_gt, alt_token)
        and alleles_only_ref_and_alt(f_gt, alt_token)
        and alleles_only_ref_and_alt(c_gt, alt_token)
    ):
        return "inherited"

    # Other alleles present: require full genotype Mendelian consistency.
    if is_mendelian_diploid(m_gt, f_gt, c_gt):
        return "inherited"
    return "mendelian_bad"

def classify_mother_son(ac: int, mac: int) -> str | None:
    """Classify a male nonPAR chrX mother–son pair."""
    if ac <= 0 or mac < 0:
        return None
    if mac > 0:
        return "inherited"
    return "denovo"


def classify_father_son(ac: int, fac: int) -> str | None:
    """Classify a chrY father–son pair."""
    if ac <= 0 or fac < 0:
        return None
    if fac > 0:
        return "inherited"
    return "denovo"


def classify_trio_genotypes(
    ac: int,
    mac: int,
    fac: int,
    m_gt: str,
    f_gt: str,
    c_gt: str,
    c_gq: str,
    alt_index: int = 1,
) -> tuple[str, str, str, str] | None:
    """Return the genotype tuple to store, or None if the site should be skipped."""
    call_class = classify_trio(ac, mac, fac, m_gt, f_gt, c_gt, alt_index)
    if call_class is None:
        return None
    return m_gt, f_gt, c_gt, c_gq


def trio_bucket(
    ac: int,
    mac: int,
    fac: int,
    m_gt: str,
    f_gt: str,
    c_gt: str,
    alt_index: int = 1,
) -> str | None:
    """Return ``inherited``, ``mendelian_bad``, ``denovo``, or None."""
    return classify_trio(ac, mac, fac, m_gt, f_gt, c_gt, alt_index)
