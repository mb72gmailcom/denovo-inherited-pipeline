from __future__ import annotations

from inherited.genotype import is_hom_ref, is_mendelian_diploid


def classify_trio(
    ac: int,
    mac: int,
    fac: int,
    m_gt: str,
    f_gt: str,
    c_gt: str,
) -> str | None:
    """Classify a diploid trio for one alternate allele.

    Returns:
        ``inherited``, ``mendelian_bad``, ``denovo``, or ``None`` to skip.

    Both parents must pass QC (``mac >= 0`` and ``fac >= 0``). Denovo requires
    homozygous-reference parental genotypes. Otherwise classification uses full
    genotype Mendelian consistency rather than allele-count dosage alone.
    """
    if ac <= 0:
        return None
    if mac < 0 or fac < 0:
        return None

    if mac == 0 and fac == 0:
        if is_hom_ref(m_gt) and is_hom_ref(f_gt):
            return "denovo"
        return None

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
) -> tuple[str, str, str, str] | None:
    """Return the genotype tuple to store, or None if the site should be skipped."""
    call_class = classify_trio(ac, mac, fac, m_gt, f_gt, c_gt)
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
) -> str | None:
    """Return ``inherited``, ``mendelian_bad``, ``denovo``, or None."""
    return classify_trio(ac, mac, fac, m_gt, f_gt, c_gt)
