from __future__ import annotations


def classify_trio(ac: int, mac: int, fac: int) -> str | None:
    """Classify a diploid trio for one alternate allele.

    Returns:
        ``inherited``, ``mendelian_bad``, ``denovo``, or ``None`` to skip.
    """
    if ac <= 0:
        return None

    if mac > 0 and fac > 0:
        if mac == 2 and fac == 2 and ac == 1:
            return "mendelian_bad"
        return "inherited"

    if mac > 0:
        return "mendelian_bad" if ac == 2 else "inherited"

    if fac > 0:
        return "mendelian_bad" if ac == 2 else "inherited"

    if mac == 0 and fac == 0:
        return "denovo"

    return None


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
    call_class = classify_trio(ac, mac, fac)
    if call_class is None:
        return None

    if mac > 0 and fac > 0:
        return m_gt, f_gt, c_gt, c_gq

    if mac > 0:
        return m_gt, "0/0", c_gt, c_gq

    if fac > 0:
        return "0/0", f_gt, c_gt, c_gq

    # denovo: both parents are good-quality non-carriers
    return m_gt, f_gt, c_gt, c_gq


def trio_bucket(
    ac: int,
    mac: int,
    fac: int,
) -> str | None:
    """Return ``inherited``, ``mendelian_bad``, ``denovo``, or None."""
    return classify_trio(ac, mac, fac)
