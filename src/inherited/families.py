from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


MALE_SEX_VALUES = frozenset({"1", "male"})
FEMALE_SEX_VALUES = frozenset({"2", "female"})
REQUIRED_FAMILY_COLUMNS = ("spid", "sfid", "father", "mother", "sex")
FAMILY_COLUMN_ALIASES: dict[str, tuple[str, ...]] = {
    "spid": ("spid", "ind_id"),
    "sfid": ("sfid", "family_id"),
    "father": ("father", "father_id"),
    "mother": ("mother", "mother_id"),
    "sex": ("sex",),
}


def load_family_column_map(path: Path) -> dict[str, str]:
    """Load a JSON object mapping internal column names to family-file headers."""
    with path.open(encoding="utf-8") as handle:
        data: Any = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"family-map must be a JSON object: {path}")
    mapping: dict[str, str] = {}
    for key, value in data.items():
        name = str(key)
        if name not in FAMILY_COLUMN_ALIASES:
            raise ValueError(
                f"Unknown family-map key {name!r}; expected one of "
                f"{list(REQUIRED_FAMILY_COLUMNS)}: {path}"
            )
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"family-map value for {name!r} must be a column name: {path}")
        mapping[name] = value.strip()
    return mapping


def resolve_family_columns(
    header: list[str],
    overrides: dict[str, str] | None = None,
) -> dict[str, int]:
    """Map internal column names to header indexes.

    Built-in aliases are used when exactly one match is present. ``overrides``
    from ``--family-map`` win for those keys. Two aliases in the same file
    without an override is an error.
    """
    present: dict[str, int] = {}
    for index, raw_name in enumerate(header):
        name = raw_name.strip()
        if name and name not in present:
            present[name] = index

    overrides = overrides or {}
    resolved: dict[str, int] = {}
    missing: list[str] = []
    for internal in REQUIRED_FAMILY_COLUMNS:
        if internal in overrides:
            target = overrides[internal]
            if target not in present:
                raise ValueError(
                    f"family-map {internal}={target!r} is not a column in the family file"
                )
            resolved[internal] = present[target]
            continue
        matches = [alias for alias in FAMILY_COLUMN_ALIASES[internal] if alias in present]
        if not matches:
            missing.append(internal)
            continue
        if len(matches) > 1:
            raise ValueError(
                f"Family file has ambiguous columns for {internal}: {matches}. "
                "Pass --family-map to choose one."
            )
        resolved[internal] = present[matches[0]]
    if missing:
        aliases = ", ".join(
            f"{name} ({' / '.join(FAMILY_COLUMN_ALIASES[name])})" for name in missing
        )
        raise ValueError(f"family file missing required columns: {aliases}")
    return resolved


def normalize_sex(value: str) -> str | None:
    """Return ``male``, ``female``, or None if sex is missing/unrecognized."""
    text = value.strip().lower()
    if not text:
        return None
    if text in MALE_SEX_VALUES:
        return "male"
    if text in FEMALE_SEX_VALUES:
        return "female"
    return None


@dataclass
class FamilyRelations:
    """Family relation tables built from the family TSV file."""

    trio: dict[str, tuple[str, str]] = field(default_factory=dict)
    trio_cl: dict[str, tuple[str, str]] = field(default_factory=dict)
    trio_all: dict[str, tuple[str, str]] = field(default_factory=dict)
    counts: dict[str, int] = field(default_factory=dict)
    family_size: dict[str, int] = field(default_factory=dict)
    trios_ids: list[list[str]] = field(default_factory=list)
    female_children: set[str] = field(default_factory=set)
    male_children: set[str] = field(default_factory=set)


def load_family_relations(
    path: Path,
    column_map: dict[str, str] | None = None,
) -> FamilyRelations:
    """Load family relations from a tab-separated file with a header row.

    Required fields (names or built-in aliases)::

        spid / ind_id
        sfid / family_id
        father / father_id
        mother / mother_id
        sex

    ``column_map`` overrides aliases (internal name → header in this file).
    Extra columns are ignored. Complete trios (``father`` and ``mother`` both
    not ``0``) with a recognized sex are retained for analysis.
    """
    relations = FamilyRelations()

    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"family file is empty: {path}") from exc

        inds = resolve_family_columns(header, column_map)

        for row in reader:
            if not row:
                continue
            if max(inds.values()) >= len(row):
                continue

            spid = row[inds["spid"]]
            family_id = row[inds["sfid"]]
            father_id = row[inds["father"]]
            mother_id = row[inds["mother"]]
            sex = normalize_sex(row[inds["sex"]])

            relations.family_size[family_id] = relations.family_size.get(family_id, 0) + 1
            relations.counts[spid] = relations.counts.get(spid, 0) + 1
            relations.trio_all[spid] = (father_id, mother_id)

            if father_id != "0" or mother_id != "0":
                relations.trio[spid] = (mother_id, father_id)

            if father_id != "0" and mother_id != "0":
                if sex is None:
                    continue
                relations.trio_cl[spid] = (mother_id, father_id)
                relations.trios_ids.append([spid, father_id, mother_id])
                if sex == "female":
                    relations.female_children.add(spid)
                else:
                    relations.male_children.add(spid)

    return relations


def build_trio_indices(
    sample_header: list[str],
    trio_cl: dict[str, tuple[str, str]],
    *,
    allowed_children: set[str] | None = None,
) -> tuple[dict[int, tuple[int, int]], list[tuple[int, int, int]]]:
    """Map VCF column indices for complete child-mother-father trios."""
    pid_to_idx = {pid: i for i, pid in enumerate(sample_header)}
    trio_ind: dict[int, tuple[int, int]] = {}
    trios_ind: list[tuple[int, int, int]] = []

    for child_idx, pid in enumerate(sample_header):
        if pid not in trio_cl:
            continue
        if allowed_children is not None and pid not in allowed_children:
            continue
        mother_id, father_id = trio_cl[pid]
        mother_idx = pid_to_idx.get(mother_id)
        father_idx = pid_to_idx.get(father_id)
        if mother_idx is None or father_idx is None:
            continue
        trio_ind[child_idx] = (mother_idx, father_idx)
        trios_ind.append((child_idx, mother_idx, father_idx))

    return trio_ind, trios_ind


def build_sexed_trio_indices(
    sample_header: list[str],
    relations: FamilyRelations,
) -> tuple[list[tuple[int, int, int]], list[tuple[int, int, int]]]:
    """Return (female_trios, male_trios) index lists present in the VCF."""
    _, female_trios = build_trio_indices(
        sample_header,
        relations.trio_cl,
        allowed_children=relations.female_children,
    )
    _, male_trios = build_trio_indices(
        sample_header,
        relations.trio_cl,
        allowed_children=relations.male_children,
    )
    return female_trios, male_trios
