"""Streamlit app for generating unique simulated STR child profiles.

Supported input layouts
-----------------------
1. Instrument matrix export (recommended): first row contains locus names,
   first column contains sample IDs, and cells contain reported genotypes.
   This is the layout used by Investigator 24plex and GlobalFiler exports.
2. Legacy structured table: family_id, parent_role, person_id plus one locus
   column per marker.

For an instrument matrix, choose the maternal and paternal sample IDs and any
known real children in the Streamlit interface. The workbook output preserves
the matrix layout: a blank first header cell, the original locus order, and
one simulated-child row per generated profile.
"""
from __future__ import annotations

import io
import itertools
import math
import re
from collections import Counter
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
import streamlit as st

REQUIRED_PARENT_COLUMNS = {"family_id", "parent_role", "person_id"}
GENOTYPE_SEPARATORS = re.compile(r"\s*[,/|;]\s*|\s+")
MAX_ENUMERATE_PROFILES = 1_000_000
# QS1 and QS2 are PCR quality-control markers in Investigator 24plex, not
# inherited DNA loci. Their expected successful-run calls are fixed.
QUALITY_MARKER_CALLS = {"qs1": "1", "qs2": "2"}
Profile = tuple[tuple[str, str] | None, ...]
ProfileKey = tuple[tuple[str, str] | None, ...]


class InputValidationError(ValueError):
    """Raised when uploaded genotype data cannot be safely simulated."""


@dataclass(frozen=True)
class FamilyParents:
    family_id: str
    mother_id: str
    father_id: str
    mother: dict[str, tuple[str, str] | None]
    father: dict[str, tuple[str, str] | None]


def read_raw_table(uploaded_file) -> pd.DataFrame:
    """Read a file as raw cells so the first unlabeled instrument column survives."""
    name = uploaded_file.name.lower()
    if name.endswith(".csv"):
        return pd.read_csv(uploaded_file, header=None, dtype=str, keep_default_na=False)
    if name.endswith((".xlsx", ".xls")):
        return pd.read_excel(uploaded_file, header=None, dtype=str, keep_default_na=False)
    raise InputValidationError("Upload a CSV, XLSX, or XLS file.")


def clean_text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def allele_sort_key(allele: str) -> tuple[int, float | str]:
    try:
        return (0, float(allele))
    except ValueError:
        return (1, allele.casefold())


def is_amel_locus(locus: str) -> bool:
    return clean_text(locus).casefold() in {"amel", "amelogenin"}


def quality_marker_call(locus: str) -> str | None:
    """Return a fixed successful-PCR call for non-genetic QS markers."""
    return QUALITY_MARKER_CALLS.get(clean_text(locus).casefold())


def is_quality_marker(locus: str) -> bool:
    return quality_marker_call(locus) is not None


def fixed_quality_genotype(locus: str) -> tuple[str, str]:
    call = quality_marker_call(locus)
    if call is None:
        raise InputValidationError(f"'{locus}' is not a configured quality marker.")
    return (call, call)


def is_y_locus(locus: str) -> bool:
    normalized = clean_text(locus).casefold().replace("-", "")
    return normalized.startswith("dys") or "yindel" in normalized


def parse_genotype(value: object, *, row_label: str, locus: str, allow_blank: bool = False) -> tuple[str, str] | None:
    """Normalize reported STR genotypes. A single AMEL X is treated as X,X internally."""
    raw = clean_text(value)
    if not raw:
        if allow_blank:
            return None
        raise InputValidationError(f"Missing genotype at locus '{locus}' for {row_label}.")

    alleles = [a.strip() for a in GENOTYPE_SEPARATORS.split(raw) if a.strip()]
    # Laboratory exports commonly render a homozygote as one value (for example
    # "9.3" rather than "9.3, 9.3"). Treat that form as a diploid homozygote.
    if len(alleles) == 1:
        alleles = [alleles[0], alleles[0]]
    if len(alleles) != 2:
        raise InputValidationError(
            f"Genotype at locus '{locus}' for {row_label} must contain one or two alleles; received '{raw}'."
        )
    return tuple(sorted(alleles, key=allele_sort_key))  # type: ignore[return-value]


def format_genotype(genotype: tuple[str, str] | None, locus: str) -> str:
    if genotype is None:
        return ""
    # QS1/QS2 are single fixed QC calls rather than diploid genotype calls.
    fixed_qc_call = quality_marker_call(locus)
    if fixed_qc_call is not None:
        return fixed_qc_call
    # Instrument files commonly report female Amel as X rather than X, X.
    if is_amel_locus(locus) and genotype == ("X", "X"):
        return "X"
    # DYS/Yindel values are normally displayed as a single haploid call.
    if is_y_locus(locus) and genotype[0] == genotype[1]:
        return genotype[0]
    return ", ".join(genotype)


def matrix_from_raw(raw: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    if raw.shape[0] < 2 or raw.shape[1] < 2:
        raise InputValidationError("The instrument matrix must have a header row, sample-ID column, and at least one sample.")

    headers = [clean_text(v) for v in raw.iloc[0].tolist()]
    loci = [h for h in headers[1:] if h]
    if not loci:
        raise InputValidationError("No locus names were found in row 1 of the uploaded instrument file.")
    if len(set(loci)) != len(loci):
        raise InputValidationError("Locus names in row 1 must be unique.")

    width = len(loci) + 1
    data = raw.iloc[1:, :width].copy()
    data.columns = ["sample_id", *loci]
    data["sample_id"] = data["sample_id"].map(clean_text)
    data = data[data["sample_id"] != ""].reset_index(drop=True)
    if data.empty:
        raise InputValidationError("No sample rows were found below the header row.")
    if data["sample_id"].duplicated().any():
        duplicates = data.loc[data["sample_id"].duplicated(), "sample_id"].head(5).tolist()
        raise InputValidationError(f"Sample IDs must be unique. Duplicate examples: {', '.join(duplicates)}")
    return data, loci


def validate_selected_matrix_family(
    matrix: pd.DataFrame,
    loci: list[str],
    family_id: str,
    mother_id: str,
    father_id: str,
) -> FamilyParents:
    if mother_id == father_id:
        raise InputValidationError("Choose two different sample IDs for the parents.")
    lookup = matrix.set_index("sample_id", drop=False)
    mother_row = lookup.loc[mother_id]
    father_row = lookup.loc[father_id]
    mother: dict[str, tuple[str, str] | None] = {}
    father: dict[str, tuple[str, str] | None] = {}
    for locus in loci:
        if is_quality_marker(locus):
            # QS1/QS2 are PCR QC markers. Do not validate or inherit their source cells.
            mother[locus] = fixed_quality_genotype(locus)
            father[locus] = fixed_quality_genotype(locus)
            continue
        y_marker = is_y_locus(locus)
        mother[locus] = parse_genotype(mother_row[locus], row_label=f"mother '{mother_id}'", locus=locus, allow_blank=y_marker)
        father[locus] = parse_genotype(father_row[locus], row_label=f"father '{father_id}'", locus=locus, allow_blank=y_marker)
        if y_marker and father[locus] is None:
            raise InputValidationError(f"The selected father has no reported genotype at Y-linked locus '{locus}'.")
    return FamilyParents(family_id, mother_id, father_id, mother, father)


def profile_from_matrix_row(row: pd.Series, loci: list[str], row_label: str) -> tuple[tuple[str, str] | None, ...]:
    profile: list[tuple[str, str] | None] = []
    for locus in loci:
        if is_quality_marker(locus):
            # QS calls cannot make two full DNA profiles different.
            profile.append(fixed_quality_genotype(locus))
        else:
            profile.append(parse_genotype(row[locus], row_label=row_label, locus=locus, allow_blank=is_y_locus(locus)))
    return tuple(profile)


def profile_comparison_key(profile: Profile, loci: list[str]) -> ProfileKey:
    """Whole-profile identity key excluding non-genetic QC marker columns."""
    return tuple(genotype for locus, genotype in zip(loci, profile) if not is_quality_marker(locus))


def locus_outcomes(maternal: tuple[str, str], paternal: tuple[str, str]) -> Counter[tuple[str, str]]:
    outcomes: Counter[tuple[str, str]] = Counter()
    for m in maternal:
        for f in paternal:
            outcomes[tuple(sorted((m, f), key=allele_sort_key))] += 1
    return outcomes


def possible_profile_count(per_locus: list[Counter[tuple[str, str] | None]]) -> int:
    return math.prod(len(outcomes) for outcomes in per_locus)


def enumerate_weighted_profiles(per_locus: list[Counter[tuple[str, str] | None]]) -> list[tuple[tuple[tuple[str, str] | None, ...], int]]:
    profiles: list[tuple[tuple[tuple[str, str] | None, ...], int]] = []
    for combination in itertools.product(*[list(outcomes.items()) for outcomes in per_locus]):
        profiles.append((tuple(item[0] for item in combination), math.prod(item[1] for item in combination)))
    return profiles


def sample_weighted_without_replacement(candidates: list[tuple[tuple[tuple[str, str] | None, ...], int]], n: int, rng: np.random.Generator) -> list[tuple[tuple[str, str] | None, ...]]:
    weights = np.asarray([weight for _, weight in candidates], dtype=float)
    weights /= weights.sum()
    idx = rng.choice(len(candidates), size=n, replace=False, p=weights)
    return [candidates[i][0] for i in idx]


def make_per_locus_outcomes(parents: FamilyParents, loci: list[str]) -> list[Counter[tuple[str, str] | None]]:
    """Construct outcomes, including Y-STR handling based on the Amel result when present."""
    amel_locus = next((locus for locus in loci if is_amel_locus(locus)), None)
    outcomes: list[Counter[tuple[str, str] | None]] = []
    for locus in loci:
        if is_quality_marker(locus):
            outcomes.append(Counter({fixed_quality_genotype(locus): 1}))
        elif is_y_locus(locus):
            paternal = parents.father[locus]
            if paternal is None:
                raise InputValidationError(f"Father genotype missing at Y-linked locus '{locus}'.")
            # Output either a paternal Y profile or blank. When Amel is present, correlation is restored below.
            outcomes.append(Counter({None: 1, paternal: 1}))
        else:
            maternal = parents.mother[locus]
            paternal = parents.father[locus]
            if maternal is None or paternal is None:
                raise InputValidationError(f"Both parents require a genotype at autosomal locus '{locus}'.")
            outcomes.append(Counter(locus_outcomes(maternal, paternal)))
    return outcomes


def reconcile_y_markers(profile: tuple[tuple[str, str] | None, ...], parents: FamilyParents, loci: list[str]) -> tuple[tuple[str, str] | None, ...]:
    amel_index = next((i for i, locus in enumerate(loci) if is_amel_locus(locus)), None)
    values = list(profile)
    if amel_index is not None:
        amel = values[amel_index]
        has_y = amel is not None and "Y" in {a.upper() for a in amel}
        for i, locus in enumerate(loci):
            if is_y_locus(locus):
                values[i] = parents.father[locus] if has_y else None
    return tuple(values)


def simulate_family(
    parents: FamilyParents,
    loci: list[str],
    blocked: set[ProfileKey],
    target: int,
    rng: np.random.Generator,
) -> tuple[list[tuple[tuple[str, str] | None, ...]], int, int]:
    per_locus = make_per_locus_outcomes(parents, loci)
    raw_total = possible_profile_count(per_locus)

    if raw_total <= MAX_ENUMERATE_PROFILES:
        merged: Counter[tuple[tuple[str, str] | None, ...]] = Counter()
        for profile, weight in enumerate_weighted_profiles(per_locus):
            merged[reconcile_y_markers(profile, parents, loci)] += weight
        candidates = [
            (profile, weight)
            for profile, weight in merged.items()
            if profile_comparison_key(profile, loci) not in blocked
        ]
        available = len(candidates)
        if available < target:
            return [], len(merged), available
        return sample_weighted_without_replacement(candidates, target, rng), len(merged), available

    selected: dict[ProfileKey, Profile] = {}
    attempts = 0
    max_attempts = max(100_000, target * 2_000)
    while len(selected) < target and attempts < max_attempts:
        attempts += 1
        sampled = []
        for outcomes in per_locus:
            states = list(outcomes)
            weights = np.asarray([outcomes[state] for state in states], dtype=float)
            weights /= weights.sum()
            sampled.append(states[int(rng.choice(len(states), p=weights))])
        profile = reconcile_y_markers(tuple(sampled), parents, loci)
        key = profile_comparison_key(profile, loci)
        if key not in blocked:
            selected[key] = profile
    return (list(selected.values()), raw_total, -1) if len(selected) == target else ([], raw_total, len(selected))


def output_matrix(profiles: list[tuple[tuple[str, str] | None, ...]], loci: list[str], family_id: str) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    for n, profile in enumerate(profiles, start=1):
        row = {"": f"{family_id}_SIM_{n:03d}"}
        row.update({locus: format_genotype(genotype, locus) for locus, genotype in zip(loci, profile)})
        rows.append(row)
    return pd.DataFrame(rows, columns=["", *loci])


def dataframe_to_xlsx(df: pd.DataFrame, sheet_name: str = "Simulated Children") -> bytes:
    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
        df.to_excel(writer, index=False, sheet_name=sheet_name)
        ws = writer.book[sheet_name]
        ws.freeze_panes = "B2"
        ws.auto_filter.ref = ws.dimensions
        for column in ws.columns:
            width = min(max(len(clean_text(cell.value)) for cell in column) + 2, 42)
            ws.column_dimensions[column[0].column_letter].width = width
    return buffer.getvalue()


def instrument_workflow(raw: pd.DataFrame) -> None:
    matrix, loci = matrix_from_raw(raw)
    st.success(f"Detected instrument matrix: {len(matrix)} sample row(s), {len(loci)} locus/loci.")
    st.dataframe(matrix, use_container_width=True, hide_index=True)

    names = matrix["sample_id"].tolist()
    c1, c2, c3 = st.columns(3)
    with c1:
        family_id = st.text_input("Case / family ID", value="CASE001")
    with c2:
        mother_id = st.selectbox("Mother sample", names)
    with c3:
        father_id = st.selectbox("Father sample", names, index=min(1, len(names) - 1))
    known_ids = st.multiselect("Known real-child sample(s) to exclude", names, help="Select any actual children present in this same export.")
    target, seed_text = generation_controls()

    if st.button("Generate simulated children", type="primary", key="matrix_generate"):
        try:
            parents = validate_selected_matrix_family(matrix, loci, clean_text(family_id) or "CASE001", mother_id, father_id)
            known = {
                profile_comparison_key(
                    profile_from_matrix_row(matrix.set_index("sample_id").loc[child_id], loci, f"known child '{child_id}'"),
                    loci,
                )
                for child_id in known_ids
            }
            seed = None if not seed_text.strip() else int(seed_text.strip())
            profiles, total, available = simulate_family(parents, loci, known, int(target), np.random.default_rng(seed))
            if not profiles:
                raise InputValidationError(
                    f"Requested {target} unique profiles, but only {available} permissible profiles are available "
                    f"from {total} possible child profiles after excluding known children."
                )
            st.session_state["matrix_output"] = output_matrix(profiles, loci, parents.family_id)
            st.session_state["matrix_summary"] = {"Generated": len(profiles), "Loci": len(loci), "Known child profiles excluded": len(known), "Theoretical profiles": total}
        except Exception as exc:
            st.error(str(exc))

    if "matrix_output" in st.session_state:
        st.subheader("Generation summary")
        st.json(st.session_state["matrix_summary"])
        output = st.session_state["matrix_output"]
        st.subheader("Simulated children — instrument matrix layout")
        st.dataframe(output, use_container_width=True, hide_index=True)
        st.download_button("Download Excel export", dataframe_to_xlsx(output), "simulated_child_genotypes.xlsx", "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
        st.download_button("Download CSV export", output.to_csv(index=False).encode("utf-8"), "simulated_child_genotypes.csv", "text/csv")


def generation_controls() -> tuple[int, str]:
    c1, c2 = st.columns([1, 2])
    with c1:
        target = st.number_input("Unique simulated children", min_value=1, max_value=10_000, value=200)
    with c2:
        seed = st.text_input("Optional random seed", value="", help="The same files and seed reproduce a run.")
    return int(target), seed


def main() -> None:
    st.set_page_config(page_title="Child Genotype Simulator", page_icon="🧬", layout="wide")
    st.title("Child Genotype Simulator")
    st.caption("Generate unique Mendelian-consistent simulated STR profiles and exclude complete profiles of known children.")
    st.info("Use the ArmedXpert allele table exports. The application detects locus names from row 1 and keeps that locus order in its export.")

    uploaded = st.file_uploader("Upload genotype export", type=["csv", "xlsx", "xls"])
    if not uploaded:
        st.stop()
    try:
        raw = read_raw_table(uploaded)
        instrument_workflow(raw)
    except Exception as exc:
        st.error(str(exc))


if __name__ == "__main__":
    main()
