from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

from rna_mod_pipeline.config import MODIFICATIONS


REGIONS = ("5UTR", "CDS", "3UTR")
DATABASE_LABELS = {
    "ornament": "oRNAment",
    "encori": "ENCORI",
    "postar3": "POSTAR3",
}
DATABASE_KEYS = {label: key for key, label in DATABASE_LABELS.items()}
EXPECTED_DATABASE_COUNTS = {
    "ornament": 133,
    "encori": 281,
    "postar3": 216,
}


@dataclass(frozen=True)
class Analysis:
    analysis_id: str
    modification: str
    display: str
    sites: str
    bedmethyl: str


ANALYSES: Mapping[str, Analysis] = {
    key: Analysis(
        analysis_id=key,
        modification="m6a" if key == "m6a_rep2" else key,
        display=definition.display,
        sites=(
            f"refactored_outputs/{key}/phase1/"
            f"filtered_{definition.short_name}_metagene.tsv"
        ),
        bedmethyl=definition.bedmethyl,
    )
    for key, definition in MODIFICATIONS.items()
}
