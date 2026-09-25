from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping


@dataclass(frozen=True)
class Modification:
    key: str
    display: str
    short_name: str
    expected_base: str
    expected_mod_codes: frozenset[str]
    folder: str
    bedmethyl: str
    filtered_sites: str
    metagene_sites: str
    drach_sites: str | None = None


MODIFICATIONS: Mapping[str, Modification] = {
    "m6a": Modification(
        key="m6a",
        display="m⁶A",
        short_name="m6A",
        expected_base="A",
        expected_mod_codes=frozenset({"a", "m6a"}),
        folder="m6a",
        bedmethyl="m6a/data/data.bedmethyl",
        filtered_sites="m6a/phase_1_results/filtered_m6A.tsv",
        drach_sites="m6a/phase_1_results/filtered_m6A_drach.tsv",
        metagene_sites="m6a/phase_1_results/filtered_m6A_metagene.tsv",
    ),
    "m5c": Modification(
        key="m5c",
        display="m⁵C",
        short_name="m5C",
        expected_base="C",
        expected_mod_codes=frozenset({"m", "m5c"}),
        folder="m5c",
        bedmethyl="m5c/data/m5c_pseU__m5C.bedmethyl",
        filtered_sites="m5c/phase_1_results/filtered_m5C.tsv",
        metagene_sites="m5c/phase_1_results/filtered_m5C_metagene.tsv",
    ),
    "pseu": Modification(
        key="pseu",
        display="Ψ",
        short_name="pseU",
        expected_base="T",
        expected_mod_codes=frozenset({"17802", "pseu", "psi"}),
        folder="pseudouridine",
        bedmethyl="pseudouridine/data/m5c_pseU__pseU_17802.bedmethyl",
        filtered_sites="pseudouridine/phase_1_results/filtered_pseU.tsv",
        metagene_sites="pseudouridine/phase_1_results/filtered_pseU_metagene.tsv",
    ),
    "m6a_rep2": Modification(
        key="m6a_rep2",
        display="m⁶A replicate 2",
        short_name="m6A",
        expected_base="A",
        expected_mod_codes=frozenset({"a", "m6a"}),
        folder="m6a_rep2",
        bedmethyl="m6a_rep2/data/data_rep2_m6a.bedmethyl",
        filtered_sites="m6a_rep2/phase_1_results/filtered_m6A.tsv",
        drach_sites="m6a_rep2/phase_1_results/filtered_m6A_drach.tsv",
        metagene_sites="m6a_rep2/phase_1_results/filtered_m6A_metagene.tsv",
    ),
}


RBP_ALIASES = {
    "A2BP1": "RBFOX1",
    "BRUNOL4": "CELF4",
    "BRUNOL5": "CELF5",
    "BRUNOL6": "CELF6",
    "DDX3": "DDX3X",
    "FMRP": "FMR1",
    "Fmr1": "FMR1",
    "HNRNPH": "HNRNPH1",
    "Khdrbs1": "KHDRBS1",
    "KSRP": "KHSRP",
    "P54NRB": "NONO",
    "Rbm38": "RBM38",
    "Rsf1": "RSF1",
    "SF2": "SRSF1",
    "STARPAP": "TUT1",
    "TDP43": "TARDBP",
    "U2AF65": "U2AF2",
    "YB1": "YBX1",
    "ZC3H14CONSTRUCT": "ZC3H14",
    "hnRNPLL": "HNRNPLL",
    "sf3b4": "SF3B4",
}

EXPRESSION_UNMAPPABLE_RBPS = frozenset({"Rbp1", "Rbp1like", "Rbp9"})


KNOWN_RELATED = {
    "m6a": {
        "METTL3",
        "METTL14",
        "WTAP",
        "VIRMA",
        "ZC3H13",
        "RBM15",
        "RBM15B",
        "HAKAI",
        "FTO",
        "ALKBH5",
        "YTHDC1",
        "YTHDC2",
        "YTHDF1",
        "YTHDF2",
        "YTHDF3",
        "IGF2BP1",
        "IGF2BP2",
        "IGF2BP3",
        "HNRNPC",
        "HNRNPA2B1",
        "ELAVL1",
        "HNRNPD",
        "FMR1",
        "FXR1",
        "FXR2",
        "SRSF3",
        "SRSF10",
    },
    "m5c": {
        "NSUN2",
        "NSUN3",
        "NSUN4",
        "NSUN5",
        "NSUN6",
        "NSUN7",
        "NOP2",
        "TRDMT1",
        "ALYREF",
        "YBX1",
    },
    "pseu": {
        "DKC1",
        "PUS1",
        "PUS7",
        "PUS7L",
        "PUS10",
        "TRUB1",
        "TRUB2",
        "RPUSD1",
        "RPUSD2",
        "RPUSD3",
        "RPUSD4",
    },
}


BEDMETHYL_COLUMNS = [
    "chrom",
    "start",
    "end",
    "mod_code",
    "score",
    "strand",
    "thick_start",
    "thick_end",
    "item_rgb",
    "Nvalid_cov",
    "fraction_modified",
    "Nmod",
    "Ncanonical",
    "Nother_mod",
    "Ndelete",
    "Nfail",
    "Ndiff",
    "Nnocall",
]


FILTERED_COLUMNS = [
    "chrom",
    "start",
    "end",
    "mod_code",
    "strand",
    "Nvalid_cov",
    "fraction_modified",
    "Nmod",
    "Ncanonical",
    "Nother_mod",
    "Ndelete",
    "Nfail",
    "Ndiff",
    "Nnocall",
]


def canonical_rbp(name: str) -> str:
    cleaned = str(name).strip()
    return RBP_ALIASES.get(cleaned, cleaned)


def expression_gene_symbol(name: str) -> str | None:
    cleaned = str(name).strip()
    if cleaned in EXPRESSION_UNMAPPABLE_RBPS:
        return None
    return canonical_rbp(cleaned).upper()


def modification(key: str) -> Modification:
    try:
        return MODIFICATIONS[key]
    except KeyError as exc:
        choices = ", ".join(MODIFICATIONS)
        raise ValueError(f"Unknown analysis '{key}'. Choose from: {choices}") from exc


def project_path(root: str | Path, relative: str) -> Path:
    return Path(root).resolve() / relative
