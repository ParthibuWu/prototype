# composition.py
from dataclasses import dataclass
from typing import Dict, Optional

CANONICAL = set("ACGT")


@dataclass(frozen=True)
class BaseComposition:
    counts: Dict[str, int]                 # A,T,G,C,N,AMB
    denom: int                             # chosen denominator
    fractions: Dict[str, Optional[float]]  # fractions for A,T,G,C


def atgc_content(seq: str, mode: str = "canonical", treat_u_as_t: bool = True) -> BaseComposition:
    s = (seq or "").strip().upper()

    counts = {"A": 0, "T": 0, "G": 0, "C": 0, "N": 0, "AMB": 0}

    for ch in s:
        if treat_u_as_t and ch == "U":
            ch = "T"
        if ch in CANONICAL:
            counts[ch] += 1
        elif ch == "N":
            counts["N"] += 1
        else:
            counts["AMB"] += 1

    if mode == "raw":
        denom = len(s)
    elif mode == "canonical":
        denom = counts["A"] + counts["T"] + counts["G"] + counts["C"]
    else:
        raise ValueError("mode must be 'raw' or 'canonical'")

    fractions = {b: (counts[b] / denom if denom > 0 else None) for b in "ATGC"}
    return BaseComposition(counts=counts, denom=denom, fractions=fractions)


def gc_fraction(seq: str, mode: str = "canonical", treat_u_as_t: bool = True) -> Optional[float]:
    comp = atgc_content(seq, mode=mode, treat_u_as_t=treat_u_as_t)
    if comp.denom == 0:
        return None
    return (comp.counts["G"] + comp.counts["C"]) / comp.denom
