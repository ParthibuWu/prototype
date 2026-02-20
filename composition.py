# composition.py
from dataclasses import dataclass
from typing import Dict, Optional

CANONICAL = set("ACGT")

#It only counts true DNA bases, It excludes N,AMB


@dataclass(frozen=True)
#Immutable after creation,Safe to share,Cannot accidentally mutate counts later,preserves invariant
#Base composition must reflect the original sequence.
class BaseComposition:
    counts: Dict[str, int]                 # A,T,G,C,N,AMB
    denom: int                             # chosen denominator
    fractions: Dict[str, Optional[float]]  # fractions for A,T,G,C


def atgc_content(seq: str, mode: str = "canonical", treat_u_as_t: bool = True) -> BaseComposition:
    s = (seq or "").strip().upper()
#biologically U is equivalent to T for composition

    counts = {"A": 0, "T": 0, "G": 0, "C": 0, "N": 0, "AMB": 0}
#AMB = ambiguous base (not A, T, G, C, and not N) it could be R/Y/etc which are partially known
#Tracking AMB helps with knowing How much uncertainty exists in this sequence?

    for ch in s:
        if treat_u_as_t and ch == "U":  #If RNA sequence contains U, Convert U -> T, This lets RNA behave like DNA
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

    fractions = {b: (counts[b] / denom if denom > 0 else None) for b in "ATGC"}  #Canonical mode,counts A/T/G/C = 0,So denom = 0, Division would fail So instead Fractions become None.
    return BaseComposition(counts=counts, denom=denom, fractions=fractions)


def gc_fraction(seq: str, mode: str = "canonical", treat_u_as_t: bool = True) -> Optional[float]:
    comp = atgc_content(seq, mode=mode, treat_u_as_t=treat_u_as_t)
    if comp.denom == 0:
        return None
    return (comp.counts["G"] + comp.counts["C"]) / comp.denom
