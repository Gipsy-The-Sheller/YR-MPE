"""Codon-based selection analysis utilities.

V1 implements the Nei-Gojobori (NG86) counting method for pairwise dN/dS and
an analytic Z-test (delta method). MLE dN/dS is delegated to IQ-TREE 3 with a
codon model (GY or MG); this module only parses its omega.

The implementation is intentionally dependency-light: Biopython is used for the
genetic-code tables, everything else is stdlib.
"""

from __future__ import annotations

import itertools
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from Bio.Data import CodonTable
from Bio.SeqRecord import SeqRecord

DNA_BASES = "ACGT"

GENETIC_CODES = {
    "Standard (1)": 1,
    "Vertebrate mitochondrial (2)": 2,
    "Mold/protozoan mitochondrial (4)": 4,
    "Invertebrate mitochondrial (5)": 5,
    "Bacterial/plastid (11)": 11,
}


@dataclass
class PairwiseSelection:
    """One pairwise NG86 result."""

    seq1: str
    seq2: str
    sites: int
    synonymous_sites: float
    nonsynonymous_sites: float
    synonymous_differences: float
    nonsynonymous_differences: float
    pS: float
    pN: float
    dS: Optional[float]
    dN: Optional[float]
    omega: Optional[float]
    z: Optional[float]
    p_two_sided: Optional[float]
    p_positive: Optional[float]
    p_negative: Optional[float]
    conclusion: str
    alignment: str = ""
    se_dS: Optional[float] = None
    se_dN: Optional[float] = None
    method: str = "NG86"


@dataclass
class AlignmentSummary:
    alignment: str
    sequences: int
    codons: int
    usable_codons: int
    skipped_stop: int
    skipped_ambiguous: int
    results: List[PairwiseSelection] = field(default_factory=list)
    names: List[str] = field(default_factory=list)
    method: str = "NG86"


def get_codon_table(code_id: int):
    try:
        return CodonTable.unambiguous_dna_by_id[int(code_id)]
    except Exception:
        return CodonTable.unambiguous_dna_by_id[1]


def _clean_sequence(sequence: str) -> str:
    return str(sequence).upper().replace("U", "T").replace(".", "-")


def _is_clean_codon(codon: str) -> bool:
    return len(codon) == 3 and all(base in DNA_BASES for base in codon)


def synonymous_sites(codon: str, table) -> Optional[float]:
    """Fraction of synonymous sites in one codon (NG86).

    Changes that create a stop codon are excluded from both the numerator and
    the denominator, following the NG86 implementation in MEGA.
    """
    if not _is_clean_codon(codon):
        return None
    aa = table.forward_table.get(codon)
    if aa is None:  # stop codon
        return None
    total = 0.0
    for position in range(3):
        synonymous = 0
        valid = 0
        for base in DNA_BASES:
            if base == codon[position]:
                continue
            alt = codon[:position] + base + codon[position + 1 :]
            alt_aa = table.forward_table.get(alt)
            if alt_aa is None:  # stop codon: not a possible substitution
                continue
            valid += 1
            if alt_aa == aa:
                synonymous += 1
        if valid:
            total += synonymous / valid
    return total


def _step_type(codon1: str, codon2: str, table) -> Optional[str]:
    aa1 = table.forward_table.get(codon1)
    aa2 = table.forward_table.get(codon2)
    if aa1 is None or aa2 is None:
        return None
    return "S" if aa1 == aa2 else "N"


def count_codon_differences(codon1: str, codon2: str, table) -> Optional[Tuple[float, float]]:
    """Return (synonymous, nonsynonymous) differences, averaged over paths.

    Paths that pass through a stop codon are discarded; the remaining
    shortest paths are averaged. MEGA's pathway tables normalise by the same
    set of valid paths.
    """
    if not (_is_clean_codon(codon1) and _is_clean_codon(codon2)):
        return None
    positions = [i for i in range(3) if codon1[i] != codon2[i]]
    if not positions:
        return (0.0, 0.0)

    if len(positions) == 1:
        kind = _step_type(codon1, codon2, table)
        if kind is None:
            return None
        return (1.0, 0.0) if kind == "S" else (0.0, 1.0)

    paths = []
    for order in itertools.permutations(positions):
        current = codon1
        syn = 0.0
        non = 0.0
        valid = True
        for position in order:
            nxt = current[:position] + codon2[position] + current[position + 1 :]
            kind = _step_type(current, nxt, table)
            if kind is None:
                valid = False
                break
            if kind == "S":
                syn += 1.0
            else:
                non += 1.0
            current = nxt
        if valid and current == codon2:
            paths.append((syn, non))

    if not paths:
        return None
    return (
        sum(path[0] for path in paths) / len(paths),
        sum(path[1] for path in paths) / len(paths),
    )


def _jc_distance(p: float) -> Optional[float]:
    if p < 0.0 or p > 0.74:
        return None
    return -0.75 * math.log(1.0 - (4.0 / 3.0) * p)


def _normal_sf(z: float) -> float:
    """Upper tail of the standard normal."""
    return 0.5 * math.erfc(z / math.sqrt(2.0))


def _z_from_counts(
    synonymous_sites: float,
    nonsynonymous_sites: float,
    synonymous_differences: float,
    nonsynonymous_differences: float,
) -> Tuple[Optional[float], Optional[float], Optional[float], Optional[float], Optional[float]]:
    # MEGA's large-sample Z-test: Var(dN - dS) = Var(dN) + Var(dS).
    # The two variances use the number of synonymous / nonsynonymous sites,
    # not the number of codons, as the binomial denominator.
    if synonymous_sites < 1.0 or nonsynonymous_sites < 1.0:
        return None, None, None, None, None
    pS = synonymous_differences / synonymous_sites
    pN = nonsynonymous_differences / nonsynonymous_sites
    if pS > 1.0 or pN > 1.0:
        return pS, pN, None, None, None
    dS = _jc_distance(pS)
    dN = _jc_distance(pN)
    if dS is None or dN is None:
        return pS, pN, dS, dN, None

    var_pS = pS * (1.0 - pS) / synonymous_sites
    var_pN = pN * (1.0 - pN) / nonsynonymous_sites
    dS_prime = 1.0 / (1.0 - (4.0 / 3.0) * pS)
    dN_prime = 1.0 / (1.0 - (4.0 / 3.0) * pN)
    var_dS = dS_prime * dS_prime * var_pS
    var_dN = dN_prime * dN_prime * var_pN
    variance = var_dN + var_dS
    if variance <= 0:
        return pS, pN, dS, dN, None
    z = (dN - dS) / math.sqrt(variance)
    return pS, pN, dS, dN, z


def pairwise_ng86(seq1: str, seq2: str, table, alignment: str = "",
                  name1: str = "", name2: str = "") -> Optional[PairwiseSelection]:
    """NG86 pairwise dN/dS + analytic Z-test."""
    s1 = _clean_sequence(seq1)
    s2 = _clean_sequence(seq2)
    if len(s1) != len(s2) or len(s1) % 3 != 0:
        return None

    S = 0.0
    N = 0.0
    Sd = 0.0
    Nd = 0.0
    used = 0
    skipped_stop = 0
    skipped_ambiguous = 0

    for offset in range(0, len(s1), 3):
        c1 = s1[offset : offset + 3]
        c2 = s2[offset : offset + 3]
        if not (_is_clean_codon(c1) and _is_clean_codon(c2)):
            skipped_ambiguous += 1
            continue
        if table.forward_table.get(c1) is None or table.forward_table.get(c2) is None:
            skipped_stop += 1
            continue
        s1_sites = synonymous_sites(c1, table)
        s2_sites = synonymous_sites(c2, table)
        if s1_sites is None or s2_sites is None:
            skipped_stop += 1
            continue
        diffs = count_codon_differences(c1, c2, table)
        if diffs is None:
            skipped_stop += 1
            continue
        site_s = (s1_sites + s2_sites) / 2.0
        S += site_s
        N += 3.0 - site_s
        Sd += diffs[0]
        Nd += diffs[1]
        used += 1

    if used == 0:
        return None

    pS, pN, dS, dN, z = _z_from_counts(S, N, Sd, Nd)
    omega = None
    if dS is not None and dN is not None and dS > 0:
        omega = dN / dS

    p_two = p_positive = p_negative = None
    conclusion = "not estimable"
    if z is not None:
        p_two = 2.0 * _normal_sf(abs(z))
        p_positive = _normal_sf(z)
        p_negative = 1.0 - _normal_sf(z)
        if p_two < 0.05:
            if dN > dS:
                conclusion = "positive selection"
            elif dN < dS:
                conclusion = "purifying selection"
            else:
                conclusion = "not significant"
        else:
            conclusion = "not significant"
    elif dS is None or dN is None:
        conclusion = "saturated/not estimable"

    return PairwiseSelection(
        seq1=name1 or "seq1",
        seq2=name2 or "seq2",
        sites=used,
        synonymous_sites=S,
        nonsynonymous_sites=N,
        synonymous_differences=Sd,
        nonsynonymous_differences=Nd,
        pS=pS if pS is not None else float("nan"),
        pN=pN if pN is not None else float("nan"),
        dS=dS,
        dN=dN,
        omega=omega,
        z=z,
        p_two_sided=p_two,
        p_positive=p_positive,
        p_negative=p_negative,
        conclusion=conclusion,
        alignment=alignment,
    )


def analyze_alignment(records: Iterable[SeqRecord], code_id: int = 1,
                      alignment: str = "") -> AlignmentSummary:
    records = list(records)
    if not records:
        return AlignmentSummary(alignment, 0, 0, 0, 0, 0)
    table = get_codon_table(code_id)
    sequences = [_clean_sequence(record.seq) for record in records]
    names = [record.id or "seq%d" % (i + 1) for i, record in enumerate(records)]
    lengths = {len(sequence) for sequence in sequences}
    if len(lengths) != 1:
        raise ValueError("Sequences are not aligned to the same length")
    full_length = lengths.pop()
    if full_length % 3 != 0:
        raise ValueError("Alignment length is not a multiple of 3")

    results = []
    for i, j in itertools.combinations(range(len(sequences)), 2):
        result = pairwise_ng86(
            sequences[i], sequences[j], table,
            alignment=alignment, name1=names[i], name2=names[j],
        )
        if result is not None:
            results.append(result)

    # Count skipped categories once for the summary.
    skipped_stop = 0
    skipped_ambiguous = 0
    usable = 0
    for offset in range(0, full_length, 3):
        codons = [sequence[offset : offset + 3] for sequence in sequences]
        if any(not _is_clean_codon(codon) for codon in codons):
            skipped_ambiguous += 1
        elif any(table.forward_table.get(codon) is None for codon in codons):
            skipped_stop += 1
        else:
            usable += 1

    return AlignmentSummary(
        alignment=alignment or "alignment",
        sequences=len(records),
        codons=full_length // 3,
        usable_codons=usable,
        skipped_stop=skipped_stop,
        skipped_ambiguous=skipped_ambiguous,
        results=results,
        names=names,
        method="NG86",
    )


# ----------------------------------------------------------------------
# Yang & Nielsen (2000), equal pathway weighting (PAML yn00 default)
# ----------------------------------------------------------------------

_PAML_BASES = "TCAG"
_BASE_INDEX = {base: index for index, base in enumerate(_PAML_BASES)}


def _codon_blocks(sequence: str) -> List[str]:
    return [sequence[i : i + 3] for i in range(0, len(sequence) - 2, 3)]


def _trim_pair_yn00(seq1: str, seq2: str, table) -> Optional[Tuple[str, str]]:
    """Keep codons that are unambiguous and sense in both sequences."""
    s1 = _clean_sequence(seq1)
    s2 = _clean_sequence(seq2)
    if len(s1) != len(s2) or len(s1) % 3 != 0:
        return None
    kept1 = []
    kept2 = []
    for offset in range(0, len(s1), 3):
        c1 = s1[offset : offset + 3]
        c2 = s2[offset : offset + 3]
        if not (_is_clean_codon(c1) and _is_clean_codon(c2)):
            continue
        if table.forward_table.get(c1) is None or table.forward_table.get(c2) is None:
            continue
        kept1.append(c1)
        kept2.append(c2)
    if not kept1:
        return None
    return "".join(kept1), "".join(kept2)


def _sense_codons(table) -> List[str]:
    return [
        a + b + c
        for a in _PAML_BASES
        for b in _PAML_BASES
        for c in _PAML_BASES
        if table.forward_table.get(a + b + c) is not None
    ]


def _f3x4_and_pi(seq1: str, seq2: str, table) -> Tuple[List[List[float]], Dict[str, float]]:
    """Average the two sequences' position-specific base frequencies."""
    codons1 = _codon_blocks(seq1)
    codons2 = _codon_blocks(seq2)
    n = len(codons1)
    counts = [[0.0] * 4 for _ in range(3)]
    for codons in (codons1, codons2):
        for codon in codons:
            for position, base in enumerate(codon):
                counts[position][_BASE_INDEX[base]] += 0.5
    for position in range(3):
        total = sum(counts[position])
        if total > 0:
            counts[position] = [value / total for value in counts[position]]
        else:
            counts[position] = [0.25] * 4
    pi = {}
    for codon in _sense_codons(table):
        value = (
            counts[0][_BASE_INDEX[codon[0]]]
            * counts[1][_BASE_INDEX[codon[1]]]
            * counts[2][_BASE_INDEX[codon[2]]]
        )
        if value > 0:
            pi[codon] = value
    total = sum(pi.values())
    if total > 0:
        for codon in pi:
            pi[codon] /= total
    return counts, pi


def _count_sites_yn00(seq: str, pi: Dict[str, float], kappa: float, table):
    """Return S, N and base frequencies at synonymous/nonsynonymous sites."""
    codons = _codon_blocks(seq)
    synonymous = 0.0
    nonsynonymous = 0.0
    fb_s = [0.0] * 4
    fb_n = [0.0] * 4
    for codon in codons:
        aa = table.forward_table.get(codon)
        if aa is None:
            continue
        for position in range(3):
            original = _BASE_INDEX[codon[position]]
            for alt_base in _PAML_BASES:
                alt_index = _BASE_INDEX[alt_base]
                if alt_index == original:
                    continue
                alt = codon[:position] + alt_base + codon[position + 1 :]
                alt_aa = table.forward_table.get(alt)
                if alt_aa is None:  # stop codon
                    continue
                rate = pi.get(alt, 0.0)
                if original + alt_index in (1, 5):  # transition
                    rate *= kappa
                if alt_aa == aa:
                    synonymous += rate
                    fb_s[original] += rate
                else:
                    nonsynonymous += rate
                    fb_n[original] += rate
    scale = 3.0 * len(codons) / (synonymous + nonsynonymous) if synonymous + nonsynonymous else 0.0
    synonymous *= scale
    nonsynonymous *= scale
    for frequencies in (fb_s, fb_n):
        total = sum(frequencies)
        if total > 0:
            for index in range(4):
                frequencies[index] /= total
    return synonymous, nonsynonymous, fb_s, fb_n


def _is_transition(base1: int, base2: int) -> bool:
    return base1 + base2 in (1, 5)


def _count_diffs_yn00(seq1: str, seq2: str, table) -> Tuple[float, float, float, float]:
    """Syn/nonsyn transitions/transversions, equal weighting of valid paths."""
    totals = [0.0, 0.0, 0.0, 0.0]  # sts, stv, nts, ntv
    for codon1, codon2 in zip(_codon_blocks(seq1), _codon_blocks(seq2)):
        if codon1 == codon2:
            continue
        aa1 = table.forward_table.get(codon1)
        aa2 = table.forward_table.get(codon2)
        if aa1 is None or aa2 is None:
            continue
        positions = [i for i in range(3) if codon1[i] != codon2[i]]
        ndiff = len(positions)
        if ndiff == 1:
            position = positions[0]
            trans = _is_transition(_BASE_INDEX[codon1[position]],
                                   _BASE_INDEX[codon2[position]])
            if aa1 == aa2:
                totals[0 if trans else 1] += 1.0
            else:
                totals[2 if trans else 3] += 1.0
            continue

        valid_paths = []
        for order in itertools.permutations(positions):
            current = list(codon1)
            path = [0.0, 0.0, 0.0, 0.0]
            valid = True
            for position in order:
                nxt = current[:]
                nxt[position] = codon2[position]
                current_codon = "".join(current)
                next_codon = "".join(nxt)
                current_aa = table.forward_table.get(current_codon)
                next_aa = table.forward_table.get(next_codon)
                if next_aa is None:
                    valid = False
                    break
                trans = _is_transition(_BASE_INDEX[codon1[position]],
                                       _BASE_INDEX[codon2[position]])
                if current_aa == next_aa:
                    path[0 if trans else 1] += 1.0
                else:
                    path[2 if trans else 3] += 1.0
                current = nxt
            if valid:
                valid_paths.append(path)
        if valid_paths:
            for index in range(4):
                totals[index] += sum(path[index] for path in valid_paths) / len(valid_paths)
        elif ndiff == 2:
            totals[2] += 0.5
            totals[3] += 1.5
        else:
            totals[2] += 0.5
            totals[3] += 2.5
    return tuple(totals)


def _distance_f84(n: float, P: float, Q: float, pi: Sequence[float]):
    """PAML's F84 distance, with K80/JC fallbacks. Returns t, SE, kappa, status."""
    if n <= 0:
        return None, None, None, 1
    P = max(0.0, float(P))
    Q = max(0.0, float(Q))
    if P + Q > 1:
        return 99.0, None, 1.0, 3
    if len(pi) != 4 or abs(sum(pi) - 1.0) > 1e-6:
        total = sum(pi)
        if total <= 0:
            pi = [0.25] * 4
        else:
            pi = [value / total for value in pi]

    q_small = min(1e-10, 0.1 / n)
    Y = pi[0] + pi[1]
    R = pi[2] + pi[3]
    tc = pi[0] * pi[1]
    ag = pi[2] * pi[3]
    fail_f84 = False
    fail_k80 = False
    if Q < q_small:
        fail_f84 = fail_k80 = True
    elif Y <= 0 or R <= 0 or (tc <= 0 and ag <= 0):
        fail_f84 = True
    else:
        A = tc / Y + ag / R
        B = tc + ag
        C = Y * R
        a = (2 * B + 2 * (tc * R / Y + ag * Y / R) * (1 - Q / (2 * C)) - P) / (2 * A)
        b = 1 - Q / (2 * C)
        if a <= 0 or b <= 0:
            fail_f84 = True
        else:
            a = -0.5 * math.log(a)
            b = -0.5 * math.log(b)
            if b <= 0:
                fail_f84 = True
            else:
                k_f84 = a / b - 1
                t = 4 * b * (tc * (1 + k_f84 / Y) + ag * (1 + k_f84 / R) + C)
                kappa = (B + (tc / Y + ag / R) * k_f84) / B
                se = None
                den = A * C - C * P / 2 - (A - B) * Q / 2
                if den != 0 and C - Q / 2 != 0:
                    sa = A * C / den
                    sb = A * (A - B) / den - (A - B - C) / (C - Q / 2)
                    var = (sa * sa * P + sb * sb * Q - (sa * P + sb * Q) ** 2) / n
                    se = math.sqrt(max(0.0, var))
                return t, se, min(kappa, 999.0), 0

    if fail_f84 and not fail_k80:
        a = 1 - 2 * P - Q
        b = 1 - 2 * Q
        if a <= 0 or b <= 0:
            fail_k80 = True
        else:
            a = -math.log(a)
            b = -math.log(b)
            if b <= 0:
                fail_k80 = True
            else:
                kappa = (0.5 * a - 0.25 * b) / (0.25 * b)
                t = 0.5 * a + 0.25 * b
                se = None
                sa = 1 / (1 - 2 * P - Q)
                sb = (sa + 1 / (1 - 2 * Q)) / 2
                var = (sa * sa * P + sb * sb * Q - (sa * P + sb * Q) ** 2) / n
                se = math.sqrt(max(0.0, var))
                return t, se, min(kappa, 999.0), 1

    p = P + Q
    if p >= 0.75:
        p = 0.75 * (n - 1.0) / n
    t = -0.75 * math.log(1 - p * 4 / 3)
    se = None
    if 3 - 4 * p > 0:
        se = math.sqrt(9 * p * (1 - p) / n) / (3 - 4 * p)
    return t, se, 1.0, 2


def _is_non_degenerate(codon: str, aa: str, position: int, table) -> bool:
    for alt in _PAML_BASES:
        if alt == codon[position]:
            continue
        alt_codon = codon[:position] + alt + codon[position + 1 :]
        if table.forward_table.get(alt_codon) == aa:
            return False
    return True


def _is_fourfold(codon: str, aa: str, table) -> bool:
    for alt in _PAML_BASES:
        if table.forward_table.get(codon[:2] + alt) != aa:
            return False
    return True


def _estimate_kappa_pair(seq1: str, seq2: str, table, code_id: int = 1) -> float:
    """PAML's kappa estimate from non-degenerate and 4-fold sites."""
    f0 = [[0.0] * 4 for _ in range(4)]
    f1 = [[0.0] * 4 for _ in range(4)]
    for codon1, codon2 in zip(_codon_blocks(seq1), _codon_blocks(seq2)):
        aa1 = table.forward_table.get(codon1)
        aa2 = table.forward_table.get(codon2)
        if aa1 is None or aa2 is None:
            continue
        b1 = [_BASE_INDEX[base] for base in codon1]
        b2 = [_BASE_INDEX[base] for base in codon2]
        for position in range(3):
            if (_is_non_degenerate(codon1, aa1, position, table)
                    and _is_non_degenerate(codon2, aa2, position, table)):
                f0[b1[position]][b2[position]] += 0.5
                f0[b2[position]][b1[position]] += 0.5
        if (aa1 == aa2 and _is_fourfold(codon1, aa1, table)
                and _is_fourfold(codon2, aa2, table)):
            f1[b1[2]][b2[2]] += 0.5
            f1[b2[2]][b1[2]] += 0.5

    estimates = []
    for matrix in (f0, f1):
        size = sum(sum(row) for row in matrix)
        if size <= 0:
            continue
        normalized = [[value / size for value in row] for row in matrix]
        P = 2 * (normalized[0][1] + normalized[2][3])
        diagonal = sum(normalized[i][i] for i in range(4))
        Q = 1 - diagonal - P
        pi = [sum(normalized[i]) for i in range(4)]
        _, _, kappa, _ = _distance_f84(size, P, Q, pi)
        if kappa is not None and kappa > 0:
            estimates.append((kappa, size))
    if not estimates:
        return 10.0 if code_id == 2 else 2.0
    return sum(kappa * weight for kappa, weight in estimates) / sum(
        weight for _, weight in estimates
    )


def pairwise_yn00(seq1: str, seq2: str, table, kappa: Optional[float] = None,
                  alignment: str = "", name1: str = "", name2: str = "",
                  code_id: int = 1) -> Optional[PairwiseSelection]:
    """Pairwise Yang & Nielsen (2000) dN/dS, equal pathway weighting."""
    trimmed = _trim_pair_yn00(seq1, seq2, table)
    if trimmed is None:
        return None
    seq1, seq2 = trimmed
    _, pi = _f3x4_and_pi(seq1, seq2, table)
    if kappa is None:
        kappa = _estimate_kappa_pair(seq1, seq2, table, code_id=code_id)

    s1, n1, fbs1, fbn1 = _count_sites_yn00(seq1, pi, kappa, table)
    s2, n2, fbs2, fbn2 = _count_sites_yn00(seq2, pi, kappa, table)
    S = (s1 + s2) / 2.0
    N = (n1 + n2) / 2.0
    fb_s = [(fbs1[i] + fbs2[i]) / 2.0 for i in range(4)]
    fb_n = [(fbn1[i] + fbn2[i]) / 2.0 for i in range(4)]
    sdts, sdtv, ndts, ndtv = _count_diffs_yn00(seq1, seq2, table)
    sd = sdts + sdtv
    nd = ndts + ndtv

    dS = dN = se_ds = se_dn = None
    omega = z = p_two = p_positive = p_negative = None
    saturated = S < 1.0 or N < 1.0
    if not saturated:
        dS, se_ds, _, status_s = _distance_f84(S, sdts / S, sdtv / S, fb_s)
        dN, se_dn, _, status_n = _distance_f84(N, ndts / N, ndtv / N, fb_n)
        saturated = (
            dS is None or dN is None or status_s >= 3 or status_n >= 3
            or dS >= 99.0 or dN >= 99.0
        )
    if not saturated and dS is not None and dN is not None:
        if dS > 1e-9:
            omega = dN / dS
        if (se_ds is not None and se_dn is not None
                and (se_ds > 0 or se_dn > 0)):
            variance = se_ds * se_ds + se_dn * se_dn
            if variance > 0:
                z = (dN - dS) / math.sqrt(variance)

    if z is not None:
        p_two = 2.0 * _normal_sf(abs(z))
        p_positive = _normal_sf(z)
        p_negative = 1.0 - _normal_sf(z)

    if saturated:
        conclusion = "saturated/not estimable"
        omega = None
        z = p_two = p_positive = p_negative = None
    elif p_two is not None and p_two < 0.05:
        if dN > dS:
            conclusion = "positive selection"
        elif dN < dS:
            conclusion = "purifying selection"
        else:
            conclusion = "not significant"
    else:
        conclusion = "not significant"

    return PairwiseSelection(
        seq1=name1 or "seq1",
        seq2=name2 or "seq2",
        sites=len(seq1) // 3,
        synonymous_sites=S,
        nonsynonymous_sites=N,
        synonymous_differences=sd,
        nonsynonymous_differences=nd,
        pS=sd / S if S > 0 else float("nan"),
        pN=nd / N if N > 0 else float("nan"),
        dS=dS,
        dN=dN,
        omega=omega,
        z=z,
        p_two_sided=p_two,
        p_positive=p_positive,
        p_negative=p_negative,
        conclusion=conclusion,
        alignment=alignment,
        se_dS=se_ds,
        se_dN=se_dn,
        method="YN00",
    )


def analyze_alignment_yn00(records: Iterable[SeqRecord], code_id: int = 1,
                           alignment: str = "") -> AlignmentSummary:
    records = list(records)
    if not records:
        return AlignmentSummary(alignment, 0, 0, 0, 0, 0)
    table = get_codon_table(code_id)
    sequences = [_clean_sequence(record.seq) for record in records]
    names = [record.id or "seq%d" % (i + 1) for i, record in enumerate(records)]
    lengths = {len(sequence) for sequence in sequences}
    if len(lengths) != 1:
        raise ValueError("Sequences are not aligned to the same length")
    full_length = lengths.pop()
    if full_length % 3 != 0:
        raise ValueError("Alignment length is not a multiple of 3")

    pairs = list(itertools.combinations(range(len(sequences)), 2))
    kappas = []
    for i, j in pairs:
        trimmed = _trim_pair_yn00(sequences[i], sequences[j], table)
        if trimmed is not None:
            kappas.append(_estimate_kappa_pair(trimmed[0], trimmed[1], table,
                                               code_id=code_id))
    kappa = sum(kappas) / len(kappas) if kappas else (10.0 if code_id == 2 else 2.0)

    results = []
    for i, j in pairs:
        result = pairwise_yn00(
            sequences[i], sequences[j], table, kappa=kappa,
            alignment=alignment, name1=names[i], name2=names[j], code_id=code_id,
        )
        if result is not None:
            results.append(result)

    skipped_stop = 0
    skipped_ambiguous = 0
    usable = 0
    for offset in range(0, full_length, 3):
        codons = [sequence[offset : offset + 3] for sequence in sequences]
        if any(not _is_clean_codon(codon) for codon in codons):
            skipped_ambiguous += 1
        elif any(table.forward_table.get(codon) is None for codon in codons):
            skipped_stop += 1
        else:
            usable += 1

    return AlignmentSummary(
        alignment=alignment or "alignment",
        sequences=len(records),
        codons=full_length // 3,
        usable_codons=usable,
        skipped_stop=skipped_stop,
        skipped_ambiguous=skipped_ambiguous,
        results=results,
        names=names,
        method="YN00",
    )


def analyze_file_yn00(path: str, code_id: int = 1) -> AlignmentSummary:
    path = Path(path)
    return analyze_alignment_yn00(read_fasta(str(path)), code_id=code_id,
                                  alignment=path.name)


def read_fasta(path: str) -> List[SeqRecord]:
    from Bio import SeqIO

    return list(SeqIO.parse(str(path), "fasta"))


def analyze_file(path: str, code_id: int = 1) -> AlignmentSummary:
    path = Path(path)
    return analyze_alignment(read_fasta(str(path)), code_id=code_id,
                             alignment=path.name)


_CONCLUSION_ORDER = {
    "positive selection": 0,
    "purifying selection": 1,
    "not significant": 2,
    "saturated/not estimable": 3,
    "not estimable": 4,
}


def selection_sort_key(result: PairwiseSelection):
    """Sort by conclusion, then by |Z| (strongest evidence first)."""
    rank = _CONCLUSION_ORDER.get(result.conclusion, 99)
    strength = -abs(result.z) if result.z is not None else 0.0
    return rank, strength, result.seq1, result.seq2


def format_results_csv(results: Sequence[PairwiseSelection]) -> str:
    columns = [
        "alignment", "seq1", "seq2", "method", "sites", "S", "N", "Sd", "Nd",
        "pS", "pN", "dS", "SE_dS", "dN", "SE_dN", "omega", "Z", "p_two_sided",
        "p_positive", "p_negative", "conclusion",
    ]
    lines = ["\t".join(columns)]
    for result in results:
        values = [
            result.alignment, result.seq1, result.seq2, result.method,
            str(result.sites),
            "%.3f" % result.synonymous_sites,
            "%.3f" % result.nonsynonymous_sites,
            "%.3f" % result.synonymous_differences,
            "%.3f" % result.nonsynonymous_differences,
            _fmt(result.pS), _fmt(result.pN), _fmt(result.dS), _fmt(result.se_dS),
            _fmt(result.dN), _fmt(result.se_dN), _fmt(result.omega),
            _fmt(result.z), _fmt(result.p_two_sided),
            _fmt(result.p_positive), _fmt(result.p_negative), result.conclusion,
        ]
        lines.append("\t".join(values))
    return "\n".join(lines) + "\n"


def _fmt(value) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float) and math.isnan(value):
        return "NA"
    if isinstance(value, float):
        if value == 0.0:
            value = 0.0  # normalise -0.0
        return "%.6g" % value
    return str(value)


_OMEGA_RE = re.compile(
    r"Nonsynonymous/synonymous ratio \(omega\):\s*([0-9.]+)", re.I
)
_LOG_RE = re.compile(r"(?:Optimal log-likelihood|BEST SCORE FOUND)\s*:?\s*(-?[0-9.]+)", re.I)


def parse_iqtree_omega(text: str) -> Tuple[Optional[float], Optional[float]]:
    """Return (omega, log_likelihood) from IQ-TREE stdout/report."""
    omega = None
    logl = None
    match = _OMEGA_RE.search(text)
    if match:
        try:
            omega = float(match.group(1))
        except ValueError:
            omega = None
    match = _LOG_RE.search(text)
    if match:
        try:
            logl = float(match.group(1))
        except ValueError:
            logl = None
    return omega, logl
