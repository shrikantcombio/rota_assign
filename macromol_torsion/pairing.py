"""
Bound/Unbound Residue Pairing Module
------------------------------------
Pairs residues of a bound structure with their counterparts in an unbound structure.

Equal-residue files (e.g. 1ASY_1EOV_C.pdb / 1ASY_1EOV_U.pdb) list the same residues in the
same order, but the unbound chains usually carry different chain IDs (5X6B: IJEF vs ABCD) and
may be renumbered. Matching residue labels ("LYS A   72") therefore drops or mis-pairs residues.
Here chains are paired by order of appearance, and residues within each chain pair by aligning
their residue-name sequences (identical to positional pairing for equal-residue files); only
residues with identical names are paired.
"""

from collections import OrderedDict
from difflib import SequenceMatcher
from typing import Dict, List, Tuple


def _split_label(label: str) -> Tuple[str, str]:
    """('LYS', 'A') from a residue label such as 'LYS A   72'."""
    parts = label.split()
    return parts[0], parts[1] if len(parts) > 1 else ""


def _by_chain(labels: List[str]) -> "OrderedDict[str, List[str]]":
    chains: "OrderedDict[str, List[str]]" = OrderedDict()
    for lbl in labels:
        chains.setdefault(_split_label(lbl)[1], []).append(lbl)
    return chains


def pair_equal_residues(bound_labels: List[str], unbound_labels: List[str]) -> Tuple[Dict[str, str], Dict[str, int]]:
    """
    Returns ({bound_label: unbound_label}, stats).

    stats keys: 'paired', 'name_mismatch', 'unpaired_bound', 'chain_length_mismatch'.
    Raises ValueError when the two structures have different numbers of chains, since the
    chain correspondence is then ambiguous.
    """
    b_chains, u_chains = _by_chain(bound_labels), _by_chain(unbound_labels)
    if len(b_chains) != len(u_chains):
        raise ValueError(
            f"bound structure has {len(b_chains)} chains ({''.join(b_chains)}) but unbound has "
            f"{len(u_chains)} ({''.join(u_chains)}); chain correspondence is ambiguous")

    mapping: Dict[str, str] = {}
    stats = {"paired": 0, "name_mismatch": 0, "unpaired_bound": 0, "chain_length_mismatch": 0}
    for b_res, u_res in zip(b_chains.values(), u_chains.values()):
        if len(b_res) != len(u_res):
            stats["chain_length_mismatch"] += 1
        b_names = [_split_label(lbl)[0] for lbl in b_res]
        u_names = [_split_label(lbl)[0] for lbl in u_res]
        matcher = SequenceMatcher(a=b_names, b=u_names, autojunk=False)
        n_paired_chain = 0
        for blk in matcher.get_matching_blocks():
            for k in range(blk.size):
                mapping[b_res[blk.a + k]] = u_res[blk.b + k]
                n_paired_chain += 1
        stats["paired"] += n_paired_chain
        stats["name_mismatch"] += min(len(b_res), len(u_res)) - n_paired_chain
    stats["unpaired_bound"] = len(bound_labels) - stats["paired"]
    return mapping, stats


def _aligned_pairs(b_res: List[str], u_res: List[str], min_block: int) -> List[Tuple[str, str]]:
    b_names = [_split_label(lbl)[0] for lbl in b_res]
    u_names = [_split_label(lbl)[0] for lbl in u_res]
    matcher = SequenceMatcher(a=b_names, b=u_names, autojunk=False)
    pairs = []
    for blk in matcher.get_matching_blocks():
        if blk.size >= min_block:
            pairs.extend((b_res[blk.a + k], u_res[blk.b + k]) for k in range(blk.size))
    return pairs


def pair_residues_across_entries(bound_labels: List[str], unbound_labels: List[str],
                                 bound_chains: str, unbound_chains: str,
                                 min_block: int = 3, min_coverage: float = 0.5) -> Tuple[Dict[str, str], Dict[str, object]]:
    """
    Pairs residues of the given bound chains with an unbound structure from a different PDB entry
    (numbering and chain IDs generally differ).

    Each bound chain is matched to the unbound chain (from `unbound_chains`, or all chains if empty)
    whose residue-name sequence aligns best, preferring chains not yet used (an unbound chain is
    shared only by homo-oligomer copies); residues are paired within aligned blocks of at
    least `min_block` identical residues. A chain pair is rejected when the paired residues cover less
    than `min_coverage` of the shorter chain (the unbound partner may be a single domain).
    Returns ({bound_label: unbound_label}, stats).
    """
    b_chains, u_chains = _by_chain(bound_labels), _by_chain(unbound_labels)
    allowed_u = [c for c in (unbound_chains or "".join(u_chains)) if c in u_chains]
    mapping: Dict[str, str] = {}
    chain_map: Dict[str, str] = {}
    used = set()
    def acceptable(bc, uc, pairs):
        return len(pairs) >= min_coverage * min(len(b_chains[bc]), len(u_chains[uc]))

    for bc in [c for c in bound_chains if c in b_chains]:
        candidates = []
        for uc in allowed_u:
            pairs = _aligned_pairs(b_chains[bc], u_chains[uc], min_block)
            if acceptable(bc, uc, pairs):
                candidates.append((uc in used, -len(pairs), uc, pairs))
        if not candidates:
            continue
        # Prefer an unused unbound chain; reuse one only for homo-oligomers whose unbound entry
        # contains fewer copies than the bound complex (e.g. bound dimer DE vs unbound monomer A)
        _, _, uc, pairs = min(candidates, key=lambda c: c[:3])
        used.add(uc)
        chain_map[bc] = uc
        mapping.update(dict(pairs))
    n_bound = sum(len(b_chains[c]) for c in bound_chains if c in b_chains)
    stats = {"paired": len(mapping), "unpaired_bound": n_bound - len(mapping),
             "chain_map": "".join(f"{b}>{u}" for b, u in chain_map.items())}
    return mapping, stats
