"""Row-aware supplied-MSA specificity, revision intended-sites-v1.

Terminal matches predict possible products; they do not prove full binding.
Only pair-local, fully observed anchored signatures authorize exemptions.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from itertools import product
from typing import TYPE_CHECKING

from primalscheme3.core.config import ALL_DNA
from primalscheme3.core.seq_functions import reverse_complement
from primalscheme3.panel.coverage_catalog import _anchored_matches
from primalscheme3.panel.coverage_types import Candidate, Catalog

if TYPE_CHECKING:
    from primalscheme3.panel.coverage_validation import ConstraintProfile


@dataclass(frozen=True, order=True)
class Hit:
    target_id: str
    row_id: str
    oligo: str
    orientation: str
    start: int
    end: int
    terminal_start: int
    terminal_end: int
    mismatches: int
    classification: str

    @property
    def signature(self) -> tuple:
        return (
            self.target_id,
            self.row_id,
            self.oligo,
            self.orientation,
            self.start,
            self.end,
        )


@dataclass(frozen=True)
class SpecificityResult:
    rejected_products: tuple[dict, ...]
    allowed_secondary_products: tuple[dict, ...]
    support: dict


class SpecificityChecker:
    """Lazy complete per-oligo hit sets; use cache=False for fresh validation.

    Non-N IUPAC uses compatible sets without combinatorial expansion. N and
    missing terminal observations never supply matches or intended support.
    Ambiguous terminal hits can predict conflicts but cannot grant exemptions.
    """

    def __init__(self, catalog: Catalog, profile: ConstraintProfile, *, cache=True):
        self.catalog = catalog
        self.profile = profile
        self.cache = cache
        self._hits = {}
        self._support = {}
        self.counters = defaultdict(int)
        self._rows = []
        self._terminal_index = None
        self._ambiguous_windows = []
        for target in catalog.targets:
            for row_id, row in zip(target.row_ids, target.rows, strict=True):
                bases = []
                missing_boundaries = set()
                for cell in row:
                    if cell == "":
                        missing_boundaries.add(len(bases))
                    elif cell != "-":
                        bases.append(cell.upper())
                self._rows.append(
                    (target.id, row_id, "".join(bases), missing_boundaries)
                )

    def _index_terminals(self):
        if self._terminal_index is not None:
            return
        self._terminal_index = defaultdict(list)
        k = self.profile.mismatch_kmersize
        for row_index, (_, _, sequence, missing) in enumerate(self._rows):
            for pos in range(len(sequence) - k + 1):
                self.counters["row_windows_indexed"] += 1
                window = sequence[pos : pos + k]
                if any(b not in ALL_DNA for b in window) or any(
                    pos < m < pos + k for m in missing
                ):
                    self.counters["unknown_terminal_windows"] += 1
                elif any(b not in "ACGT" for b in window):
                    self._ambiguous_windows.append((window, row_index, pos))
                else:
                    self._terminal_index[window].append((row_index, pos))

    def hits(self, oligo: str) -> tuple[Hit, ...]:
        if self.cache and oligo in self._hits:
            self.counters["hit_cache_hits"] += 1
            return self._hits[oligo]
        self.counters["hit_evaluations"] += 1
        k = self.profile.mismatch_kmersize
        if len(oligo) < k or any(b not in "ACGT" for b in oligo):
            return ()
        self._index_terminals()
        terminal = oligo[-k:]
        queries = (("+", terminal), ("-", reverse_complement(terminal)))
        hits = []
        allowance = int(self.profile.mismatch_fuzzy)
        for orientation, query in queries:
            # At most 3*k substitutions; never expand ambiguous row bases.
            variants = {query: 0}
            if allowance:
                for i, base in enumerate(query):
                    for alternate in "ACGT":
                        if alternate != base:
                            variants[query[:i] + alternate + query[i + 1 :]] = 1
            observations = []
            for variant, mismatches in variants.items():
                observations.extend(
                    (ri, pos, mismatches, False)
                    for ri, pos in self._terminal_index.get(variant, ())
                )
            for window, ri, pos in self._ambiguous_windows:
                mismatches = sum(
                    base not in ALL_DNA[observed]
                    for base, observed in zip(query, window, strict=True)
                )
                if mismatches <= allowance:
                    observations.append((ri, pos, mismatches, True))
            for ri, pos, mismatches, ambiguous in observations:
                target_id, row_id, sequence, missing = self._rows[ri]
                start, end = (
                    (pos + k - len(oligo), pos + k)
                    if orientation == "+"
                    else (pos, pos + len(oligo))
                )
                unavailable = (
                    start < 0
                    or end > len(sequence)
                    or any(start < m < end for m in missing)
                )
                classification = (
                    "uncertain-footprint"
                    if unavailable
                    else "ambiguous"
                    if ambiguous
                    else "single-mismatch"
                    if mismatches
                    else "exact-terminal"
                )
                hits.append(
                    Hit(
                        target_id,
                        row_id,
                        oligo,
                        orientation,
                        start,
                        end,
                        pos,
                        pos + k,
                        mismatches,
                        classification,
                    )
                )
        result = tuple(sorted(hits))
        if self.cache:
            self._hits[oligo] = result
        return result

    def intended(self, candidate: Candidate) -> tuple[frozenset, dict]:
        """Reconstruct all exact concrete pair signatures from anchor cells."""
        if self.cache and candidate.id in self._support:
            return self._support[candidate.id]
        self.counters["support_evaluations"] += 1
        target = self.catalog.target_by_id.get(candidate.target_id)
        signatures = set()
        joint, unknown, spans = [], [], []
        if target is not None:
            for row_id, row in zip(target.row_ids, target.rows, strict=True):
                forward, fu = _anchored_matches(
                    target,
                    row,
                    candidate.forward_oligos,
                    candidate.interior_interval[0],
                    reverse=False,
                )
                reverse, ru = _anchored_matches(
                    target,
                    row,
                    candidate.reverse_oligos,
                    candidate.interior_interval[1],
                    reverse=True,
                )
                if forward and reverse:
                    joint.append(row_id)
                    row_spans = set()
                    for f, r in product(forward, reverse):
                        if f[2] <= r[1]:
                            signatures.add(
                                (
                                    (target.id, row_id, f[0], "+", f[1], f[2]),
                                    (target.id, row_id, r[0], "-", r[1], r[2]),
                                )
                            )
                            row_spans.add((f[1], r[2]))
                    spans.extend(
                        [row_id, start, end] for start, end in sorted(row_spans)
                    )
                if fu or ru:
                    unknown.append(row_id)
        result = (
            frozenset(signatures),
            {"joint_rows": joint, "unknown_rows": unknown, "row_product_spans": spans},
        )
        if self.cache:
            self._support[candidate.id] = result
        return result

    def intrinsic(self, candidate: Candidate) -> SpecificityResult:
        signatures, support = self.intended(candidate)
        return self._evaluate((candidate,), signatures, frozenset(), support)

    def pair(self, a: Candidate, b: Candidate) -> SpecificityResult:
        a, b = sorted((a, b), key=lambda c: c.id)
        sa, _ = self.intended(a)
        sb, _ = self.intended(b)
        secondary = set()
        if a.target_id == b.target_id:
            for left, right in ((sa, sb), (sb, sa)):
                for lp, rp in product(left, right):
                    if lp[0][:2] == rp[0][:2] and lp[1][5] <= rp[0][4]:
                        secondary.add((lp[0], rp[1]))
        return self._evaluate((a, b), sa | sb, frozenset(secondary), {})

    def _evaluate(self, candidates, declared, secondary, support):
        groups = []
        for candidate in candidates:
            grouped = defaultdict(lambda: {"+": set(), "-": set()})
            for oligo in sorted(
                set(candidate.forward_oligos + candidate.reverse_oligos)
            ):
                for hit in self.hits(oligo):
                    grouped[(hit.target_id, hit.row_id)][hit.orientation].add(hit)
            groups.append(grouped)
        physical = set()
        left, right = groups if len(groups) == 2 else (groups[0], groups[0])
        for row in sorted(left.keys() & right.keys()):
            directions = ((left[row]["+"], right[row]["-"]),)
            if len(groups) == 2:
                directions += ((right[row]["+"], left[row]["-"]),)
            for plus, minus in directions:
                for f, r in product(plus, minus):
                    if (
                        f.end <= r.start
                        and 0 < r.end - f.start <= self.profile.mismatch_product_size
                    ):
                        physical.add((f, r))
        rejected, allowed = [], []
        for f, r in sorted(physical):
            self.counters["products_evaluated"] += 1
            signature = (f.signature, r.signature)
            if signature in declared:
                continue
            tolerated = signature in secondary
            witness = {
                "target_id": f.target_id,
                "row_id": f.row_id,
                "owners": [c.id for c in candidates],
                "start": f.start,
                "end": r.end,
                "length": r.end - f.start,
                "plus": self._hit_record(f, candidates),
                "minus": self._hit_record(r, candidates),
                "reason": (
                    "ordered-disjoint-intended-sites"
                    if tolerated
                    else "unexpected-product"
                ),
            }
            (allowed if tolerated else rejected).append(witness)
        return SpecificityResult(tuple(rejected), tuple(allowed), support)

    @staticmethod
    def _hit_record(hit, candidates):
        return {
            "oligo": hit.oligo,
            "orientation": hit.orientation,
            "start": hit.start,
            "end": hit.end,
            "terminal_interval": [hit.terminal_start, hit.terminal_end],
            "mismatches": hit.mismatches,
            "classification": hit.classification,
            "owners": [
                c.id
                for c in candidates
                if hit.oligo in c.forward_oligos + c.reverse_oligos
            ],
        }
