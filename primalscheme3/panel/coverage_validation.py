"""Versioned constraints, boolean search oracle, and fresh final validation.

Search uses candidate_valid(id) and conflict(a,b); reasons and witnesses are
separate diagnostic APIs. Final validation never reads an oracle instance.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from itertools import combinations
from types import SimpleNamespace

from primalschemers import do_pool_interact

from primalscheme3.core.config import Config
from primalscheme3.core.thermo import THERMO_RESULT, thermo_check
from primalscheme3.panel.coverage_specificity import SpecificityChecker
from primalscheme3.panel.coverage_types import Catalog

_THERMO_FIELDS = (
    "primer_size_min",
    "primer_size_max",
    "primer_gc_min",
    "primer_gc_max",
    "primer_tm_min",
    "primer_tm_max",
    "primer_hairpin_th_max",
    "primer_homopolymer_max",
    "use_annealing",
    "primer_annealing_prop",
    "primer_annealing_tempc",
    "mv_conc",
    "dv_conc",
    "dntp_conc",
    "dna_conc",
    "dimer_score",
)


@dataclass(frozen=True)
class ConstraintProfile:
    """Immutable resolved scientific settings, including caller-owned caps.

    Construct with from_config after the CLI resolves coverage defaults. Zero
    specificity distance and disabled MatchDB are errors, never fallbacks.
    """

    amplicon_size: int
    amplicon_size_min: int
    amplicon_size_max: int
    n_pools: int
    mismatch_product_size: int
    mismatch_kmersize: int
    mismatch_fuzzy: bool
    thermochemistry_json: str
    max_amplicons: int | None = None
    max_amplicons_msa: int | None = None

    def __post_init__(self):
        for field in (
            "amplicon_size",
            "amplicon_size_min",
            "amplicon_size_max",
            "n_pools",
            "mismatch_product_size",
            "mismatch_kmersize",
        ):
            value = getattr(self, field)
            if type(value) is not int or value <= 0:
                raise ValueError(f"{field} must be a positive integer")
        if not self.amplicon_size_min <= self.amplicon_size <= self.amplicon_size_max:
            raise ValueError("amplicon bounds must contain the nominal size")
        for field in ("max_amplicons", "max_amplicons_msa"):
            value = getattr(self, field)
            if value is not None and (type(value) is not int or value < 0):
                raise ValueError(f"{field} must be a nonnegative integer or None")
        if type(self.mismatch_fuzzy) is not bool:
            raise ValueError("mismatch_fuzzy must be boolean")
        chemistry = json.loads(self.thermochemistry_json)
        if set(chemistry) != set(_THERMO_FIELDS):
            raise ValueError(
                "thermochemistry must contain all resolved kernel settings"
            )
        if not 0 < chemistry["primer_size_min"] <= chemistry["primer_size_max"]:
            raise ValueError("invalid oligo length bounds")
        if any(
            isinstance(v, int | float) and not math.isfinite(v)
            for v in chemistry.values()
        ):
            raise ValueError("thermochemistry must be finite")

    @classmethod
    def from_config(cls, config: Config, *, max_amplicons=None, max_amplicons_msa=None):
        if not config.use_matchdb:
            raise ValueError("panel-v1 requires supplied-MSA specificity (use_matchdb)")
        if config.editdist_max != 1:
            raise ValueError(
                "panel-v1 supports the existing single-mismatch policy only"
            )
        return cls(
            amplicon_size=config.amplicon_size,
            amplicon_size_min=config.amplicon_size_min,
            amplicon_size_max=config.amplicon_size_max,
            n_pools=config.n_pools,
            mismatch_product_size=config.mismatch_product_size,
            mismatch_kmersize=config.mismatch_kmersize,
            mismatch_fuzzy=config.mismatch_fuzzy,
            thermochemistry_json=json.dumps(
                {k: getattr(config, k) for k in _THERMO_FIELDS},
                sort_keys=True,
                separators=(",", ":"),
            ),
            max_amplicons=max_amplicons,
            max_amplicons_msa=max_amplicons_msa,
        )

    def to_dict(self) -> dict:
        result = asdict(self)
        result["thermochemistry"] = json.loads(result.pop("thermochemistry_json"))
        result["thermochemistry"]["effective_tm_upper_offset"] = 2
        result.update(
            {
                "name": "panel-v1",
                "specificity_revision": "intended-sites-v1",
                "allowed_secondary_product_policy": "ordered-disjoint-intended-sites",
                "reference_coordinates": "zero-based-half-open",
                "product_coordinates": "ungapped-row-full-footprint-half-open",
                "product_bound": "0 < length <= mismatch_product_size; plus_end <= minus_start",
                "mismatch_policy": "exact-or-single-substitution"
                if self.mismatch_fuzzy
                else "exact",
                "ambiguity_policy": "non-N-IUPAC-compatible-uncertain; N/missing-unknown-no-support",
                "specificity_scope": "all-supplied-MSA-rows-and-occurrences",
            }
        )
        return result

    @property
    def cache_key(self) -> str:
        return hashlib.sha256(
            json.dumps(self.to_dict(), sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()


class _Evaluator:
    def __init__(self, catalog, profile, *, cache):
        self.catalog, self.profile = catalog, profile
        self.specificity = SpecificityChecker(catalog, profile, cache=cache)
        self.chemistry = SimpleNamespace(**json.loads(profile.thermochemistry_json))
        self.counters = defaultdict(
            int,
            {
                key: 0
                for key in (
                    "candidate_evaluations",
                    "pair_evaluations",
                    "thermo_calls",
                    "interaction_calls",
                )
            },
        )
        self.cache = cache
        self._thermo = {}
        self._interactions = {}

    def _oligo_reason(self, oligo):
        if self.cache and oligo in self._thermo:
            return self._thermo[oligo]
        if not oligo or any(b not in "ACGT" for b in oligo):
            reason = "nonconcrete_oligo"
        else:
            self.counters["thermo_calls"] += 1
            result = thermo_check(oligo, self.chemistry)
            reason = None if result == THERMO_RESULT.PASS else f"thermo:{result.name}"
        if self.cache:
            self._thermo[oligo] = reason
        return reason

    def _interacts(self, a, b):
        key = tuple(sorted((tuple(sorted(set(a))), tuple(sorted(set(b))))))
        if self.cache and key in self._interactions:
            return self._interactions[key]
        if any(
            not seq or any(base not in "ACGT" for base in seq)
            for group in key
            for seq in group
        ):
            return False  # Intrinsic nonconcrete_oligo rejects these before selection.
        self.counters["interaction_calls"] += 1
        result = bool(
            do_pool_interact(
                [s.encode() for s in key[0]],
                [s.encode() for s in key[1]],
                self.chemistry.dimer_score,
            )
        )
        if self.cache:
            self._interactions[key] = result
        return result

    def candidate(self, candidate_id):
        self.counters["candidate_evaluations"] += 1
        candidate = self.catalog.candidate_by_id.get(candidate_id)
        if candidate is None:
            return {
                "reasons": ["unknown_candidate"],
                "support": {},
                "rejected_products": [],
                "allowed_secondary_products": [],
            }
        reasons = []
        target = self.catalog.target_by_id.get(candidate.target_id)
        start, end = candidate.full_interval
        fe, rs = candidate.interior_interval
        oligos = candidate.forward_oligos + candidate.reverse_oligos
        if (
            target is None
            or not 0 <= start < fe <= rs < end <= target.reference_length
            or not candidate.forward_oligos
            or not candidate.reverse_oligos
        ):
            reasons.append("invalid_geometry")
        elif start != fe - max(map(len, candidate.forward_oligos)) or end != rs + max(
            map(len, candidate.reverse_oligos)
        ):
            reasons.append("inconsistent_oligo_envelope")
        if (
            not self.profile.amplicon_size_min
            <= end - start
            <= self.profile.amplicon_size_max
        ):
            reasons.append("amplicon_size_bounds")
        for oligo in sorted(set(oligos)):
            if len(oligo) < self.chemistry.primer_size_min:
                reasons.append("oligo_too_short")
            if len(oligo) < self.profile.mismatch_kmersize:
                reasons.append("oligo_shorter_than_terminal_kmer")
            reason = self._oligo_reason(oligo)
            if reason:
                reasons.append(reason)
        if self._interacts(oligos, oligos):
            reasons.append("oligo_interaction")
        specificity = self.specificity.intrinsic(candidate)
        if not specificity.support["joint_rows"]:
            reasons.append("no_joint_support")
        if specificity.rejected_products:
            reasons.append("specificity")
        return {
            "reasons": sorted(set(reasons)),
            "support": specificity.support,
            "rejected_products": list(specificity.rejected_products),
            "allowed_secondary_products": list(specificity.allowed_secondary_products),
        }

    def pair(self, a_id, b_id):
        self.counters["pair_evaluations"] += 1
        a, b = (self.catalog.candidate_by_id.get(i) for i in sorted((a_id, b_id)))
        if a is None or b is None:
            return {
                "reasons": ["unknown_candidate"],
                "rejected_products": [],
                "allowed_secondary_products": [],
            }
        reasons = []
        if a.id == b.id:
            reasons.append("duplicate_candidate")
        if a.target_id == b.target_id and max(
            a.full_interval[0], b.full_interval[0]
        ) < min(a.full_interval[1], b.full_interval[1]):
            reasons.append("same_target_overlap")
        if self._interacts(
            a.forward_oligos + a.reverse_oligos, b.forward_oligos + b.reverse_oligos
        ):
            reasons.append("oligo_interaction")
        specificity = self.specificity.pair(a, b)
        if specificity.rejected_products:
            reasons.append("specificity")
        return {
            "reasons": reasons,
            "rejected_products": list(specificity.rejected_products),
            "allowed_secondary_products": list(specificity.allowed_secondary_products),
        }


class CompatibilityOracle:
    """Pure memoized predicates scoped to one immutable catalogue and profile.

    conflict is pair-only; callers must separately require candidate_valid for
    both members. Count/pool limits apply to assignments, not these predicates.
    No verdict depends on selected, removed, or other catalogue candidates.
    """

    def __init__(self, catalog: Catalog, profile: ConstraintProfile):
        self.catalog, self.profile = catalog, profile
        self.cache_identity = (catalog.semantic_digest, profile.cache_key)
        self._evaluator = _Evaluator(catalog, profile, cache=True)
        self._candidate_results, self._pair_results = {}, {}

    def candidate_diagnostics(self, candidate_id) -> dict:
        if candidate_id not in self._candidate_results:
            self._candidate_results[candidate_id] = self._evaluator.candidate(
                candidate_id
            )
        # Do not expose mutable cache entries to diagnostic consumers.
        return json.loads(json.dumps(self._candidate_results[candidate_id]))

    def pair_diagnostics(self, a, b) -> dict:
        key = tuple(sorted((a, b)))
        if key not in self._pair_results:
            self._pair_results[key] = self._evaluator.pair(*key)
        return json.loads(json.dumps(self._pair_results[key]))

    def candidate_valid(self, candidate_id) -> bool:
        if candidate_id not in self._candidate_results:
            self._candidate_results[candidate_id] = self._evaluator.candidate(
                candidate_id
            )
        return not self._candidate_results[candidate_id]["reasons"]

    def conflict(self, a, b) -> bool:
        key = tuple(sorted((a, b)))
        if key not in self._pair_results:
            self._pair_results[key] = self._evaluator.pair(*key)
        return bool(self._pair_results[key]["reasons"])

    def candidate_reasons(self, candidate_id) -> tuple[str, ...]:
        self.candidate_valid(candidate_id)
        return tuple(self._candidate_results[candidate_id]["reasons"])

    def conflict_reasons(self, a, b) -> tuple[str, ...]:
        self.conflict(a, b)
        return tuple(self._pair_results[tuple(sorted((a, b)))]["reasons"])

    @property
    def counters(self) -> dict:
        return dict(self._evaluator.counters) | {
            f"specificity_{k}": v
            for k, v in self._evaluator.specificity.counters.items()
        }


def _coverage(intervals, reference_length):
    merged = []
    for start, end in sorted(intervals):
        start, end = max(0, start), min(reference_length, end)
        if start >= end:
            continue
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(end, merged[-1][1])
        else:
            merged.append([start, end])
    covered = sum(end - start for start, end in merged)
    return {
        "intervals": merged,
        "covered_bases": covered,
        "fraction": covered / reference_length if reference_length else 0.0,
    }


def validate_assignments(
    catalog: Catalog,
    assignments,
    profile: ConstraintProfile,
    metric: str,
    coverage_target: float,
) -> dict:
    """Fresh uncached I/P, reconstructed support, and fresh interval unions.

    Coverage target is an objective: falling short is diagnostic, not an invalid
    panel. All failures remain explicit; empty assignments are feasible.
    """
    if metric not in ("full-span", "primer-trimmed"):
        raise ValueError("unknown coverage metric")
    if not math.isfinite(coverage_target) or not 0 <= coverage_target <= 1:
        raise ValueError("coverage_target must be between zero and one")
    assignments = tuple(assignments)
    evaluator = _Evaluator(catalog, profile, cache=False)
    violations, allowed = [], []
    support = {}
    counts = Counter(a.candidate_id for a in assignments)
    for candidate_id, count in sorted(counts.items()):
        if count > 1:
            violations.append(
                {"reason": "duplicate_candidate", "candidate_id": candidate_id}
            )
    if profile.max_amplicons is not None and len(assignments) > profile.max_amplicons:
        violations.append({"reason": "max_amplicons"})
    target_counts = Counter()
    selected = []
    for assignment in assignments:
        if (
            type(assignment.pool) is not int
            or not 0 <= assignment.pool < profile.n_pools
        ):
            violations.append(
                {
                    "reason": "pool_bounds",
                    "candidate_id": assignment.candidate_id,
                    "pool": assignment.pool,
                }
            )
        result = evaluator.candidate(assignment.candidate_id)
        support[assignment.candidate_id] = result["support"]
        for reason in result["reasons"]:
            violations.append(
                {
                    "reason": reason,
                    "candidate_id": assignment.candidate_id,
                    "rejected_products": result["rejected_products"]
                    if reason == "specificity"
                    else [],
                }
            )
        candidate = catalog.candidate_by_id.get(assignment.candidate_id)
        if candidate is not None:
            selected.append((assignment, candidate))
            target_counts[candidate.target_id] += 1
    for target_id, count in sorted(target_counts.items()):
        if profile.max_amplicons_msa is not None and count > profile.max_amplicons_msa:
            violations.append({"reason": "max_amplicons_msa", "target_id": target_id})
    for (aa, a), (ab, b) in combinations(selected, 2):
        if aa.pool != ab.pool:
            continue
        result = evaluator.pair(a.id, b.id)
        for reason in result["reasons"]:
            violations.append(
                {
                    "reason": reason,
                    "candidate_ids": sorted((a.id, b.id)),
                    "pool": aa.pool,
                    "rejected_products": result["rejected_products"]
                    if reason == "specificity"
                    else [],
                }
            )
        allowed.extend(
            dict(p, pool=aa.pool) for p in result["allowed_secondary_products"]
        )
    per_target = {}
    for target in catalog.targets:
        candidates = [c for _, c in selected if c.target_id == target.id]
        full = _coverage([c.full_interval for c in candidates], target.reference_length)
        interior = _coverage(
            [c.interior_interval for c in candidates], target.reference_length
        )
        chosen = full if metric == "full-span" else interior
        per_target[target.id] = {
            "reference_length": target.reference_length,
            "full": full,
            "interior": interior,
            "selected_count": len(candidates),
            "coverage_fraction": chosen["fraction"],
            "shortfall": max(0.0, coverage_target - chosen["fraction"]),
            "normalized_shortfall": max(0.0, coverage_target - chosen["fraction"])
            / coverage_target
            if coverage_target
            else 0.0,
            "target_met": chosen["fraction"] >= coverage_target,
        }
    return {
        "schemaVersion": "primalscheme3.panel-validation/v1",
        "valid": not violations,
        "metric": metric,
        "coverage_target": coverage_target,
        "profile": profile.to_dict(),
        "catalog_semantic_digest": catalog.semantic_digest,
        "assignments": [asdict(a) for a in assignments],
        "violations": violations,
        "per_target": per_target,
        "support_diagnostics": support,
        "allowed_secondary_products": allowed,
        "counters": dict(evaluator.counters)
        | {f"specificity_{k}": v for k, v in evaluator.specificity.counters.items()},
        "limitations": [
            "Specificity evaluates supplied MSA rows only, not whole-genome background.",
            "N and missing terminal windows are unknown and do not establish hits or support.",
            "Non-N IUPAC-compatible hits are uncertain and cannot authorize intended exemptions.",
            "Terminal-only hits predict possible products, not laboratory amplification.",
            "Ordered disjoint intended-site secondary products are tolerated and never add coverage.",
        ],
    }
