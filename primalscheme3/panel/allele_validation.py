"""Selected-sites-v2 specificity, numerical pool constraints and fresh validation.

The supplied-row uncertainty scope intentionally does not invent terminal seeds
outside observed row coordinates. Seeded full footprints extending beyond those
coordinates are conservative uncertainty blocks, not confirmed binding evidence.
"""

from __future__ import annotations

import json
import math
from collections import Counter, OrderedDict, defaultdict
from copy import deepcopy
from dataclasses import asdict, dataclass, replace
from importlib.metadata import version
from itertools import combinations, combinations_with_replacement, product

from primalschemers import calc_at_offset_py

from primalscheme3.core import thermo
from primalscheme3.core.config import ALL_DNA_WITH_N
from primalscheme3.core.seq_functions import reverse_complement
from primalscheme3.panel.allele_coverage import (
    allele_summary,
    binding_support,
    canonical_observations,
    configuration_products,
)
from primalscheme3.panel.coverage_specificity import (
    Hit,
    SpecificityChecker,
    SpecificityResult,
)
from primalscheme3.panel.coverage_types import (
    ConfigurationLedger,
    semantic_id,
)
from primalscheme3.panel.coverage_variants import (
    ConfigurationIneligible,
    make_configuration,
)


@dataclass(frozen=True)
class AlleleConstraintProfile:
    amplicon_size_min: int
    amplicon_size_max: int
    n_pools: int
    mismatch_kmersize: int = 17
    mismatch_product_size: int = 2000
    mismatch_fuzzy: bool = True
    secondary_product_policy: str = "ordered-disjoint-intended-sites"
    max_amplicons: int | None = None
    max_amplicons_msa: int | None = None

    def __post_init__(self):
        for name in (
            "amplicon_size_min",
            "amplicon_size_max",
            "n_pools",
            "mismatch_kmersize",
            "mismatch_product_size",
        ):
            value = getattr(self, name)
            if type(value) is not int or value <= 0:
                raise ValueError(f"{name} must be a positive integer")
        if self.amplicon_size_min > self.amplicon_size_max:
            raise ValueError("invalid amplicon bounds")
        if self.mismatch_fuzzy is not True:
            raise ValueError("selected-sites-v2 requires at most one substitution")
        if self.secondary_product_policy not in (
            "ordered-disjoint-intended-sites",
            "reject-secondary-products/v1",
        ):
            raise ValueError("unknown secondary product policy")
        for name in ("max_amplicons", "max_amplicons_msa"):
            value = getattr(self, name)
            if value is not None and (type(value) is not int or value < 0):
                raise ValueError(f"{name} must be a nonnegative integer")

    @classmethod
    def from_config(cls, config, **kwargs):
        if (
            not config.use_matchdb
            or config.editdist_max != 1
            or not config.mismatch_fuzzy
        ):
            raise ValueError(
                "selected-sites-v2 requires supplied-row single-substitution specificity"
            )
        if config.dimer_score != -26:
            raise ValueError("strict dimer baseline must be -26")
        return cls(
            config.amplicon_size_min,
            config.amplicon_size_max,
            config.n_pools,
            config.mismatch_kmersize,
            config.mismatch_product_size,
            **kwargs,
        )

    def to_dict(self):
        return asdict(self) | {
            "name": "allele-panel-v1",
            "specificity_revision": "selected-sites-v2",
            "uncertainty_scope": "supplied-row-seeds-and-seeded-full-footprints/v1",
            "specificity_scope": "all-supplied-MSA-rows-and-occurrences",
            "mismatch_policy": "at-most-one-substitution-no-indels",
            "product_coordinates": "ungapped-row-full-footprint-half-open",
            "product_bound": "0 < length <= D; plus_end <= minus_start",
            "missing_padding_policy": "no-invented-terminal-seeds; seeded-unavailable-full-footprints-block",
            "chemistry_policy": "any-complete-declared-profile; effective-tm-upper-offset-2",
        }

    @property
    def cache_key(self):
        return semantic_id("allele-constraint-profile", self.to_dict())


@dataclass(frozen=True)
class StagePolicy:
    stage_id: str = "strict"
    active_cutoff: float = -26
    strict_cutoff: float = -26
    max_violating_edges: int = 0
    max_incident_species: int = 0

    def __post_init__(self):
        if (
            self.strict_cutoff != -26
            or not math.isfinite(self.active_cutoff)
            or self.active_cutoff > -26
        ):
            raise ValueError(
                "strict baseline must be -26; active cutoff must be finite and <= -26"
            )
        if not self.stage_id or (
            self.stage_id == "strict" and self.active_cutoff != -26
        ):
            raise ValueError("strict stage cannot relax the cutoff")
        for name in ("max_violating_edges", "max_incident_species"):
            if type(getattr(self, name)) is not int or getattr(self, name) < 0:
                raise ValueError("exposure limits must be nonnegative integers")

    def rejects(self, score):
        return score <= self.active_cutoff

    @property
    def cache_key(self):
        return semantic_id("allele-stage-policy", asdict(self))


DEFAULT_STAGE_POLICY = StagePolicy()


@dataclass(frozen=True)
class DimerScore:
    sequences: tuple[str, str]
    orientation_scores: tuple[float, float]
    orientation_offsets: tuple[int, int]
    kernel: str

    @property
    def score(self):
        return min(self.orientation_scores)

    def to_dict(self):
        return asdict(self) | {"score": self.score}


def numerical_dimer_score(a, b):
    """Exact installed boolean-kernel offset range, in both ordered directions."""
    if any(len(s) < 3 or any(c not in "ACGT" for c in s) for s in (a, b)):
        raise ValueError("native dimer score requires concrete oligos of length >= 3")
    ordered = tuple(sorted((a, b)))
    measured = []
    for x, y in (ordered, ordered[::-1]):
        measured.append(
            min(
                (float(calc_at_offset_py(x, y, offset)), offset)
                for offset in range(-(len(x) - 2), len(y) - len(x))
            )
        )
    return DimerScore(
        ordered,
        tuple(x[0] for x in measured),
        tuple(x[1] for x in measured),
        "primalschemers/"
        + version("primalschemers")
        + "/native-offsets-bidirectional/v1",
    )


@dataclass(frozen=True)
class _SelectedView:
    id: str
    target_id: str
    forward_oligos: tuple[str, ...]
    reverse_oligos: tuple[str, ...]
    selected_site_ids: tuple[str, ...]


class SelectedSiteSpecificityChecker(SpecificityChecker):
    """Pair-local complete certificates; compatible N/IUPAC never exemptions."""

    def _view(self, config):
        return _SelectedView(
            config.id,
            config.target_id,
            tuple(self.catalog.site_by_id[x].sequence for x in config.forward_site_ids),
            tuple(self.catalog.site_by_id[x].sequence for x in config.reverse_site_ids),
            config.forward_site_ids + config.reverse_site_ids,
        )

    def _index_terminals(self):
        if self._terminal_index is not None:
            return
        self._terminal_index = defaultdict(list)
        k = self.profile.mismatch_kmersize
        for ri, (_, _, sequence, missing) in enumerate(self._rows):
            for pos in range(len(sequence) - k + 1):
                self.counters["row_windows_indexed"] += 1
                window = sequence[pos : pos + k]
                if any(c not in "ACGT" for c in window) or any(
                    pos < m < pos + k for m in missing
                ):
                    self._ambiguous_windows.append((window, ri, pos))
                else:
                    self._terminal_index[window].append((ri, pos))

    def hits(self, oligo):
        if self.cache and oligo in self._hits:
            return self._hits[oligo]
        self.counters["hit_evaluations"] += 1
        k = self.profile.mismatch_kmersize
        if len(oligo) < k:
            return ()  # Oracle independently rejects unsupported oligo lengths.
        self._index_terminals()
        found = set()
        terminal = oligo[-k:]
        for orientation, query in (
            ("+", terminal),
            ("-", reverse_complement(terminal)),
        ):
            variants = {query: 0}
            for i, base in enumerate(query):
                for alt in "ACGT":
                    if alt != base:
                        variants[query[:i] + alt + query[i + 1 :]] = 1
            observations = []
            for word, mismatches in variants.items():
                observations.extend(
                    (ri, pos, mismatches, False)
                    for ri, pos in self._terminal_index.get(word, ())
                )
            for window, ri, pos in self._ambiguous_windows:
                mismatches = sum(
                    base not in ALL_DNA_WITH_N.get(cell, "ACGT")
                    for base, cell in zip(query, window, strict=True)
                )
                if mismatches <= 1:
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
                ambiguous = ambiguous or any(
                    c not in "ACGT"
                    for c in sequence[max(0, start) : min(len(sequence), end)]
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
                found.add(
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
        result = tuple(sorted(found))
        if self.cache:
            self._hits[oligo] = result
        return result

    def intended(self, config):
        if self.cache and config.id in self._support:
            return self._support[config.id]
        signatures, joint, unknown, spans = set(), set(), set(), set()
        products = configuration_products(self.catalog, config)
        for p in products:
            allele = self.catalog.observation_by_id[p.allele_id]
            f, r = (
                self.catalog.site_by_id[p.forward_site_id],
                self.catalog.site_by_id[p.reverse_site_id],
            )
            for row_id in allele.row_ids:
                signatures.add(
                    (
                        (
                            config.target_id,
                            row_id,
                            f.sequence,
                            "+",
                            *p.forward_footprint,
                        ),
                        (
                            config.target_id,
                            row_id,
                            r.sequence,
                            "-",
                            *p.reverse_footprint,
                        ),
                    )
                )
                joint.add(row_id)
                spans.add((row_id, p.forward_footprint[0], p.reverse_footprint[1]))
        for allele in self.catalog.observations:
            if allele.target_id == config.target_id and any(
                binding_support(self.catalog.site_by_id[s], allele).status == "unknown"
                for s in config.forward_site_ids + config.reverse_site_ids
            ):
                unknown.update(allele.row_ids)
        result = (
            frozenset(signatures),
            {
                "joint_rows": sorted(joint),
                "unknown_rows": sorted(unknown),
                "row_product_spans": sorted(spans),
            },
        )
        if self.cache:
            self._support[config.id] = result
        return result

    def intrinsic(self, config):
        declared, support = self.intended(config)
        return self._evaluate((self._view(config),), declared, frozenset(), support)

    def pair(self, a, b):
        a, b = sorted((a, b), key=lambda c: c.id)
        sa, _ = self.intended(a)
        sb, _ = self.intended(b)
        secondary = set()
        if (
            self.profile.secondary_product_policy == "ordered-disjoint-intended-sites"
            and a.target_id == b.target_id
        ):
            for left, right in ((sa, sb), (sb, sa)):
                for lp, rp in product(left, right):
                    if lp[0][:2] == rp[0][:2] and lp[1][5] <= rp[0][4]:
                        secondary.add((lp[0], rp[1]))
        return self._evaluate(
            (self._view(a), self._view(b)), sa | sb, frozenset(secondary), {}
        )

    def _evaluate(self, candidates, declared, secondary, support):
        result = super()._evaluate(candidates, declared, secondary, support)
        rejected, allowed = [], []

        def identified(witness):
            sequences = {witness[side]["oligo"] for side in ("plus", "minus")}
            ids = tuple(
                sorted(
                    {
                        sid
                        for c in candidates
                        for sid in c.selected_site_ids
                        if self.catalog.site_by_id[sid].sequence in sequences
                    }
                )
            )
            return witness | {
                "id": semantic_id("selected-site-product-witness", witness),
                "site_ids": ids,
            }

        for witness in result.rejected_products:
            witness = identified(witness)
            uncertain = any(
                witness[side]["classification"] in ("ambiguous", "uncertain-footprint")
                for side in ("plus", "minus")
            )
            rejected.append(
                witness
                | {
                    "uncertain": uncertain,
                    "reason": "uncertain-potential-product"
                    if uncertain
                    else witness["reason"],
                    "coverage_credit": 0,
                }
            )
        for witness in result.allowed_secondary_products:
            witness = identified(witness)
            allowed.append(witness | {"uncertain": False, "coverage_credit": 0})
        return SpecificityResult(tuple(rejected), tuple(allowed), support)


class AlleleCompatibilityOracle:
    """Memoized exact subset predicates; pool exposure is an aggregate guard."""

    def __init__(
        self,
        catalog,
        ledger,
        profile,
        policy=DEFAULT_STAGE_POLICY,
        *,
        cache=True,
        score_cache=None,
    ):
        if ledger.catalog_digest != catalog.semantic_digest:
            raise ValueError("ledger belongs to another catalog")
        self.catalog, self.profile, self.policy, self.cache = (
            catalog,
            profile,
            policy,
            cache,
        )
        self.configurations = dict(ledger.configuration_by_id)
        self.specificity = SelectedSiteSpecificityChecker(catalog, profile, cache=cache)
        self.score_cache = {} if score_cache is None else score_cache
        self._scores = self.score_cache
        self._score_kernel = (
            "primalschemers/"
            + version("primalschemers")
            + "/native-offsets-bidirectional/v1"
        )
        self._candidates, self._pairs = {}, {}
        self._pool_states = OrderedDict()
        self.cache_identity = (
            catalog.semantic_digest,
            profile.cache_key,
            policy.cache_key,
        )

    def register(self, *configurations):
        for c in configurations:
            if c.id in self.configurations and self.configurations[c.id] != c:
                raise ValueError("configuration identity collision")
            self.configurations[c.id] = c

    def _score(self, a, b):
        sequences = tuple(sorted((a, b)))
        key = (self._score_kernel, *sequences)
        if not self.cache or key not in self._scores:
            measured = numerical_dimer_score(*sequences)
            if self.cache:
                self._scores[key] = measured
            return measured
        return self._scores[key]

    def pool_diagnostics(self, configuration_ids):
        owners = defaultdict(set)
        for cid in sorted(set(configuration_ids)):
            c = self.configurations[cid]
            for sid in c.forward_site_ids + c.reverse_site_ids:
                owners[self.catalog.site_by_id[sid].sequence].add((cid, sid))
        edges, violating, incident, active = [], [], set(), []
        for a, b in combinations_with_replacement(sorted(owners), 2):
            measured = self._score(a, b)
            edge = measured.to_dict() | {
                "owners": {s: sorted(owners[s]) for s in sorted({a, b})}
            }
            edge["id"] = semantic_id("physical-dimer-edge", measured.to_dict())
            edge["site_ids"] = tuple(
                sorted({sid for seq in (a, b) for _, sid in owners[seq]})
            )
            edges.append(edge)
            if measured.score <= self.policy.strict_cutoff:
                violating.append(edge)
                incident.update((a, b))
            if self.policy.rejects(measured.score):
                active.append(edge)
        reasons = []
        if active:
            reasons.append("dimer-threshold")
        if (
            len(violating) > self.policy.max_violating_edges
            or len(incident) > self.policy.max_incident_species
        ):
            reasons.append("salvage-exposure-budget")
        return {
            "reasons": reasons,
            "dimer_edges": edges,
            "active_violations": active,
            "strict_violating_edges": violating,
            "violating_edge_count": len(violating),
            "incident_species": sorted(incident),
            "incident_species_count": len(incident),
        }

    def pool_valid(self, configuration_ids):
        # Exact incremental aggregate guard. No diagnostic dictionaries on hot fits.
        key = frozenset(configuration_ids)
        if self.cache and key in self._pool_states:
            self._pool_states.move_to_end(key)
            return True
        base_species, violating, incident = frozenset(), frozenset(), frozenset()
        if self.cache:
            for cid in sorted(key):
                base = self._pool_states.get(key - {cid})
                if base is not None:
                    base_species, violating, incident = base
                    break
        species = frozenset(
            self.catalog.site_by_id[sid].sequence
            for cid in key
            for sid in self.configurations[cid].forward_site_ids
            + self.configurations[cid].reverse_site_ids
        )
        added = species - base_species
        violations = set(violating)
        incidents = set(incident)
        edges = combinations_with_replacement(sorted(added), 2)
        cross = product(sorted(added), sorted(base_species))
        for group in (edges, cross):
            for a, b in group:
                score = self._score(a, b).score
                if self.policy.rejects(score):
                    return False
                if score <= self.policy.strict_cutoff:
                    violations.add(tuple(sorted((a, b))))
                    incidents.update((a, b))
                    if (
                        len(violations) > self.policy.max_violating_edges
                        or len(incidents) > self.policy.max_incident_species
                    ):
                        return False
        if self.cache:
            self._pool_states[key] = (
                species,
                frozenset(violations),
                frozenset(incidents),
            )
            if len(self._pool_states) > 2048:
                self._pool_states.popitem(last=False)
        return True

    def pool_exposure_counts(self, configuration_ids):
        """Exact (strict-violating edges, incident species) for a valid pool.

        Invalid pools raise; callers must never turn unknown/failed guards into
        a zero-exposure objective. Full witnesses remain pool_diagnostics work.
        """
        key = frozenset(configuration_ids)
        if not self.pool_valid(key):
            raise ValueError("exposure objective requested for invalid pool")
        if self.cache:
            _, edges, species = self._pool_states[key]
            return len(edges), len(species)
        result = self.pool_diagnostics(key)
        return result["violating_edge_count"], result["incident_species_count"]

    def _chemistry(self, site):
        profiles = json.loads(self.catalog.resolved_config_json).get("profiles", {})
        measurements = {}
        passing = set()
        for name in sorted(set(site.accepting_profile_ids)):
            if name not in profiles:
                measurements[name] = {"outcome": "fail", "reason": "undeclared-profile"}
                continue
            try:
                checks = chemistry_measurements(site.sequence, profiles[name])
                measurements[name] = checks
                if all(v["outcome"] == "pass" for v in checks.values()):
                    passing.add(name)
            except (ValueError, KeyError, TypeError) as exc:
                measurements[name] = {
                    "outcome": "fail",
                    "reason": "invalid-complete-profile",
                    "detail": str(exc),
                }
        reasons = []
        if not passing:
            reasons.append("chemistry")
        if set(site.accepting_profile_ids) - passing:
            reasons.append("false-profile-membership")
        return {
            "reasons": reasons,
            "accepting_profile_ids": sorted(passing),
            "measurements": measurements,
        }

    def candidate_diagnostics(self, cid):
        if self.cache and cid in self._candidates:
            return deepcopy(self._candidates[cid])
        reasons, chemistry = [], {}
        checks = {
            name: "not-evaluated"
            for name in ("geometry-support", "chemistry", "dimer", "specificity")
        }
        result = {
            "reasons": reasons,
            "checks": checks,
            "chemistry": chemistry,
            "support": {},
            "rejected_products": [],
            "uncertainty_blocks": [],
            "allowed_secondary_products": [],
            "dimer_edges": [],
        }
        c = self.configurations.get(cid)
        try:
            if c is None:
                raise ConfigurationIneligible("unknown-configuration")
            rebuilt = make_configuration(
                self.catalog,
                c.family_id,
                c.forward_site_ids,
                c.reverse_site_ids,
                min_size=self.profile.amplicon_size_min,
                max_size=self.profile.amplicon_size_max,
            )
            if (
                c.target_id != rebuilt.target_id
                or c.anchor_pair != rebuilt.anchor_pair
                or c.full_interval != rebuilt.full_interval
            ):
                raise ConfigurationIneligible("inconsistent-configuration-geometry")
        except (ConfigurationIneligible, KeyError, ValueError) as exc:
            reasons.append(str(exc))
            checks["geometry-support"] = "fail"
        else:
            checks["geometry-support"] = "pass"
            sequences = []
            for sid in c.forward_site_ids + c.reverse_site_ids:
                site = self.catalog.site_by_id[sid]
                chemistry[sid] = self._chemistry(site)
                reasons.extend(chemistry[sid]["reasons"])
                sequences.append(site.sequence)
                if len(site.sequence) < self.profile.mismatch_kmersize:
                    reasons.append("oligo-shorter-than-terminal-kmer")
            checks["chemistry"] = (
                "fail" if any(x["reasons"] for x in chemistry.values()) else "pass"
            )
            if all(len(s) >= 3 for s in sequences):
                exposure = self.pool_diagnostics((cid,))
                result["dimer_edges"] = exposure["dimer_edges"]
                reasons.extend(exposure["reasons"])
                checks["dimer"] = "fail" if exposure["reasons"] else "pass"
            else:
                reasons.append("unsupported-native-dimer-length")
            if all(len(s) >= self.profile.mismatch_kmersize for s in sequences):
                specificity = self.specificity.intrinsic(rebuilt)
                result.update(
                    support=specificity.support,
                    rejected_products=list(specificity.rejected_products),
                    allowed_secondary_products=list(
                        specificity.allowed_secondary_products
                    ),
                    uncertainty_blocks=[
                        p for p in specificity.rejected_products if p["uncertain"]
                    ],
                )
                checks["specificity"] = (
                    "fail" if specificity.rejected_products else "pass"
                )
                if specificity.rejected_products:
                    reasons.append("specificity")
        result["reasons"] = sorted(set(reasons))
        if self.cache:
            self._candidates[cid] = deepcopy(result)
        return result

    def candidate_valid(self, cid):
        if self.cache and cid in self._candidates:
            return not self._candidates[cid]["reasons"]
        return not self.candidate_diagnostics(cid)["reasons"]

    def candidate_reasons(self, cid):
        return tuple(self.candidate_diagnostics(cid)["reasons"])

    def pair_diagnostics(self, a_id, b_id):
        key = tuple(sorted((a_id, b_id)))
        if self.cache and key in self._pairs:
            return deepcopy(self._pairs[key])
        reasons = []
        result = {
            "reasons": reasons,
            "rejected_products": [],
            "allowed_secondary_products": [],
            "uncertainty_blocks": [],
            "dimer_edges": [],
            "checks": {
                "overlap": "not-evaluated",
                "dimer": "not-evaluated",
                "specificity": "not-evaluated",
            },
        }
        a, b = (self.configurations.get(x) for x in key)
        if a is None or b is None:
            reasons.append("unknown-configuration")
        elif not self.candidate_valid(a_id) or not self.candidate_valid(b_id):
            reasons.append("invalid-configuration")
        else:
            overlap = a.target_id == b.target_id and max(
                a.full_interval[0], b.full_interval[0]
            ) < min(a.full_interval[1], b.full_interval[1])
            if overlap:
                reasons.append("same-target-overlap")
            result["checks"]["overlap"] = "fail" if overlap else "pass"
            exposure = self.pool_diagnostics(key)
            reasons.extend(exposure["reasons"])
            result["dimer_edges"] = exposure["dimer_edges"]
            result["checks"]["dimer"] = "fail" if exposure["reasons"] else "pass"
            specificity = self.specificity.pair(a, b)
            result.update(
                rejected_products=list(specificity.rejected_products),
                allowed_secondary_products=list(specificity.allowed_secondary_products),
                uncertainty_blocks=[
                    p for p in specificity.rejected_products if p["uncertain"]
                ],
            )
            result["checks"]["specificity"] = (
                "fail" if specificity.rejected_products else "pass"
            )
            if specificity.rejected_products:
                reasons.append("specificity")
        result["reasons"] = sorted(set(reasons))
        if self.cache:
            self._pairs[key] = deepcopy(result)
        return result

    def conflict(self, a, b):
        key = tuple(sorted((a, b)))
        if self.cache and key in self._pairs:
            return bool(self._pairs[key]["reasons"])
        return bool(self.pair_diagnostics(a, b)["reasons"])

    def conflict_reasons(self, a, b):
        return tuple(self.pair_diagnostics(a, b)["reasons"])


def chemistry_measurements(sequence, settings):
    """Uncached independent native thermochemistry, with all checks measured."""
    from primalscheme3.core.variant_thermo import _FIELDS

    if set(settings) != set(_FIELDS):
        raise ValueError("complete profile must contain all resolved chemistry fields")
    if any(
        isinstance(x, float | int) and not math.isfinite(x) for x in settings.values()
    ):
        raise ValueError("chemistry settings must be finite")
    c = settings
    salts = {k: c[k] for k in ("mv_conc", "dv_conc", "dntp_conc", "dna_conc")}
    checks = {}

    def interval(name, value, lo, hi):
        if lo > hi or not math.isfinite(value):
            raise ValueError("invalid chemistry measurement/bounds")
        checks[name] = {
            "value": value,
            "minimum": lo,
            "maximum": hi,
            "outcome": "pass" if lo <= value <= hi else "fail",
        }

    interval("length", len(sequence), c["primer_size_min"], c["primer_size_max"])
    interval("gc", thermo.gc(sequence), c["primer_gc_min"], c["primer_gc_max"])
    if c["use_annealing"] and c["primer_annealing_prop"] is not None:
        value = thermo.calc_annealing(
            sequence, **salts, temp_c=c["primer_annealing_tempc"]
        )
        interval(
            "annealing",
            value,
            c["primer_annealing_prop"] - thermo.ANNEALING_DIFF,
            c["primer_annealing_prop"] + thermo.ANNEALING_DIFF,
        )
    else:
        interval(
            "tm",
            thermo.calc_tm(sequence, **salts),
            c["primer_tm_min"],
            c["primer_tm_max"] + 2,
        )
    interval("homopolymer", thermo.max_homo(sequence), 0, c["primer_homopolymer_max"])
    interval(
        "hairpin",
        thermo.calc_hairpin_tm(sequence, **salts),
        -1e9,
        c["primer_hairpin_th_max"],
    )
    return checks


def _rebuild_target(target):
    rows = tuple(tuple(c.upper() for c in row) for row in target.rows)
    if (
        not rows
        or len(rows) != len(target.row_ids)
        or len({len(row) for row in rows}) != 1
    ):
        raise ValueError("invalid authoritative alignment rows")
    indexes = tuple(i for i, c in enumerate(rows[0]) if c not in ("", "-"))
    inverse = {col: pos for pos, col in enumerate(indexes)}
    return replace(
        target,
        rows=rows,
        reference_sequence="".join(rows[0][i] for i in indexes),
        reference_length=len(indexes),
        mapping=tuple(inverse.get(i) for i in range(len(rows[0]))),
        ref_to_alignment=indexes + ((indexes[-1] + 1 if indexes else 0),),
    )


def _selected_records(catalog, assignments, configurations):
    records = []
    for assignment in sorted(assignments, key=lambda a: (a.configuration_id, a.pool)):
        c = configurations[assignment.configuration_id]
        for sid in sorted(c.forward_site_ids + c.reverse_site_ids):
            s = catalog.site_by_id[sid]
            records.append(
                {
                    "configuration_id": c.id,
                    "pool": assignment.pool,
                    "site_id": sid,
                    "sequence": s.sequence,
                    "strand": s.strand,
                    "alignment_anchor": s.alignment_anchor,
                    "reference_footprint": list(s.reference_footprint),
                }
            )
    return records


def validate_allele_assignments(
    catalog,
    assignments,
    ledger,
    profile,
    *,
    policy=DEFAULT_STAGE_POLICY,
    goal=0.95,
    expected_summary=None,
    selected_sites=None,
    authoritative_targets=None,
    history=None,
):
    """Fresh reconstruction, never reading optimizer scientific caches.

    Publication supplies separately loaded authoritative_targets and selected_sites
    parsed from final artifacts, and compares expected_summary. With them omitted,
    this validates the in-memory selection against the catalog's authoritative rows.
    """
    assignments = tuple(assignments)
    violations = []
    if ledger.catalog_digest != catalog.semantic_digest:
        violations.append({"reason": "ledger-catalog-mismatch"})
    supplied_targets = tuple(
        catalog.targets if authoritative_targets is None else authoritative_targets
    )
    targets = tuple(_rebuild_target(t) for t in supplied_targets)
    if tuple(sorted(targets, key=lambda t: t.id)) != tuple(
        sorted(catalog.targets, key=lambda t: t.id)
    ):
        violations.append({"reason": "authoritative-target-mismatch"})
    observations = tuple(a for t in targets for a in canonical_observations(t))
    if tuple(sorted(observations, key=lambda a: a.id)) != tuple(
        sorted(catalog.observations, key=lambda a: a.id)
    ):
        violations.append({"reason": "observation-mismatch"})
    fresh = replace(catalog, targets=targets, observations=observations)
    rebuilt = {}
    for cid in sorted({a.configuration_id for a in assignments}):
        c = ledger.configuration_by_id.get(cid)
        try:
            if c is None:
                raise ConfigurationIneligible("unknown-configuration")
            actual = make_configuration(
                fresh,
                c.family_id,
                c.forward_site_ids,
                c.reverse_site_ids,
                min_size=profile.amplicon_size_min,
                max_size=profile.amplicon_size_max,
            )
            if (
                c.target_id != actual.target_id
                or c.anchor_pair != actual.anchor_pair
                or c.full_interval != actual.full_interval
            ):
                violations.append(
                    {
                        "reason": "configuration-geometry-mismatch",
                        "configuration_id": cid,
                    }
                )
            if c.supported_products != actual.supported_products:
                violations.append(
                    {
                        "reason": "configuration-support-mismatch",
                        "configuration_id": cid,
                    }
                )
            rebuilt[cid] = actual
        except (ValueError, KeyError) as exc:
            violations.append(
                {
                    "reason": "configuration-reconstruction",
                    "configuration_id": cid,
                    "detail": str(exc),
                }
            )
    fresh_ledger = ConfigurationLedger(fresh.semantic_digest, tuple(rebuilt.values()))
    evaluator = AlleleCompatibilityOracle(
        fresh, fresh_ledger, profile, policy, cache=False
    )
    counts = Counter((a.configuration_id, a.pool) for a in assignments)
    if any(n > 1 for n in counts.values()):
        violations.append({"reason": "duplicate-assignment"})
    if profile.max_amplicons is not None and len(assignments) > profile.max_amplicons:
        violations.append({"reason": "max-amplicons"})
    target_counts = Counter(
        rebuilt[a.configuration_id].target_id
        for a in assignments
        if a.configuration_id in rebuilt
    )
    if profile.max_amplicons_msa is not None:
        violations.extend(
            {"reason": "max-amplicons-msa", "target_id": tid}
            for tid, n in target_counts.items()
            if n > profile.max_amplicons_msa
        )
    accepted = []
    support = {}
    allowed = []
    uncertainty = []
    pools = {}
    kernels = {name: version(name) for name in ("primalschemers", "primer3-py")}

    def record_checks(result, entity_ids, pool):
        for reason in result["reasons"]:
            violations.append(
                {
                    "reason": reason,
                    "configuration_ids": list(entity_ids),
                    "pool": pool,
                    "rejected_products": result.get("rejected_products", [])
                    if reason == "specificity"
                    else [],
                }
            )
        uncertainty.extend(result.get("uncertainty_blocks", []))
        allowed.extend(
            dict(p, pool=pool) for p in result.get("allowed_secondary_products", [])
        )
        if history is not None:
            evidence = history.record_evidence(
                entity_ids=entity_ids,
                measurement="fresh-validation",
                dependency_key={
                    "catalog": fresh.semantic_digest,
                    "profile": profile.cache_key,
                    "policy": policy.cache_key,
                    "kernels": kernels,
                },
                values=result,
                status="measured",
            )
            for name, outcome in result.get("checks", {}).items():
                history.assess(
                    stage_id=policy.stage_id,
                    entity_ids=entity_ids,
                    pool=pool,
                    context_digest=profile.cache_key,
                    profile_id="allele-panel-v1",
                    kernel_versions=kernels,
                    thresholds=profile.to_dict() | asdict(policy),
                    check_name=name,
                    outcome=outcome,
                    reason=result.get("check_reasons", {}).get(
                        name,
                        "upstream-check-failed"
                        if outcome == "not-evaluated"
                        else "freshly-recomputed",
                    ),
                    evidence_ids=(evidence.id,),
                )

    for assignment in assignments:
        if (
            type(assignment.pool) is not int
            or not 0 <= assignment.pool < profile.n_pools
        ):
            violations.append(
                {
                    "reason": "pool-bounds",
                    "configuration_id": assignment.configuration_id,
                }
            )
        if assignment.stage_id != policy.stage_id:
            violations.append(
                {
                    "reason": "stage-mismatch",
                    "configuration_id": assignment.configuration_id,
                }
            )
        if assignment.configuration_id not in rebuilt:
            if history is not None:
                record_checks(
                    {
                        "reasons": [],
                        "checks": dict.fromkeys(
                            ("chemistry", "dimer", "specificity"), "not-evaluated"
                        )
                        | {"geometry-support": "fail"},
                    },
                    (assignment.configuration_id,),
                    assignment.pool,
                )
            continue
        result = evaluator.candidate_diagnostics(assignment.configuration_id)
        support[assignment.configuration_id] = result["support"]
        record_checks(result, (assignment.configuration_id,), assignment.pool)
        accepted.append(assignment)
    for a, b in combinations(accepted, 2):
        if a.pool != b.pool:
            continue
        ca, cb = rebuilt[a.configuration_id], rebuilt[b.configuration_id]
        overlap = ca.target_id == cb.target_id and max(
            ca.full_interval[0], cb.full_interval[0]
        ) < min(ca.full_interval[1], cb.full_interval[1])
        result = {
            "reasons": ["same-target-overlap"] if overlap else [],
            "checks": {
                "overlap": "fail" if overlap else "pass",
                "specificity": "not-evaluated",
                "dimer": "not-evaluated",
            },
            "check_reasons": {"dimer": "deferred-to-aggregate-pool-check"},
        }
        site_ids = (
            ca.forward_site_ids
            + ca.reverse_site_ids
            + cb.forward_site_ids
            + cb.reverse_site_ids
        )
        if all(
            len(fresh.site_by_id[s].sequence) >= profile.mismatch_kmersize
            for s in site_ids
        ):
            specificity = evaluator.specificity.pair(ca, cb)
            result.update(
                rejected_products=list(specificity.rejected_products),
                uncertainty_blocks=[
                    p for p in specificity.rejected_products if p["uncertain"]
                ],
                allowed_secondary_products=list(specificity.allowed_secondary_products),
            )
            result["checks"]["specificity"] = (
                "fail" if specificity.rejected_products else "pass"
            )
            if specificity.rejected_products:
                result["reasons"].append("specificity")
        record_checks(
            result, tuple(sorted((a.configuration_id, b.configuration_id))), a.pool
        )
    for pool in sorted({a.pool for a in accepted}):
        cids = tuple(a.configuration_id for a in accepted if a.pool == pool)
        try:
            result = evaluator.pool_diagnostics(cids)
        except ValueError as exc:
            result = {
                "reasons": ["unsupported-native-dimer-length"],
                "detail": str(exc),
            }
        pools[pool] = result
        record_checks(
            result
            | {
                "checks": {
                    "aggregate-pool-exposure": "fail" if result["reasons"] else "pass"
                }
            },
            tuple(sorted(set(cids))),
            pool,
        )
    unique = tuple({(a.configuration_id, a.pool): a for a in accepted}.values())
    summary = allele_summary(fresh, unique, fresh_ledger, goal=goal)
    if expected_summary is not None:
        expected = (
            expected_summary.to_dict()
            if hasattr(expected_summary, "to_dict")
            else expected_summary
        )
        if expected != summary.to_dict():
            violations.append({"reason": "coverage-summary-mismatch"})
    output = _selected_records(fresh, unique, rebuilt)
    if selected_sites is not None:

        def canonical(records):
            return sorted(
                json.dumps(r, sort_keys=True, separators=(",", ":")) for r in records
            )

        if canonical(selected_sites) != canonical(output):
            violations.append({"reason": "selected-output-mismatch"})
    return {
        "schemaVersion": "primalscheme3.allele-panel-validation/v2",
        "valid": not violations,
        "violations": violations,
        "profile": profile.to_dict(),
        "stage_policy": asdict(policy),
        "coverage": summary.to_dict(),
        "support_diagnostics": support,
        "selected_sites": output,
        "allowed_secondary_products": allowed,
        "uncertainty_blocks": uncertainty,
        "pools": pools,
        "catalog_semantic_digest": fresh.semantic_digest,
        "assignments": [asdict(a) for a in assignments],
        "validation_scope": "external-authoritative-inputs-and-selected-artifacts"
        if authoritative_targets is not None and selected_sites is not None
        else "catalog-authoritative-rows-and-in-memory-selection",
        "limitations": [
            "Supplied MSA rows only; no whole-genome specificity claim.",
            "No terminal seeds invented outside supplied row coordinates.",
            "Terminal screening predicts potential products; secondary products have zero coverage credit.",
        ],
    }
