"""Versioned allele preset resolution shared by CLI and programmatic callers."""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, fields
from enum import Enum
from pathlib import Path

PRESET = "allele-balanced-v1"
WORK_DEFAULTS = {
    "frontier_candidates": 64,
    "construction_candidate_attempts": 2048,
    "repair_candidate_probes_per_round": 128,
    "repair_neighborhoods_per_round": 256,
    "repair_trials_per_round": 256,
    "pool_lookahead_candidates": 4,
    "cleanup_moves_per_round": 64,
    "families_per_refresh": 16,
}
NEW_OPTION_NAMES = frozenset(
    (
        "preset",
        "candidate_profiles",
        "variant_selection",
        "allele_weighting",
        "discovery_length_mode",
        "specificity_terminal_k",
        "secondary_product_policy",
        "subset_beam_width",
        "subset_expansion_limit",
        "exchange_width",
        "salvage",
        "salvage_thresholds",
        "salvage_max_stages",
        "salvage_max_edges_per_pool",
        "salvage_max_oligos_per_pool",
        "salvage_time_limit",
        "primary_tier",
        *("work_" + name for name in WORK_DEFAULTS),
    )
)


def _plain(value):
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {k: _plain(v) for k, v in value.items()}
    if isinstance(value, list | tuple):
        return [_plain(v) for v in value]
    return value


@dataclass(frozen=True)
class AlleleOptions:
    preset: str = PRESET
    candidate_profiles: str = "union"
    variant_selection: str = "subsets"
    allele_weighting: str = "distinct-observed"
    discovery_length_mode: str = "first-compatible"
    coverage_metric: str = "observed-allele-primer-trimmed"
    coverage_target: float = 0.95
    amplicon_size: int = 400
    amplicon_size_min: int = 360
    amplicon_size_max: int = 440
    n_pools: int = 2
    ncores: int = 1
    min_base_freq: float = 0.0
    max_amplicons: int | None = None
    max_amplicons_msa: int | None = None
    dimer_score: float = -26.0
    specificity_terminal_k: int = 17
    mismatch_product_size: int = 2000
    secondary_product_policy: str = "ordered-disjoint-intended-sites"
    optimizer_seed: int = 0
    optimizer_starts: int = 4
    optimizer_repair_rounds: int = 2
    optimizer_time_limit: float = 120.0
    subset_beam_width: int = 16
    subset_expansion_limit: int = 256
    exchange_width: int = 2
    work_frontier_candidates: int = 64
    work_construction_candidate_attempts: int = 2048
    work_repair_candidate_probes_per_round: int = 128
    work_repair_neighborhoods_per_round: int = 256
    work_repair_trials_per_round: int = 256
    work_pool_lookahead_candidates: int = 4
    work_cleanup_moves_per_round: int = 64
    work_families_per_refresh: int = 16
    salvage: str = "off"
    salvage_thresholds: tuple[float, ...] = (-28.0, -30.0, -32.0)
    salvage_max_stages: int = 3
    salvage_max_edges_per_pool: int = 8
    salvage_max_oligos_per_pool: int = 4
    salvage_time_limit: float = 60.0
    primary_tier: str = "strict"
    requested_options_json: str = "{}"

    def __post_init__(self):
        choices = {
            "preset": (PRESET,),
            "candidate_profiles": ("union", "normal", "high-gc"),
            "variant_selection": ("full-cloud", "subsets"),
            "allele_weighting": ("distinct-observed",),
            "discovery_length_mode": ("first-compatible", "all"),
            "coverage_metric": ("observed-allele-primer-trimmed",),
            "salvage": ("off", "bounded"),
            "secondary_product_policy": (
                "ordered-disjoint-intended-sites",
                "reject-secondary-products/v1",
            ),
        }
        for name, allowed in choices.items():
            if getattr(self, name) not in allowed:
                raise ValueError(f"{name} must be one of {allowed}")
        nonnegative = {
            "optimizer_repair_rounds",
            "exchange_width",
            "salvage_max_edges_per_pool",
            "salvage_max_oligos_per_pool",
        }
        for f in fields(self):
            value = getattr(self, f.name)
            if type(f.default) is int:
                if type(value) is not int:
                    raise ValueError(f"{f.name} must be an integer")
                if f.name != "optimizer_seed" and value < (
                    0 if f.name in nonnegative else 1
                ):
                    raise ValueError(
                        f"{f.name} is outside its nonnegative/positive bounds"
                    )
            elif type(f.default) is float:
                if type(value) not in (int, float) or not math.isfinite(value):
                    raise ValueError(f"{f.name} must be finite numeric")
                object.__setattr__(self, f.name, float(value))
        if not 0 < self.coverage_target <= 1 or not 0 <= self.min_base_freq <= 1:
            raise ValueError("coverage_target/min_base_freq outside unit bounds")
        if self.optimizer_time_limit <= 0 or self.salvage_time_limit <= 0:
            raise ValueError("time limits must be positive")
        if not self.amplicon_size_min <= self.amplicon_size <= self.amplicon_size_max:
            raise ValueError("amplicon bounds must contain nominal size")
        if self.dimer_score != -26:
            raise ValueError("strict dimer_score is fixed at -26; use salvage")
        if self.exchange_width > 2:
            raise ValueError("exchange_width must be <=2")
        if self.salvage_max_stages > 3:
            raise ValueError("salvage_max_stages must be <=3")
        for name in ("max_amplicons", "max_amplicons_msa"):
            value = getattr(self, name)
            if value is not None and (type(value) is not int or value < 0):
                raise ValueError(f"{name} must be nonnegative integer")
        ladder = self.salvage_thresholds
        if not isinstance(ladder, tuple | list) or not ladder:
            raise ValueError("salvage_thresholds must be a nonempty sequence")
        prior = -26
        for cutoff in ladder:
            if (
                type(cutoff) not in (int, float)
                or not math.isfinite(cutoff)
                or cutoff >= prior
            ):
                raise ValueError(
                    "salvage thresholds must be finite and strictly decreasing from -26"
                )
            prior = cutoff
        if len(ladder) > self.salvage_max_stages:
            raise ValueError("salvage ladder exceeds explicit stage limit")
        object.__setattr__(self, "salvage_thresholds", tuple(float(x) for x in ladder))
        tiers = (
            ("strict",) + tuple(f"salvage-{i + 1}" for i in range(len(ladder)))
            if self.salvage == "bounded"
            else ("strict",)
        )
        if self.primary_tier not in tiers:
            raise ValueError("primary_tier must name an enabled requested stage")
        requested = json.loads(self.requested_options_json)
        if self.salvage == "off" and any(k.startswith("salvage_") for k in requested):
            raise ValueError("salvage controls require --salvage bounded")

    def to_dict(self):
        result = asdict(self)
        result["requested_options"] = json.loads(result.pop("requested_options_json"))
        result["salvage_thresholds"] = list(self.salvage_thresholds)
        return result

    def config_values(self):
        result = asdict(self)
        result.pop("requested_options_json")
        return result | {
            "selection_algorithm": "allele-coverage",
            "terminal_gap_policy": "observed-only",
            "amplicon_size_metric": "reference-span",
            "mapping": "first",
            "mismatch_kmersize": self.specificity_terminal_k,
        }

    @classmethod
    def from_config(cls, config):
        data = json.loads(config.allele_options_json)
        data["requested_options_json"] = json.dumps(
            data.pop("requested_options"), sort_keys=True
        )
        return cls(**data)

    def search_options(self):
        from .allele_search import AlleleSearchOptions, AlleleWorkLimits

        return AlleleSearchOptions(
            coverage_target=self.coverage_target,
            seed=self.optimizer_seed,
            starts=self.optimizer_starts,
            repair_rounds=self.optimizer_repair_rounds,
            time_limit=self.optimizer_time_limit,
            subset_beam_width=self.subset_beam_width,
            subset_expansion_limit=self.subset_expansion_limit,
            exchange_width=self.exchange_width,
            requested_amplicon_size=self.amplicon_size,
            work_limits=AlleleWorkLimits(
                **{name: getattr(self, "work_" + name) for name in WORK_DEFAULTS}
            ),
            variant_selection=self.variant_selection,
        )

    def salvage_options(self):
        from .coverage_salvage import SalvageOptions

        return SalvageOptions(
            mode=self.salvage,
            thresholds=self.salvage_thresholds,
            max_stages=self.salvage_max_stages,
            max_edges_per_pool=self.salvage_max_edges_per_pool,
            max_oligos_per_pool=self.salvage_max_oligos_per_pool,
            time_limit=self.salvage_time_limit,
        )

    def stage_policies(self):
        from .allele_validation import StagePolicy

        return (StagePolicy(),) + tuple(
            StagePolicy(
                stage_id=f"salvage-{i + 1}",
                active_cutoff=cutoff,
                max_violating_edges=self.salvage_max_edges_per_pool,
                max_incident_species=self.salvage_max_oligos_per_pool,
            )
            for i, cutoff in enumerate(self.salvage_thresholds)
            if self.salvage == "bounded"
        )


def resolve_allele_options(params):
    from primalscheme3.core.config import Config

    raw = {k: _plain(v) for k, v in params.items() if v is not None}
    known = (
        set(dir(Config))
        | {f.name for f in fields(AlleleOptions)}
        | {
            "allele_options_json",
            "allele_requested_options_json",
            "_allele_requested_options",
            "discovery_workers_by_msa",
            "mismatch_kmersize",
            "primer_size_min",
            "primer_size_max",
            "primer_gc_min",
            "primer_gc_max",
            "mode",
            "msa",
            "region_bedfile",
            "bedfile",
            "max_amplicons_region_group",
            "offline_plots",
        }
    )
    unknown = set(raw) - known
    if unknown:
        raise ValueError(
            "unknown allele configuration options: " + ", ".join(sorted(unknown))
        )
    saved = raw.get("allele_options_json")
    if saved:
        saved_values = json.loads(saved)
        original = dict(saved_values["requested_options"])
        for name in {f.name for f in fields(AlleleOptions)} - {
            "requested_options_json"
        }:
            if name in raw and raw[name] != saved_values.get(name):
                original[name] = raw[name]
    else:
        original = raw.get("_allele_requested_options", raw)
    original = {
        k: v
        for k, v in original.items()
        if not k.startswith("_")
        and k not in ("allele_options_json", "allele_requested_options_json")
    }
    ignored = (
        "downsample_target",
        "downsample_always_add_prop",
        "in_memory_db",
        "mismatch_in_memory",
        "min_overlap",
    )
    if any(name in original for name in ignored):
        raise ValueError("ignored legacy control is unsupported for allele-coverage")
    incompatible = {
        "mode": "equal",
        "terminal_gap_policy": "observed-only",
        "mapping": "first",
        "amplicon_size_metric": "reference-span",
        "use_matchdb": True,
        "mismatch_fuzzy": True,
        "editdist_max": 1,
    }
    for name, required in incompatible.items():
        if name in raw and raw[name] != required:
            raise ValueError(f"allele-coverage requires {name}={required}")
    for name in (
        "force",
        "high_gc",
        "circular",
        "backtrack",
        "downsample",
        "use_annealing",
        "ignore_n",
    ):
        if raw.get(name):
            raise ValueError(f"allele-coverage does not support {name}")
    for name in (
        "input_bedfile",
        "bedfile",
        "region_bedfile",
        "max_amplicons_region_group",
    ):
        if raw.get(name) is not None:
            raise ValueError(f"allele-coverage does not support {name}")
    chemistry = (
        "primer_tm_min",
        "primer_tm_max",
        "primer_hairpin_th_max",
        "primer_homopolymer_max",
        "primer_size_min",
        "primer_size_max",
        "primer_gc_min",
        "primer_gc_max",
        "mv_conc",
        "dv_conc",
        "dntp_conc",
        "dna_conc",
        "primer_annealing_tempc",
        "primer_annealing_prop",
    )
    if saved:
        derived = {
            "primer_size_min": 19,
            "primer_size_max": 36,
            "primer_gc_min": 30,
            "primer_gc_max": 55,
        }
        for name in chemistry:
            if name in raw and raw[name] != derived.get(
                name, getattr(Config, name, None)
            ):
                raise ValueError(
                    "saved configuration contains unsupported chemistry override"
                )
    if any(name in original for name in chemistry):
        raise ValueError(
            "unsupported chemistry override; use complete candidate profiles"
        )
    if not any(
        raw.get(name) is not None for name in ("amplicon_size_min", "amplicon_size_max")
    ):
        raise ValueError(
            "allele-coverage requires explicit reference-span amplicon bounds"
        )
    names = {f.name for f in fields(AlleleOptions)} - {"requested_options_json"}
    values = {name: raw[name] for name in names if name in raw}
    nominal = values.get("amplicon_size", 400)
    if type(nominal) is not int:
        raise ValueError("amplicon_size must be an integer")
    values.setdefault("amplicon_size_min", int(nominal * 0.9))
    values.setdefault("amplicon_size_max", int(nominal * 1.1))
    stage_limit = values.get("salvage_max_stages", 3)
    if (
        "salvage_thresholds" not in original
        and type(stage_limit) is int
        and 1 <= stage_limit <= 3
    ):
        values["salvage_thresholds"] = (-28.0, -30.0, -32.0)[:stage_limit]
    if values.get("variant_selection") == "full-cloud" and any(
        name in original for name in ("subset_beam_width", "subset_expansion_limit")
    ):
        raise ValueError(
            "subset proposal controls are ignored by full-cloud variant selection"
        )
    if "mismatch_kmersize" in raw:
        if (
            "specificity_terminal_k" in values
            and values["specificity_terminal_k"] != raw["mismatch_kmersize"]
        ):
            raise ValueError("specificity terminal k aliases disagree")
        values["specificity_terminal_k"] = raw["mismatch_kmersize"]
    return AlleleOptions(
        **values,
        requested_options_json=json.dumps(
            original, sort_keys=True, separators=(",", ":")
        ),
    )
