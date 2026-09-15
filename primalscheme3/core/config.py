import math
import pathlib
from dataclasses import dataclass
from enum import Enum
from importlib.metadata import version
from typing import Any


# Written by Andy Smith, modified by: Chris Kent
class MappingType(Enum):
    """
    Enum for the mapping type
    """

    FIRST = "first"
    CONSENSUS = "consensus"


class TerminalGapPolicy(str, Enum):
    """Explicit CLI discovery policy; legacy retains the upstream Rust backend."""

    LEGACY = "legacy"
    OBSERVED_ONLY = "observed-only"


class AmpliconSizeMetric(str, Enum):
    """Persisted size interpretation, independent of discovery backend."""

    LEGACY_PAIRING = "legacy-pairing"
    REFERENCE_SPAN = "reference-span"


PRIMER_COUNT_ATTR_STRING = "pc"


@dataclass(frozen=True)
class _CoverageNumericOption:
    expected_type: type[int] | type[float]
    coverage_only: bool
    constraint: str = "finite"


# Caller-set numeric values consumed by coverage discovery, its immutable
# scientific profile, or the selector. Derived primer size/GC bounds are
# validated by ConstraintProfile after Config resolves the high-GC preset.
_COVERAGE_NUMERIC_OPTIONS = {
    "coverage_target": _CoverageNumericOption(float, False, "unit_interval"),
    "optimizer_seed": _CoverageNumericOption(int, False),
    "optimizer_starts": _CoverageNumericOption(int, False, "positive"),
    "optimizer_repair_rounds": _CoverageNumericOption(int, False, "nonnegative"),
    "optimizer_time_limit": _CoverageNumericOption(float, False, "positive"),
    "ncores": _CoverageNumericOption(int, True, "positive"),
    "n_pools": _CoverageNumericOption(int, True, "positive"),
    "min_overlap": _CoverageNumericOption(int, True, "nonnegative"),
    "min_base_freq": _CoverageNumericOption(float, True, "unit_interval"),
    "amplicon_size": _CoverageNumericOption(int, True, "positive"),
    "amplicon_size_min": _CoverageNumericOption(int, True, "positive"),
    "amplicon_size_max": _CoverageNumericOption(int, True, "positive"),
    "primer_tm_min": _CoverageNumericOption(float, True),
    "primer_tm_max": _CoverageNumericOption(float, True),
    "primer_annealing_tempc": _CoverageNumericOption(int, True),
    "primer_hairpin_th_max": _CoverageNumericOption(float, True),
    "primer_homopolymer_max": _CoverageNumericOption(int, True),
    "primer_max_walk": _CoverageNumericOption(int, True),
    "editdist_max": _CoverageNumericOption(int, True, "single_mismatch"),
    "mismatch_product_size": _CoverageNumericOption(int, True, "positive"),
    "mv_conc": _CoverageNumericOption(float, True),
    "dv_conc": _CoverageNumericOption(float, True),
    "dntp_conc": _CoverageNumericOption(float, True),
    "dna_conc": _CoverageNumericOption(float, True),
    "dimer_score": _CoverageNumericOption(float, True),
}


def _validate_raw_coverage_numerics(
    kwargs: dict[str, Any], *, coverage_requested: bool
) -> None:
    for field, option in _COVERAGE_NUMERIC_OPTIONS.items():
        if option.coverage_only and not coverage_requested:
            continue
        value = kwargs.get(field)
        if value is None:
            continue
        if option.expected_type is int:
            if type(value) is not int:
                raise ValueError(f"{field} must be an integer")
        elif type(value) not in {int, float}:
            raise ValueError(f"{field} must be numeric")


def _validate_effective_coverage_numerics(config: "Config") -> None:
    for field, option in _COVERAGE_NUMERIC_OPTIONS.items():
        if option.coverage_only and config.selection_algorithm not in {
            "coverage",
            "allele-coverage",
        }:
            continue
        value = getattr(config, field)
        if option.expected_type is int:
            if type(value) is not int:
                raise ValueError(f"{field} must be an integer")
        elif type(value) is not float or not math.isfinite(value):
            raise ValueError(f"{field} must be finite")

        if option.constraint == "positive" and value <= 0:
            raise ValueError(f"{field} must be positive")
        if option.constraint == "nonnegative" and value < 0:
            raise ValueError(f"{field} must be nonnegative")
        if option.constraint == "unit_interval" and not 0 <= value <= 1:
            raise ValueError(f"{field} must be between zero and one")
        if option.constraint == "single_mismatch" and value != 1:
            raise ValueError(
                "coverage supports the existing single-mismatch policy only"
            )


class Config:
    """
    PrimalScheme3 configuration.
    Class properties are defaults, can be overridden
    on instantiation (and will shadow class defaults)
    """

    terminal_gap_policy: TerminalGapPolicy = TerminalGapPolicy.LEGACY

    @property
    def discovery_backend(self) -> str:
        return (
            "python-observed-only"
            if self.terminal_gap_policy == TerminalGapPolicy.OBSERVED_ONLY
            else "rust-legacy"
        )

    @property
    def discovery_core_count(self) -> int:
        return (
            max(self.discovery_workers_by_msa.values(), default=0)
            if self.terminal_gap_policy == TerminalGapPolicy.OBSERVED_ONLY
            else self.ncores
        )

    # Run Settings
    output: pathlib.Path = pathlib.Path("./output")
    force: bool = False
    high_gc: bool = False
    input_bedfile: pathlib.Path | None = None
    version: str = version("primalscheme3")
    ncores = 1
    # Panel selector settings. Legacy remains the historical default.
    selection_algorithm: str = "legacy"
    coverage_metric: str = "full-span"
    coverage_target: float = 0.90
    optimizer_seed: int = 0
    optimizer_starts: int = 4
    optimizer_repair_rounds: int = 2
    optimizer_time_limit: float = 120.0
    # Scheme Settings
    n_pools: int = 2
    min_overlap: int = 10
    mapping: MappingType = MappingType.FIRST
    circular: bool = False
    backtrack: bool = False
    ignore_n: bool = False
    # PrimerCloud settings
    min_base_freq: float = 0.0
    downsample: bool = False
    downsample_target = 0.99
    downsample_always_add_prop = 0.25
    # Amplicon Settings
    amplicon_size: int = 400
    amplicon_size_min: int = 0
    amplicon_size_max: int = 0
    amplicon_size_metric: AmpliconSizeMetric = AmpliconSizeMetric.LEGACY_PAIRING
    # Primer Settings
    _primer_size_default_min: int = 19
    _primer_size_default_max: int = 36
    _primer_size_hgc_min: int = 17
    _primer_size_hgc_max: int = 30
    _primer_gc_default_min: int = 30
    _primer_gc_default_max: int = 55
    _primer_gc_hgc_min: int = 40
    _primer_gc_hgc_max: int = 65
    # Thermo
    primer_tm_min: float = 59.5
    primer_tm_max: float = 62.5
    primer_annealing_tempc = 65
    primer_annealing_prop: float | None = None
    _primer_annealing_prop_default = 20
    use_annealing: bool = False

    primer_hairpin_th_max: float = 51
    primer_homopolymer_max: int = 5
    primer_max_walk: int = 80
    # MatchDB Settings
    use_matchdb: bool = True
    in_memory_db: bool = True
    editdist_max: int = 1
    mismatch_fuzzy: bool = True
    mismatch_kmersize: int  # Same as primer_size_min
    mismatch_product_size: int = 0
    mismatch_in_memory: bool = False  # Use an in memory dict rather than a dbm file
    # Thermodynamic Parameters
    mv_conc: float = 100.0
    dv_conc: float = 2.0
    dntp_conc: float = 0.8
    dna_conc: float = 15.0
    dimer_score: float = -26.0

    def __init__(self, **kwargs: Any) -> None:
        self.discovery_workers_by_msa = {}
        allele_options = None
        if kwargs.get("selection_algorithm", "legacy") != "allele-coverage":
            from primalscheme3.panel.allele_options import NEW_OPTION_NAMES

            if any(kwargs.get(name) is not None for name in NEW_OPTION_NAMES):
                raise ValueError(
                    "allele options require selection_algorithm=allele-coverage"
                )
        if kwargs.get("selection_algorithm") == "allele-coverage":
            from primalscheme3.panel.allele_options import resolve_allele_options

            allele_options = resolve_allele_options(kwargs)
            kwargs = kwargs | allele_options.config_values()
        coverage_requested = kwargs.get("selection_algorithm", "legacy") in {
            "coverage",
            "allele-coverage",
        }
        _validate_raw_coverage_numerics(kwargs, coverage_requested=coverage_requested)
        if coverage_requested:
            # assign_kwargs historically infers conversion from the runtime
            # default. Normalize coverage real defaults so an integer-valued
            # default cannot truncate a valid fractional scientific threshold.
            for field, option in _COVERAGE_NUMERIC_OPTIONS.items():
                if option.expected_type is float:
                    setattr(self, field, float(getattr(self, field)))
        self.assign_kwargs(**kwargs)
        if self.selection_algorithm not in {"legacy", "coverage", "allele-coverage"}:
            raise ValueError(
                "selection_algorithm must be legacy, coverage or allele-coverage"
            )
        allowed_metrics = (
            {"observed-allele-primer-trimmed"}
            if allele_options is not None
            else {"full-span", "primer-trimmed"}
        )
        if self.coverage_metric not in allowed_metrics:
            raise ValueError("coverage_metric must be full-span or primer-trimmed")
        if self.terminal_gap_policy == TerminalGapPolicy.OBSERVED_ONLY and (
            self.downsample or self.use_annealing
        ):
            raise ValueError(
                "observed-only Python discovery does not support experimental downsampling or annealing mode"
            )
        if self.amplicon_size_metric == AmpliconSizeMetric.REFERENCE_SPAN:
            for bound in ("amplicon_size_min", "amplicon_size_max"):
                if kwargs.get(bound) is not None and getattr(self, bound) <= 0:
                    raise ValueError("Explicit amplicon size bounds must be positive")
        # Set amplicon size
        if self.amplicon_size_min == 0:
            self.amplicon_size_min = int(self.amplicon_size * 0.9)
        if self.amplicon_size_max == 0:
            self.amplicon_size_max = int(self.amplicon_size * 1.1)
        if self.amplicon_size_metric == AmpliconSizeMetric.REFERENCE_SPAN:
            if (
                not 0
                < self.amplicon_size_min
                <= self.amplicon_size
                <= self.amplicon_size_max
            ):
                raise ValueError(
                    "Amplicon sizes must satisfy 0 < minimum <= nominal target <= maximum"
                )
            if self.mapping != MappingType.FIRST:
                raise ValueError(
                    "reference-span amplicon bounds require --mapping first; consensus mapping uses alignment columns"
                )
            if self.circular or self.input_bedfile is not None:
                raise ValueError(
                    "reference-span amplicon bounds require fresh linear design without imported primer pairs"
                )
        _validate_effective_coverage_numerics(self)
        if self.selection_algorithm == "coverage":
            if self.amplicon_size_metric != AmpliconSizeMetric.REFERENCE_SPAN:
                raise ValueError(
                    "coverage selection requires explicit reference-span amplicon bounds"
                )
            if self.mapping != MappingType.FIRST:
                raise ValueError("coverage selection requires --mapping first")
            if not self.use_matchdb:
                raise ValueError(
                    "coverage selection requires supplied-MSA specificity (--use-matchdb)"
                )
            if self.mismatch_product_size <= 0:
                raise ValueError(
                    "coverage selection requires a positive mispriming product size"
                )
            if self.downsample or self.use_annealing:
                raise ValueError(
                    "coverage selection does not support experimental downsampling or annealing"
                )
        if self.high_gc:
            self.primer_size_min = self._primer_size_hgc_min
            self.primer_size_max = self._primer_size_hgc_max
            self.primer_gc_min = self._primer_gc_hgc_min
            self.primer_gc_max = self._primer_gc_hgc_max
        else:
            self.primer_size_min = self._primer_size_default_min
            self.primer_size_max = self._primer_size_default_max
            self.primer_gc_min = self._primer_gc_default_min
            self.primer_gc_max = self._primer_gc_default_max
        # Set MisMatch Kmer Size
        self.mismatch_kmersize = self.primer_size_min
        if allele_options is not None:
            import json

            for name, value in allele_options.config_values().items():
                if not hasattr(self, name):
                    setattr(self, name, value)
            self.mismatch_kmersize = allele_options.specificity_terminal_k
            self.allele_options_json = json.dumps(
                allele_options.to_dict(), sort_keys=True, separators=(",", ":")
            )
            self.allele_requested_options_json = allele_options.requested_options_json
        # Set annealing
        if self.use_annealing:
            self.primer_annealing_prop = self._primer_annealing_prop_default

    def items(self) -> dict[str, Any]:
        """
        Return a dict (key, val) for non-private, non-callable members
        """
        items = {}
        for key in [x for x in dir(self) if not x.startswith("_")]:
            if not callable(getattr(self, key)):  # prevent functions
                value = getattr(self, key)
                # Convert Enum and Path objects to their values
                if isinstance(value, Enum):
                    items[key] = value.value
                if isinstance(value, pathlib.Path):
                    items[key] = str(value)
                else:
                    items[key] = value

        return items

    def to_dict(self) -> dict[str, Any]:
        """
        Return a dict (key, val) for non-private, non-callable members
        """
        dict = {}
        for key, val in self.items().items():
            if isinstance(val, Enum):
                dict[key] = val.value
            elif isinstance(val, pathlib.Path):
                dict[key] = str(val)
            else:
                dict[key] = val
        return dict

    def __str__(self) -> str:
        return "\n".join(f"{key}: {val}" for key, val in self.items())

    def assign_kwargs(self, **kwargs: Any) -> None:
        for key, value in kwargs.items():
            # Check if key is valid
            if value is None:
                continue

            if key in {
                "discovery_backend",
                "discovery_core_count",
                "discovery_workers_by_msa",
            }:
                continue  # Computed from effective policy, never caller-controlled.

            if key == "input_bedfile":
                setattr(self, key, pathlib.Path(value))
                continue

            if hasattr(self, key):
                # Convert to expected type
                if isinstance(getattr(self, key), TerminalGapPolicy):
                    setattr(self, key, TerminalGapPolicy(value))
                elif isinstance(getattr(self, key), AmpliconSizeMetric):
                    setattr(self, key, AmpliconSizeMetric(value))
                elif isinstance(getattr(self, key), MappingType):
                    setattr(self, key, MappingType(value))
                elif isinstance(getattr(self, key), pathlib.Path):
                    setattr(self, key, pathlib.Path(value))
                elif isinstance(
                    getattr(self, key), bool
                ):  # Need to check bool before int
                    parsed_bool = str(value).lower()
                    setattr(self, key, parsed_bool == "true")
                elif isinstance(getattr(self, key), int):
                    setattr(self, key, int(value))
                elif isinstance(getattr(self, key), float):
                    setattr(self, key, float(value))
                elif isinstance(getattr(self, key), str):
                    setattr(self, key, str(value))
                else:
                    print(f"Could not parse {key} with value {value} ({type(value)}")


# All bases allowed in the input MSA
IUPAC_ALL_ALLOWED_DNA = {
    "A",
    "G",
    "K",
    "Y",
    "B",
    "S",
    "N",
    "H",
    "C",
    "W",
    "D",
    "R",
    "M",
    "T",
    "V",
    "-",
}

SIMPLE_BASES = {"A", "C", "G", "T"}

AMBIGUOUS_DNA = {
    "M": "AC",
    "R": "AG",
    "W": "AT",
    "S": "CG",
    "Y": "CT",
    "K": "GT",
    "V": "ACG",
    "H": "ACT",
    "D": "AGT",
    "B": "CGT",
}
ALL_DNA: dict[str, str] = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "M": "AC",
    "R": "AG",
    "W": "AT",
    "S": "CG",
    "Y": "CT",
    "K": "GT",
    "V": "ACG",
    "H": "ACT",
    "D": "AGT",
    "B": "CGT",
}
ALL_DNA_WITH_N: dict[str, str] = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "M": "AC",
    "R": "AG",
    "W": "AT",
    "S": "CG",
    "Y": "CT",
    "K": "GT",
    "V": "ACG",
    "H": "ACT",
    "D": "AGT",
    "B": "CGT",
    "N": "ACGT",
}
ALL_BASES: set[str] = {
    "A",
    "C",
    "G",
    "T",
    "M",
    "R",
    "W",
    "S",
    "Y",
    "K",
    "V",
    "H",
    "D",
    "B",
}
ALL_BASES_WITH_N: set[str] = {
    "A",
    "C",
    "G",
    "T",
    "M",
    "R",
    "W",
    "S",
    "Y",
    "K",
    "V",
    "H",
    "D",
    "B",
    "N",
}
AMB_BASES = {"Y", "W", "R", "B", "H", "V", "D", "K", "M", "S"}
AMBIGUOUS_DNA_COMPLEMENT = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "M": "K",
    "R": "Y",
    "W": "W",
    "S": "S",
    "Y": "R",
    "K": "M",
    "V": "B",
    "H": "D",
    "D": "H",
    "B": "V",
    "X": "X",
    "N": "N",
    "-": "-",
}
