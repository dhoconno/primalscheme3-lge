"""Immutable scientific records used by the coverage panel selector."""

from __future__ import annotations

from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Mapping


@dataclass(frozen=True)
class Target:
    """One input-MSA occurrence in stable reference coordinates.

    ``mapping[i]`` is the zero-based reference coordinate represented by
    alignment column ``i``, or ``None`` when the reference has no base there.
    ``ref_to_alignment[r]`` is the alignment column containing reference base
    ``r``; its final sentinel is one column past the last reference base.
    Rows retain individual alignment cells so terminal missing values (``""``)
    remain distinguishable from internal gaps (``"-"``).
    """

    id: str
    source_msa_index: int
    occurrence: int
    row_ids: tuple[str, ...]
    rows: tuple[tuple[str, ...], ...]
    reference_sequence: str
    reference_length: int
    mapping: tuple[int | None, ...]
    ref_to_alignment: tuple[int, ...]


@dataclass(frozen=True)
class Candidate:
    """Canonical primer-pair candidate, independent of its native object.

    ``interior_interval`` is the pair of shared reference anchors: the
    exclusive forward end and inclusive reverse start. Individual cloud oligo
    lengths recover their exact anchored sites from those coordinates.
    """

    id: str
    target_id: str
    full_interval: tuple[int, int]
    interior_interval: tuple[int, int]
    forward_oligos: tuple[str, ...]
    reverse_oligos: tuple[str, ...]
    joint_rows: tuple[str, ...]
    unknown_rows: tuple[str, ...]
    row_product_spans: tuple[tuple[str, int, int], ...]


@dataclass(frozen=True)
class Catalog:
    """Ordered canonical catalogue plus non-semantic native source mapping."""

    targets: tuple[Target, ...]
    candidates: tuple[Candidate, ...]
    semantic_digest: str
    resolved_config_json: str
    source_mapping: tuple[tuple[int, str], ...]
    _native_pair_items: tuple[tuple[str, object], ...] = field(
        repr=False, compare=False, hash=False
    )
    _target_by_id: Mapping[str, Target] = field(
        init=False, repr=False, compare=False, hash=False
    )
    _candidate_by_id: Mapping[str, Candidate] = field(
        init=False, repr=False, compare=False, hash=False
    )
    _native_pair_by_candidate_id: Mapping[str, object] = field(
        init=False, repr=False, compare=False, hash=False
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "_target_by_id",
            MappingProxyType({target.id: target for target in self.targets}),
        )
        object.__setattr__(
            self,
            "_candidate_by_id",
            MappingProxyType(
                {candidate.id: candidate for candidate in self.candidates}
            ),
        )
        object.__setattr__(
            self,
            "_native_pair_by_candidate_id",
            MappingProxyType(dict(self._native_pair_items)),
        )

    @property
    def target_by_id(self) -> Mapping[str, Target]:
        return self._target_by_id

    @property
    def candidate_by_id(self) -> Mapping[str, Candidate]:
        return self._candidate_by_id

    @property
    def native_pair_by_candidate_id(self) -> Mapping[str, object]:
        """Native objects are intentionally excluded from equality and digests."""

        return self._native_pair_by_candidate_id


@dataclass(frozen=True)
class Assignment:
    """Place a canonical candidate into a zero-based pool."""

    candidate_id: str
    pool: int


# Additive v2 contracts. Historical records above retain their v1 meanings.
from dataclasses import fields, is_dataclass
from functools import cached_property
import hashlib
import json
import types as _types
from collections.abc import Mapping as _Mapping
from typing import Any, ClassVar, get_args, get_origin, get_type_hints


def record_value(value: Any) -> Any:
    """JSON-compatible payload (including frozen mappings), without loss of nulls."""
    if is_dataclass(value):
        return {f.name: record_value(getattr(value, f.name)) for f in fields(value)}
    if isinstance(value, _Mapping):
        return {str(k): record_value(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [record_value(v) for v in value]
    return value


def canonical_json(value: Any) -> str:
    return json.dumps(record_value(value), sort_keys=True, separators=(',', ':'), allow_nan=False)


def semantic_id(kind: str, value: Any) -> str:
    payload = {'kind': kind, 'canonicalVersion': 1, 'semanticFields': record_value(value)}
    return kind + '-' + hashlib.sha256(canonical_json(payload).encode()).hexdigest()


def _freeze(value: Any) -> Any:
    if isinstance(value, _Mapping):
        return MappingProxyType({k: _freeze(v) for k, v in sorted(value.items())})
    if isinstance(value, (tuple, list)):
        return tuple(_freeze(v) for v in value)
    return value


def _decode(annotation: Any, value: Any) -> Any:
    if value is None:
        return None
    origin, args = get_origin(annotation), get_args(annotation)
    if origin is _types.UnionType:
        return _decode(next(t for t in args if t is not type(None)), value)
    if origin is tuple:
        return tuple(_decode(args[0] if len(args) == 2 and args[1] is Ellipsis else args[i], v)
                     for i, v in enumerate(value))
    if isinstance(annotation, type) and is_dataclass(annotation):
        hints = get_type_hints(annotation)
        result = annotation(**{f.name: _decode(hints[f.name], value[f.name])
                               for f in fields(annotation) if f.init and f.name in value})
        if isinstance(result, ScientificRecord) and 'id' in value and result.id != value['id']:
            raise ValueError('record identity does not match payload')
        return result
    return _freeze(value)


class ScientificRecord:
    _identity_fields: ClassVar[tuple[str, ...]] = ()

    def __post_init__(self) -> None:
        for f in fields(self):
            if f.init:
                object.__setattr__(self, f.name, _freeze(getattr(self, f.name)))
        if self._identity_fields:
            object.__setattr__(self, 'id', semantic_id(type(self).__name__, {
                name: getattr(self, name) for name in self._identity_fields}))

    def to_dict(self) -> dict:
        return record_value(self)

    @classmethod
    def from_dict(cls, value: dict):
        result = _decode(cls, value)
        if 'id' in value and getattr(result, 'id', None) != value['id']:
            raise ValueError('record identity does not match payload')
        return result


@dataclass(frozen=True)
class ObservedAllele(ScientificRecord):
    target_id: str
    cells: tuple[str, ...]
    alignment_to_row: tuple[int | None, ...]
    row_to_alignment: tuple[int, ...]
    observed_positions: tuple[int, ...]
    row_ids: tuple[str, ...]
    multiplicity: int
    id: str = field(init=False)
    _identity_fields = ('target_id', 'cells')

    def __post_init__(self):
        cells = tuple(cell.upper() for cell in self.cells)
        if any(cell not in ('', '-', *'ACGTRYSWKMBDHVN') for cell in cells):
            raise ValueError('invalid aligned observation cell')
        row_to_alignment = tuple(i for i, cell in enumerate(cells) if cell not in ('', '-'))
        positions = {column: index for index, column in enumerate(row_to_alignment)}
        alignment_to_row = tuple(positions.get(i) for i in range(len(cells)))
        observed = tuple(index for index, column in enumerate(row_to_alignment) if cells[column] in 'ACGT')
        for name, expected in (('row_to_alignment', row_to_alignment),
                               ('alignment_to_row', alignment_to_row),
                               ('observed_positions', observed)):
            if tuple(getattr(self, name)) != expected:
                raise ValueError('observation ' + name + ' is inconsistent with cells')
        if (not self.row_ids or any(not isinstance(x, str) or not x for x in self.row_ids)
                or len(set(self.row_ids)) != len(self.row_ids)
                or type(self.multiplicity) is not int or self.multiplicity != len(self.row_ids)):
            raise ValueError('observation aliases/multiplicity are inconsistent')
        object.__setattr__(self, 'cells', cells)
        object.__setattr__(self, 'row_ids', tuple(sorted(self.row_ids)))
        super().__post_init__()


@dataclass(frozen=True)
class OligoSite(ScientificRecord):
    """Alignment anchor: F exclusive end; R inclusive start. Footprint is reference half-open."""
    target_id: str
    sequence: str
    strand: str
    alignment_anchor: int
    reference_footprint: tuple[int, int] | None
    mapping_failure: str | None = None
    accepting_profile_ids: tuple[str, ...] = ()
    generated_profile_ids: tuple[str, ...] = ()
    intrinsic_evidence_ids: tuple[str, ...] = ()
    id: str = field(init=False)
    _identity_fields = ('target_id', 'sequence', 'strand', 'alignment_anchor', 'reference_footprint')

    def __post_init__(self):
        object.__setattr__(self, 'sequence', self.sequence.upper())
        if self.strand not in ('+', '-') or not self.sequence or set(self.sequence) - set('ACGT'):
            raise ValueError('site requires a nonempty concrete sequence and +/- strand')
        for name in ('accepting_profile_ids', 'generated_profile_ids', 'intrinsic_evidence_ids'):
            object.__setattr__(self, name, tuple(sorted(set(getattr(self, name)))))
        super().__post_init__()

    @property
    def accepting_profiles(self) -> tuple[str, ...]:
        return self.accepting_profile_ids

    @property
    def species_id(self) -> str:
        return semantic_id('physical-species', self.sequence)


@dataclass(frozen=True)
class BindingSupport(ScientificRecord):
    site_id: str
    allele_id: str
    row_id: str
    status: str
    reason: str
    alignment_footprint: tuple[int, ...]
    row_footprint: tuple[int, int] | None
    evidence_id: str | None = None
    id: str = field(init=False)
    _identity_fields = ('site_id', 'allele_id', 'row_id', 'status', 'reason', 'alignment_footprint', 'row_footprint')


@dataclass(frozen=True)
class CandidateFamily(ScientificRecord):
    target_id: str
    anchor_pair: tuple[int, int]
    forward_site_ids: tuple[str, ...]
    reverse_site_ids: tuple[str, ...]
    discovery_profile_combinations: tuple[tuple[str, str], ...] = ()
    geometry_evidence_ids: tuple[str, ...] = ()
    id: str = field(init=False)
    _identity_fields = ('target_id', 'anchor_pair')

    def __post_init__(self):
        for name in ('forward_site_ids', 'reverse_site_ids', 'discovery_profile_combinations', 'geometry_evidence_ids'):
            object.__setattr__(self, name, tuple(sorted(set(getattr(self, name)))))
        super().__post_init__()


@dataclass(frozen=True)
class SupportedProduct(ScientificRecord):
    allele_id: str
    forward_site_id: str
    reverse_site_id: str
    forward_footprint: tuple[int, int]
    reverse_footprint: tuple[int, int]

    @property
    def interior(self) -> tuple[int, int]:
        return self.forward_footprint[1], self.reverse_footprint[0]


@dataclass(frozen=True)
class SelectionConfiguration(ScientificRecord):
    target_id: str
    family_id: str
    forward_site_ids: tuple[str, ...]
    reverse_site_ids: tuple[str, ...]
    full_interval: tuple[int, int]
    anchor_pair: tuple[int, int]
    supported_products: tuple[SupportedProduct, ...] = ()
    geometry_evidence_ids: tuple[str, ...] = ()
    support_evidence_ids: tuple[str, ...] = ()
    id: str = field(init=False)
    _identity_fields = ('target_id', 'family_id', 'forward_site_ids', 'reverse_site_ids')

    def __post_init__(self):
        for name in ('forward_site_ids', 'reverse_site_ids'):
            object.__setattr__(self, name, tuple(sorted(set(getattr(self, name)))))
            if not getattr(self, name):
                raise ValueError('configuration requires nonempty forward and reverse sites')
        super().__post_init__()


@dataclass(frozen=True)
class AlleleAssignment(ScientificRecord):
    configuration_id: str
    pool: int
    stage_id: str

    def __post_init__(self):
        if self.pool < 0:
            raise ValueError('pool must be zero-based nonnegative')
        super().__post_init__()


@dataclass(frozen=True)
class ConfigurationLedger(ScientificRecord):
    catalog_digest: str
    configurations: tuple[SelectionConfiguration, ...]
    schema_version: str = 'primalscheme3.configuration-ledger/v2'

    def __post_init__(self):
        if len({c.id for c in self.configurations}) != len(self.configurations):
            raise ValueError('duplicate configuration IDs')
        object.__setattr__(self, 'configurations', tuple(sorted(self.configurations, key=lambda x: x.id)))
        super().__post_init__()

    @cached_property
    def configuration_by_id(self) -> Mapping[str, SelectionConfiguration]:
        return MappingProxyType({c.id: c for c in self.configurations})

    @cached_property
    def semantic_digest(self) -> str:
        return semantic_id('configuration-ledger', self.to_dict())

    def appended(self, *configurations: SelectionConfiguration) -> ConfigurationLedger:
        records = dict(self.configuration_by_id)
        for c in configurations:
            if c.id in records and records[c.id] != c:
                raise ValueError('immutable configuration identity collision')
            records[c.id] = c
        return ConfigurationLedger(self.catalog_digest, tuple(records.values()))


@dataclass(frozen=True)
class VariantCatalog(ScientificRecord):
    targets: tuple[Target, ...]
    observations: tuple[ObservedAllele, ...]
    sites: tuple[OligoSite, ...]
    families: tuple[CandidateFamily, ...]
    resolved_config_json: str = '{}'
    source_mapping: tuple[tuple[int, str], ...] = ()
    schema_version: str = 'primalscheme3.variant-catalog/v2'

    def __post_init__(self):
        for name in ('targets', 'observations', 'sites', 'families'):
            values = getattr(self, name)
            if len({x.id for x in values}) != len(values):
                raise ValueError('duplicate ' + name + ' IDs')
            object.__setattr__(self, name, tuple(sorted(values, key=lambda x: x.id)))
        super().__post_init__()

    @cached_property
    def target_by_id(self):
        return MappingProxyType({x.id: x for x in self.targets})

    @cached_property
    def selectable_site_ids(self) -> tuple[str, ...]:
        return tuple(s.id for s in self.sites if s.accepting_profile_ids
                     and s.reference_footprint is not None and s.mapping_failure is None)

    @cached_property
    def site_by_id(self):
        return MappingProxyType({x.id: x for x in self.sites})

    @cached_property
    def family_by_id(self):
        return MappingProxyType({x.id: x for x in self.families})

    @cached_property
    def observation_by_id(self):
        return MappingProxyType({x.id: x for x in self.observations})

    @cached_property
    def semantic_digest(self):
        return semantic_id('variant-catalog', {
            'targets': [(x.id, x.rows, x.mapping) for x in self.targets],
            'observations': [x.id for x in self.observations],
            'sites': [x.to_dict() for x in self.sites],
            'families': [x.to_dict() for x in self.families],
            'resolved_config_json': self.resolved_config_json})


@dataclass(frozen=True)
class ClassCoverage(ScientificRecord):
    target_id: str
    allele_id: str
    observed_count: int
    covered_count: int
    unknown_count: int
    missing_count: int
    gap_count: int
    fraction: float | None
    covered_intervals: tuple[tuple[int, int], ...]


@dataclass(frozen=True)
class TargetCoverage(ScientificRecord):
    target_id: str
    fraction: float | None
    status: str
    assessable_classes: int
    unassessable_classes: int


@dataclass(frozen=True)
class CoverageSummary(ScientificRecord):
    classes: tuple[ClassCoverage, ...]
    targets: tuple[TargetCoverage, ...]
    mean_utility: float | None
    mean_coverage: float | None
    worst_fraction: float | None
    strictly_above_goal: int
    at_least_goal: int
    goal: float
    status: str
    metric_id: str = 'observed-allele-primer-trimmed/v1'
