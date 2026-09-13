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
