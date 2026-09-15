# Modules
from collections import Counter
from collections.abc import Callable
from enum import Enum

import numpy as np

# Submodules
from primalschemers import FKmer, RKmer, do_pool_interact  # type: ignore

from primalscheme3.core.classes import PrimerPair
from primalscheme3.core.config import ALL_BASES, AMB_BASES, AmpliconSizeMetric, Config
from primalscheme3.core.errors import (
    ERROR_SET,
    ContainsInvalidBase,
    CustomErrors,
    CustomRecursionError,
    EndOfSequence,
    GapOnSetBase,
    WalksOut,
    WalksTooFar,
)
from primalscheme3.core.get_window import get_r_window_FAST2
from primalscheme3.core.progress_tracker import ProgressManager
from primalscheme3.core.seq_functions import (
    expand_ambs,
    reverse_complement,
)
from primalscheme3.core.thermo import (
    THERMO_RESULT,
    calc_tm,
    thermo_check_kmers,
)

EARLY_RETURN_FREQ = -1


class DIGESTION_ERROR(Enum):
    """
    Enum for the different types of errors that can occur during digestion
    """

    WALKS_OUT = "WalksOut"
    CONTAINS_INVALID_BASE = "ContainsInvalidBase"
    CUSTOM_RECURSION_ERROR = "CustomRecursionError"
    CUSTOM_ERRORS = "CustomErrors"
    GAP_ON_SET_BASE = "GapOnSetBase"
    HAIRPIN_FAIL = "HairpinFail"
    DIMER_FAIL = "DimerFail"  # Interaction within the kmer
    WALK_TO_FAR = "WalkToFar"  # When indels causes the walk to go to far
    AMB_FAIL = "AmbFail"  # Generic error for when the error is unknown
    NO_SEQUENCES = "NoSequences"
    END_OF_SEQUENCE = "EndOfSequence"


class DIGESTION_RESULT:
    seq: str | DIGESTION_ERROR
    count: float
    status: THERMO_RESULT | None

    def __init__(
        self,
        seq: str | DIGESTION_ERROR,
        count: float,
        status: THERMO_RESULT | None = None,
    ):
        self.seq = seq
        self.count = count
        self.status = status

    def thermo_check(self, config: Config) -> THERMO_RESULT | DIGESTION_ERROR:
        # If the seq is an error return the error
        if isinstance(self.seq, DIGESTION_ERROR):
            return self.seq

        # TODO Add cache for thermo check

        # Thermo check the sequence
        self.status = thermo_check_kmers({self.seq}, config)
        return self.status


def parse_error(results: set[CustomErrors | str]) -> DIGESTION_ERROR:
    """
    Parses the error set for the error that occurred
    As only one error is returned, there is an arbitrary hierarchy of errors
    - CONTAINS_INVALID_BASE > END_OF_SEQUENCE > GAP_ON_SET_BASE > WALKS_OUT > CUSTOM_RECURSION_ERROR > WALK_TO_FAR > CUSTOM_ERRORS
    """
    if ContainsInvalidBase() in results:
        return DIGESTION_ERROR.CONTAINS_INVALID_BASE
    elif EndOfSequence() in results:
        return DIGESTION_ERROR.END_OF_SEQUENCE
    elif GapOnSetBase() in results:
        return DIGESTION_ERROR.GAP_ON_SET_BASE
    elif WalksOut() in results:
        return DIGESTION_ERROR.WALKS_OUT
    elif CustomRecursionError() in results:
        return DIGESTION_ERROR.CUSTOM_RECURSION_ERROR
    elif CustomErrors() in results:
        return DIGESTION_ERROR.CUSTOM_ERRORS
    elif WalksTooFar() in results:
        return DIGESTION_ERROR.WALK_TO_FAR
    else:  # Return a generic error
        return DIGESTION_ERROR.AMB_FAIL


def parse_error_list(
    error_list: list[str | CustomErrors],
) -> list[str | DIGESTION_ERROR]:
    """
    Parses a list of errors and returns a list of DIGESTION_ERROR
    """
    return_list = []
    for result in error_list:
        if isinstance(result, str):
            return_list.append(result)
        elif isinstance(result, CustomErrors):
            return_list.append(parse_error({result}))
    return return_list


def generate_valid_primerpairs(
    fkmers: list[FKmer],
    rkmers: list[RKmer],
    amplicon_size_min: int,
    amplicon_size_max: int,
    dimerscore: float,
    msa_index: int,
    progress_manager: ProgressManager,
    chrom: str = "",
    amplicon_size_metric: AmpliconSizeMetric = AmpliconSizeMetric.LEGACY_PAIRING,
) -> list[PrimerPair]:
    """Generates valid primer pairs for a given set of forward and reverse kmers.

    Args:
        fkmers: A list of forward kmers.
        rkmers: A list of reverse kmers.
        cfg: A dictionary containing configuration parameters.
        msa_index: An integer representing the index of the multiple sequence alignment.
        disable_progress_bar: A boolean indicating whether to disable the progress bar.

    Returns:
        A list of valid primer pairs.
    """
    ## Generate all primerpairs without checking
    checked_pp = []
    metric = AmpliconSizeMetric(amplicon_size_metric)
    # Reverse starts remain sorted; reverse outer ends need not be sorted when
    # alternative oligos differ in length. Widen first, then check the BED span.
    reverse_extension = (
        max((rkmer.region()[1] - rkmer.start for rkmer in rkmers), default=0)
        if metric == AmpliconSizeMetric.REFERENCE_SPAN else 0
    )
    pt = progress_manager.create_sub_progress(
        iter=fkmers, process="Generating primer pairs", chrom=chrom
    )
    for fkmer in pt:
        fkmer_start = min(fkmer.starts())
        # Get all rkmers that would make a valid amplicon
        pos_rkmer = get_r_window_FAST2(
            kmers=rkmers,
            start=fkmer_start + amplicon_size_min - reverse_extension,
            end=fkmer_start + amplicon_size_max,
        )
        for rkmer in pos_rkmer:
            if metric == AmpliconSizeMetric.REFERENCE_SPAN:
                span = rkmer.region()[1] - fkmer_start
                if not amplicon_size_min <= span <= amplicon_size_max:
                    continue
            # Check for interactions
            if not do_pool_interact(fkmer.seqs_bytes(), rkmer.seqs_bytes(), dimerscore):
                checked_pp.append(PrimerPair(fkmer, rkmer, msa_index))

        # Update the count
        pt.manual_update(count=len(checked_pp))

    checked_pp.sort(key=lambda pp: (pp.fprimer.end, -pp.rprimer.start))
    return checked_pp


def walk_right(
    array: np.ndarray,
    col_index_right: int,
    col_index_left: int,
    row_index: int,
    seq_str: str,
    config: Config,
) -> set[str] | Exception:
    """
    Walks to the right of the array and returns a set of valid sequences.

    Args:
        array: A numpy array of DNA sequences.
        col_index_right: The current column index to the right.
        col_index_left: The current column index to the left.
        row_index: The current row index.
        seq_str: The current sequence string.
        config: The configuration object.

    Returns:
        A set of valid DNA sequences or an exception if an error occurs.

    Raises:
        WalksOut: If the function walks out of the array size.
        ContainsInvalidBase: If the sequence contains an invalid base.
    """
    # Guard for correct tm
    if (
        calc_tm(
            seq_str,
            mv_conc=config.mv_conc,
            dv_conc=config.dv_conc,
            dna_conc=config.dna_conc,
            dntp_conc=config.dntp_conc,
        )
        >= config.primer_tm_min
    ):
        return {seq_str}

    # Guard prevents walking out of array size
    if col_index_right >= array.shape[1] - 1 or col_index_left >= array.shape[1] - 1:
        raise WalksOut()

    # Guard for walking too far
    if col_index_right - col_index_left >= config.primer_max_walk:
        raise WalksTooFar()

    new_base = array[row_index, col_index_right]

    # Terminal padding is missing coverage, not a base to impute.
    if new_base == "":
        raise EndOfSequence()
    if new_base == "-":
        raise GapOnSetBase()
    new_string = seq_str + new_base

    # Prevent Ns from being added
    if "N" in new_string:
        raise ContainsInvalidBase()

    # Guard for invalid bases in the sequence
    exp_new_string: set[str] | None = expand_ambs([new_string])
    if exp_new_string is None:
        raise ContainsInvalidBase()

    passing_str = []
    for exp_str in exp_new_string:
        results = wrap_walk(
            walk_right,
            array,
            col_index_right + 1,
            col_index_left,
            row_index,
            exp_str,
            config,
        )
        passing_str.extend(results)

    return passing_str  # type: ignore


def walk_left(
    array: np.ndarray,
    col_index_right: int,
    col_index_left: int,
    row_index: int,
    seq_str: str,
    config: Config,
) -> set[str] | Exception:
    """
    Recursively walks left from a given starting position in a 2D numpy array of DNA bases,
    constructing a set of valid DNA sequences that meet certain criteria.

    Args:
        array: A 2D numpy array of DNA bases.
        col_index_right: The rightmost column index of the region of interest.
        col_index_left: The current leftmost column index of the region of interest.
        row_index: The current row index of the region of interest.
        seq_str: The current DNA sequence being constructed.
        cfg: A dictionary of configuration parameters.

    Returns:
        A set of valid DNA sequences that meet the criteria specified in the function body.

    Raises:
        WalksOut: If the function attempts to walk out of the array.
        ContainsInvalidBase: If the constructed sequence contains an invalid DNA base.
    """

    # Guard prevents walking out of array size
    if col_index_left <= 0 or col_index_right <= 0:
        raise WalksOut()

    # Guard for correct tm
    if (
        calc_tm(
            seq_str,
            mv_conc=config.mv_conc,
            dv_conc=config.dv_conc,
            dna_conc=config.dna_conc,
            dntp_conc=config.dntp_conc,
        )
        >= config.primer_tm_min
    ):
        return {seq_str}

    # Guard for walking too far
    if col_index_right - col_index_left >= config.primer_max_walk:
        raise WalksTooFar()

    new_base = array[row_index, col_index_left - 1]

    # Terminal padding is missing coverage, not a base to impute.
    if new_base == "":
        raise EndOfSequence()
    if new_base == "-":
        raise GapOnSetBase()
    new_string = new_base + seq_str

    # Guard prevents seqs with an N
    if "N" in new_string:
        raise ContainsInvalidBase()

    # If invalid bases return None
    exp_new_string: set[str] | None = expand_ambs([new_string])
    if exp_new_string is None:
        raise ContainsInvalidBase()

    passing_str = []
    for exp_str in exp_new_string:
        results = wrap_walk(
            walk_left,
            array=array,
            col_index_right=col_index_right,
            col_index_left=col_index_left - 1,
            row_index=row_index,
            seq_str=exp_str,
            config=config,
        )
        passing_str.extend(results)

    return passing_str  # type: ignore


def wrap_walk(
    walkfunction: Callable,
    array: np.ndarray,
    col_index_right: int,
    col_index_left: int,
    row_index: int,
    seq_str: str,
    config: Config,
) -> list[str | CustomErrors]:
    return_list = []
    try:
        seqs = walkfunction(
            array=array,
            col_index_right=col_index_right,
            col_index_left=col_index_left,
            row_index=row_index,
            seq_str=seq_str,
            config=config,
        )
    except CustomErrors as e:
        return_list.append(e)
    except Exception as e:
        raise e
    else:
        return_list.extend(seqs)

    return return_list


def r_digest_to_result(
    align_array: np.ndarray, config: Config, start_col: int, min_freq: float
) -> tuple[int, list[DIGESTION_RESULT]]:
    """
    Returns the count of each sequence / error at a given index
    A value of -1 in the return dict means the function returned early, and not all seqs were counted. Only used for WALKS_OUT and GAP_ON_SET_BASE
    """

    ### Process early return conditions
    # If the initial slice is outside the range of the array
    if start_col + config.primer_size_min > align_array.shape[1]:
        return (
            start_col,
            [DIGESTION_RESULT(DIGESTION_ERROR.WALKS_OUT, EARLY_RETURN_FREQ, None)],
        )

    # Check for gap frequency on first base
    base, counts = np.unique(align_array[:, start_col], return_counts=True)
    first_base_counter = dict(zip(base, counts, strict=False))
    first_base_counter.pop("", None)

    num_seqs = sum(first_base_counter.values())
    first_base_freq = {k: v / num_seqs for k, v in first_base_counter.items()}

    # If the freq of gap is above minfreq
    gap_freq = first_base_freq.get("-")
    if gap_freq is not None and gap_freq >= min_freq:
        return (
            start_col,
            [
                DIGESTION_RESULT(
                    DIGESTION_ERROR.GAP_ON_SET_BASE, EARLY_RETURN_FREQ, None
                )
            ],
        )

    ### Calculate the total number of sequences
    # Create a counter
    total_col_seqs: Counter[str | DIGESTION_ERROR] = Counter()
    for row_index in range(0, align_array.shape[0]):
        if align_array[row_index, start_col] == "":
            total_col_seqs.update([DIGESTION_ERROR.END_OF_SEQUENCE])
            continue

        # Check if this row starts on a gap, and if so update the counter and skip
        if align_array[row_index, start_col] == "-":
            total_col_seqs.update([DIGESTION_ERROR.GAP_ON_SET_BASE])
            continue

        start_array = align_array[
            row_index, start_col : start_col + config.primer_size_min
        ]
        if "" in start_array:
            total_col_seqs.update([DIGESTION_ERROR.END_OF_SEQUENCE])
            continue
        if "-" in start_array:
            total_col_seqs.update([DIGESTION_ERROR.GAP_ON_SET_BASE])
            continue
        start_seq = "".join(start_array)

        if not start_seq:  # If the start seq is empty go to the next row
            continue

        start_seq_bases = set(start_seq)

        # Prevent Ns from being added
        if "N" in start_seq_bases:
            total_col_seqs.update([DIGESTION_ERROR.CONTAINS_INVALID_BASE])
            continue

        # Check for Non DNA bases
        if start_seq_bases - ALL_BASES:
            total_col_seqs.update([DIGESTION_ERROR.CONTAINS_INVALID_BASE])
            continue

        # expand any ambs
        if AMB_BASES & start_seq_bases:
            expanded_start_seq = expand_ambs([start_seq])
            assert expanded_start_seq is not None
        else:
            expanded_start_seq = [start_seq]

        results = []
        for start_seq in expanded_start_seq:
            # Get all sequences
            results.extend(
                wrap_walk(
                    walk_right,
                    array=align_array,
                    col_index_right=start_col + config.primer_size_min,
                    col_index_left=start_col,
                    row_index=row_index,
                    seq_str=start_seq,
                    config=config,
                )
            )

        # If all mutations matter, return on any Error
        if min_freq == 0 and set(results) & ERROR_SET:
            return (
                start_col,
                [DIGESTION_RESULT(parse_error(set(results)), EARLY_RETURN_FREQ, None)],
            )

        # Add the results to the Counter
        total_col_seqs.update(
            {seq: 1 / len(results) for seq in parse_error_list(results)}
        )
    # parse total_col_seqs into DIGESTION_RESULT
    digestion_result_list = [
        DIGESTION_RESULT(seq=k, count=v) for k, v in total_col_seqs.items()
    ]

    return (start_col, digestion_result_list)


def process_results(
    digestion_results: list[DIGESTION_RESULT],
    min_freq,
    ignore_n: bool = False,
) -> DIGESTION_ERROR | list[DIGESTION_RESULT]:
    """
    Takes the output from *_digest_to_count and returns a set of valid sequences. Or the error that occurred.

    Args:
        col (int): The column number.
        seq_counts (dict[str | DIGESTION_ERROR, int]): A dictionary containing sequence counts.
        min_freq: The minimum frequency threshold.

    Returns:
        DIGESTION_ERROR | list[DIGESTION_RESULT]: either an error or a dictionary of parsed sequences.
    """
    # Check for early return conditions
    for dr in digestion_results:
        if dr.count == EARLY_RETURN_FREQ and isinstance(dr.seq, DIGESTION_ERROR):
            return dr.seq

    # Remove Ns if asked
    if ignore_n:
        digestion_results = [
            dr
            for dr in digestion_results
            if dr.seq != DIGESTION_ERROR.CONTAINS_INVALID_BASE
        ]

    # Terminal padding is missing coverage. It neither contributes to the
    # frequency denominator nor participates in primer acceptance/rejection.
    covered_results = [
        dr for dr in digestion_results if dr.seq != DIGESTION_ERROR.END_OF_SEQUENCE
    ]
    total_values = sum(dr.count for dr in covered_results)

    if total_values == 0:
        return []

    results_above_freq: list[DIGESTION_RESULT] = [
        dr for dr in covered_results if dr.count / total_values >= min_freq
    ]

    # Check for Digestion Errors above the threshold
    for dr in results_above_freq:
        if isinstance(dr.seq, DIGESTION_ERROR):
            return dr.seq

    return results_above_freq


def r_digest_index(
    align_array: np.ndarray,
    config: Config,
    start_col: int,
    min_freq: float,
    early_return=False,
) -> RKmer | tuple[int, DIGESTION_ERROR | THERMO_RESULT]:
    """
    This will try and create a RKmer started at the given index
    :align_array: The alignment array
    :config: The configuration object
    :start_col: The column index to start the RKmer
    :min_freq: The minimum frequency threshold

    :return: A RKmer object or a tuple of (start_col, error)
    """
    # Count how many times each sequence / error occurs
    _start_col, digestion_results = r_digest_to_result(
        align_array, config, start_col, min_freq
    )
    parsed_digestion_results = process_results(
        digestion_results, min_freq, ignore_n=config.ignore_n
    )
    if isinstance(parsed_digestion_results, DIGESTION_ERROR):
        return (start_col, parsed_digestion_results)

    if not parsed_digestion_results:
        return (start_col, DIGESTION_ERROR.NO_SEQUENCES)

    # Rc the sequences
    for dr in parsed_digestion_results:
        dr.seq = reverse_complement(dr.seq)  # type: ignore

    # Thermo check the results
    for dr in parsed_digestion_results:
        if dr.thermo_check(config) is not THERMO_RESULT.PASS and not early_return:
            return (start_col, dr.status)  # type: ignore

    # Check for dimer
    seqs = [dr.seq.encode() for dr in parsed_digestion_results]  # type: ignore
    if do_pool_interact(seqs, seqs, config.dimer_score):
        return (start_col, DIGESTION_ERROR.DIMER_FAIL)

    # All checks pass return the kmer
    return RKmer(seqs, start_col, [dr.count for dr in parsed_digestion_results])


def f_digest_to_result(
    align_array: np.ndarray,
    config: Config,
    end_col: int,
    min_freq: float,
) -> tuple[int, list[DIGESTION_RESULT]]:
    """
    This will try and create a FKmer ended at the given index
    :return: A FKmer object or a tuple of (end_col, error)
    """

    # Check for gap frequency on first base
    base, counts = np.unique(align_array[:, end_col - 1], return_counts=True)
    first_base_counter = dict(zip(base, counts, strict=False))
    first_base_counter.pop("", None)

    num_seqs = sum(first_base_counter.values())

    first_base_freq = {k: v / num_seqs for k, v in first_base_counter.items()}

    # If the freq of gap is above minfreq
    gap_freq = first_base_freq.get("-")
    if gap_freq is not None and gap_freq >= min_freq:
        return (
            end_col,
            [
                DIGESTION_RESULT(
                    DIGESTION_ERROR.GAP_ON_SET_BASE, EARLY_RETURN_FREQ, None
                )
            ],
        )

    # If the initial slice is outside the range of the array
    if end_col - config.primer_size_min < 0:
        return (
            end_col,
            [DIGESTION_RESULT(DIGESTION_ERROR.WALKS_OUT, EARLY_RETURN_FREQ, None)],
        )

    total_col_seqs: Counter[str | DIGESTION_ERROR] = Counter()
    for row_index in range(0, align_array.shape[0]):
        if align_array[row_index, end_col - 1] == "":
            total_col_seqs.update([DIGESTION_ERROR.END_OF_SEQUENCE])
            continue

        # Check if this row starts on a gap, and if so update the counter and skip
        if align_array[row_index, end_col - 1] == "-":
            total_col_seqs.update([DIGESTION_ERROR.GAP_ON_SET_BASE])
            # Skip to next row
            continue

        start_array = align_array[row_index, end_col - config.primer_size_min : end_col]
        if "" in start_array:
            total_col_seqs.update([DIGESTION_ERROR.END_OF_SEQUENCE])
            continue
        if "-" in start_array:
            total_col_seqs.update([DIGESTION_ERROR.GAP_ON_SET_BASE])
            continue
        start_seq = "".join(start_array)

        if not start_seq:  # If the start seq is empty go to the next row
            continue

        start_seq_bases = set(start_seq)

        # Prevent Ns from being added
        if "N" in start_seq_bases:
            total_col_seqs.update([DIGESTION_ERROR.CONTAINS_INVALID_BASE])
            continue

        # Check for Non DNA bases
        if start_seq_bases - ALL_BASES:
            total_col_seqs.update([DIGESTION_ERROR.CONTAINS_INVALID_BASE])
            continue

        # expand any ambs
        if AMB_BASES & start_seq_bases:
            expanded_start_seq = expand_ambs([start_seq])
            assert expanded_start_seq is not None
        else:
            expanded_start_seq = [start_seq]

        results = []
        for start_seq in expanded_start_seq:
            results.extend(
                wrap_walk(
                    walk_left,
                    array=align_array,
                    col_index_right=end_col,
                    col_index_left=end_col - config.primer_size_min,
                    row_index=row_index,
                    seq_str=start_seq,
                    config=config,
                )
            )
        # Early return if all errors matter
        if min_freq == 0 and set(results) & ERROR_SET:
            return (
                end_col,
                [DIGESTION_RESULT(parse_error(set(results)), EARLY_RETURN_FREQ, None)],
            )

        # Add the results to the Counter
        total_col_seqs.update(
            {seq: 1 / len(results) for seq in parse_error_list(results)}
        )

    # Parse total_col_seqs into DIGESTION_RESULT
    digestion_result_list = [
        DIGESTION_RESULT(seq=k, count=v) for k, v in total_col_seqs.items()
    ]

    return (end_col, digestion_result_list)


def f_digest_index(
    align_array: np.ndarray,
    config: Config,
    end_col: int,
    min_freq: float,
    early_return=False,
) -> FKmer | tuple[int, DIGESTION_ERROR | THERMO_RESULT]:
    """
    This will try and create a FKmer ended at the given index
    :align_array: The alignment array
    :config: The configuration object
    :end_col: The column index to end the FKmer
    :min_freq: The minimum frequency threshold

    :return: A FKmer object or a tuple of (end_col, error)
    """

    # Count how many times each sequence / error occurs
    _end_col, digestion_results = f_digest_to_result(
        align_array, config, end_col, min_freq
    )

    parsed_digestion_results = process_results(
        digestion_results, min_freq, ignore_n=config.ignore_n
    )

    if isinstance(parsed_digestion_results, DIGESTION_ERROR):
        return (end_col, parsed_digestion_results)

    # Thermo check the results
    for dr in parsed_digestion_results:
        if dr.thermo_check(config) is not THERMO_RESULT.PASS and not early_return:
            return (end_col, dr.status)  # type: ignore

    # Check for dimer
    seqs = [dr.seq.encode() for dr in parsed_digestion_results]  # type: ignore
    if do_pool_interact(seqs, seqs, config.dimer_score):
        return (end_col, DIGESTION_ERROR.DIMER_FAIL)

    if not parsed_digestion_results:
        return (end_col, DIGESTION_ERROR.NO_SEQUENCES)

    return FKmer(seqs, end_col, [dr.count for dr in parsed_digestion_results])


def hamming_dist(s1, s2) -> int:
    """
    Return the number of substitutions, starting from the 3p end
    """
    return sum((x != y for x, y in zip(s1[::-1], s2[::-1], strict=False)))


def f_digest(
    msa_array: np.ndarray, config: Config, findexes: list[int], logger
) -> list[FKmer]:
    fkmers = []
    for findex in findexes:
        fkmer = f_digest_index(msa_array, config, findex, config.min_base_freq)

        # Append valid FKmers
        if isinstance(fkmer, FKmer) and fkmer.all_seqs():  # type: ignore
            fkmers.append(fkmer)

        # Log the Digestion
        if logger is not None:
            if isinstance(fkmer, tuple):
                logger.debug(f"FKmer: [red]{fkmer[0]}[/red]\t{fkmer[1].value}")
            else:
                logger.debug(f"FKmer: [green]{fkmer.end}[/green]: AllPass")

    return fkmers


def r_digest(
    msa_array: np.ndarray, config: Config, rindexes: list[int], logger
) -> list[RKmer]:
    rkmers = []
    for rindex in rindexes:
        rkmer = r_digest_index(msa_array, config, rindex, config.min_base_freq)

        # Append valid RKmers
        if isinstance(rkmer, RKmer) and rkmer.all_seqs():  # type: ignore
            rkmers.append(rkmer)

        # Log the Digestion
        if logger is not None:
            if isinstance(rkmer, tuple):
                logger.debug(f"RKmer: [red]{rkmer[0]}[/red]\t{rkmer[1].value}")
            else:
                logger.debug(f"RKmer: [green]{rkmer.start}[/green]: AllPass")
    return rkmers


def digest(
    msa_array: np.ndarray,
    config: Config,
    progress_manager: ProgressManager,
    indexes: tuple[list[int], list[int]] | None = None,
    logger=None,
    chrom: str = "",
) -> tuple[list[FKmer], list[RKmer]]:
    """
    Digest the given MSA array and return the FKmers and RKmers.

    :param msa_array: The input MSA array.
    :param cfg: A dictionary containing configuration parameters.
    :param indexes: A tuple of MSA indexes for (FKmers, RKmers), or None to use all indexes.
    :param logger: None or the logger object.
    :return: A tuple containing lists of sorted FKmers and RKmers.
    """
    # Guard for invalid indexes
    if indexes is not None:
        if min(indexes[0]) < 0 or max(indexes[0]) > msa_array.shape[1]:
            raise IndexError("FIndexes are out of range")
        if min(indexes[1]) < 0 or max(indexes[1]) >= msa_array.shape[1]:
            raise IndexError("RIndexes are out of range")

    # Get the indexes to digest
    findexes = (
        indexes[0]
        if indexes is not None
        else range(config.primer_size_min, msa_array.shape[1] + 1)
    )
    rindexes = (
        indexes[1]
        if indexes is not None
        else range(msa_array.shape[1] - config.primer_size_min + 1)
    )

    # Digest the findexes
    fkmers = []
    pt = progress_manager.create_sub_progress(
        iter=findexes, process="Creating forward primers", chrom=chrom
    )
    for findex in pt:
        fkmer = f_digest_index(msa_array, config, findex, config.min_base_freq)

        if logger is not None:
            if isinstance(fkmer, tuple):
                logger.debug(f"{chrom}:FKmer: {fkmer[0]}\t{fkmer[1]}")
            else:
                logger.debug(f"{chrom}:FKmer: {fkmer.end}\t{THERMO_RESULT.PASS}")

        # Append valid FKmers
        if isinstance(fkmer, FKmer) and fkmer.seqs_bytes():  # type: ignore
            fkmers.append(fkmer)

        # Update the count
        pt.manual_update(count=len(fkmers))

    # Digest the rindexes
    rkmers = []
    pt = progress_manager.create_sub_progress(
        iter=rindexes, process="Creating reverse primers", chrom=chrom
    )
    for rindex in pt:
        rkmer = r_digest_index(msa_array, config, rindex, config.min_base_freq)

        if logger is not None:
            if isinstance(rkmer, tuple):
                logger.debug(f"{chrom}:RKmer: {rkmer[0]}\t{rkmer[1]}")
            else:
                logger.debug(f"{chrom}:RKmer: {rkmer.start}\t{THERMO_RESULT.PASS}")

        # Append valid RKmers
        if isinstance(rkmer, RKmer) and rkmer.seqs_bytes():  # type: ignore
            rkmers.append(rkmer)

        # Update the count
        pt.manual_update(count=len(rkmers))

    return (fkmers, rkmers)


# Variant discovery deliberately does not call process_results or the cloud
# digesters: those routines discard row identity and chemically useful siblings.
def variant_digest_index(array, config, direction, index, *, expansion_limit=256, length_mode="first-compatible"):
    """Return primitive row-local evidence; no distinct-member dimer gate.

    Enumerate a measured length prefix (or all lengths), walking across internal
    gaps. Missing cells stop the walk. IUPAC expansions are bounded explicitly;
    omitted branches are recorded, never silently treated as observed sequence.
    """
    from itertools import product
    from math import prod
    from primalscheme3.core.config import ALL_DNA_WITH_N
    from primalscheme3.core.variant_thermo import variant_measurements

    if length_mode not in ('first-compatible', 'all'):
        raise ValueError('discovery length mode must be first-compatible or all')
    if direction not in ('f', 'r'):
        raise ValueError('Invalid discovery task direction')
    records = []
    for row_index, row in enumerate(array):
        selected = []
        cursor = index - 1 if direction == 'f' else index
        step = -1 if direction == 'f' else 1
        walked = 0
        for length in range(1, config.primer_size_max + 1):
            failure = None
            while True:
                if not 0 <= cursor < len(row) or row[cursor] == '':
                    failure = 'unavailable-footprint'
                    break
                if walked >= config.primer_max_walk:
                    failure = 'maximum-alignment-walk'
                    break
                cell = str(row[cursor]).upper()
                walked += 1
                column = cursor
                cursor += step
                if cell == '-':
                    continue
                if cell not in ALL_DNA_WITH_N:
                    failure = 'invalid-base'
                    break
                selected.append(column)
                break
            base = dict(direction=direction, index=index, row_index=row_index,
                        alignment_footprint=tuple(sorted(selected)), length=length)
            if failure:
                records.append(dict(base, sequence=None, accepted=False, reason=failure, checks={}))
                break
            if length < config.primer_size_min:
                continue
            cells = [str(row[i]).upper() for i in sorted(selected)]
            expansions = prod(len(ALL_DNA_WITH_N[x]) for x in cells)
            if expansions > expansion_limit:
                records.append(dict(base, sequence=None, accepted=False,
                                    reason='ambiguous-expansion-limit', checks={},
                                    expansion_count=expansions, expansion_limit=expansion_limit))
                continue
            compatible = False
            for bases in product(*(ALL_DNA_WITH_N[x] for x in cells)):
                seq = ''.join(bases)
                if direction == 'r':
                    seq = reverse_complement(seq)
                checks = variant_measurements(seq, config)
                accepted = all(x['outcome'] == 'pass' for k, x in checks.items() if k != 'self-dimer')
                compatible |= accepted
                records.append(dict(base, sequence=seq,
                                    accepted=accepted,
                                    reason='profile-chemistry', checks=checks,
                                    support='unknown' if expansions > 1 else 'confirmed'))
            if compatible and length_mode == 'first-compatible':
                records.append(dict(base, sequence=None, accepted=False, checks={},
                                    reason='first-compatible-length-stop',
                                    unexamined_length_interval=(length + 1, config.primer_size_max + 1)))
                break
    # Frequency is a per-sequence enumeration policy, never an error quorum.
    counts = Counter((r['sequence'], r['length']) for r in records if r['sequence'])
    for record in records:
        if record['sequence']:
            frequency = counts[record['sequence'], record['length']] / len(array)
            passes = frequency >= config.min_base_freq
            record['frequency'] = frequency
            record['checks']['minimum-frequency'] = dict(value=frequency, minimum=config.min_base_freq,
                                                         outcome='pass' if passes else 'fail')
            record['accepted'] &= passes
    return records
