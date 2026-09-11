"""Process scheduling for LGE Python discovery; biological calculations stay in digestion.

Each spawned worker owns an alignment copy. Only primitive records cross process
boundaries; native kmer objects are reconstructed in stable input position order.
"""
import multiprocessing
import os
import signal
import threading
from contextlib import contextmanager

from primalschemers import FKmer, RKmer

from primalscheme3.core.digestion import f_digest_index, r_digest_index
from primalscheme3.core.progress_tracker import ProgressManager

CHUNK_SIZE = 128
_worker_array = None
_worker_config = None


def resolve_worker_count(requested, chunk_count, available=None):
    if requested < 1:
        raise ValueError('ncores must be positive')
    available = available if available is not None else (os.cpu_count() or 1)
    return max(1, min(requested, max(1, available), max(1, chunk_count)))


def _initialize_worker(array, config, inherited_mask=None):
    global _worker_array, _worker_config
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    signal.signal(signal.SIGTERM, signal.SIG_DFL)
    if inherited_mask is not None:
        signal.pthread_sigmask(signal.SIG_SETMASK, inherited_mask)
    array.flags.writeable = False
    _worker_array, _worker_config = array, config


def evaluate_task(array, config, task):
    direction, index = task
    if direction == 'f':
        result = f_digest_index(array, config, index, config.min_base_freq)
    elif direction == 'r':
        result = r_digest_index(array, config, index, config.min_base_freq)
    else:
        raise ValueError('Invalid discovery task direction')
    if isinstance(result, tuple):
        error = result[1]
        return ('error', direction, index, type(error).__name__, error.value)
    # Python set/hash order differs between spawned interpreters. Canonicalize
    # alternatives in both serial and parallel execution, retaining their counts.
    alternatives = sorted(zip(result.seqs_bytes(), result.counts(), strict=True))
    return ('candidate', direction, index, tuple(x[0] for x in alternatives), tuple(x[1] for x in alternatives))


def _evaluate_chunk(chunk):
    return [evaluate_task(_worker_array, _worker_config, task) for task in chunk]


@contextmanager
def _termination_cleanup():
    # The app may terminate the parent process rather than its process group.
    # Unwind the Pool context so children are terminated and joined in either case.
    can_handle = threading.current_thread() is threading.main_thread()
    previous = signal.getsignal(signal.SIGTERM) if can_handle else None
    def terminate(signum, frame):
        raise KeyboardInterrupt('Discovery terminated')
    if can_handle:
        signal.signal(signal.SIGTERM, terminate)
    try:
        yield
    finally:
        if can_handle:
            signal.signal(signal.SIGTERM, previous)


def pool_results(chunks, array, config, workers):
    with _termination_cleanup():
        pool = None
        try:
            # Do not interrupt Pool.__init__ between spawning its first child
            # and assigning ownership. Workers restore the inherited mask only
            # after installing their own signal handlers.
            blocked = {signal.SIGTERM, signal.SIGINT}
            old_mask = signal.pthread_sigmask(signal.SIG_BLOCK, blocked) if hasattr(signal, 'pthread_sigmask') else None
            try:
                pool = multiprocessing.get_context('spawn').Pool(
                    workers, initializer=_initialize_worker,
                    initargs=(array, config, old_mask)
                )
            finally:
                if old_mask is not None:
                    signal.pthread_sigmask(signal.SIG_SETMASK, old_mask)
            yield from pool.imap(_evaluate_chunk, chunks, chunksize=1)
        finally:
            if pool is not None:
                pool.terminate()
                pool.join()


def discover(array, config, progress_manager=None, indexes=None, logger=None, chrom=''):
    width = array.shape[1]
    if indexes is None:
        forward = range(config.primer_size_min, width + 1)
        reverse = range(max(0, width - config.primer_size_min + 1))
    else:
        forward, reverse = indexes
        if any(i < config.primer_size_min or i > width for i in forward):
            raise IndexError('FIndexes are out of range')
        if any(i < 0 or i + config.primer_size_min > width for i in reverse):
            raise IndexError('RIndexes are out of range')
    tasks = [('f', int(i)) for i in forward] + [('r', int(i)) for i in reverse]
    if not tasks:
        return ([], []), 0
    chunks = [tasks[i:i + CHUNK_SIZE] for i in range(0, len(tasks), CHUNK_SIZE)]
    workers = resolve_worker_count(config.ncores, len(chunks))
    records = (
        ([evaluate_task(array, config, task) for task in chunk] for chunk in chunks)
        if workers == 1 else pool_results(chunks, array, config, workers)
    )
    manager = progress_manager or ProgressManager()
    tracker = manager.create_sub_progress(iter=records, process='Python discovery', chrom=chrom, total=len(chunks))
    fkmers, rkmers = [], []
    try:
        for batch in tracker:
            for kind, direction, index, data, values in batch:
                if kind == 'error':
                    if logger:
                        logger.debug(f'{chrom}:{direction.upper()}Kmer: {index}\t{data}.{values}')
                    continue
                if not data:
                    continue
                if direction == 'f':
                    fkmers.append(FKmer(list(data), index, list(values)))
                else:
                    rkmers.append(RKmer(list(data), index, list(values)))
                if logger:
                    logger.debug(f'{chrom}:{direction.upper()}Kmer: {index}\tPASS')
            tracker.manual_update(count=len(fkmers) + len(rkmers))
    finally:
        # Close generator promptly if consumer/progress/logging raises.
        records.close()
        tracker.close()
    return (fkmers, rkmers), workers
