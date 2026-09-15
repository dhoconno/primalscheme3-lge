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
    if len(task) >= 3 and task[2] == 'variants':
        from primalscheme3.core.digestion import variant_digest_index
        return variant_digest_index(array, config, task[0], task[1], length_mode=task[3] if len(task) > 3 else 'first-compatible')
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
    # Python dispatches process signals on the main thread even when another
    # unblocked thread received them. Defer exceptions until Pool ownership is
    # established; a main-thread pthread mask alone cannot protect construction.
    can_handle = threading.current_thread() is threading.main_thread()
    signals = (signal.SIGTERM, signal.SIGINT)
    previous = {sig: signal.getsignal(sig) for sig in signals} if can_handle else {}
    deferred = True
    pending = set()

    def interrupt(signum, frame):
        nonlocal deferred
        if signum == signal.SIGINT and previous[signum] == signal.SIG_IGN:
            return
        if deferred:
            pending.add(signum)
            return
        # Further signals must not interrupt termination/join after this raises.
        deferred = True
        if signum == signal.SIGTERM:
            raise KeyboardInterrupt('Discovery terminated')
        if callable(previous[signum]):
            previous[signum](signum, frame)
        else:
            raise KeyboardInterrupt('Discovery interrupted')
        deferred = False

    def defer_interruptions(value):
        nonlocal deferred
        deferred = value
        if not deferred:
            queued = sorted(pending, reverse=True)  # Termination precedes SIGINT.
            pending.clear()
            for signum in queued:
                interrupt(signum, None)

    if can_handle:
        for sig in signals:
            signal.signal(sig, interrupt)
    try:
        yield defer_interruptions
    finally:
        if can_handle:
            for sig, handler in previous.items():
                signal.signal(sig, handler)


def pool_results(chunks, array, config, workers):
    with _termination_cleanup() as defer_interruptions:
        pool = None
        try:
            # Keep new workers protected until they install their own signal
            # handlers. The Python handler also defers parent interruptions.
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
            defer_interruptions(False)
            yield from pool.imap(_evaluate_chunk, chunks, chunksize=1)
        finally:
            defer_interruptions(True)
            if pool is not None:
                pool.terminate()
                pool.join()


def discover(array, config, progress_manager=None, indexes=None, logger=None, chrom='', *, variant_mode=False, discovery_length_mode='first-compatible'):
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
    if variant_mode:
        if discovery_length_mode not in ('first-compatible', 'all'):
            raise ValueError('discovery length mode must be first-compatible or all')
        tasks = [(direction, index, 'variants', discovery_length_mode) for direction, index in tasks]
    if not tasks:
        if variant_mode:
            return [], 0
        return ([], []), 0
    chunks = [tasks[i:i + CHUNK_SIZE] for i in range(0, len(tasks), CHUNK_SIZE)]
    workers = resolve_worker_count(config.ncores, len(chunks))
    records = (
        ([evaluate_task(array, config, task) for task in chunk] for chunk in chunks)
        if workers == 1 else pool_results(chunks, array, config, workers)
    )
    if variant_mode:
        try:
            return [record for batch in records for task_records in batch for record in task_records], workers
        finally:
            records.close()
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
