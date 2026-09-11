import multiprocessing
import unittest
import numpy as np
from primalscheme3.core.config import Config
from primalscheme3.core.progress_tracker import ProgressManager
from primalscheme3.core.parallel_discovery import discover, resolve_worker_count, pool_results


class ParallelDiscoveryTests(unittest.TestCase):
    def test_worker_bounds(self):
        self.assertEqual(resolve_worker_count(99, 3, 8), 3)
        self.assertEqual(resolve_worker_count(99, 20, 4), 4)
        self.assertEqual(resolve_worker_count(1, 20, 4), 1)
        with self.assertRaises(ValueError):
            resolve_worker_count(0, 20, 4)

    def test_one_two_four_workers_match_nonempty_candidates(self):
        # Deterministic synthetic array, not a biological fixture.
        rng = np.random.default_rng(20260911)
        row = rng.choice(list('ACGT'), 380)
        array = np.array([row, row, row])
        array[2, :60] = ''
        array[2, -60:] = ''
        expected = None
        for cores in [1, 2, 4]:
            config = Config(terminal_gap_policy='observed-only', ncores=cores)
            result, workers = discover(array, config, ProgressManager(), chrom='synthetic')
            values = tuple(tuple((x.region(), tuple(x.seqs()), tuple(x.counts())) for x in direction) for direction in result)
            self.assertTrue(values[0])
            self.assertTrue(values[1])
            self.assertLessEqual(workers, cores)
            if expected is None:
                expected = values
            else:
                self.assertEqual(values, expected)

    def test_worker_exception_cleans_up_children(self):
        before = {x.pid for x in multiprocessing.active_children()}
        config = Config(terminal_gap_policy='observed-only', ncores=2)
        with self.assertRaises(ValueError):
            list(pool_results([[('invalid', 0)]], np.array([list('ACGT')]), config, 2))
        self.assertEqual({x.pid for x in multiprocessing.active_children()}, before)
