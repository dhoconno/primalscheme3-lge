"""Tiny synthetic arrays verify missing-data accounting, not assay performance."""
import unittest
from unittest.mock import patch
import numpy as np
from primalscheme3.core.config import Config
from primalscheme3.core.digestion import (
    DIGESTION_ERROR, DIGESTION_RESULT, process_results,
    r_digest_to_result, f_digest_to_result, walk_left, walk_right,
)
from primalscheme3.core.errors import EndOfSequence, GapOnSetBase


class TerminalObservationTests(unittest.TestCase):
    def setUp(self):
        self.config = Config(terminal_gap_policy='observed-only')
        self.config.primer_size_min = 4
        self.config.primer_tm_min = -100

    def test_terminal_absence_excluded_from_denominator_both_directions(self):
        complete = list('ACGTACGTACGT')
        for function, index, partial in [
            (r_digest_to_result, 4, complete[:6] + [''] * 6),
            (f_digest_to_result, 7, [''] * 5 + complete[5:]),
        ]:
            _, observations = function(np.array([complete, complete, partial]), self.config, index, .8)
            processed = process_results(observations, .8)
            self.assertIsInstance(processed, list)
            self.assertEqual([(d.seq, d.count) for d in processed], [('ACGT' if index == 4 else 'TACG', 2)])
            self.assertEqual(sum(d.count for d in observations), 3)

    def test_internal_gap_remains_an_observation(self):
        for function, index in [(r_digest_to_result, 4), (f_digest_to_result, 8)]:
            rows = [list('ACGTACGTACGT'), list('ACGTA-GTACGT'), [''] * 12]
            _, observations = function(np.array(rows), self.config, index, .5)
            self.assertEqual(process_results(observations, .5), DIGESTION_ERROR.GAP_ON_SET_BASE)

    def test_all_missing_has_no_candidates_and_no_division_by_zero(self):
        rows = np.array([[''] * 12, [''] * 12])
        for function, index in [(r_digest_to_result, 4), (f_digest_to_result, 8)]:
            _, observations = function(rows, self.config, index, .8)
            self.assertEqual(process_results(observations, .8), [])

    def test_observed_variant_counts_and_threshold_equality(self):
        values = [DIGESTION_RESULT('ACGT', 1), DIGESTION_RESULT('TCGT', 1), DIGESTION_RESULT(DIGESTION_ERROR.END_OF_SEQUENCE, 8)]
        self.assertEqual([(x.seq, x.count) for x in process_results(values, .5)], [('ACGT', 1), ('TCGT', 1)])

    def test_walk_never_imputes_terminal_missing_or_erases_internal_gap(self):
        for symbol, error in [('', EndOfSequence), ('-', GapOnSetBase)]:
            rows = np.array([list('ACGTACGTACGT')])
            rows[0, 5] = symbol
            with patch('primalscheme3.core.digestion.calc_tm', return_value=-200):
                with self.assertRaises(error):
                    walk_right(rows, 5, 1, 0, 'CGTA', self.config)
                with self.assertRaises(error):
                    walk_left(rows, 10, 6, 0, 'GTAC', self.config)
