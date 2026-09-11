"""Python discovery integration with the installed native kmer API."""
import unittest
from unittest.mock import patch
import numpy as np
from primalschemers import FKmer, RKmer
from primalscheme3.core.config import Config
from primalscheme3.core.msa import MSA


class MappingCompatibilityTests(unittest.TestCase):
    def test_current_native_remap_api_and_reference_boundaries(self):
        msa = object.__new__(MSA)
        msa.array = np.array([list('ACGTACGT')])
        msa._mapping_array = np.array([None, None, 0, 1, 2, 3, None, None], dtype=object)
        msa.name = 'synthetic'
        msa.msa_index = 0
        msa.logger = None
        msa.progress_manager = None
        forward = [FKmer([b'AC'], 4, [2.0]), FKmer([b'ACGT'], 3, [2.0])]
        reverse = [RKmer([b'AC'], 4, [2.0]), RKmer([b'ACGT'], 5, [2.0])]
        with patch('primalscheme3.core.parallel_discovery.discover', return_value=((forward, reverse), 1)):
            msa.digest(Config(terminal_gap_policy='observed-only'))
        self.assertEqual([(x.end, x.counts()) for x in msa.fkmers], [(2, [2.0])])
        self.assertEqual([(x.start, x.counts()) for x in msa.rkmers], [(2, [2.0])])
