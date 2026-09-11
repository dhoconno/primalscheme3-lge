"""CLI/configuration contracts; no biological engine execution."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from typer.testing import CliRunner
from primalscheme3.cli import app
from primalscheme3.core.config import Config
from primalscheme3.core.msa import MSA


class PolicyTests(unittest.TestCase):
    def test_default_and_explicit_config_roundtrip(self):
        for policy, backend in [('legacy', 'rust-legacy'), ('observed-only', 'python-observed-only')]:
            config = Config(terminal_gap_policy=policy)
            data = json.loads(json.dumps(config.to_dict()))
            self.assertEqual(data['terminal_gap_policy'], policy)
            self.assertEqual(data['discovery_backend'], backend)
            self.assertEqual(Config(**data).to_dict(), data)
        self.assertEqual(Config().to_dict()['terminal_gap_policy'], 'legacy')

    def test_invalid_policy_rejected(self):
        with self.assertRaises(ValueError):
            Config(terminal_gap_policy='guess')

    def test_cli_passes_policy_to_both_native_workflows(self):
        with tempfile.TemporaryDirectory() as temp:
            msa = Path(temp) / 'opaque.fasta'
            msa.write_text('opaque unused file\n')
            for command, function in [('scheme-create', 'schemecreate'), ('panel-create', 'panelcreate')]:
                for policy in ['legacy', 'observed-only']:
                    with patch('primalscheme3.cli.' + function) as run:
                        result = CliRunner().invoke(app, [command, '--msa', str(msa), '--output', str(Path(temp) / 'output'), '--terminal-gap-policy', policy])
                    self.assertEqual(result.exit_code, 0, result.output + str(result.exception))
                    self.assertEqual(run.call_args.kwargs['config'].to_dict()['terminal_gap_policy'], policy)

    def test_observed_only_dispatches_python_without_rust(self):
        instance = object.__new__(MSA)
        instance.logger = None
        with patch.object(MSA, 'digest') as python:
            instance.digest_rs(Config(terminal_gap_policy='observed-only'))
            python.assert_called_once()

    def test_observed_only_rejects_unimplemented_experimental_modes(self):
        for setting in ['downsample', 'use_annealing']:
            with self.assertRaisesRegex(ValueError, 'observed-only'):
                Config(terminal_gap_policy='observed-only', **{setting: True})
