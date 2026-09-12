"""Size contracts on synthetic coordinates; no discovery or biological datasets."""
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

from click import UsageError
from primalschemers import FKmer, RKmer
from typer.testing import CliRunner

from primalscheme3.cli import app
from primalscheme3.core.bedfiles import create_amplicon_str
from primalscheme3.core.config import Config
from primalscheme3.core.digestion import generate_valid_primerpairs
from primalscheme3.core.msa import MSA
from primalscheme3.panel.panel_main import PanelRunModes, panelcreate
from primalscheme3.replace.replace import ReplaceRunModes, replace
from primalscheme3.scheme.scheme_main import schemecreate


class QuietProgress:
    def create_sub_progress(self, iter, **kwargs):
        return self.Tracker(iter)

    class Tracker:
        def __init__(self, items):
            self.items = items

        def __iter__(self):
            return iter(self.items)

        def manual_update(self, **kwargs):
            pass


def geometry():
    # Forward envelope [100, 125); reverse envelope length 30 despite alternatives.
    forward = FKmer([b'A' * 20, b'C' * 25], 125)
    reverse = [RKmer([b'A' * 20, b'C' * 30], start) for start in (219, 220, 320, 321)]
    return [forward], reverse


def bed_lengths(pairs):
    return [int(row.split('\t')[2]) - int(row.split('\t')[1])
            for row in create_amplicon_str(pairs).splitlines() if row.strip()]


class AmpliconConfigTests(unittest.TestCase):
    def test_metric_survives_json_and_old_resolved_ranges_stay_legacy(self):
        for metric in ('legacy-pairing', 'reference-span'):
            config = Config(amplicon_size=200, amplicon_size_min=150,
                            amplicon_size_max=250, amplicon_size_metric=metric)
            data = json.loads(json.dumps(config.to_dict()))
            self.assertEqual(data.get('amplicon_size_metric'), metric)
            self.assertEqual(Config(**data).to_dict(), data)
        old = Config(amplicon_size=200, amplicon_size_min=150, amplicon_size_max=250)
        self.assertEqual(old.to_dict().get('amplicon_size_metric'), 'legacy-pairing')

    def test_reference_span_resolves_each_omitted_bound(self):
        cases = [({}, (180, 220)), ({'amplicon_size_min': 150}, (150, 220)),
                 ({'amplicon_size_max': 250}, (180, 250))]
        for bounds, expected in cases:
            config = Config(amplicon_size=200, amplicon_size_metric='reference-span', **bounds)
            self.assertEqual((config.amplicon_size_min, config.amplicon_size_max), expected)

    def test_invalid_reference_ranges_and_unsupported_inputs_rejected(self):
        cases = [dict(amplicon_size_min=0), dict(amplicon_size_max=-1),
                 dict(amplicon_size_min=201), dict(amplicon_size_max=199),
                 dict(amplicon_size_min=250, amplicon_size_max=150),
                 dict(circular=True), dict(input_bedfile='existing.bed'), dict(mapping='consensus')]
        for options in cases:
            with self.subTest(options=options), self.assertRaises(ValueError):
                Config(amplicon_size=200, amplicon_size_metric='reference-span', **options)
        with self.assertRaises(ValueError):
            Config(amplicon_size_metric='unknown')


class AmpliconGeometryTests(unittest.TestCase):
    def test_reference_span_accepts_exact_boundaries_and_excludes_neighbors(self):
        forward, reverse = geometry()
        with patch('primalscheme3.core.digestion.do_pool_interact', return_value=False):
            pairs = generate_valid_primerpairs(forward, reverse, 150, 250, -26, 0,
                                               QuietProgress(), amplicon_size_metric='reference-span')
        # Existing ordering is descending reverse start for equal forward ends.
        self.assertEqual(bed_lengths(pairs), [250, 150])

    def test_reference_span_handles_variable_reverse_lengths(self):
        forward = [FKmer([b'A' * 20], 120)]  # Starts at 100.
        reverse = [RKmer([b'C' * length], start) for start, length in
                   ((210, 40), (230, 20), (231, 20), (300, 50), (301, 50))]
        with patch('primalscheme3.core.digestion.do_pool_interact', return_value=False):
            pairs = generate_valid_primerpairs(forward, reverse, 150, 250, -26, 0,
                                               QuietProgress(), amplicon_size_metric='reference-span')
        self.assertEqual(bed_lengths(pairs), [250, 151, 150, 150])

    def test_legacy_default_retains_previous_start_based_window(self):
        forward, reverse = geometry()
        with patch('primalscheme3.core.digestion.do_pool_interact', return_value=False):
            pairs = generate_valid_primerpairs(forward, reverse, 150, 250, -26, 0,
                                               QuietProgress())
        self.assertEqual(bed_lengths(pairs), [251, 250])

    def test_interaction_rejection_still_applies(self):
        forward, reverse = geometry()
        with patch('primalscheme3.core.digestion.do_pool_interact', return_value=True):
            pairs = generate_valid_primerpairs(forward, reverse, 150, 250, -26, 0,
                                               QuietProgress(), amplicon_size_metric='reference-span')
        self.assertEqual(pairs, [])


class AmpliconCLITests(unittest.TestCase):
    def test_create_flags_select_metric_and_resolved_bounds_for_both_backends(self):
        cases = [([], 'legacy-pairing', (180, 220)),
                 (['--amplicon-size-min', '150'], 'reference-span', (150, 220)),
                 (['--amplicon-size-max', '250'], 'reference-span', (180, 250)),
                 (['--amplicon-size-min', '150', '--amplicon-size-max', '250'],
                  'reference-span', (150, 250))]
        with tempfile.TemporaryDirectory() as temp:
            msa = Path(temp) / 'opaque.fasta'
            msa.write_text('unused opaque input\n')
            for command, function, mode in [('scheme-create', 'schemecreate', []),
                                            ('panel-create', 'panelcreate', ['--mode', 'equal']),
                                            ('panel-create', 'panelcreate', ['--mode', 'entropy'])]:
                for policy in ('legacy', 'observed-only'):
                    for flags, metric, bounds in cases:
                        with self.subTest(command=command, policy=policy, flags=flags), \
                                patch('primalscheme3.cli.' + function) as run:
                            args = [command, '--msa', str(msa), '--output', str(Path(temp) / 'out'),
                                    '--amplicon-size', '200', '--terminal-gap-policy', policy, *mode, *flags]
                            result = CliRunner().invoke(app, args)
                            self.assertEqual(result.exit_code, 0, result.output + str(result.exception))
                            config = run.call_args.kwargs['config']
                            self.assertEqual(config.to_dict().get('amplicon_size_metric'), metric)
                            self.assertEqual((config.amplicon_size_min, config.amplicon_size_max), bounds)

    def test_invalid_flags_fail_before_workflow_and_output_creation(self):
        with tempfile.TemporaryDirectory() as temp:
            msa = Path(temp) / 'opaque.fasta'
            msa.write_text('unused opaque input\n')
            existing = Path(temp) / 'existing.bed'
            existing.write_text('unused opaque input\n')
            for command in ('scheme-create', 'panel-create'):
                cases = [['--amplicon-size-min', '0'], ['--amplicon-size-max', '-1'],
                         ['--amplicon-size-min', '201'], ['--amplicon-size-max', '199'],
                         ['--amplicon-size-min', '250', '--amplicon-size-max', '150'],
                         ['--amplicon-size-min', '150', '--input-bedfile', str(existing)],
                         ['--amplicon-size-min', '150', '--mapping', 'consensus']]
                if command == 'scheme-create':
                    cases += [['--amplicon-size-min', '150', '--circular'],
                              ['--amplicon-size-min', '150', '--bedfile', str(existing)]]
                else:
                    cases += [['--amplicon-size-min', '150', '--mode', 'region-only'],
                              ['--amplicon-size-min', '150', '--region-bedfile', str(existing)]]
                for index, flags in enumerate(cases):
                    out = Path(temp) / f'{command}-{index}'
                    mode = ['--mode', 'equal'] if command == 'panel-create' else []
                    with self.subTest(command=command, flags=flags), \
                            patch('primalscheme3.cli.schemecreate', side_effect=AssertionError('validation reached workflow')), \
                            patch('primalscheme3.cli.panelcreate', side_effect=AssertionError('validation reached workflow')):
                        result = CliRunner().invoke(app, [command, '--msa', str(msa), '--output', str(out),
                                                       '--amplicon-size', '200', *mode, *flags])
                        self.assertEqual(result.exit_code, 2, result.output + str(result.exception))
                        self.assertFalse(out.exists())


class StopAfterPairs(Exception):
    pass


class GeometryMSA:
    """Replace only input parsing/discovery, then execute real shared pairing."""
    def __init__(self):
        self.fkmers, self.rkmers = geometry()
        self.primerpairs = []
        self.msa_index = 0
        self.name = self._chrom_name = self._uuid = 'synthetic'
        self.progress_manager = QuietProgress()
        self.array = SimpleNamespace(shape=(1, 500))
        self.regions = None

    def write_msa_to_file(self, path):
        path.write_text('unused opaque input\n')

    def digest_rs(self, *args, **kwargs):
        pass

    def add_regions(self, *args):
        pass

    def create_score_array(self, *args, **kwargs):
        pass

    def generate_primerpairs(self, **kwargs):
        MSA.generate_primerpairs(self, **kwargs)
        raise StopAfterPairs


class AmpliconWorkflowTests(unittest.TestCase):
    def test_scheme_and_panel_use_saved_metric_and_real_minimum(self):
        for workflow, module, msa_type, extra in [
            (schemecreate, 'primalscheme3.scheme.scheme_main', 'MSA', {}),
            (panelcreate, 'primalscheme3.panel.panel_main', 'PanelMSA', {'mode': PanelRunModes.EQUAL}),
            (panelcreate, 'primalscheme3.panel.panel_main', 'PanelMSA', {'mode': PanelRunModes.ENTROPY}),
        ]:
            for metric, expected in [('reference-span', [250, 150]), ('legacy-pairing', [251, 250])]:
                with self.subTest(workflow=workflow.__name__, extra=extra, metric=metric), tempfile.TemporaryDirectory() as temp:
                    instance = GeometryMSA()
                    config = Config(amplicon_size=200, amplicon_size_min=150, amplicon_size_max=250,
                                    amplicon_size_metric=metric, use_matchdb=False)
                    with patch(module + '.' + msa_type, return_value=instance), \
                            patch(module + '.MatchDB'), patch(module + '.setup_rich_logger', return_value=Mock()), \
                            patch('primalscheme3.core.digestion.do_pool_interact', return_value=False):
                        with self.assertRaises(StopAfterPairs):
                            workflow(msa=[Path(temp) / 'opaque.fasta'], output_dir=Path(temp) / 'out',
                                     config=config, pm=QuietProgress(), **extra)
                    self.assertEqual(bed_lengths(instance.primerpairs), expected)

    def test_direct_workflows_reject_unsupported_scope_before_outputs(self):
        config = Config(amplicon_size=200, amplicon_size_min=150,
                        amplicon_size_max=250, amplicon_size_metric='reference-span')
        cases = [(schemecreate, {'input_bedfile': Path('existing.bed')}),
                 (panelcreate, {'mode': PanelRunModes.REGION_ONLY}),
                 (panelcreate, {'mode': PanelRunModes.EQUAL, 'region_bedfile': Path('regions.bed')}),
                 (panelcreate, {'mode': PanelRunModes.EQUAL, 'input_bedfile': Path('existing.bed')})]
        for workflow, extra in cases:
            with self.subTest(extra=extra), tempfile.TemporaryDirectory() as temp:
                out = Path(temp) / 'out'
                module = workflow.__module__
                with patch(module + '.MatchDB'), \
                        patch(module + '.setup_rich_logger', return_value=Mock()), \
                        self.assertRaises(UsageError):
                    workflow(msa=[], output_dir=out, config=config, pm=QuietProgress(), **extra)
                self.assertFalse(out.exists())

    def test_replacement_cannot_silently_reinterpret_saved_reference_span(self):
        config = Config(amplicon_size=200, amplicon_size_min=150,
                        amplicon_size_max=250, amplicon_size_metric='reference-span')
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp) / 'out'
            with patch('primalscheme3.replace.replace.setup_rich_logger', return_value=Mock()), \
                    self.assertRaises(UsageError):
                replace(config, Path('existing.bed'), 'opaque', Path('opaque.fasta'),
                        QuietProgress(), out, False, ReplaceRunModes.ListAll)
            self.assertFalse(out.exists())
