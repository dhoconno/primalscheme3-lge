"""Independent measurements using the existing thermodynamic kernels."""
from copy import deepcopy
from functools import lru_cache
from primalschemers import calc_at_offset_py
from primalscheme3.core import thermo

_FIELDS = ('primer_size_min', 'primer_size_max', 'primer_gc_min', 'primer_gc_max',
           'primer_tm_min', 'primer_tm_max', 'primer_homopolymer_max', 'primer_hairpin_th_max',
           'mv_conc', 'dv_conc', 'dntp_conc', 'dna_conc', 'use_annealing',
           'primer_annealing_prop', 'primer_annealing_tempc')


def variant_measurements(sequence, config):
    return deepcopy(_measure(sequence, tuple(getattr(config, k) for k in _FIELDS)))


@lru_cache(maxsize=65536)
def _measure(sequence, values):
    c = dict(zip(_FIELDS, values, strict=True))
    salts = {k: c[k] for k in ('mv_conc', 'dv_conc', 'dntp_conc', 'dna_conc')}
    checks = {}
    def interval(name, value, lower, upper):
        checks[name] = dict(value=value, minimum=lower, maximum=upper,
                            outcome='pass' if lower <= value <= upper else 'fail')
    interval('length', len(sequence), c['primer_size_min'], c['primer_size_max'])
    interval('gc', thermo.gc(sequence), c['primer_gc_min'], c['primer_gc_max'])
    if c['use_annealing'] and c['primer_annealing_prop'] is not None:
        value = thermo.calc_annealing(sequence, **salts, temp_c=c['primer_annealing_tempc'])
        interval('annealing', value, c['primer_annealing_prop'] - thermo.ANNEALING_DIFF,
                 c['primer_annealing_prop'] + thermo.ANNEALING_DIFF)
    else:
        interval('tm', thermo.calc_tm(sequence, **salts), c['primer_tm_min'], c['primer_tm_max'] + 2)
    interval('homopolymer', thermo.max_homo(sequence), 0, c['primer_homopolymer_max'])
    interval('hairpin', thermo.calc_hairpin_tm(sequence, **salts), -1e9, c['primer_hairpin_th_max'])
    # Exactly the native boolean kernel offset range (self orientations coincide).
    # This is evidence only: stage policy
    # owns the cutoff and may reconsider a self edge during salvage.
    score = min(calc_at_offset_py(sequence, sequence, offset)
                for offset in range(-(len(sequence)-2), 0))
    checks['self-dimer'] = dict(value=score, outcome='not-evaluated', reason='deferred-to-stage-policy')
    return checks
