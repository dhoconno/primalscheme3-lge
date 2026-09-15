from primalscheme3.panel.coverage_types import Target


def target(*rows, id='target-test', aliases=None):
    cells = tuple(tuple(row) for row in rows)
    ref = cells[0]
    indexes = tuple(i for i, x in enumerate(ref) if x not in ('', '-'))
    mapping = tuple(indexes.index(i) if i in indexes else None for i in range(len(ref)))
    return Target(id, 0, 0, aliases or tuple(f'row-{i}' for i in range(len(rows))), cells,
                  ''.join(ref[i] for i in indexes), len(indexes), mapping,
                  indexes + ((indexes[-1] + 1) if indexes else 0,))


def literal_coverage(observed, interiors):
    recovered = set().union(*(set(range(a, b)) for a, b in interiors))
    return None if not observed else len(set(observed) & recovered) / len(observed)
