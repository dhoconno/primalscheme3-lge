# Independent legacy gap completion

## Goal

`panel-create --gap-completion-parent PATH` creates a new, independent follow-up
scheme for uncovered trimmed intervals while preserving the completed parent
scheme byte-for-byte. The follow-up uses the existing legacy candidate
generation, chemistry, specificity, and dimer admission rules. Follow-up pool
count is the existing `--n-pools` option (default two).

## Contract

The command accepts whole, linear, first-row mapped MSAs. Circular inputs,
imported primer BEDs, regions, salvage mode, and allele-coverage mode are
rejected explicitly. The parent must contain `primer.bed`,
`primertrim.amplicon.bed`, and `reference.fasta`; its target identities and
reference sequence are checked against the supplied MSAs. Parent files are
copied under `parent/` and their hashes are recorded.

Parent pairs seed only the coverage union. They are not inserted into
follow-up pools. Each follow-up candidate is admitted independently into one
of the configured pools using the native legacy overlap, dimer, and MatchDB
checks. A candidate is accepted only when it adds positive trimmed coverage to
the combined parent and follow-up union. Selection is deterministic and
retains full MSA coordinates, including primer binding flanks.

The normal `primer.bed`, `amplicon.bed`, `primertrim.amplicon.bed`, and
`reference.fasta` outputs contain follow-up pairs only. `gap-completion.json`
contains the parent/follow-up namespace, candidate status manifest, and fresh
validation result. `gap-completion-coverage.json` reports primary, follow-up,
combined, candidate-union ceiling, and remaining-gap coverage.

Every run writes ordinary scientific provenance with exact argv, resolved
options, runtime, copied raw inputs, hashes, sizes, exit status, wall time, and
stderr. Existing legacy panel creation is unchanged when the option is absent.
