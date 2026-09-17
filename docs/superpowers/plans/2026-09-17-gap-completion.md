# Gap-completion implementation plan

1. Add a narrow `gap_completion_parent` hook to `panelcreate`; reject unsupported
   mode combinations before ordinary legacy selection starts.
2. Snapshot source MSAs and parent files, copy unique raw inputs, and validate
   parent reference/primer/trimmed BED identity against the current MSAs.
3. Reuse legacy digestion and `Multiplex` admission to enumerate the retained
   candidate pool once. Rank candidates by deterministic positive trimmed gain
   against the parent plus accepted follow-up union, then place each candidate
   in an independent follow-up pool.
4. Publish follow-up-only standard BED/reference outputs plus compact candidate,
   coverage, namespace, and parent manifests. Keep the resolved config small by
   referencing reports rather than embedding the candidate manifest.
5. Replay the published follow-up from fresh parser/MatchDB state, compare
   candidate IDs and coverage to the selection snapshot, and verify parent and
   source hashes did not change during the run.
6. Exercise pure selector geometry, tiny successful parent-to-follow-up workflow,
   corrupt-parent failure provenance, stale-output rejection, and full focused
   tests before any biological benchmark.
