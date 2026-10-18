# Result feature slices

`ResultsView.vue` composes the canonical public summary into mode-specific tabs
and invokes request lifecycles exposed by `features/results/`. Tabs use roving
keyboard focus and keep inactive panels mounted but hidden so history selections
survive navigation. Result tabs do not modify the recovery URL. Workflow warnings
remain visible across tabs. Filtering, provenance, likelihood comparisons, and
history controls are optional disclosures; model scores are opt-in. Child directories own narrow
scientific or publication views:

- `overview/` — synopsis, detection, preflight, and evidence chain
- `orientation/` — likelihood comparison and orientation decision
- `reconstruction/` — selected history, events, model, and spacer inventory
- `history/` — tree layout, ancestral nodes, exact events, and interaction
- `publication/` — provenance and capability-protected downloads

Keep imports directed from the composition root into slices. A slice may use
shared utilities and common UI, but should not import another slice's internals
except the explicit reconstruction-to-history composition.
