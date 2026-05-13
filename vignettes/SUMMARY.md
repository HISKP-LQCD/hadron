# Folder summary: `vignettes/`

## Folder purpose

Original R package vignettes and bibliography.

## Porting relevance

Map each workflow to a Python example, notebook, or documentation page, preserving reference-output checkpoints.

## Inventory table

| Category | Count | Notes |
|---|---|---|
| Files excluding summary | 9 | 8 Rmd plus 1 bibliography |
| R Markdown vignettes | 8 | Workflow docs |
| Bibliography files | 1 | Reference metadata |

## Mapping table

| Status | Source path | Workflow topic | Proposed Python target | Owner | Notes |
|---|---|---|---|---|---|
| `[ ]` | vignettes/Two_Amplitudes_Model.Rmd | R vignette: Two Amplitudes Model | docs/examples/two_amplitudes_model.md or notebooks/two_amplitudes_model.ipynb | unassigned | convert after owning APIs are ported |
| `[ ]` | vignettes/gevp.Rmd | R vignette: gevp | docs/examples/gevp.md or notebooks/gevp.ipynb | unassigned | convert after owning APIs are ported |
| `[ ]` | vignettes/hankel.Rmd | R vignette: hankel | docs/examples/hankel.md or notebooks/hankel.ipynb | unassigned | convert after owning APIs are ported |
| `[ ]` | vignettes/hankel.bib | bibliography for Hankel/PGEVM docs | docs/examples/hankel.md or notebooks/hankel.ipynb | unassigned | convert after owning APIs are ported |
| `[ ]` | vignettes/jackknife_cov_and_missing_values.Rmd | R vignette: jackknife cov and missing values | docs/examples/jackknife_cov_and_missing_values.md or notebooks/jackknife_cov_and_missing_values.ipynb | unassigned | convert after owning APIs are ported |
| `[ ]` | vignettes/jackknife_error_normalization.Rmd | R vignette: jackknife error normalization | docs/examples/jackknife_error_normalization.md or notebooks/jackknife_error_normalization.ipynb | unassigned | convert after owning APIs are ported |
| `[ ]` | vignettes/multi_particle_fit_test.Rmd | R vignette: multi particle fit test | docs/examples/multi_particle_fit_test.md or notebooks/multi_particle_fit_test.ipynb | unassigned | convert after owning APIs are ported |
| `[ ]` | vignettes/pgevm.Rmd | R vignette: pgevm | docs/examples/pgevm.md or notebooks/pgevm.ipynb | unassigned | convert after owning APIs are ported |
| `[ ]` | vignettes/truncated_pgevm.Rmd | R vignette: truncated pgevm | docs/examples/truncated_pgevm.md or notebooks/truncated_pgevm.ipynb | unassigned | convert after owning APIs are ported |

## Checklist

- [ ] Convert vignettes only after dependent APIs are stable.
- [ ] Use vignette examples as integration-test candidates.
- [ ] Keep bibliography references for migrated spectroscopy documentation.
- [ ] Record any workflow intentionally left R-only with `[X]` and reason.

## Testing requirements

| Source item | Required Python test or validation | Reference data needed | Owner/status | Notes |
|---|---|---|---|---|
| Rmd vignettes | example execution or notebook smoke test with fixture parity | yes where numerical | `[ ]` / unassigned | 8 workflows |
| hankel.bib | documentation build/reference check | no | `[ ]` / unassigned | bibliography only |

## Known issues and ambiguities

| Source path | Issue | Proposed resolution | Owner/status |
|---|---|---|---|
