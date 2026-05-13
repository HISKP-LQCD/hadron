# Folder summary: `tests/`

## Folder purpose

Original R tests and R Markdown validation workflows.

## Porting relevance

Map each parsed R test to pytest and identify which tests should emit frozen R reference outputs.

## Inventory table

| Category | Count | Notes |
|---|---|---|
| Test files | 13 | All parsed with R or Rmd chunk parsing |
| Parse errors | 0 | None |
| testthat R files | 11 | Includes `tests/testthat.R` runner |
| R Markdown test workflows | 2 | Code chunks parsed |

## Mapping table

| Status | Source path | R parse status | Current scope | Detected tested symbols | Proposed Python target | Owner | Notes |
|---|---|---|---|---|---|---|---|
| `[ ]` | tests/new_matrixfit_all_points.Rmd | ok | root documentation/workflow file | bootstrap.cf; matrixfit; new_matrixfit; residual_plot | tests/reference_workflows/test_new_matrixfit_all_points.py | unassigned | convert or use as reference fixture driver |
| `[ ]` | tests/single_constant_model.Rmd | ok | root documentation/workflow file | bootstrap.cf; bootstrap.effectivemass; new_matrixfit; residual_plot | tests/reference_workflows/test_single_constant_model.py | unassigned | convert or use as reference fixture driver |
| `[ ]` | tests/testthat.R | ok | root documentation/workflow file | - | tests/test_testthat.py | unassigned | convert or use as reference fixture driver |
| `[ ]` | tests/testthat/test_bootstrapfit.R | ok | y errors; y errors cov; xy errors; xy errors cov; y errors with priors; y errors cov with priors; xy errors with priors; xy errors cov with priors | parametric.bootstrap; bootstrap.nlsfit | tests/test_bootstrapfit.py | unassigned | convert or use as reference fixture driver |
| `[ ]` | tests/testthat/test_computeDisc.R | ok | cross_vs_diagonal | computeDisc | tests/test_compute_disc.py | unassigned | convert or use as reference fixture driver |
| `[ ]` | tests/testthat/test_dummy.R | ok | root documentation/workflow file | - | tests/test_dummy.py | unassigned | convert or use as reference fixture driver |
| `[ ]` | tests/testthat/test_extractSingleCor_cf.R | ok | extract_symmetrized | bootstrap.cf; mul.cf; extractSingleCor.cf | tests/test_extract_single_cor_cf.py | unassigned | convert or use as reference fixture driver |
| `[ ]` | tests/testthat/test_new_matrixfit_higher_states.R | ok | two states with cov; two states | new_matrixfit | tests/test_new_matrixfit_higher_states.py | unassigned | convert or use as reference fixture driver |
| `[ ]` | tests/testthat/test_new_matrixfit.R | ok | SingleModelJacobian; SingleModel; ShiftedModelPrediction; ShiftedModel; TwoStateModel | make_parlist; dmatrixChi; make_parind; make_sign_vec; make_ov_sign_vec; bootstrap.cf; matrixChi.shifted; dmatrixChi.shifted; takeTimeDiff.cf; bootstrap.gevp; gevp2cf | tests/test_new_matrixfit.py | unassigned | convert or use as reference fixture driver |
| `[ ]` | tests/testthat/test_parlist.R | ok | parlist_1; parlist_4; parind_1; parind_4 | make_parlist; make_parind | tests/test_parlist.py | unassigned | convert or use as reference fixture driver |
| `[ ]` | tests/testthat/test_removeTemporal.R | ok | equality | bootstrap.cf; matrixfit; old_removeTemporal.cf; removeTemporal.cf | tests/test_remove_temporal.py | unassigned | convert or use as reference fixture driver |
| `[ ]` | tests/testthat/test_string2error.R | ok | 1 | string2error | tests/test_string2error.py | unassigned | convert or use as reference fixture driver |
| `[ ]` | tests/testthat/test_tex_catwitherror.R | ok | small_error; borderline_error; another_borderline_error; even_nastier_borderline_error; scientific_notation; very_small_number; same_error; intermediate_error; large_error; similar_error; no_error; zero_error; zero_val_zero_err; NA; vector | tex.catwitherror | tests/test_tex_catwitherror.py | unassigned | convert or use as reference fixture driver |

## Checklist

- [ ] Convert each `testthat` file into a focused pytest module.
- [ ] Extract deterministic R outputs before implementing equivalent Python algorithms.
- [ ] Keep R Markdown test workflows as higher-level regression examples.
- [ ] Add fixture ownership when datasets or `inst/extdata` are needed.

## Testing requirements

| Source item | Required Python test or validation | Reference data needed | Owner/status | Notes |
|---|---|---|---|---|
| testthat files | pytest assertions plus R fixture parity where numerical | yes | `[ ]` / unassigned | R parser detected test names and called symbols |
| Rmd workflows | workflow smoke test or converted notebook/example | yes where numerical | `[ ]` / unassigned | 2 files |

## Known issues and ambiguities

| Source path | Issue | Proposed resolution | Owner/status |
|---|---|---|---|
