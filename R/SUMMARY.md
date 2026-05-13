# Folder summary: `R/`

## Folder purpose

Original R implementation and algorithmic reference for the Python port.

## Porting relevance

Every R file and R-parser-backed top-level symbol maps to a Python target. Public exports and S3 methods are cross-checked against `NAMESPACE`.

## Inventory table

| Category | Count | Notes |
|---|---|---|
| R source files | 64 | All parsed by R `parse()` |
| Top-level functions | 325 | R parser, top-level assignments only |
| R6 classes | 8 | All in `R/new_matrixfit.R` |
| NAMESPACE exports | 157 | All matched to parsed definitions |
| S3 methods | 48 | All matched to parsed definitions |
| UseMethod generics | 2 | extract_mass, residual_plot |
| Parse errors | 0 | None |

## Mapping table

| Status | Source path | Proposed Python target | Top-level functions | R6 classes | Exports | S3 methods | Owner | Notes |
|---|---|---|---|---|---|---|---|---|
| `[ ]` | R/alpha_s.R | src/pydron/alpha_s.py | 2 | 0 | 1 | 0 | unassigned | running coupling wrapper and R fallback |
| `[ ]` | R/analysis_gradient_flow.R | src/pydron/analysis_gradient_flow.py | 1 | 0 | 1 | 0 | unassigned | analysis gradient flow utilities |
| `[ ]` | R/analysis_online.R | src/pydron/analysis_online.py | 3 | 0 | 1 | 0 | unassigned | analysis online utilities |
| `[ ]` | R/block.R | src/pydron/block.py | 1 | 0 | 0 | 0 | unassigned | block utilities |
| `[ ]` | R/boot_ts_array.R | src/pydron/boot_ts_array.py | 1 | 0 | 1 | 0 | unassigned | boot ts array utilities |
| `[ ]` | R/bootstrap.nlsfit.R | src/pydron/bootstrap_nlsfit.py | 18 | 0 | 7 | 5 | unassigned | parametric/simple/bootstrap nonlinear fits |
| `[ ]` | R/bootstrapnumber.R | src/pydron/bootstrapnumber.py | 3 | 0 | 1 | 0 | unassigned | bootstrapnumber utilities |
| `[ ]` | R/cdh.R | src/pydron/cdh.py | 4 | 0 | 2 | 0 | unassigned | finite-size corrections with C fallback |
| `[ ]` | R/CExp.R | src/pydron/cexp.py | 2 | 0 | 1 | 0 | unassigned | CExp utilities |
| `[ ]` | R/cf.R | src/pydron/cf.py | 41 | 0 | 29 | 8 | unassigned | correlation-function container, resampling, arithmetic, summaries |
| `[ ]` | R/ChiSqr.R | src/pydron/chi_sqr.py | 8 | 0 | 0 | 0 | unassigned | ChiSqr utilities |
| `[ ]` | R/computeacf.R | src/pydron/computeacf.py | 3 | 0 | 1 | 2 | unassigned | autocorrelation computation and presentation |
| `[ ]` | R/computeDisc.R | src/pydron/compute_disc.py | 1 | 0 | 1 | 0 | unassigned | computeDisc utilities |
| `[ ]` | R/computefps.R | src/pydron/computefps.py | 2 | 0 | 2 | 0 | unassigned | computefps utilities |
| `[ ]` | R/correlatedRNG.R | src/pydron/correlated_rng.py | 2 | 0 | 0 | 0 | unassigned | correlatedRNG utilities |
| `[ ]` | R/cosh_nlsfit.R | src/pydron/cosh_nlsfit.py | 5 | 0 | 1 | 2 | unassigned | cosh nlsfit utilities |
| `[ ]` | R/cvc_readutils.R | src/pydron/cvc_readutils.py | 8 | 0 | 8 | 0 | unassigned | CVC correlator and loop readers/key builders |
| `[ ]` | R/cyprus_readutils.R | src/pydron/cyprus_readutils.py | 3 | 0 | 3 | 0 | unassigned | Cyprus-format loop readers/key builders |
| `[ ]` | R/deriv_utils.R | src/pydron/deriv_utils.py | 1 | 0 | 0 | 0 | unassigned | deriv utils utilities |
| `[ ]` | R/effectivemass.R | src/pydron/effectivemass.py | 8 | 0 | 3 | 4 | unassigned | effective-mass extraction/fits/plots |
| `[ ]` | R/fit.plateau2cf.R | src/pydron/fit_plateau2cf.py | 1 | 0 | 1 | 0 | unassigned | fit plateau2cf utilities |
| `[ ]` | R/fitmass.R | src/pydron/fitmass.py | 1 | 0 | 0 | 0 | unassigned | fitmass utilities |
| `[ ]` | R/fs.mpia0.R | src/pydron/fs_mpia0.py | 3 | 0 | 3 | 0 | unassigned | fs mpia0 utilities |
| `[ ]` | R/fs.R | src/pydron/fs.py | 6 | 0 | 0 | 0 | unassigned | fs utilities |
| `[ ]` | R/functional.R | src/pydron/functional.py | 1 | 0 | 1 | 0 | unassigned | functional utilities |
| `[ ]` | R/g1.R | src/pydron/g1.py | 1 | 0 | 1 | 0 | unassigned | g1 utilities |
| `[ ]` | R/gamma.R | src/pydron/gamma.py | 2 | 0 | 1 | 0 | unassigned | gamma utilities |
| `[ ]` | R/getCor.R | src/pydron/get_cor.py | 1 | 0 | 0 | 0 | unassigned | getCor utilities |
| `[ ]` | R/getfit.boot.R | src/pydron/getfit_boot.py | 1 | 0 | 0 | 0 | unassigned | getfit boot utilities |
| `[ ]` | R/getNxNmatrix.R | src/pydron/get_nx_nmatrix.py | 1 | 0 | 0 | 0 | unassigned | getNxNmatrix utilities |
| `[ ]` | R/gevp.R | src/pydron/gevp.py | 7 | 0 | 4 | 2 | unassigned | generalized eigenvalue analysis and amplitudes |
| `[ ]` | R/h5utils.R | src/pydron/h5utils.py | 2 | 0 | 1 | 0 | unassigned | h5utils utilities |
| `[ ]` | R/hadron-package.R | src/pydron/hadron_package.py | 0 | 0 | 0 | 0 | unassigned | package-level roxygen documentation |
| `[ ]` | R/hankel.R | src/pydron/hankel.py | 13 | 0 | 7 | 1 | unassigned | Hankel and PGEVM analysis |
| `[ ]` | R/hankel.truncated.R | src/pydron/hankel_truncated.py | 8 | 0 | 3 | 1 | unassigned | truncated Hankel/PGEVM analysis |
| `[ ]` | R/inv_cosh.R | src/pydron/inv_cosh.py | 1 | 0 | 1 | 0 | unassigned | inverse cosh native wrapper |
| `[ ]` | R/invertCovMatrix.R | src/pydron/invert_cov_matrix.py | 1 | 0 | 1 | 0 | unassigned | invertCovMatrix utilities |
| `[ ]` | R/jackknifeafterboot.R | src/pydron/jackknifeafterboot.py | 6 | 0 | 0 | 0 | unassigned | jackknifeafterboot utilities |
| `[ ]` | R/kappa.R | src/pydron/kappa.py | 2 | 0 | 0 | 0 | unassigned | kappa utilities |
| `[ ]` | R/lanczos.R | src/pydron/lanczos.py | 2 | 0 | 1 | 0 | unassigned | lanczos utilities |
| `[ ]` | R/legacy_functions.R | src/pydron/legacy_functions.py | 3 | 0 | 0 | 0 | unassigned | legacy functions utilities |
| `[ ]` | R/looptools.R | src/pydron/looptools.py | 5 | 0 | 5 | 0 | unassigned | looptools utilities |
| `[ ]` | R/LuescherMethod.R | src/pydron/luescher_method.py | 2 | 0 | 0 | 0 | unassigned | LuescherMethod utilities |
| `[ ]` | R/matrixfit.R | src/pydron/matrixfit.py | 23 | 0 | 3 | 2 | unassigned | legacy matrix-fit machinery |
| `[ ]` | R/momentum_utils.R | src/pydron/momentum_utils.py | 1 | 0 | 1 | 0 | unassigned | momentum utils utilities |
| `[ ]` | R/new_matrixfit.R | src/pydron/new_matrixfit.py | 5 | 8 | 3 | 0 | unassigned | R6 matrix-fit models and new matrix-fit driver |
| `[ ]` | R/nucleonfs.R | src/pydron/nucleonfs.py | 1 | 0 | 0 | 0 | unassigned | nucleonfs utilities |
| `[ ]` | R/onlinemeas.R | src/pydron/onlinemeas.py | 6 | 0 | 1 | 0 | unassigned | onlinemeas utilities |
| `[ ]` | R/pcac.R | src/pydron/pcac.py | 3 | 0 | 2 | 0 | unassigned | pcac utilities |
| `[ ]` | R/plotutils.R | src/pydron/plotutils.py | 16 | 0 | 3 | 7 | unassigned | plotting helpers |
| `[ ]` | R/prop_error.R | src/pydron/prop_error.py | 6 | 0 | 0 | 0 | unassigned | prop error utilities |
| `[ ]` | R/raw_cf.R | src/pydron/raw_cf.py | 27 | 0 | 17 | 8 | unassigned | raw correlator container, conversion, blocking, plotting |
| `[ ]` | R/RcppExports.R | src/pydron/rcpp_exports.py | 1 | 0 | 0 | 0 | unassigned | generated R wrapper for Rcpp export |
| `[ ]` | R/readutils.R | src/pydron/readutils.py | 19 | 0 | 17 | 0 | unassigned | text/binary/CMI/NISSA/gradient-flow readers |
| `[ ]` | R/removeTemporal.cf.R | src/pydron/remove_temporal_cf.py | 10 | 0 | 7 | 2 | unassigned | temporal-pollution removal and weighting |
| `[ ]` | R/seed.R | src/pydron/seed.py | 2 | 0 | 0 | 0 | unassigned | seed utilities |
| `[ ]` | R/string2error.R | src/pydron/string2error.py | 1 | 0 | 1 | 0 | unassigned | string2error utilities |
| `[ ]` | R/summary.ofit.R | src/pydron/summary_ofit.py | 2 | 0 | 0 | 2 | unassigned | summary ofit utilities |
| `[ ]` | R/tex-utils.R | src/pydron/tex_utils.py | 2 | 0 | 2 | 0 | unassigned | tex-utils utilities |
| `[ ]` | R/tflops.R | src/pydron/tflops.py | 1 | 0 | 0 | 0 | unassigned | tflops utilities |
| `[ ]` | R/tikzutils.R | src/pydron/tikzutils.py | 2 | 0 | 2 | 0 | unassigned | tikzutils utilities |
| `[ ]` | R/timeseries.R | src/pydron/timeseries.py | 2 | 0 | 2 | 0 | unassigned | timeseries utilities |
| `[ ]` | R/UWerr.R | src/pydron/uwerr.py | 7 | 0 | 3 | 2 | unassigned | autocorrelation/error analysis |
| `[ ]` | R/zeta_zp.R | src/pydron/zeta_zp.py | 1 | 0 | 1 | 0 | unassigned | zeta zp utilities |

## Checklist

- [ ] Keep original R files unchanged while Python parity is developed.
- [ ] Use `PORTING_INDEX.md` symbol rows as the authoritative implementation backlog.
- [ ] Rerun the R-backed parser whenever R source changes.
- [ ] For each public symbol, decide Pythonic API and compatibility wrapper before implementation.
- [ ] Promote local helper functions only when needed by the object model or tests.

## Testing requirements

| Source item | Required Python test or validation | Reference data needed | Owner/status | Notes |
|---|---|---|---|---|
| Public exports and S3 methods | R reference-output equivalence tests for numerical behavior | yes | `[ ]` / unassigned | 157 exports and 48 S3 methods |
| Internal helpers | Unit tests through public behavior unless promoted | usually no | `[ ]` / unassigned | R parser inventory is top-level only |
| R6 model classes | Constructor, prediction, Jacobian, and fit workflow parity | yes | `[ ]` / unassigned | 8 classes in `R/new_matrixfit.R` |

## Known issues and ambiguities

| Source path | Issue | Proposed resolution | Owner/status |
|---|---|---|---|
| Local/nested functions | R-backed symbol table intentionally tracks top-level declarations only | Treat local closures as implementation details unless promoted | `[?]` / unassigned |
