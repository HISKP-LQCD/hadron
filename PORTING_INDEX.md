# pydron porting index

Baseline commit: `42f8841c88e8ecb94a01a4abc3efc7a0a569fdc1`. R-backed inventory was generated with `R version 4.3.3 (2024-02-29)`. This remains a planning-only index; no Python algorithms are implemented here.

## Status labels

| Label | Meaning |
|---|---|
| [ ] | not started |
| [~] | in progress |
| [T] | translated, not fully tested |
| [E] | equivalence-tested against R output |
| [D] | documented |
| [X] | intentionally not ported, with reason |
| [?] | needs decision |

## Inventory summary

| Category | Count | Notes |
|---|---|---|
| Baseline commit | 42f8841c88e8ecb94a01a4abc3efc7a0a569fdc1 | Branch `port-plan` at inspection time |
| R/Rscript version | R version 4.3.3 (2024-02-29) | `Rscript --version` confirmed |
| R source files | 64 | All parsed successfully with `parse()` |
| Top-level R function definitions | 325 | R parser, top-level assignments only |
| Parsed R symbol rows | 333 | Top-level functions plus R6 class declarations |
| NAMESPACE exports | 157 | All matched to parsed top-level definitions |
| NAMESPACE S3 methods | 48 | All matched to parsed top-level definitions |
| Roxygen export tags | 208 | Roxygen `@export` tags mapped to following declarations where possible |
| UseMethod generics | 2 | extract_mass, residual_plot |
| R6 classes | 8 | All in `R/new_matrixfit.R` |
| Native source/header files | 8 | C/C++/header files |
| Compiled native artifacts in `src/` | 7 | Classified as build artifacts, not port sources |
| Registered native routines | 5 | From `src/RcppExports.cpp` |
| R-to-native calls | 5 | All `.Call`; no `.C` or `.Fortran` found |
| Dataset files inspected | 6 | All `.RData` files loaded in clean R environments |
| Dataset objects inspected | 6 | Object names/classes/dimensions/lengths/sizes recorded |
| Test files inspected | 13 | All parsed with R or Rmd chunk parsing |
| Vignettes/workflows indexed | 124 | `vignettes/`, `exec/`, selected `inst/`, `notes/`, root workflows |
| Man pages indexed | 231 | All parsed with `tools::parse_Rd()` |
| Files still not parsed reliably | 10 | Listed explicitly |

## Top-level tree

| Path | Kind | Porting note |
|---|---|---|
| .Rbuildignore | file | source/workflow/reference |
| .agents | dir | hidden local tooling - preserved |
| .clang-format | file | source/workflow/reference |
| .codex | dir | hidden local tooling - preserved |
| .gitignore | file | source/workflow/reference |
| .vimrc | file | source/workflow/reference |
| CODEX_PROMPT.md | file | source/workflow/reference |
| CONTRIBUTING.md | file | source/workflow/reference |
| DESCRIPTION | file | source/workflow/reference |
| NAMESPACE | file | source/workflow/reference |
| NEWS.md | file | source/workflow/reference |
| PLAN.md | file | source/workflow/reference |
| PORTING_INDEX.md | file | source/workflow/reference |
| R | dir | source/workflow/reference |
| README.md | file | source/workflow/reference |
| Weighted_Model.nb | file | source/workflow/reference |
| autom4te.cache | dir | build/cache - no summary |
| check | file | source/workflow/reference |
| cleanup | file | source/workflow/reference |
| config.log | file | build/cache - no summary |
| config.status | file | build/cache - no summary |
| configure | file | source/workflow/reference |
| configure.ac | file | source/workflow/reference |
| data | dir | source/workflow/reference |
| document | file | source/workflow/reference |
| exec | dir | source/workflow/reference |
| hadron.Rproj | file | source/workflow/reference |
| hooks | dir | source/workflow/reference |
| inst | dir | source/workflow/reference |
| install | file | source/workflow/reference |
| man | dir | source/workflow/reference |
| notes | dir | source/workflow/reference |
| src | dir | source/workflow/reference |
| test | file | source/workflow/reference |
| tests | dir | source/workflow/reference |
| verify-exports | file | source/workflow/reference |
| vignettes | dir | source/workflow/reference |

## R source files

| Status | Source path | Proposed Python target | Top-level functions | R6 classes | NAMESPACE exports | S3 methods | Owner | Main responsibility |
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

## R symbols

| Source path | R symbol | Category | Proposed Python path | Proposed Python API name | Compatibility wrapper name | Status label | Owner | Test requirement | Notes |
|---|---|---|---|---|---|---|---|---|---|
| R/CExp.R | CExp | exported function | src/pydron/cexp.py | cexp | cexp | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 23 |
| R/CExp.R | dCExpdm | internal helper | src/pydron/cexp.py | d_cexpdm | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 27 |
| R/ChiSqr.R | ChiSqr.singleCor | internal helper | src/pydron/chi_sqr.py | chi_sqr_single_cor | chi_sqr_single_cor | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 1 |
| R/ChiSqr.R | ChiSqr.pcac | internal helper | src/pydron/chi_sqr.py | chi_sqr_pcac | chi_sqr_pcac | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 6 |
| R/ChiSqr.R | ChiSqr.cst | internal helper | src/pydron/chi_sqr.py | chi_sqr_cst | chi_sqr_cst | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 12 |
| R/ChiSqr.R | ChiSqr.smeared | internal helper | src/pydron/chi_sqr.py | chi_sqr_smeared | chi_sqr_smeared | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 16 |
| R/ChiSqr.R | dChiSqrdpar.smeared | internal helper | src/pydron/chi_sqr.py | d_chi_sqrdpar_smeared | d_chi_sqrdpar_smeared | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 24 |
| R/ChiSqr.R | ChiSqr.1mass | internal helper | src/pydron/chi_sqr.py | chi_sqr_1mass | chi_sqr_1mass | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 38 |
| R/ChiSqr.R | ChiSqr.2mass | internal helper | src/pydron/chi_sqr.py | chi_sqr_2mass | chi_sqr_2mass | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 109 |
| R/ChiSqr.R | ChiSqr.3mass | internal helper | src/pydron/chi_sqr.py | chi_sqr_3mass | chi_sqr_3mass | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 208 |
| R/LuescherMethod.R | compute.qtildesq | internal helper | src/pydron/luescher_method.py | compute_qtildesq | compute_qtildesq | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 8 |
| R/LuescherMethod.R | compute.qtildesq.contdisp | internal helper | src/pydron/luescher_method.py | compute_qtildesq_contdisp | compute_qtildesq_contdisp | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 20 |
| R/RcppExports.R | read_nissa_textcf_kernel | I/O helper | src/pydron/rcpp_exports.py | read_nissa_textcf_kernel | - | `[ ]` | unassigned | fixture read/write parity or explicit non-port decision | line 4 |
| R/UWerr.R | uwerr | exported function | src/pydron/uwerr.py | uwerr | uwerr | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 63 |
| R/UWerr.R | uwerrprimary | exported function | src/pydron/uwerr.py | uwerrprimary | uwerrprimary | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 85 |
| R/UWerr.R | uwerrderived | exported function | src/pydron/uwerr.py | uwerrderived | uwerrderived | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 212 |
| R/UWerr.R | summary.uwerr | S3 method | src/pydron/uwerr.py | summary_uwerr | summary_uwerr | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 412 |
| R/UWerr.R | plot.uwerr | S3 method | src/pydron/uwerr.py | plot_uwerr | plot_uwerr | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 467 |
| R/UWerr.R | gammaerror | internal helper | src/pydron/uwerr.py | gammaerror | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 507 |
| R/UWerr.R | tauintplot | internal helper | src/pydron/uwerr.py | tauintplot | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 519 |
| R/alpha_s.R | alphas | exported function | src/pydron/alpha_s.py | alphas | alphas | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 23 |
| R/alpha_s.R | alphas.R | internal helper | src/pydron/alpha_s.py | alphas_r | alphas_r | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 31 |
| R/analysis_gradient_flow.R | analysis_gradient_flow | exported function | src/pydron/analysis_gradient_flow.py | analysis_gradient_flow | analysis_gradient_flow | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 52 |
| R/analysis_online.R | append_pdf_filename | internal helper | src/pydron/analysis_online.py | append_pdf_filename | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 1 |
| R/analysis_online.R | analysis_online | exported function | src/pydron/analysis_online.py | analysis_online | analysis_online | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 88 |
| R/analysis_online.R | construct_onlinemeas_rundir | internal helper | src/pydron/analysis_online.py | construct_onlinemeas_rundir | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 594 |
| R/block.R | block.ts | method | src/pydron/block.py | block_ts | block_ts | `[ ]` | unassigned | covered through owning object/generic tests | line 5 |
| R/boot_ts_array.R | boot_ts_array | exported function | src/pydron/boot_ts_array.py | boot_ts_array | boot_ts_array | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 42 |
| R/bootstrap.nlsfit.R | parametric.bootstrap | exported function | src/pydron/bootstrap_nlsfit.py | parametric_bootstrap | parametric_bootstrap | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 24 |
| R/bootstrap.nlsfit.R | parametric.bootstrap.cov | exported function | src/pydron/bootstrap_nlsfit.py | parametric_bootstrap_cov | parametric_bootstrap_cov | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 58 |
| R/bootstrap.nlsfit.R | parametric.nlsfit | exported function | src/pydron/bootstrap_nlsfit.py | parametric_nlsfit | parametric_nlsfit | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 109 |
| R/bootstrap.nlsfit.R | parametric.nlsfit.cov | exported function | src/pydron/bootstrap_nlsfit.py | parametric_nlsfit_cov | parametric_nlsfit_cov | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 155 |
| R/bootstrap.nlsfit.R | get.errors | internal helper | src/pydron/bootstrap_nlsfit.py | get_errors | get_errors | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 185 |
| R/bootstrap.nlsfit.R | get.errors.wo.bootstrap | internal helper | src/pydron/bootstrap_nlsfit.py | get_errors_wo_bootstrap | get_errors_wo_bootstrap | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 252 |
| R/bootstrap.nlsfit.R | set.fitchi | internal helper | src/pydron/bootstrap_nlsfit.py | set_fitchi | set_fitchi | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 293 |
| R/bootstrap.nlsfit.R | set.dfitchi | internal helper | src/pydron/bootstrap_nlsfit.py | set_dfitchi | set_dfitchi | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 324 |
| R/bootstrap.nlsfit.R | set.dfitchisqr | internal helper | src/pydron/bootstrap_nlsfit.py | set_dfitchisqr | set_dfitchisqr | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 376 |
| R/bootstrap.nlsfit.R | set.wrapper | internal helper | src/pydron/bootstrap_nlsfit.py | set_wrapper | set_wrapper | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 385 |
| R/bootstrap.nlsfit.R | simple.nlsfit | exported function | src/pydron/bootstrap_nlsfit.py | simple_nlsfit | simple_nlsfit | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 462 |
| R/bootstrap.nlsfit.R | bootstrap.nlsfit | exported function | src/pydron/bootstrap_nlsfit.py | bootstrap_nlsfit | bootstrap_nlsfit | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 825 |
| R/bootstrap.nlsfit.R | summary.bootstrapfit | S3 method | src/pydron/bootstrap_nlsfit.py | summary_bootstrapfit | summary_bootstrapfit | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 1126 |
| R/bootstrap.nlsfit.R | print.bootstrapfit | S3 method | src/pydron/bootstrap_nlsfit.py | print_bootstrapfit | print_bootstrapfit | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 1202 |
| R/bootstrap.nlsfit.R | plot.bootstrapfit | S3 method | src/pydron/bootstrap_nlsfit.py | plot_bootstrapfit | plot_bootstrapfit | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 1234 |
| R/bootstrap.nlsfit.R | residual_plot | generic | src/pydron/bootstrap_nlsfit.py | residual_plot | residual_plot | `[ ]` | unassigned | dispatch parity tests for registered methods | NAMESPACE export; UseMethod: residual_plot; line 1296 |
| R/bootstrap.nlsfit.R | residual_plot.bootstrapfit | S3 method | src/pydron/bootstrap_nlsfit.py | residual_plot_bootstrapfit | residual_plot_bootstrapfit | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 1301 |
| R/bootstrap.nlsfit.R | predict.bootstrapfit | S3 method | src/pydron/bootstrap_nlsfit.py | predict_bootstrapfit | predict_bootstrapfit | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 1365 |
| R/bootstrapnumber.R | meanindexed | internal helper | src/pydron/bootstrapnumber.py | meanindexed | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 1 |
| R/bootstrapnumber.R | sd.index | internal helper | src/pydron/bootstrapnumber.py | sd_index | sd_index | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 5 |
| R/bootstrapnumber.R | bootstrap.analysis | exported function | src/pydron/bootstrapnumber.py | bootstrap_analysis | bootstrap_analysis | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 41 |
| R/cdh.R | cdh | exported function | src/pydron/cdh.py | cdh | cdh | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 47 |
| R/cdh.R | cdhnew | exported function | src/pydron/cdh.py | cdhnew | cdhnew | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 120 |
| R/cdh.R | cdh.R | internal helper | src/pydron/cdh.py | cdh_r | cdh_r | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 136 |
| R/cdh.R | cdhnew.R | internal helper | src/pydron/cdh.py | cdhnew_r | cdhnew_r | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 257 |
| R/cf.R | cf | exported function | src/pydron/cf.py | cf | cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 27 |
| R/cf.R | cf_meta | exported function | src/pydron/cf.py | cf_meta | cf_meta | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 51 |
| R/cf.R | cf_boot | exported function | src/pydron/cf.py | cf_boot | cf_boot | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 90 |
| R/cf.R | jackknife_error | exported function | src/pydron/cf.py | jackknife_error | jackknife_error | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 162 |
| R/cf.R | jackknife_cov | exported function | src/pydron/cf.py | jackknife_cov | jackknife_cov | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 196 |
| R/cf.R | cf_orig | exported function | src/pydron/cf.py | cf_orig | cf_orig | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 239 |
| R/cf.R | cf_principal_correlator | exported function | src/pydron/cf.py | cf_principal_correlator | cf_principal_correlator | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 263 |
| R/cf.R | cf_shifted | exported function | src/pydron/cf.py | cf_shifted | cf_shifted | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 291 |
| R/cf.R | cf_smeared | exported function | src/pydron/cf.py | cf_smeared | cf_smeared | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 323 |
| R/cf.R | cf_subtracted | exported function | src/pydron/cf.py | cf_subtracted | cf_subtracted | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 349 |
| R/cf.R | cf_weighted | exported function | src/pydron/cf.py | cf_weighted | cf_weighted | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 377 |
| R/cf.R | is_empty.cf | exported function | src/pydron/cf.py | is_empty_cf | is_empty_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 403 |
| R/cf.R | resampling_is_compatible | internal helper | src/pydron/cf.py | resampling_is_compatible | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 422 |
| R/cf.R | resampling_is_concatenable | internal helper | src/pydron/cf.py | resampling_is_concatenable | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 456 |
| R/cf.R | has_icf | exported function | src/pydron/cf.py | has_icf | has_icf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 483 |
| R/cf.R | gen.block.array | internal helper | src/pydron/cf.py | gen_block_array | gen_block_array | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 488 |
| R/cf.R | bootstrap.cf | exported function | src/pydron/cf.py | bootstrap_cf | bootstrap_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 535 |
| R/cf.R | double_bootstrap.cf | exported function | src/pydron/cf.py | double_bootstrap_cf | double_bootstrap_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 601 |
| R/cf.R | jackknife.cf | exported function | src/pydron/cf.py | jackknife_cf | jackknife_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 683 |
| R/cf.R | uwerr.cf | exported function | src/pydron/cf.py | uwerr_cf | uwerr_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 768 |
| R/cf.R | addConfIndex2cf | exported function | src/pydron/cf.py | add_conf_index2cf | add_conf_index2cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 812 |
| R/cf.R | addStat.cf | exported function | src/pydron/cf.py | add_stat_cf | add_stat_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 844 |
| R/cf.R | avg.cbt.cf | exported function | src/pydron/cf.py | avg_cbt_cf | avg_cbt_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 887 |
| R/cf.R | add.cf | exported function | src/pydron/cf.py | add_cf | add_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 916 |
| R/cf.R | +.cf | S3 method | src/pydron/cf.py | add_cf | add_cf | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 983 |
| R/cf.R | -.cf | S3 method | src/pydron/cf.py | sub_cf | sub_cf | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 996 |
| R/cf.R | *.cf | S3 method | src/pydron/cf.py | mul_cf | mul_cf | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 1002 |
| R/cf.R | /.cf | S3 method | src/pydron/cf.py | div_cf | div_cf | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 1023 |
| R/cf.R | apply_elementwise.cf | internal helper | src/pydron/cf.py | apply_elementwise_cf | apply_elementwise_cf | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 1027 |
| R/cf.R | mul.cf | exported function | src/pydron/cf.py | mul_cf | mul_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 1069 |
| R/cf.R | extractSingleCor.cf | exported function | src/pydron/cf.py | extract_single_cor_cf | extract_single_cor_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 1124 |
| R/cf.R | is.cf | exported function | src/pydron/cf.py | is_cf | is_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 1168 |
| R/cf.R | c.cf | S3 method | src/pydron/cf.py | c_cf | c_cf | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 1181 |
| R/cf.R | concat.cf | exported function | src/pydron/cf.py | concat_cf | concat_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 1195 |
| R/cf.R | plot.cf | S3 method | src/pydron/cf.py | plot_cf | plot_cf | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 1298 |
| R/cf.R | shift.cf | exported function | src/pydron/cf.py | shift_cf | shift_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 1333 |
| R/cf.R | invalidate.samples.cf | exported function | src/pydron/cf.py | invalidate_samples_cf | invalidate_samples_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 1407 |
| R/cf.R | symmetrise.cf | exported function | src/pydron/cf.py | symmetrise_cf | symmetrise_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 1450 |
| R/cf.R | unsymmetrise.cf | exported function | src/pydron/cf.py | unsymmetrise_cf | unsymmetrise_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 1548 |
| R/cf.R | summary.cf | S3 method | src/pydron/cf.py | summary_cf | summary_cf | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 1649 |
| R/cf.R | print.cf | S3 method | src/pydron/cf.py | print_cf | print_cf | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 1691 |
| R/computeDisc.R | computeDisc | exported function | src/pydron/compute_disc.py | compute_disc | compute_disc | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 62 |
| R/computeacf.R | computeacf | exported function | src/pydron/computeacf.py | computeacf | computeacf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 41 |
| R/computeacf.R | plot.hadronacf | S3 method | src/pydron/computeacf.py | plot_hadronacf | plot_hadronacf | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 88 |
| R/computeacf.R | summary.hadronacf | S3 method | src/pydron/computeacf.py | summary_hadronacf | summary_hadronacf | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 118 |
| R/computefps.R | computefps | exported function | src/pydron/computefps.py | computefps | computefps | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 64 |
| R/computefps.R | computefpsOS | exported function | src/pydron/computefps.py | computefps_os | computefps_os | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 174 |
| R/correlatedRNG.R | corrnorm | internal helper | src/pydron/correlated_rng.py | corrnorm | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 1 |
| R/correlatedRNG.R | corrnorm2 | internal helper | src/pydron/correlated_rng.py | corrnorm2 | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 9 |
| R/cosh_nlsfit.R | sum_cosh | internal helper | src/pydron/cosh_nlsfit.py | sum_cosh | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 4 |
| R/cosh_nlsfit.R | cosh_to_effmass | internal helper | src/pydron/cosh_nlsfit.py | cosh_to_effmass | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 10 |
| R/cosh_nlsfit.R | fit.cosh | exported function | src/pydron/cosh_nlsfit.py | fit_cosh | fit_cosh | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 93 |
| R/cosh_nlsfit.R | plot.coshfit | S3 method | src/pydron/cosh_nlsfit.py | plot_coshfit | plot_coshfit | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 251 |
| R/cosh_nlsfit.R | summary.coshfit | S3 method | src/pydron/cosh_nlsfit.py | summary_coshfit | summary_coshfit | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 317 |
| R/cvc_readutils.R | cvc_local_loop_key | exported function | src/pydron/cvc_readutils.py | cvc_local_loop_key | cvc_local_loop_key | `[ ]` | unassigned | I/O fixture plus R reference object comparison | NAMESPACE export; line 11 |
| R/cvc_readutils.R | cvc_local_loop_key | exported function | src/pydron/cvc_readutils.py | cvc_local_loop_key | cvc_local_loop_key | `[ ]` | unassigned | I/O fixture plus R reference object comparison | NAMESPACE export; line 11 |
| R/cvc_readutils.R | correlators_key_meson_2pt | exported function | src/pydron/cvc_readutils.py | correlators_key_meson_2pt | correlators_key_meson_2pt | `[ ]` | unassigned | I/O fixture plus R reference object comparison | NAMESPACE export; line 57 |
| R/cvc_readutils.R | cf_key_meson_2pt | exported function | src/pydron/cvc_readutils.py | cf_key_meson_2pt | cf_key_meson_2pt | `[ ]` | unassigned | I/O fixture plus R reference object comparison | NAMESPACE export; line 90 |
| R/cvc_readutils.R | correlators_key_meson_3pt | exported function | src/pydron/cvc_readutils.py | correlators_key_meson_3pt | correlators_key_meson_3pt | `[ ]` | unassigned | I/O fixture plus R reference object comparison | NAMESPACE export; line 161 |
| R/cvc_readutils.R | cf_key_meson_3pt | exported function | src/pydron/cvc_readutils.py | cf_key_meson_3pt | cf_key_meson_3pt | `[ ]` | unassigned | I/O fixture plus R reference object comparison | NAMESPACE export; line 262 |
| R/cvc_readutils.R | cvc_to_raw_cf | exported function | src/pydron/cvc_readutils.py | cvc_to_raw_cf | cvc_to_raw_cf | `[ ]` | unassigned | I/O fixture plus R reference object comparison | NAMESPACE export; line 324 |
| R/cvc_readutils.R | cvc_read_loops | exported function | src/pydron/cvc_readutils.py | cvc_read_loops | cvc_read_loops | `[ ]` | unassigned | I/O fixture plus R reference object comparison | NAMESPACE export; line 385 |
| R/cyprus_readutils.R | cyprus_make_key_scalar | exported function | src/pydron/cyprus_readutils.py | cyprus_make_key_scalar | cyprus_make_key_scalar | `[ ]` | unassigned | I/O fixture plus R reference object comparison | NAMESPACE export; line 18 |
| R/cyprus_readutils.R | cyprus_make_key_vector | exported function | src/pydron/cyprus_readutils.py | cyprus_make_key_vector | cyprus_make_key_vector | `[ ]` | unassigned | I/O fixture plus R reference object comparison | NAMESPACE export; line 50 |
| R/cyprus_readutils.R | cyprus_read_loops | exported function | src/pydron/cyprus_readutils.py | cyprus_read_loops | cyprus_read_loops | `[ ]` | unassigned | I/O fixture plus R reference object comparison | NAMESPACE export; line 130 |
| R/deriv_utils.R | create_displ_chains | internal helper | src/pydron/deriv_utils.py | create_displ_chains | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 10 |
| R/effectivemass.R | effectivemass.cf | exported function | src/pydron/effectivemass.py | effectivemass_cf | effectivemass_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 71 |
| R/effectivemass.R | bootstrap.effectivemass | exported function | src/pydron/effectivemass.py | bootstrap_effectivemass | bootstrap_effectivemass | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 226 |
| R/effectivemass.R | fit.constant | internal helper | src/pydron/effectivemass.py | fit_constant | fit_constant | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 264 |
| R/effectivemass.R | fit.effectivemass | exported function | src/pydron/effectivemass.py | fit_effectivemass | fit_effectivemass | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 350 |
| R/effectivemass.R | summary.effectivemass | S3 method | src/pydron/effectivemass.py | summary_effectivemass | summary_effectivemass | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 507 |
| R/effectivemass.R | summary.effectivemassfit | S3 method | src/pydron/effectivemass.py | summary_effectivemassfit | summary_effectivemassfit | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 529 |
| R/effectivemass.R | print.effectivemassfit | S3 method | src/pydron/effectivemass.py | print_effectivemassfit | print_effectivemassfit | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 566 |
| R/effectivemass.R | plot.effectivemass | S3 method | src/pydron/effectivemass.py | plot_effectivemass | plot_effectivemass | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 584 |
| R/fit.plateau2cf.R | fit.plateau2cf | exported function | src/pydron/fit_plateau2cf.py | fit_plateau2cf | fit_plateau2cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 23 |
| R/fitmass.R | fitmass | internal helper | src/pydron/fitmass.py | fitmass | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 1 |
| R/fs.R | fs | internal helper | src/pydron/fs.py | fs | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 1 |
| R/fs.R | ChiSqr.fs | internal helper | src/pydron/fs.py | chi_sqr_fs | chi_sqr_fs | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 87 |
| R/fs.R | ChiSqr.fs.fps | internal helper | src/pydron/fs.py | chi_sqr_fs_fps | chi_sqr_fs_fps | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 93 |
| R/fs.R | ChiSqr.fs.comb | internal helper | src/pydron/fs.py | chi_sqr_fs_comb | chi_sqr_fs_comb | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 99 |
| R/fs.R | ChiSqr.fs.comb2 | internal helper | src/pydron/fs.py | chi_sqr_fs_comb2 | chi_sqr_fs_comb2 | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 107 |
| R/fs.R | ChiSqr.pow | internal helper | src/pydron/fs.py | chi_sqr_pow | chi_sqr_pow | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 115 |
| R/fs.mpia0.R | fs.qcotdelta | exported function | src/pydron/fs_mpia0.py | fs_qcotdelta | fs_qcotdelta | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 34 |
| R/fs.mpia0.R | fs.a0 | exported function | src/pydron/fs_mpia0.py | fs_a0 | fs_a0 | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 60 |
| R/fs.mpia0.R | fs.mpia0 | exported function | src/pydron/fs_mpia0.py | fs_mpia0 | fs_mpia0 | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 87 |
| R/functional.R | foldr1 | exported function | src/pydron/functional.py | foldr1 | foldr1 | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 32 |
| R/g1.R | g1 | exported function | src/pydron/g1.py | g1 | g1 | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 10 |
| R/gamma.R | gm_mu | exported function | src/pydron/gamma.py | gm_mu | gm_mu | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 33 |
| R/gamma.R | Tr | internal helper | src/pydron/gamma.py | tr | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 180 |
| R/getCor.R | getCor | internal helper | src/pydron/get_cor.py | get_cor | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 1 |
| R/getNxNmatrix.R | getNxNmatrix | internal helper | src/pydron/get_nx_nmatrix.py | get_nx_nmatrix | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 1 |
| R/getfit.boot.R | getfit.boot | internal helper | src/pydron/getfit_boot.py | getfit_boot | getfit_boot | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 1 |
| R/gevp.R | permutations | internal helper | src/pydron/gevp.py | permutations | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 2 |
| R/gevp.R | gevp | exported function | src/pydron/gevp.py | gevp | gevp | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 80 |
| R/gevp.R | bootstrap.gevp | exported function | src/pydron/gevp.py | bootstrap_gevp | bootstrap_gevp | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 302 |
| R/gevp.R | gevp2cf | exported function | src/pydron/gevp.py | gevp2cf | gevp2cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 368 |
| R/gevp.R | gevp2amplitude | exported function | src/pydron/gevp.py | gevp2amplitude | gevp2amplitude | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 482 |
| R/gevp.R | summary.gevp.amplitude | S3 method | src/pydron/gevp.py | summary_gevp_amplitude | summary_gevp_amplitude | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 598 |
| R/gevp.R | plot.gevp.amplitude | S3 method | src/pydron/gevp.py | plot_gevp_amplitude | plot_gevp_amplitude | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 638 |
| R/h5utils.R | h5_get_dataset | exported function | src/pydron/h5utils.py | h5_get_dataset | h5_get_dataset | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 10 |
| R/h5utils.R | h5_names_exist | I/O helper | src/pydron/h5utils.py | h5_names_exist | - | `[ ]` | unassigned | fixture read/write parity or explicit non-port decision | line 37 |
| R/hankel.R | hankel.matrix | internal helper | src/pydron/hankel.py | hankel_matrix | hankel_matrix | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 2 |
| R/hankel.R | gevp.hankel_summed | internal helper | src/pydron/hankel.py | gevp_hankel_summed | gevp_hankel_summed | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 28 |
| R/hankel.R | bootstrap.hankel_summed | exported function | src/pydron/hankel.py | bootstrap_hankel_summed | bootstrap_hankel_summed | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 90 |
| R/hankel.R | summary.hankel_summed | S3 method | src/pydron/hankel.py | summary_hankel_summed | summary_hankel_summed | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 134 |
| R/hankel.R | gevp.hankel | internal helper | src/pydron/hankel.py | gevp_hankel | gevp_hankel | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 186 |
| R/hankel.R | bootstrap.hankel | exported function | src/pydron/hankel.py | bootstrap_hankel | bootstrap_hankel | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 325 |
| R/hankel.R | bootstrap.pgevm | exported function | src/pydron/hankel.py | bootstrap_pgevm | bootstrap_pgevm | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 437 |
| R/hankel.R | pgevm2effectivemass | exported function | src/pydron/hankel.py | pgevm2effectivemass | pgevm2effectivemass | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 579 |
| R/hankel.R | plot_hankel_spectrum | plotting helper | src/pydron/hankel.py | plot_hankel_spectrum | - | `[ ]` | unassigned | plot smoke test or extracted plotting-data parity | line 802 |
| R/hankel.R | hankel2cf | exported function | src/pydron/hankel.py | hankel2cf | hankel2cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 862 |
| R/hankel.R | hankel2effectivemass | exported function | src/pydron/hankel.py | hankel2effectivemass | hankel2effectivemass | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 1046 |
| R/hankel.R | hankeldensity2effectivemass | internal helper | src/pydron/hankel.py | hankeldensity2effectivemass | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 1136 |
| R/hankel.R | resample_hankel | exported function | src/pydron/hankel.py | resample_hankel | resample_hankel | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 1229 |
| R/hankel.truncated.R | spectrum.truncated.gevp | internal helper | src/pydron/hankel_truncated.py | spectrum_truncated_gevp | spectrum_truncated_gevp | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 3 |
| R/hankel.truncated.R | coeffs.truncated.gevp | internal helper | src/pydron/hankel_truncated.py | coeffs_truncated_gevp | coeffs_truncated_gevp | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 40 |
| R/hankel.truncated.R | vec.coeffs.truncated.gevp | internal helper | src/pydron/hankel_truncated.py | vec_coeffs_truncated_gevp | vec_coeffs_truncated_gevp | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 66 |
| R/hankel.truncated.R | reconstruct.correlators | internal helper | src/pydron/hankel_truncated.py | reconstruct_correlators | reconstruct_correlators | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 82 |
| R/hankel.truncated.R | gevp.truncated.hankel | exported function | src/pydron/hankel_truncated.py | gevp_truncated_hankel | gevp_truncated_hankel | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 146 |
| R/hankel.truncated.R | bootstrap.truncated.pgevm | exported function | src/pydron/hankel_truncated.py | bootstrap_truncated_pgevm | bootstrap_truncated_pgevm | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 280 |
| R/hankel.truncated.R | pgevm2bootstrapfit | exported function | src/pydron/hankel_truncated.py | pgevm2bootstrapfit | pgevm2bootstrapfit | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 413 |
| R/hankel.truncated.R | plot.truncated.pgevm | S3 method | src/pydron/hankel_truncated.py | plot_truncated_pgevm | plot_truncated_pgevm | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 471 |
| R/inv_cosh.R | invcosh | exported function | src/pydron/inv_cosh.py | invcosh | invcosh | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 19 |
| R/invertCovMatrix.R | invertCovMatrix | exported function | src/pydron/invert_cov_matrix.py | invert_cov_matrix | invert_cov_matrix | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 41 |
| R/jackknifeafterboot.R | jab | internal helper | src/pydron/jackknifeafterboot.py | jab | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 1 |
| R/jackknifeafterboot.R | jab.cf | internal helper | src/pydron/jackknifeafterboot.py | jab_cf | jab_cf | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 45 |
| R/jackknifeafterboot.R | jab.cf.derived | internal helper | src/pydron/jackknifeafterboot.py | jab_cf_derived | jab_cf_derived | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 79 |
| R/jackknifeafterboot.R | jab.matrixfit | internal helper | src/pydron/jackknifeafterboot.py | jab_matrixfit | jab_matrixfit | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 96 |
| R/jackknifeafterboot.R | jab.effectivemass | internal helper | src/pydron/jackknifeafterboot.py | jab_effectivemass | jab_effectivemass | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 105 |
| R/jackknifeafterboot.R | jab.effectivemassfit | internal helper | src/pydron/jackknifeafterboot.py | jab_effectivemassfit | jab_effectivemassfit | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 113 |
| R/kappa.R | kappa | internal helper | src/pydron/kappa.py | kappa | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 1 |
| R/kappa.R | m0 | internal helper | src/pydron/kappa.py | m0 | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 5 |
| R/lanczos.R | bootstrap.lanczos | exported function | src/pydron/lanczos.py | bootstrap_lanczos | bootstrap_lanczos | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 48 |
| R/lanczos.R | lanczos.solve | internal helper | src/pydron/lanczos.py | lanczos_solve | lanczos_solve | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 134 |
| R/legacy_functions.R | effmass | internal helper | src/pydron/legacy_functions.py | effmass | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 17 |
| R/legacy_functions.R | effmass2 | internal helper | src/pydron/legacy_functions.py | effmass2 | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 32 |
| R/legacy_functions.R | effectivemass | internal helper | src/pydron/legacy_functions.py | effectivemass | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 53 |
| R/looptools.R | disc_3pt | exported function | src/pydron/looptools.py | disc_3pt | disc_3pt | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 25 |
| R/looptools.R | loop_2pt | exported function | src/pydron/looptools.py | loop_2pt | loop_2pt | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 136 |
| R/looptools.R | loop_stochav | exported function | src/pydron/looptools.py | loop_stochav | loop_stochav | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 233 |
| R/looptools.R | loop_vev_subtract | exported function | src/pydron/looptools.py | loop_vev_subtract | loop_vev_subtract | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 273 |
| R/looptools.R | loop_spin_project | exported function | src/pydron/looptools.py | loop_spin_project | loop_spin_project | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 326 |
| R/matrixfit.R | bootstrap.meanerror | exported function | src/pydron/matrixfit.py | bootstrap_meanerror | bootstrap_meanerror | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 15 |
| R/matrixfit.R | matrixModel | internal helper | src/pydron/matrixfit.py | matrix_model | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 38 |
| R/matrixfit.R | pcModel | internal helper | src/pydron/matrixfit.py | pc_model | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 59 |
| R/matrixfit.R | matrixChisqr | internal helper | src/pydron/matrixfit.py | matrix_chisqr | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 63 |
| R/matrixfit.R | matrixChi | internal helper | src/pydron/matrixfit.py | matrix_chi | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 68 |
| R/matrixfit.R | dmatrixChi | internal helper | src/pydron/matrixfit.py | dmatrix_chi | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 74 |
| R/matrixfit.R | dmatrixChisqr | internal helper | src/pydron/matrixfit.py | dmatrix_chisqr | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 90 |
| R/matrixfit.R | pcChi | internal helper | src/pydron/matrixfit.py | pc_chi | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 117 |
| R/matrixfit.R | pcChisqr | internal helper | src/pydron/matrixfit.py | pc_chisqr | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 122 |
| R/matrixfit.R | dpcChi | internal helper | src/pydron/matrixfit.py | dpc_chi | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 128 |
| R/matrixfit.R | dpcChisqr | internal helper | src/pydron/matrixfit.py | dpc_chisqr | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 139 |
| R/matrixfit.R | matrixChisqr.shifted | internal helper | src/pydron/matrixfit.py | matrix_chisqr_shifted | matrix_chisqr_shifted | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 152 |
| R/matrixfit.R | dmatrixChisqr.shifted | internal helper | src/pydron/matrixfit.py | dmatrix_chisqr_shifted | dmatrix_chisqr_shifted | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 158 |
| R/matrixfit.R | matrixChi.shifted | internal helper | src/pydron/matrixfit.py | matrix_chi_shifted | matrix_chi_shifted | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 177 |
| R/matrixfit.R | dmatrixChi.shifted | internal helper | src/pydron/matrixfit.py | dmatrix_chi_shifted | dmatrix_chi_shifted | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 183 |
| R/matrixfit.R | deriv.CExp | method | src/pydron/matrixfit.py | deriv_cexp | deriv_cexp | `[ ]` | unassigned | covered through owning object/generic tests | line 201 |
| R/matrixfit.R | deriv.CExp.shifted | method | src/pydron/matrixfit.py | deriv_cexp_shifted | deriv_cexp_shifted | `[ ]` | unassigned | covered through owning object/generic tests | line 212 |
| R/matrixfit.R | deriv.pcModel | method | src/pydron/matrixfit.py | deriv_pc_model | deriv_pc_model | `[ ]` | unassigned | covered through owning object/generic tests | line 223 |
| R/matrixfit.R | matrixfit | exported function | src/pydron/matrixfit.py | matrixfit | matrixfit | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 347 |
| R/matrixfit.R | plot.matrixfit | S3 method | src/pydron/matrixfit.py | plot_matrixfit | plot_matrixfit | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 649 |
| R/matrixfit.R | summary.matrixfit | S3 method | src/pydron/matrixfit.py | summary_matrixfit | summary_matrixfit | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 781 |
| R/matrixfit.R | fit.formatrixboot | internal helper | src/pydron/matrixfit.py | fit_formatrixboot | fit_formatrixboot | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 844 |
| R/matrixfit.R | subtract.excitedstates | exported function | src/pydron/matrixfit.py | subtract_excitedstates | subtract_excitedstates | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 894 |
| R/momentum_utils.R | mom_combinations | exported function | src/pydron/momentum_utils.py | mom_combinations | mom_combinations | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 8 |
| R/new_matrixfit.R | make_sign_vec | internal helper | src/pydron/new_matrixfit.py | make_sign_vec | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 2 |
| R/new_matrixfit.R | make_ov_sign_vec | internal helper | src/pydron/new_matrixfit.py | make_ov_sign_vec | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 16 |
| R/new_matrixfit.R | MatrixModel | R6 class | src/pydron/new_matrixfit.py | MatrixModel | MatrixModel | `[ ]` | unassigned | constructor/API tests plus R reference predictions when numerical | R6Class declaration: MatrixModel; line 27 |
| R/new_matrixfit.R | SingleModel | R6 class | src/pydron/new_matrixfit.py | SingleModel | SingleModel | `[ ]` | unassigned | constructor/API tests plus R reference predictions when numerical | R6Class declaration: SingleModel; line 79 |
| R/new_matrixfit.R | TwoAmplitudesModel | R6 class | src/pydron/new_matrixfit.py | TwoAmplitudesModel | TwoAmplitudesModel | `[ ]` | unassigned | constructor/API tests plus R reference predictions when numerical | R6Class declaration: TwoAmplitudesModel; line 127 |
| R/new_matrixfit.R | ShiftedModel | R6 class | src/pydron/new_matrixfit.py | ShiftedModel | ShiftedModel | `[ ]` | unassigned | constructor/API tests plus R reference predictions when numerical | R6Class declaration: ShiftedModel; line 159 |
| R/new_matrixfit.R | WeightedModel | R6 class | src/pydron/new_matrixfit.py | WeightedModel | WeightedModel | `[ ]` | unassigned | constructor/API tests plus R reference predictions when numerical | R6Class declaration: WeightedModel; line 219 |
| R/new_matrixfit.R | TwoStateModel | R6 class | src/pydron/new_matrixfit.py | TwoStateModel | TwoStateModel | `[ ]` | unassigned | constructor/API tests plus R reference predictions when numerical | R6Class declaration: TwoStateModel; line 276 |
| R/new_matrixfit.R | NParticleModel | R6 class | src/pydron/new_matrixfit.py | NParticleModel | NParticleModel | `[ ]` | unassigned | constructor/API tests plus R reference predictions when numerical | R6Class declaration: NParticleModel; line 330 |
| R/new_matrixfit.R | SingleConstantModel | R6 class | src/pydron/new_matrixfit.py | SingleConstantModel | SingleConstantModel | `[ ]` | unassigned | constructor/API tests plus R reference predictions when numerical | R6Class declaration: SingleConstantModel; line 376 |
| R/new_matrixfit.R | new_matrixfit | exported function | src/pydron/new_matrixfit.py | new_matrixfit | new_matrixfit | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 480 |
| R/new_matrixfit.R | make_parlist | exported function | src/pydron/new_matrixfit.py | make_parlist | make_parlist | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 648 |
| R/new_matrixfit.R | make_parind | exported function | src/pydron/new_matrixfit.py | make_parind | make_parind | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 669 |
| R/nucleonfs.R | nucleonfs | internal helper | src/pydron/nucleonfs.py | nucleonfs | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 4 |
| R/onlinemeas.R | onlinemeas | exported function | src/pydron/onlinemeas.py | onlinemeas | onlinemeas | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 95 |
| R/onlinemeas.R | fitmpcac.online | internal helper | src/pydron/onlinemeas.py | fitmpcac_online | fitmpcac_online | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 302 |
| R/onlinemeas.R | fitmass.online | internal helper | src/pydron/onlinemeas.py | fitmass_online | fitmass_online | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 320 |
| R/onlinemeas.R | fitf.online | internal helper | src/pydron/onlinemeas.py | fitf_online | fitf_online | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 338 |
| R/onlinemeas.R | pcacsym.online | internal helper | src/pydron/onlinemeas.py | pcacsym_online | pcacsym_online | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 357 |
| R/onlinemeas.R | fit.online.boot | internal helper | src/pydron/onlinemeas.py | fit_online_boot | fit_online_boot | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 362 |
| R/pcac.R | pcacsym | internal helper | src/pydron/pcac.py | pcacsym | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 1 |
| R/pcac.R | pcacfit | exported function | src/pydron/pcac.py | pcacfit | pcacfit | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 39 |
| R/pcac.R | pcac | exported function | src/pydron/pcac.py | pcac | pcac | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 79 |
| R/plotutils.R | compute.plotlims | plotting helper | src/pydron/plotutils.py | compute_plotlims | compute_plotlims | `[ ]` | unassigned | plot smoke test or extracted plotting-data parity | line 14 |
| R/plotutils.R | is.vectorial | plotting helper | src/pydron/plotutils.py | is_vectorial | is_vectorial | `[ ]` | unassigned | plot smoke test or extracted plotting-data parity | line 35 |
| R/plotutils.R | errorpos | plotting helper | src/pydron/plotutils.py | errorpos | - | `[ ]` | unassigned | plot smoke test or extracted plotting-data parity | line 39 |
| R/plotutils.R | plotwitherror | exported function | src/pydron/plotutils.py | plotwitherror | plotwitherror | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 169 |
| R/plotutils.R | plothlinewitherror | exported function | src/pydron/plotutils.py | plothlinewitherror | plothlinewitherror | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 401 |
| R/plotutils.R | plot.massfit | S3 method | src/pydron/plotutils.py | plot_massfit | plot_massfit | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 424 |
| R/plotutils.R | plot.cfit | S3 method | src/pydron/plotutils.py | plot_cfit | plot_cfit | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 441 |
| R/plotutils.R | plot.ofit | S3 method | src/pydron/plotutils.py | plot_ofit | plot_ofit | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 518 |
| R/plotutils.R | plot.effmass | S3 method | src/pydron/plotutils.py | plot_effmass | plot_effmass | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 535 |
| R/plotutils.R | plot.averx | S3 method | src/pydron/plotutils.py | plot_averx | plot_averx | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 574 |
| R/plotutils.R | plot.pionff | S3 method | src/pydron/plotutils.py | plot_pionff | plot_pionff | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 611 |
| R/plotutils.R | plot.outputdata | S3 method | src/pydron/plotutils.py | plot_outputdata | plot_outputdata | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 660 |
| R/plotutils.R | drawYbars | plotting helper | src/pydron/plotutils.py | draw_ybars | - | `[ ]` | unassigned | plot smoke test or extracted plotting-data parity | line 679 |
| R/plotutils.R | drawXbars | plotting helper | src/pydron/plotutils.py | draw_xbars | - | `[ ]` | unassigned | plot smoke test or extracted plotting-data parity | line 689 |
| R/plotutils.R | new_window_if_appropriate | plotting helper | src/pydron/plotutils.py | new_window_if_appropriate | - | `[ ]` | unassigned | plot smoke test or extracted plotting-data parity | line 699 |
| R/plotutils.R | pointswithslantederror | exported function | src/pydron/plotutils.py | pointswithslantederror | pointswithslantederror | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 740 |
| R/prop_error.R | compute_square | internal helper | src/pydron/prop_error.py | compute_square | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 3 |
| R/prop_error.R | compute_sqrt | internal helper | src/pydron/prop_error.py | compute_sqrt | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 14 |
| R/prop_error.R | compute_ratio | internal helper | src/pydron/prop_error.py | compute_ratio | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 30 |
| R/prop_error.R | compute_product | internal helper | src/pydron/prop_error.py | compute_product | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 41 |
| R/prop_error.R | compute_sum | internal helper | src/pydron/prop_error.py | compute_sum | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 52 |
| R/prop_error.R | compute_difference | internal helper | src/pydron/prop_error.py | compute_difference | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 63 |
| R/raw_cf.R | raw_cf | exported function | src/pydron/raw_cf.py | raw_cf | raw_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 17 |
| R/raw_cf.R | raw_cf_meta | exported function | src/pydron/raw_cf.py | raw_cf_meta | raw_cf_meta | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 39 |
| R/raw_cf.R | raw_cf_data | exported function | src/pydron/raw_cf.py | raw_cf_data | raw_cf_data | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 67 |
| R/raw_cf.R | raw_cf_to_cf | exported function | src/pydron/raw_cf.py | raw_cf_to_cf | raw_cf_to_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 92 |
| R/raw_cf.R | uwerr.raw_cf | exported function | src/pydron/raw_cf.py | uwerr_raw_cf | uwerr_raw_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 136 |
| R/raw_cf.R | block.raw_cf | exported function | src/pydron/raw_cf.py | block_raw_cf | block_raw_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 216 |
| R/raw_cf.R | addStat.raw_cf | exported function | src/pydron/raw_cf.py | add_stat_raw_cf | add_stat_raw_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 274 |
| R/raw_cf.R | add.raw_cf | exported function | src/pydron/raw_cf.py | add_raw_cf | add_raw_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 312 |
| R/raw_cf.R | +.raw_cf | S3 method | src/pydron/raw_cf.py | add_raw_cf | add_raw_cf | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 335 |
| R/raw_cf.R | -.raw_cf | S3 method | src/pydron/raw_cf.py | sub_raw_cf | sub_raw_cf | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 345 |
| R/raw_cf.R | /.raw_cf | S3 method | src/pydron/raw_cf.py | div_raw_cf | div_raw_cf | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 355 |
| R/raw_cf.R | *.raw_cf | S3 method | src/pydron/raw_cf.py | mul_raw_cf | mul_raw_cf | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 375 |
| R/raw_cf.R | conj_raw_cf | internal helper | src/pydron/raw_cf.py | conj_raw_cf | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 392 |
| R/raw_cf.R | mul.raw_cf | exported function | src/pydron/raw_cf.py | mul_raw_cf | mul_raw_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 405 |
| R/raw_cf.R | is.raw_cf | exported function | src/pydron/raw_cf.py | is_raw_cf | is_raw_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 420 |
| R/raw_cf.R | is_empty.raw_cf | exported function | src/pydron/raw_cf.py | is_empty_raw_cf | is_empty_raw_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 431 |
| R/raw_cf.R | c.raw_cf | S3 method | src/pydron/raw_cf.py | c_raw_cf | c_raw_cf | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 446 |
| R/raw_cf.R | concat.raw_cf | exported function | src/pydron/raw_cf.py | concat_raw_cf | concat_raw_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 464 |
| R/raw_cf.R | get_plotdata_raw_cf | exported function | src/pydron/raw_cf.py | get_plotdata_raw_cf | get_plotdata_raw_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 531 |
| R/raw_cf.R | plot.raw_cf | S3 method | src/pydron/raw_cf.py | plot_raw_cf | plot_raw_cf | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 627 |
| R/raw_cf.R | overview_plot_raw_cf | exported function | src/pydron/raw_cf.py | overview_plot_raw_cf | overview_plot_raw_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 698 |
| R/raw_cf.R | shift.raw_cf | exported function | src/pydron/raw_cf.py | shift_raw_cf | shift_raw_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 839 |
| R/raw_cf.R | idx_matrix.raw_cf | exported function | src/pydron/raw_cf.py | idx_matrix_raw_cf | idx_matrix_raw_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 910 |
| R/raw_cf.R | int_idx_matrix.raw_cf | exported function | src/pydron/raw_cf.py | int_idx_matrix_raw_cf | int_idx_matrix_raw_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 940 |
| R/raw_cf.R | summary.raw_cf | S3 method | src/pydron/raw_cf.py | summary_raw_cf | summary_raw_cf | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 961 |
| R/raw_cf.R | print.raw_cf | S3 method | src/pydron/raw_cf.py | print_raw_cf | print_raw_cf | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 1002 |
| R/raw_cf.R | store_correl | internal helper | src/pydron/raw_cf.py | store_correl | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 1020 |
| R/readutils.R | readcmicor | exported function | src/pydron/readutils.py | readcmicor | readcmicor | `[ ]` | unassigned | I/O fixture plus R reference object comparison | NAMESPACE export; line 2 |
| R/readutils.R | getorderedfilelist | exported function | src/pydron/readutils.py | getorderedfilelist | getorderedfilelist | `[ ]` | unassigned | I/O fixture plus R reference object comparison | NAMESPACE export; line 35 |
| R/readutils.R | getconfignumbers | I/O helper | src/pydron/readutils.py | getconfignumbers | - | `[ ]` | unassigned | fixture read/write parity or explicit non-port decision | line 41 |
| R/readutils.R | getorderedconfigindices | I/O helper | src/pydron/readutils.py | getorderedconfigindices | - | `[ ]` | unassigned | fixture read/write parity or explicit non-port decision | line 54 |
| R/readutils.R | getorderedconfignumbers | exported function | src/pydron/readutils.py | getorderedconfignumbers | getorderedconfignumbers | `[ ]` | unassigned | I/O fixture plus R reference object comparison | NAMESPACE export; line 91 |
| R/readutils.R | readcmifiles | exported function | src/pydron/readutils.py | readcmifiles | readcmifiles | `[ ]` | unassigned | I/O fixture plus R reference object comparison | NAMESPACE export; line 212 |
| R/readutils.R | readcmidatafiles | exported function | src/pydron/readutils.py | readcmidatafiles | readcmidatafiles | `[ ]` | unassigned | I/O fixture plus R reference object comparison | NAMESPACE export; line 305 |
| R/readutils.R | readcmiloopfiles | exported function | src/pydron/readutils.py | readcmiloopfiles | readcmiloopfiles | `[ ]` | unassigned | I/O fixture plus R reference object comparison | NAMESPACE export; line 316 |
| R/readutils.R | extract.loop | exported function | src/pydron/readutils.py | extract_loop | extract_loop | `[ ]` | unassigned | I/O fixture plus R reference object comparison | NAMESPACE export; line 361 |
| R/readutils.R | extract.obs | exported function | src/pydron/readutils.py | extract_obs | extract_obs | `[ ]` | unassigned | I/O fixture plus R reference object comparison | NAMESPACE export; line 440 |
| R/readutils.R | readhlcor | exported function | src/pydron/readutils.py | readhlcor | readhlcor | `[ ]` | unassigned | I/O fixture plus R reference object comparison | NAMESPACE export; line 525 |
| R/readutils.R | readoutputdata | exported function | src/pydron/readutils.py | readoutputdata | readoutputdata | `[ ]` | unassigned | I/O fixture plus R reference object comparison | NAMESPACE export; line 553 |
| R/readutils.R | readtextcf | exported function | src/pydron/readutils.py | readtextcf | readtextcf | `[ ]` | unassigned | I/O fixture plus R reference object comparison | NAMESPACE export; line 601 |
| R/readutils.R | readnissatextcf | exported function | src/pydron/readutils.py | readnissatextcf | readnissatextcf | `[ ]` | unassigned | I/O fixture plus R reference object comparison | NAMESPACE export; line 699 |
| R/readutils.R | readbinarycf | exported function | src/pydron/readutils.py | readbinarycf | readbinarycf | `[ ]` | unassigned | I/O fixture plus R reference object comparison | NAMESPACE export; line 789 |
| R/readutils.R | readbinarysamples | exported function | src/pydron/readutils.py | readbinarysamples | readbinarysamples | `[ ]` | unassigned | I/O fixture plus R reference object comparison | NAMESPACE export; line 905 |
| R/readutils.R | readbinarydisc | exported function | src/pydron/readutils.py | readbinarydisc | readbinarydisc | `[ ]` | unassigned | I/O fixture plus R reference object comparison | NAMESPACE export; line 1007 |
| R/readutils.R | readcmidisc | exported function | src/pydron/readutils.py | readcmidisc | readcmidisc | `[ ]` | unassigned | I/O fixture plus R reference object comparison | NAMESPACE export; line 1084 |
| R/readutils.R | readgradflow | exported function | src/pydron/readutils.py | readgradflow | readgradflow | `[ ]` | unassigned | I/O fixture plus R reference object comparison | NAMESPACE export; line 1172 |
| R/removeTemporal.cf.R | old_removeTemporal.cf | exported function | src/pydron/remove_temporal_cf.py | old_remove_temporal_cf | old_remove_temporal_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 33 |
| R/removeTemporal.cf.R | takeTimeDiff.cf | exported function | src/pydron/remove_temporal_cf.py | take_time_diff_cf | take_time_diff_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 188 |
| R/removeTemporal.cf.R | dispersion_relation | exported function | src/pydron/remove_temporal_cf.py | dispersion_relation | dispersion_relation | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 280 |
| R/removeTemporal.cf.R | extract_mass | generic | src/pydron/remove_temporal_cf.py | extract_mass | extract_mass | `[ ]` | unassigned | dispatch parity tests for registered methods | NAMESPACE export; UseMethod: extract_mass; line 308 |
| R/removeTemporal.cf.R | extract_mass.effectivemassfit | S3 method | src/pydron/remove_temporal_cf.py | extract_mass_effectivemassfit | extract_mass_effectivemassfit | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 320 |
| R/removeTemporal.cf.R | extract_mass.matrixfit | S3 method | src/pydron/remove_temporal_cf.py | extract_mass_matrixfit | extract_mass_matrixfit | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 333 |
| R/removeTemporal.cf.R | make_weight_factor | internal helper | src/pydron/remove_temporal_cf.py | make_weight_factor | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 338 |
| R/removeTemporal.cf.R | weight.cf | exported function | src/pydron/remove_temporal_cf.py | weight_cf | weight_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 367 |
| R/removeTemporal.cf.R | weight_shift_reweight.cf | exported function | src/pydron/remove_temporal_cf.py | weight_shift_reweight_cf | weight_shift_reweight_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 401 |
| R/removeTemporal.cf.R | removeTemporal.cf | exported function | src/pydron/remove_temporal_cf.py | remove_temporal_cf | remove_temporal_cf | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 448 |
| R/seed.R | swap_seed | internal helper | src/pydron/seed.py | swap_seed | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 13 |
| R/seed.R | restore_seed | internal helper | src/pydron/seed.py | restore_seed | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 32 |
| R/string2error.R | string2error | exported function | src/pydron/string2error.py | string2error | string2error | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 23 |
| R/summary.ofit.R | print.ofit | S3 method | src/pydron/summary_ofit.py | print_ofit | print_ofit | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 10 |
| R/summary.ofit.R | summary.ofit | S3 method | src/pydron/summary_ofit.py | summary_ofit | summary_ofit | `[ ]` | unassigned | R dispatch and R reference-output equivalence where numerical | NAMESPACE S3 method; line 23 |
| R/tex-utils.R | tex.catwitherror | exported function | src/pydron/tex_utils.py | tex_catwitherror | tex_catwitherror | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 37 |
| R/tex-utils.R | escapeLatexSpecials | exported function | src/pydron/tex_utils.py | escape_latex_specials | escape_latex_specials | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 98 |
| R/tflops.R | tflops | internal helper | src/pydron/tflops.py | tflops | - | `[ ]` | unassigned | unit test if promoted to public/internal Python helper | line 1 |
| R/tikzutils.R | tikz.init | exported function | src/pydron/tikzutils.py | tikz_init | tikz_init | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 32 |
| R/tikzutils.R | tikz.finalize | exported function | src/pydron/tikzutils.py | tikz_finalize | tikz_finalize | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 76 |
| R/timeseries.R | plot_timeseries | exported function | src/pydron/timeseries.py | plot_timeseries | plot_timeseries | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 35 |
| R/timeseries.R | plot_eigenvalue_timeseries | exported function | src/pydron/timeseries.py | plot_eigenvalue_timeseries | plot_eigenvalue_timeseries | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 184 |
| R/zeta_zp.R | zetazp | exported function | src/pydron/zeta_zp.py | zetazp | zetazp | `[ ]` | unassigned | R reference-output equivalence test required if numerical | NAMESPACE export; line 28 |

## Documentation pages (`man/*.Rd`)

| Status | Source path | Rd name | Aliases | Matched R/data symbols | Proposed Python doc target | Owner | Title |
|---|---|---|---|---|---|---|---|
| `[ ]` | man/add.cf.Rd | add.cf | add.cf | add.cf | docs/api/add_cf.md | unassigned | Arithmetically adds two correlation functions |
| `[ ]` | man/add.raw_cf.Rd | add.raw_cf | add.raw_cf | add.raw_cf | docs/api/add_raw_cf.md | unassigned | add two raw_cf objects |
| `[ ]` | man/addConfIndex2cf.Rd | addConfIndex2cf | addConfIndex2cf | addConfIndex2cf | docs/api/add_conf_index2cf.md | unassigned | add a configuration index to an cf object |
| `[ ]` | man/addStat.cf.Rd | addStat.cf | addStat.cf | addStat.cf | docs/api/add_stat_cf.md | unassigned | Combine statistics of two cf objects |
| `[ ]` | man/addStat.raw_cf.Rd | addStat.raw_cf | addStat.raw_cf | addStat.raw_cf | docs/api/add_stat_raw_cf.md | unassigned | Extend statistics of an existing raw_cf container |
| `[ ]` | man/alphas.Rd | alphas | alphas | alphas | docs/api/alphas.md | unassigned | compute alpha strong at given scale |
| `[ ]` | man/analysis_gradient_flow.Rd | analysis_gradient_flow | analysis_gradient_flow | analysis_gradient_flow | docs/api/analysis_gradient_flow.md | unassigned | analysis_gradient_flow |
| `[ ]` | man/analysis_online.Rd | analysis_online | analysis_online | analysis_online | docs/api/analysis_online.md | unassigned | analysis_online |
| `[ ]` | man/avg.cbt.cf.Rd | avg.cbt.cf | avg.cbt.cf | avg.cbt.cf | docs/api/avg_cbt_cf.md | unassigned | average close-by-times in a correlation function |
| `[ ]` | man/block.raw_cf.Rd | block.raw_cf | block.raw_cf | block.raw_cf | docs/api/block_raw_cf.md | unassigned | Block average correlation function data |
| `[ ]` | man/boot_ts_array.Rd | boot_ts_array | boot_ts_array | boot_ts_array | docs/api/boot_ts_array.md | unassigned | boot_ts_array |
| `[ ]` | man/bootstrap.analysis.Rd | bootstrap.analysis | bootstrap.analysis | bootstrap.analysis | docs/api/bootstrap_analysis.md | unassigned | Performs a Bootstrap with Blocking Analysis of a Timeseries |
| `[ ]` | man/bootstrap.cf.Rd | bootstrap.cf | bootstrap.cf | bootstrap.cf | docs/api/bootstrap_cf.md | unassigned | bootstrap a set of correlation functions |
| `[ ]` | man/bootstrap.effectivemass.Rd | bootstrap.effectivemass | bootstrap.effectivemass | bootstrap.effectivemass | docs/api/bootstrap_effectivemass.md | unassigned | Computes effective masses with bootstrapping errors |
| `[ ]` | man/bootstrap.gevp.Rd | bootstrap.gevp | bootstrap.gevp | bootstrap.gevp | docs/api/bootstrap_gevp.md | unassigned | perform a bootstrap analysis of a GEVP |
| `[ ]` | man/bootstrap.hankel_summed.Rd | bootstrap.hankel_summed | bootstrap.hankel_summed | bootstrap.hankel_summed | docs/api/bootstrap_hankel_summed.md | unassigned | GEVP method based on Hankel matrices. |
| `[ ]` | man/bootstrap.hankel.Rd | bootstrap.hankel | bootstrap.hankel | bootstrap.hankel | docs/api/bootstrap_hankel.md | unassigned | GEVP method based on Hankel matrices. |
| `[ ]` | man/bootstrap.lanczos.Rd | bootstrap.lanczos | bootstrap.lanczos | bootstrap.lanczos | docs/api/bootstrap_lanczos.md | unassigned | Lanczos method for LQCD correlators |
| `[ ]` | man/bootstrap.meanerror.Rd | bootstrap.meanerror | bootstrap.meanerror | bootstrap.meanerror | docs/api/bootstrap_meanerror.md | unassigned | Compute the bootstrap error of the mean |
| `[ ]` | man/bootstrap.nlsfit.Rd | bootstrap.nlsfit | bootstrap.nlsfit | bootstrap.nlsfit | docs/api/bootstrap_nlsfit.md | unassigned | Bootstrap a non-linear least-squares fit |
| `[ ]` | man/bootstrap.pgevm.Rd | bootstrap.pgevm | bootstrap.pgevm | bootstrap.pgevm | docs/api/bootstrap_pgevm.md | unassigned | PGEVM |
| `[ ]` | man/bootstrap.truncated.pgevm.Rd | bootstrap.truncated.pgevm | bootstrap.truncated.pgevm | bootstrap.truncated.pgevm | docs/api/bootstrap_truncated_pgevm.md | unassigned | Truncated PGEVM |
| `[ ]` | man/c.cf.Rd | c.cf | c.cf | c.cf | docs/api/c_cf.md | unassigned | Concatenate correlation function objects |
| `[ ]` | man/c.raw_cf.Rd | c.raw_cf | c.raw_cf | c.raw_cf | docs/api/c_raw_cf.md | unassigned | Concatenate raw_cf correlation function objects |
| `[ ]` | man/cA2.09.48_3pi_I3_0_A1u_1_pc.Rd | cA2.09.48_3pi_I3_0_A1u_1_pc | cA2.09.48_3pi_I3_0_A1u_1_pc | cA2.09.48_3pi_I3_0_A1u_1_pc | docs/api/c_a2_09_48_3pi_i3_0_a1u_1_pc.md | unassigned | A three pion correlator with significant thermal states. |
| `[ ]` | man/cdh.Rd | cdh | cdh | cdh | docs/api/cdh.md | unassigned | finite size corrections a la Colangelo, Duerr, Haefeli |
| `[ ]` | man/cdhnew.Rd | cdhnew | cdhnew | cdhnew | docs/api/cdhnew.md | unassigned | finite size corrections a la Colangelo, Duerr, Haefeli, but re-expanded as<br>series in the quark mass |
| `[ ]` | man/CExp.Rd | CExp | CExp | CExp | docs/api/cexp.md | unassigned | Cosh Or Sinh Build Out Of Two Exps |
| `[ ]` | man/cf_boot.Rd | cf_boot | cf_boot | cf_boot | docs/api/cf_boot.md | unassigned | Bootstrapped CF mixin constructor |
| `[ ]` | man/cf_key_meson_2pt.Rd | cf_key_meson_2pt | cf_key_meson_2pt | cf_key_meson_2pt | docs/api/cf_key_meson_2pt.md | unassigned | Generate key string to identify a meson 2pt function |
| `[ ]` | man/cf_key_meson_3pt.Rd | cf_key_meson_3pt | cf_key_meson_3pt | cf_key_meson_3pt | docs/api/cf_key_meson_3pt.md | unassigned | Generate HDF5 key for CVC 'correlators' meson 3pt function with a local or derivative insertion |
| `[ ]` | man/cf_meta.Rd | cf_meta | cf_meta | cf_meta | docs/api/cf_meta.md | unassigned | CF metadata mixin constructor |
| `[ ]` | man/cf_orig.Rd | cf_orig | cf_orig | cf_orig | docs/api/cf_orig.md | unassigned | Original data CF mixin constructor |
| `[ ]` | man/cf_principal_correlator.Rd | cf_principal_correlator | cf_principal_correlator | cf_principal_correlator | docs/api/cf_principal_correlator.md | unassigned | Principal correlator CF mixin constructor |
| `[ ]` | man/cf_shifted.Rd | cf_shifted | cf_shifted | cf_shifted | docs/api/cf_shifted.md | unassigned | Shifted CF mixin constructor |
| `[ ]` | man/cf_smeared.Rd | cf_smeared | cf_smeared | cf_smeared | docs/api/cf_smeared.md | unassigned | Smeared CF mixin constructor |
| `[ ]` | man/cf_subtracted.Rd | cf_subtracted | cf_subtracted | cf_subtracted | docs/api/cf_subtracted.md | unassigned | Subtracted CF mixin constructor |
| `[ ]` | man/cf_weighted.Rd | cf_weighted | cf_weighted | cf_weighted | docs/api/cf_weighted.md | unassigned | Weighted CF mixin constructor |
| `[ ]` | man/cf.Rd | cf | cf | cf | docs/api/cf.md | unassigned | Correlation function container |
| `[ ]` | man/compute.plotlims.Rd | compute.plotlims | compute.plotlims | compute.plotlims | docs/api/compute_plotlims.md | unassigned | compute.plotlims |
| `[ ]` | man/computeacf.Rd | computeacf | computeacf | computeacf | docs/api/computeacf.md | unassigned | Computes The ACF and Integrated AC Time |
| `[ ]` | man/computeDisc.Rd | computeDisc | computeDisc | computeDisc | docs/api/compute_disc.md | unassigned | computes a disconnected correlation function from loops |
| `[ ]` | man/computefps.Rd | computefps | computefps | computefps | docs/api/computefps.md | unassigned | Computes the pseudoscalar decay constant for the twisted mass case from the<br>pseudoscalar amplitude and mass |
| `[ ]` | man/computefpsOS.Rd | computefpsOS | computefpsOS | computefpsOS | docs/api/computefps_os.md | unassigned | Computes the pseudoscalar decay constant for the Osterwalder Seiler case<br>from the pseudoscalar amplitude and mass |
| `[ ]` | man/concat.cf.Rd | concat.cf | concat.cf | concat.cf | docs/api/concat_cf.md | unassigned | Concatenate two correlation function objects |
| `[ ]` | man/concat.raw_cf.Rd | concat.raw_cf | concat.raw_cf | concat.raw_cf | docs/api/concat_raw_cf.md | unassigned | Concatenate two raw_cf correlation function objects |
| `[ ]` | man/conj_raw_cf.Rd | conj_raw_cf | conj_raw_cf | conj_raw_cf | docs/api/conj_raw_cf.md | unassigned | Take the complex conjugate of a raw_cf object |
| `[ ]` | man/construct_onlinemeas_rundir.Rd | construct_onlinemeas_rundir | construct_onlinemeas_rundir | construct_onlinemeas_rundir | docs/api/construct_onlinemeas_rundir.md | unassigned | Construct a run directory string for analysis_online |
| `[ ]` | man/correlatormatrix.Rd | correlatormatrix | correlatormatrix | correlatormatrix | docs/api/correlatormatrix.md | unassigned | Sample correlator matrix |
| `[ ]` | man/correlators_key_meson_2pt.Rd | correlators_key_meson_2pt | correlators_key_meson_2pt | correlators_key_meson_2pt | docs/api/correlators_key_meson_2pt.md | unassigned | Generate HDF5 key for CVC 'correlators' meson 2pt function |
| `[ ]` | man/correlators_key_meson_3pt.Rd | correlators_key_meson_3pt | correlators_key_meson_3pt | correlators_key_meson_3pt | docs/api/correlators_key_meson_3pt.md | unassigned | Generate HDF5 key for CVC 'correlators' meson 3pt function with a local or derivative insertion |
| `[ ]` | man/cosh_to_effmass.Rd | cosh_to_effmass | cosh_to_effmass | cosh_to_effmass | docs/api/cosh_to_effmass.md | unassigned | extract the effective mass from the sum of coshs |
| `[ ]` | man/create_displ_chains.Rd | create_displ_chains | create_displ_chains | create_displ_chains | docs/api/create_displ_chains.md | unassigned | create list of chains of displacements<br>Multilpe covariant displacements, when applied in order, form<br>a list of displacments. Each consists of a direction and a dimension. |
| `[ ]` | man/cvc_local_loop_key.Rd | cvc_local_loop_key | cvc_local_loop_key | cvc_local_loop_key | docs/api/cvc_local_loop_key.md | unassigned | Generate HDF5 key for a momentum and spin-projected CVC loop |
| `[ ]` | man/cvc_read_loops.Rd | cvc_read_loops | cvc_read_loops | cvc_read_loops | docs/api/cvc_read_loops.md | unassigned | read HDF5 loop files in the CVC loop format |
| `[ ]` | man/cvc_to_raw_cf.Rd | cvc_to_raw_cf | cvc_to_raw_cf | cvc_to_raw_cf | docs/api/cvc_to_raw_cf.md | unassigned | Convert correlation function read from CVC HDF5 or AFF format to 'raw_cf' |
| `[ ]` | man/cyprus_make_key_scalar.Rd | cyprus_make_key_scalar | cyprus_make_key_scalar | cyprus_make_key_scalar | docs/api/cyprus_make_key_scalar.md | unassigned | HDF5 key for Cyprus CalcLoops scalar-type loops |
| `[ ]` | man/cyprus_make_key_vector.Rd | cyprus_make_key_vector | cyprus_make_key_vector | cyprus_make_key_vector | docs/api/cyprus_make_key_vector.md | unassigned | HDF5 key for Cyprus CalcLoops derivative-type loops |
| `[ ]` | man/cyprus_read_loops.Rd | cyprus_read_loops | cyprus_read_loops | cyprus_read_loops | docs/api/cyprus_read_loops.md | unassigned | read HDF5 loop files in the Cyprus CalcLoops format |
| `[ ]` | man/deriv.CExp.Rd | deriv.CExp | deriv.CExp | deriv.CExp | docs/api/deriv_cexp.md | unassigned | the calling code must supply the correct three parameters |
| `[ ]` | man/disc_3pt.Rd | disc_3pt | disc_3pt | disc_3pt | docs/api/disc_3pt.md | unassigned | disconnected contribution to current insertion three-point function |
| `[ ]` | man/dispersion_relation.Rd | dispersion_relation | dispersion_relation | dispersion_relation | docs/api/dispersion_relation.md | unassigned | Continuum dispersion relation for CM to lattice frame |
| `[ ]` | man/dot-cf.Rd | -.cf | -.cf | -.cf | docs/api/sub_cf.md | unassigned | Arithmetically subtract correlators |
| `[ ]` | man/dot-raw_cf.Rd | -.raw_cf | -.raw_cf | -.raw_cf | docs/api/sub_raw_cf.md | unassigned | add two raw_cf objects |
| `[ ]` | man/double_bootstrap.cf.Rd | double_bootstrap.cf | double_bootstrap.cf | double_bootstrap.cf | docs/api/double_bootstrap_cf.md | unassigned | double bootstrap function for cf |
| `[ ]` | man/effectivemass.cf.Rd | effectivemass.cf | effectivemass.cf | effectivemass.cf | docs/api/effectivemass_cf.md | unassigned | Computes effective mass values for a correlation function |
| `[ ]` | man/effectivemass.Rd | effectivemass | effectivemass | effectivemass | docs/api/effectivemass.md | unassigned | effectivemass |
| `[ ]` | man/effmass.Rd | effmass | effmass | effmass | docs/api/effmass.md | unassigned | effmass |
| `[ ]` | man/effmass2.Rd | effmass2 | effmass2 | effmass2 | docs/api/effmass2.md | unassigned | effmass2 |
| `[ ]` | man/escapeLatexSpecials.Rd | escapeLatexSpecials | escapeLatexSpecials | escapeLatexSpecials | docs/api/escape_latex_specials.md | unassigned | Escape special LaTeX characters for use in LaTeX labels |
| `[ ]` | man/extract_mass.effectivemassfit.Rd | extract_mass.effectivemassfit | extract_mass.effectivemassfit | extract_mass.effectivemassfit | docs/api/extract_mass_effectivemassfit.md | unassigned | specialisation of extract_mass to objects of type<br>effectivemassfit |
| `[ ]` | man/extract_mass.matrixfit.Rd | extract_mass.matrixfit | extract_mass.matrixfit | extract_mass.matrixfit | docs/api/extract_mass_matrixfit.md | unassigned | specialisation of extract_mass to objects of type<br>matrixfit |
| `[ ]` | man/extract_mass.Rd | extract_mass | extract_mass | extract_mass | docs/api/extract_mass.md | unassigned | generic function to extract a fitted mass |
| `[ ]` | man/extract.loop.Rd | extract.loop | extract.loop | extract.loop | docs/api/extract_loop.md | unassigned | Extract a single loop from an object of class cmiloop |
| `[ ]` | man/extract.obs.Rd | extract.obs | extract.obs | extract.obs | docs/api/extract_obs.md | unassigned | Extract One or More Gamma Combinations from am CMI Correlator |
| `[ ]` | man/extractSingleCor.cf.Rd | extractSingleCor.cf | extractSingleCor.cf | extractSingleCor.cf | docs/api/extract_single_cor_cf.md | unassigned | extract one single correlator object as cf object from a large<br>cf object. |
| `[ ]` | man/fit.cosh.Rd | fit.cosh | fit.cosh | fit.cosh | docs/api/fit_cosh.md | unassigned | Fits a sum of several cosh-functions |
| `[ ]` | man/fit.effectivemass.Rd | fit.effectivemass | fit.effectivemass | fit.effectivemass | docs/api/fit_effectivemass.md | unassigned | Fits a constant to effective mass data |
| `[ ]` | man/fit.plateau2cf.Rd | fit.plateau2cf | fit.plateau2cf | fit.plateau2cf | docs/api/fit_plateau2cf.md | unassigned | fits a plateau to an object of class cf |
| `[ ]` | man/foldr1.Rd | foldr1 | foldr1 | foldr1 | docs/api/foldr1.md | unassigned | Folds the non-empty list with the binary function |
| `[ ]` | man/fs.a0.Rd | fs.a0 | fs.a0 | fs.a0 | docs/api/fs_a0.md | unassigned | Finite Size Corrections to q\cot\deltaqcotdelta for I=2<br>\pi\pipipi near threshold |
| `[ ]` | man/fs.mpia0.Rd | fs.mpia0 | fs.mpia0 | fs.mpia0 | docs/api/fs_mpia0.md | unassigned | Finite Size Corrections to q\cot\deltaqcotdelta for I=2<br>\pi\pipipi near threshold |
| `[ ]` | man/fs.qcotdelta.Rd | fs.qcotdelta | fs.qcotdelta | fs.qcotdelta | docs/api/fs_qcotdelta.md | unassigned | Finite Size Corrections to q\cot\deltaqcotdelta for I=2<br>\pi\pipipi near threshold |
| `[ ]` | man/g1.Rd | g1 | g1 | g1 | docs/api/g1.md | unassigned | g1 |
| `[ ]` | man/get_plotdata_raw_cf.Rd | get_plotdata_raw_cf | get_plotdata_raw_cf | get_plotdata_raw_cf | docs/api/get_plotdata_raw_cf.md | unassigned | extract data from 'raw_cf' in format convenient to plot |
| `[ ]` | man/getorderedconfignumbers.Rd | getorderedconfignumbers | getorderedconfignumbers | getorderedconfignumbers | docs/api/getorderedconfignumbers.md | unassigned | Creates an ordered vector of gauge config file numbers |
| `[ ]` | man/getorderedfilelist.Rd | getorderedfilelist | getorderedfilelist | getorderedfilelist | docs/api/getorderedfilelist.md | unassigned | Creates an ordered filelist from a basename and a path |
| `[ ]` | man/gevp.hankel_summed.Rd | gevp.hankel_summed | gevp.hankel_summed | gevp.hankel_summed | docs/api/gevp_hankel_summed.md | unassigned | GEVP method based on Hankel matrices. |
| `[ ]` | man/gevp.hankel.Rd | gevp.hankel | gevp.hankel | gevp.hankel | docs/api/gevp_hankel.md | unassigned | GEVP method based on Hankel matrices. |
| `[ ]` | man/gevp.Rd | gevp | gevp | gevp | docs/api/gevp.md | unassigned | solve GEVP for correlator matrix |
| `[ ]` | man/gevp.truncated.hankel.Rd | gevp.truncated.hankel | gevp.truncated.hankel | gevp.truncated.hankel | docs/api/gevp_truncated_hankel.md | unassigned | GEVP method based on truncated Hankel matrices. |
| `[ ]` | man/gevp2amplitude.Rd | gevp2amplitude | gevp2amplitude | gevp2amplitude | docs/api/gevp2amplitude.md | unassigned | Extracts physical amplitudes from a GEVP |
| `[ ]` | man/gevp2cf.Rd | gevp2cf | gevp2cf | gevp2cf | docs/api/gevp2cf.md | unassigned | Extracts a principle correlator from a GEVEP |
| `[ ]` | man/gm_mu.Rd | gm_mu | gm_mu | gm_mu | docs/api/gm_mu.md | unassigned | Accessor function for gm |
| `[ ]` | man/gm.Rd | gm | gm | - | docs/api/gm.md | unassigned | List of arrays of gamma structures |
| `[ ]` | man/h5_get_dataset.Rd | h5_get_dataset | h5_get_dataset | h5_get_dataset | docs/api/h5_get_dataset.md | unassigned | get dataset from HDF5 file |
| `[ ]` | man/h5_names_exist.Rd | h5_names_exist | h5_names_exist | h5_names_exist | docs/api/h5_names_exist.md | unassigned | check if group names exist in HDF5 file |
| `[ ]` | man/hadron.Rd | hadron | hadron; hadron-package | - | docs/api/hadron.md | unassigned | The Hadron Package |
| `[ ]` | man/hankel2cf.Rd | hankel2cf | hankel2cf | hankel2cf | docs/api/hankel2cf.md | unassigned | hankel2cf |
| `[ ]` | man/hankel2effectivemass.Rd | hankel2effectivemass | hankel2effectivemass | hankel2effectivemass | docs/api/hankel2effectivemass.md | unassigned | hankel2effectivemass |
| `[ ]` | man/hankeldensity2effectivemass.Rd | hankeldensity2effectivemass | hankeldensity2effectivemass | hankeldensity2effectivemass | docs/api/hankeldensity2effectivemass.md | unassigned | hankeldensity2effectivemass |
| `[ ]` | man/has_icf.Rd | has_icf | has_icf | has_icf | docs/api/has_icf.md | unassigned | Checks whether the cf object contains an imaginary part |
| `[ ]` | man/idx_matrix.raw_cf.Rd | idx_matrix.raw_cf | idx_matrix.raw_cf | idx_matrix.raw_cf | docs/api/idx_matrix_raw_cf.md | unassigned | Construct the tensor index set for the entire raw correlator |
| `[ ]` | man/int_idx_matrix.raw_cf.Rd | int_idx_matrix.raw_cf | int_idx_matrix.raw_cf | int_idx_matrix.raw_cf | docs/api/int_idx_matrix_raw_cf.md | unassigned | Construct tensor index set for the internal degrees of freedom |
| `[ ]` | man/InternalHadronFunctions.Rd | InternalHadronFunctions | InternalHadronFunctions; arrangeCor.vector; arrangeCor.pion; arrangeCor.b1; arrangeCor.a0; getNxNmatrix; ChiSqr.1mass; ChiSqr.2mass; ChiSqr.3mass; fitmasses.vector; fitmasses.vector.boot; fitmasses.pion; fitmasses.pion.boot; fitmasses.b1; fitmasses.b1.boot; fitmasses.a0; fitmasses.a0.boot; fitf.pion; mean.index; fitmpcac.pion; read_nissa_textcf_kernel; compute.qtildesq; compute.qtildesq.contdisp | getNxNmatrix; ChiSqr.1mass; ChiSqr.2mass; ChiSqr.3mass; read_nissa_textcf_kernel; compute.qtildesq; compute.qtildesq.contdisp | docs/api/internal_hadron_functions.md | unassigned | Internal Hadron Functions |
| `[ ]` | man/invalidate.samples.cf.Rd | invalidate.samples.cf | invalidate.samples.cf | invalidate.samples.cf | docs/api/invalidate_samples_cf.md | unassigned | Invalidate samples |
| `[ ]` | man/invcosh.Rd | invcosh | invcosh | invcosh | docs/api/invcosh.md | unassigned | numerically invert the cosh function for the mass |
| `[ ]` | man/invertCovMatrix.Rd | invertCovMatrix | invertCovMatrix | invertCovMatrix | docs/api/invert_cov_matrix.md | unassigned | Inverts the covariance matrix for noisy data |
| `[ ]` | man/is_empty.cf.Rd | is_empty.cf | is_empty.cf | is_empty.cf | docs/api/is_empty_cf.md | unassigned | Checks whether the cf object contains no data |
| `[ ]` | man/is_empty.raw_cf.Rd | is_empty.raw_cf | is_empty.raw_cf | is_empty.raw_cf | docs/api/is_empty_raw_cf.md | unassigned | check if an obect is of class raw_cf and empty otherwise |
| `[ ]` | man/is.cf.Rd | is.cf | is.cf | is.cf | docs/api/is_cf.md | unassigned | Checks whether an object is a cf |
| `[ ]` | man/is.raw_cf.Rd | is.raw_cf | is.raw_cf | is.raw_cf | docs/api/is_raw_cf.md | unassigned | check if an object is of class raw_cf |
| `[ ]` | man/jackknife_cov.Rd | jackknife_cov | jackknife_cov | jackknife_cov | docs/api/jackknife_cov.md | unassigned | jackknife_cov |
| `[ ]` | man/jackknife_error.Rd | jackknife_error | jackknife_error | jackknife_error | docs/api/jackknife_error.md | unassigned | Estimates error from jackknife samples |
| `[ ]` | man/jackknife-after-bootstrap.Rd | jackknife-after-bootstrap | jackknife-after-bootstrap; jackknifeafterboot; jab.cf; jab.cf.derived; jab.effectivemass; jab.effectivemassfit; jab.matrixfit | jab.cf; jab.cf.derived; jab.effectivemass; jab.effectivemassfit; jab.matrixfit | docs/api/jackknife_after_bootstrap.md | unassigned | jackknife-after-bootstrap analysis |
| `[ ]` | man/jackknife.cf.Rd | jackknife.cf | jackknife.cf | jackknife.cf | docs/api/jackknife_cf.md | unassigned | jackknife a set of correlation functions |
| `[ ]` | man/lanczos.solve.Rd | lanczos.solve | lanczos.solve | lanczos.solve | docs/api/lanczos_solve.md | unassigned | Lanczos solver |
| `[ ]` | man/loop_2pt.Rd | loop_2pt | loop_2pt | loop_2pt | docs/api/loop_2pt.md | unassigned | compute two-point correlation function between quark loops |
| `[ ]` | man/loop_spin_project.Rd | loop_spin_project | loop_spin_project | loop_spin_project | docs/api/loop_spin_project.md | unassigned | spin projection of quark loop data |
| `[ ]` | man/loop_stochav.Rd | loop_stochav | loop_stochav | loop_stochav | docs/api/loop_stochav.md | unassigned | average over stochastic samples of loop |
| `[ ]` | man/loop_vev_subtract.Rd | loop_vev_subtract | loop_vev_subtract | loop_vev_subtract | docs/api/loop_vev_subtract.md | unassigned | subtract vev from loop data |
| `[ ]` | man/loopdata.Rd | loopdata | loopdata | loopdata | docs/api/loopdata.md | unassigned | Sample loop data |
| `[ ]` | man/make_parind.Rd | make_parind | make_parind | make_parind | docs/api/make_parind.md | unassigned | Create a parameter index matrix for matrixfit |
| `[ ]` | man/make_parlist.Rd | make_parlist | make_parlist | make_parlist | docs/api/make_parlist.md | unassigned | Create a parameter list for matrixfit |
| `[ ]` | man/matrixfit.Rd | matrixfit | matrixfit | matrixfit | docs/api/matrixfit.md | unassigned | Routine For A Factorising Matrix Fit |
| `[ ]` | man/matrixModel.Rd | matrixModel | matrixModel | matrixModel | docs/api/matrix_model.md | unassigned | Correlator matrix model. |
| `[ ]` | man/mom_combinations.Rd | mom_combinations | mom_combinations | mom_combinations | docs/api/mom_combinations.md | unassigned | Generate table of momentum component combinations |
| `[ ]` | man/mul.cf.Rd | mul.cf | mul.cf | mul.cf | docs/api/mul_cf.md | unassigned | Arithmetically scale a correlator by a scalar a |
| `[ ]` | man/mul.raw_cf.Rd | mul.raw_cf | mul.raw_cf | mul.raw_cf | docs/api/mul_raw_cf.md | unassigned | scale raw_cf data |
| `[ ]` | man/new_matrixfit.Rd | new_matrixfit | new_matrixfit | new_matrixfit | docs/api/new_matrixfit.md | unassigned | perform a factorising fit of a matrix of correlation functions |
| `[ ]` | man/old_removeTemporal.cf.Rd | old_removeTemporal.cf | old_removeTemporal.cf | old_removeTemporal.cf | docs/api/old_remove_temporal_cf.md | unassigned | Remove temporal states |
| `[ ]` | man/onlinemeas.Rd | onlinemeas | onlinemeas | onlinemeas | docs/api/onlinemeas.md | unassigned | determines pion mass and pcac mass from online measured correlator of the<br>HMC code |
| `[ ]` | man/overview_plot_raw_cf.Rd | overview_plot_raw_cf | overview_plot_raw_cf | overview_plot_raw_cf | docs/api/overview_plot_raw_cf.md | unassigned | create convenient overview plots for a raw_cf object |
| `[ ]` | man/parametric.bootstrap.cov.Rd | parametric.bootstrap.cov | parametric.bootstrap.cov | parametric.bootstrap.cov | docs/api/parametric_bootstrap_cov.md | unassigned | Parametric bootstrap with covariance |
| `[ ]` | man/parametric.bootstrap.Rd | parametric.bootstrap | parametric.bootstrap | parametric.bootstrap | docs/api/parametric_bootstrap.md | unassigned | Parametric bootstrap |
| `[ ]` | man/parametric.nlsfit.cov.Rd | parametric.nlsfit.cov | parametric.nlsfit.cov | parametric.nlsfit.cov | docs/api/parametric_nlsfit_cov.md | unassigned | parametric.nlsfit.cov |
| `[ ]` | man/parametric.nlsfit.Rd | parametric.nlsfit | parametric.nlsfit | parametric.nlsfit | docs/api/parametric_nlsfit.md | unassigned | NLS fit with parametric bootstrap |
| `[ ]` | man/pcac.Rd | pcac | pcac | pcac | docs/api/pcac.md | unassigned | Computes the pcac mass |
| `[ ]` | man/pcacfit.Rd | pcacfit | pcacfit | pcacfit | docs/api/pcacfit.md | unassigned | pcacfit |
| `[ ]` | man/pcModel.Rd | pcModel | pcModel | pcModel | docs/api/pc_model.md | unassigned | Principal correlator two state model. |
| `[ ]` | man/pgevm2bootstrapfit.Rd | pgevm2bootstrapfit | pgevm2bootstrapfit | pgevm2bootstrapfit | docs/api/pgevm2bootstrapfit.md | unassigned | pgevm2bootstrapfit |
| `[ ]` | man/pgevm2effectivemass.Rd | pgevm2effectivemass | pgevm2effectivemass | pgevm2effectivemass | docs/api/pgevm2effectivemass.md | unassigned | pgevm2effectivemass |
| `[ ]` | man/plaq.sample.Rd | plaq.sample | plaq.sample | plaq.sample | docs/api/plaq_sample.md | unassigned | Sample plaquette time series |
| `[ ]` | man/plot_eigenvalue_timeseries.Rd | plot_eigenvalue_timeseries | plot_eigenvalue_timeseries | plot_eigenvalue_timeseries | docs/api/plot_eigenvalue_timeseries.md | unassigned | plot_eigenvalue_timeseries |
| `[ ]` | man/plot_hankel_spectrum.Rd | plot_hankel_spectrum | plot_hankel_spectrum | plot_hankel_spectrum | docs/api/plot_hankel_spectrum.md | unassigned | plot_hankel_spectrum |
| `[ ]` | man/plot_timeseries.Rd | plot_timeseries | plot_timeseries | plot_timeseries | docs/api/plot_timeseries.md | unassigned | plot_timeseries |
| `[ ]` | man/plot.averx.Rd | plot.averx | plot.averx | plot.averx | docs/api/plot_averx.md | unassigned | Plots averx data |
| `[ ]` | man/plot.bootstrapfit.Rd | plot.bootstrapfit | plot.bootstrapfit | plot.bootstrapfit | docs/api/plot_bootstrapfit.md | unassigned | Plot a bootstrap NLS fit |
| `[ ]` | man/plot.cf.Rd | plot.cf | plot.cf | plot.cf | docs/api/plot_cf.md | unassigned | Plot a correlation function |
| `[ ]` | man/plot.cfit.Rd | plot.cfit | plot.cfit | plot.cfit | docs/api/plot_cfit.md | unassigned | plot.c1fit |
| `[ ]` | man/plot.coshfit.Rd | plot.coshfit | plot.coshfit | plot.coshfit | docs/api/plot_coshfit.md | unassigned | Plot a cosh-fit |
| `[ ]` | man/plot.effectivemass.Rd | plot.effectivemass | plot.effectivemass | plot.effectivemass | docs/api/plot_effectivemass.md | unassigned | plot.effectivemass |
| `[ ]` | man/plot.effmass.Rd | plot.effmass | plot.effmass | plot.effmass | docs/api/plot_effmass.md | unassigned | plot.effmass |
| `[ ]` | man/plot.gevp.amplitude.Rd | plot.gevp.amplitude | plot.gevp.amplitude | plot.gevp.amplitude | docs/api/plot_gevp_amplitude.md | unassigned | plot.gevp.amplitude |
| `[ ]` | man/plot.hadronacf.Rd | plot.hadronacf | plot.hadronacf | plot.hadronacf | docs/api/plot_hadronacf.md | unassigned | plot.hadronacf |
| `[ ]` | man/plot.massfit.Rd | plot.massfit | plot.massfit | plot.massfit | docs/api/plot_massfit.md | unassigned | plot.massfit |
| `[ ]` | man/plot.matrixfit.Rd | plot.matrixfit | plot.matrixfit | plot.matrixfit | docs/api/plot_matrixfit.md | unassigned | Plot a matrixfit |
| `[ ]` | man/plot.ofit.Rd | plot.ofit | plot.ofit | plot.ofit | docs/api/plot_ofit.md | unassigned | plot.ofit |
| `[ ]` | man/plot.outputdata.Rd | plot.outputdata | plot.outputdata | plot.outputdata | docs/api/plot_outputdata.md | unassigned | Plot Command For Class Ouputdata |
| `[ ]` | man/plot.pionff.Rd | plot.pionff | plot.pionff | plot.pionff | docs/api/plot_pionff.md | unassigned | plot.pionff |
| `[ ]` | man/plot.raw_cf.Rd | plot.raw_cf | plot.raw_cf | plot.raw_cf | docs/api/plot_raw_cf.md | unassigned | plot all correlators in raw_cf object |
| `[ ]` | man/plot.truncated.pgevm.Rd | plot.truncated.pgevm | plot.truncated.pgevm | plot.truncated.pgevm | docs/api/plot_truncated_pgevm.md | unassigned | plot.truncated.pgevm |
| `[ ]` | man/plot.uwerr.Rd | plot.uwerr | plot.uwerr | plot.uwerr | docs/api/plot_uwerr.md | unassigned | Plot Command For Class UWerr |
| `[ ]` | man/plothlinewitherror.Rd | plothlinewitherror | plothlinewitherror | plothlinewitherror | docs/api/plothlinewitherror.md | unassigned | plothlinewitherror |
| `[ ]` | man/plotwitherror.Rd | plotwitherror | plotwitherror | plotwitherror | docs/api/plotwitherror.md | unassigned | Plot Command For XY Plots With Error Bars |
| `[ ]` | man/plus-.cf.Rd | +.cf | +.cf | +.cf | docs/api/add_cf.md | unassigned | Arithmetically add correlators |
| `[ ]` | man/plus-.raw_cf.Rd | +.raw_cf | +.raw_cf | +.raw_cf | docs/api/add_raw_cf.md | unassigned | add two raw_cf objects |
| `[ ]` | man/pointswithslantederror.Rd | pointswithslantederror | pointswithslantederror | pointswithslantederror | docs/api/pointswithslantederror.md | unassigned | pointswithslantederror |
| `[ ]` | man/predict.bootstrapfit.Rd | predict.bootstrapfit | predict.bootstrapfit | predict.bootstrapfit | docs/api/predict_bootstrapfit.md | unassigned | Predict values for bootstrapfit |
| `[ ]` | man/print.bootstrapfit.Rd | print.bootstrapfit | print.bootstrapfit | print.bootstrapfit | docs/api/print_bootstrapfit.md | unassigned | Print a bootstrap NLS fit |
| `[ ]` | man/print.cf.Rd | print.cf | print.cf | print.cf | docs/api/print_cf.md | unassigned | print.cf |
| `[ ]` | man/print.effectivemassfit.Rd | print.effectivemassfit | print.effectivemassfit | print.effectivemassfit | docs/api/print_effectivemassfit.md | unassigned | print.effectivemassfit |
| `[ ]` | man/print.ofit.Rd | print.ofit | print.ofit | print.ofit | docs/api/print_ofit.md | unassigned | print.ofit |
| `[ ]` | man/print.raw_cf.Rd | print.raw_cf | print.raw_cf | print.raw_cf | docs/api/print_raw_cf.md | unassigned | Print summary of data contained in raw_cf container |
| `[ ]` | man/pscor.sample.Rd | pscor.sample | pscor.sample | pscor.sample | docs/api/pscor_sample.md | unassigned | Sample pseudoscalar correlator |
| `[ ]` | man/raw_cf_data.Rd | raw_cf_data | raw_cf_data | raw_cf_data | docs/api/raw_cf_data.md | unassigned | Original data mixin constructor for raw_cf |
| `[ ]` | man/raw_cf_meta.Rd | raw_cf_meta | raw_cf_meta | raw_cf_meta | docs/api/raw_cf_meta.md | unassigned | raw_cf metadata mixin constructor |
| `[ ]` | man/raw_cf_to_cf.Rd | raw_cf_to_cf | raw_cf_to_cf | raw_cf_to_cf | docs/api/raw_cf_to_cf.md | unassigned | Extract a particular internal component of a 'raw_cf' into a 'cf' |
| `[ ]` | man/raw_cf.Rd | raw_cf | raw_cf | raw_cf | docs/api/raw_cf.md | unassigned | Container for raw correlation functions |
| `[ ]` | man/readbinarycf.Rd | readbinarycf | readbinarycf | readbinarycf | docs/api/readbinarycf.md | unassigned | read correlation function from binary files |
| `[ ]` | man/readbinarydisc.Rd | readbinarydisc | readbinarydisc | readbinarydisc | docs/api/readbinarydisc.md | unassigned | read disconnected loops from binary files |
| `[ ]` | man/readbinarysamples.Rd | readbinarysamples | readbinarysamples | readbinarysamples | docs/api/readbinarysamples.md | unassigned | Read binary correlation function by sample |
| `[ ]` | man/readcmidisc.Rd | readcmidisc | readcmidisc | readcmidisc | docs/api/readcmidisc.md | unassigned | reads disconnected loops in cmi format |
| `[ ]` | man/readcmifiles.Rd | readcmifiles | readcmifiles; readcmicor; readcmidatafiles; readcmiloopfiles | readcmifiles; readcmicor; readcmidatafiles; readcmiloopfiles | docs/api/readcmifiles.md | unassigned | Read Single Data Files in Chris Michael Format |
| `[ ]` | man/readgradflow.Rd | readgradflow | readgradflow | readgradflow | docs/api/readgradflow.md | unassigned | Read Gradient Flow Output Files in tmLQCD format |
| `[ ]` | man/readhlcor.Rd | readhlcor | readhlcor | readhlcor | docs/api/readhlcor.md | unassigned | readhlcor |
| `[ ]` | man/readnissatextcf.Rd | readnissatextcf | readnissatextcf | readnissatextcf | docs/api/readnissatextcf.md | unassigned | reader for Nissa text format correlation functions |
| `[ ]` | man/readoutputdata.Rd | readoutputdata | readoutputdata | readoutputdata | docs/api/readoutputdata.md | unassigned | Read Data In output.data Format of tmLQCD |
| `[ ]` | man/readtextcf.Rd | readtextcf | readtextcf | readtextcf | docs/api/readtextcf.md | unassigned | Read correlator data from single file |
| `[ ]` | man/removeTemporal.cf.Rd | removeTemporal.cf | removeTemporal.cf | removeTemporal.cf | docs/api/remove_temporal_cf.md | unassigned | Remove Thermal States by Weighting and Shifting |
| `[ ]` | man/resample_hankel.Rd | resample_hankel | resample_hankel | resample_hankel | docs/api/resample_hankel.md | unassigned | Resample bootstrap samples in Hankel effmass |
| `[ ]` | man/resampling_is_compatible.Rd | resampling_is_compatible | resampling_is_compatible | resampling_is_compatible | docs/api/resampling_is_compatible.md | unassigned | Checks whether the resampling of two cf objects is compatible |
| `[ ]` | man/resampling_is_concatenable.Rd | resampling_is_concatenable | resampling_is_concatenable | resampling_is_concatenable | docs/api/resampling_is_concatenable.md | unassigned | Checks whether the resampling of two cf objects is concatenable |
| `[ ]` | man/residual_plot.Rd | residual_plot | residual_plot | residual_plot | docs/api/residual_plot.md | unassigned | residual_plot |
| `[ ]` | man/restore_seed.Rd | restore_seed | restore_seed | restore_seed | docs/api/restore_seed.md | unassigned | Restore random number generator state |
| `[ ]` | man/samplecf.Rd | samplecf | samplecf | samplecf | docs/api/samplecf.md | unassigned | Sample cf data |
| `[ ]` | man/shift.cf.Rd | shift.cf | shift.cf | shift.cf | docs/api/shift_cf.md | unassigned | shift a correlation function by 'places' time-slices |
| `[ ]` | man/shift.raw_cf.Rd | shift.raw_cf | shift.raw_cf | shift.raw_cf | docs/api/shift_raw_cf.md | unassigned | shift a raw_cf correlation function by 'places' time-slices |
| `[ ]` | man/simple.nlsfit.Rd | simple.nlsfit | simple.nlsfit | simple.nlsfit | docs/api/simple_nlsfit.md | unassigned | NLS fit with without bootstrap |
| `[ ]` | man/slash-.cf.Rd | *.cf | *.cf; /.cf | *.cf; /.cf | docs/api/mul_cf.md | unassigned | Divide two cf objects by each other measurement by measurement |
| `[ ]` | man/slash-.raw_cf.Rd | /.raw_cf | /.raw_cf | /.raw_cf | docs/api/div_raw_cf.md | unassigned | divide two raw_cf objects |
| `[ ]` | man/store_correl.Rd | store_correl | store_correl | store_correl | docs/api/store_correl.md | unassigned | Store a 'raw_cf' correlator in an associative array together with a description<br>The object cf will be stored as an element of cmap under key out_key<br>in the member obj of cmap. The data frame passed via desc will be<br>appended as a row to cmap[[out_key]]$map. If out_key does not exist<br>as a key in cmap, a new element will be created. If it already exists,<br>addStat.raw_cf is called to add statistics to the existing raw_cf. Requires<br>the 'hash' package. |
| `[ ]` | man/string2error.Rd | string2error | string2error | string2error | docs/api/string2error.md | unassigned | string2error |
| `[ ]` | man/subtract.excitedstates.Rd | subtract.excitedstates | subtract.excitedstates | subtract.excitedstates | docs/api/subtract_excitedstates.md | unassigned | Substract excited states. |
| `[ ]` | man/sum_cosh.Rd | sum_cosh | sum_cosh | sum_cosh | docs/api/sum_cosh.md | unassigned | sum_i (a_i cosh(m_i*t))<br>a_i are the amplitudes, m_i the masses and t is a vector of times |
| `[ ]` | man/summary.bootstrapfit.Rd | summary.bootstrapfit | summary.bootstrapfit | summary.bootstrapfit | docs/api/summary_bootstrapfit.md | unassigned | Summarize a bootstrap NLS fit |
| `[ ]` | man/summary.cf.Rd | summary.cf | summary.cf | summary.cf | docs/api/summary_cf.md | unassigned | summary.cf |
| `[ ]` | man/summary.coshfit.Rd | summary.coshfit | summary.coshfit | summary.coshfit | docs/api/summary_coshfit.md | unassigned | Summarize a cosh-fit |
| `[ ]` | man/summary.effectivemass.Rd | summary.effectivemass | summary.effectivemass | summary.effectivemass | docs/api/summary_effectivemass.md | unassigned | summary.effectivemass |
| `[ ]` | man/summary.effectivemassfit.Rd | summary.effectivemassfit | summary.effectivemassfit | summary.effectivemassfit | docs/api/summary_effectivemassfit.md | unassigned | summary.effectivemassfit |
| `[ ]` | man/summary.gevp.amplitude.Rd | summary.gevp.amplitude | summary.gevp.amplitude | summary.gevp.amplitude | docs/api/summary_gevp_amplitude.md | unassigned | summary.gevp.amplitude |
| `[ ]` | man/summary.hadronacf.Rd | summary.hadronacf | summary.hadronacf | summary.hadronacf | docs/api/summary_hadronacf.md | unassigned | summary.hadronacf |
| `[ ]` | man/summary.hankel_summed.Rd | summary.hankel_summed | summary.hankel_summed | summary.hankel_summed | docs/api/summary_hankel_summed.md | unassigned | summary.hankel_summed |
| `[ ]` | man/summary.matrixfit.Rd | summary.matrixfit | summary.matrixfit | summary.matrixfit | docs/api/summary_matrixfit.md | unassigned | summary.matrixfit |
| `[ ]` | man/summary.ofit.Rd | summary.ofit | summary.ofit | summary.ofit | docs/api/summary_ofit.md | unassigned | summary.ofit |
| `[ ]` | man/summary.raw_cf.Rd | summary.raw_cf | summary.raw_cf | summary.raw_cf | docs/api/summary_raw_cf.md | unassigned | Print summary of data contained in raw_cf container |
| `[ ]` | man/summary.uwerr.Rd | summary.uwerr | summary.uwerr | summary.uwerr | docs/api/summary_uwerr.md | unassigned | summary.uwerr |
| `[ ]` | man/swap_seed.Rd | swap_seed | swap_seed | swap_seed | docs/api/swap_seed.md | unassigned | Set seed and store a seed which can be used to<br>reset the random number generator |
| `[ ]` | man/symmetrise.cf.Rd | symmetrise.cf | symmetrise.cf | symmetrise.cf | docs/api/symmetrise_cf.md | unassigned | Average backward and forward-dominated parts of the correlation function |
| `[ ]` | man/takeTimeDiff.cf.Rd | takeTimeDiff.cf | takeTimeDiff.cf | takeTimeDiff.cf | docs/api/take_time_diff_cf.md | unassigned | Take time difference |
| `[ ]` | man/tex.catwitherror.Rd | tex.catwitherror | tex.catwitherror | tex.catwitherror | docs/api/tex_catwitherror.md | unassigned | paste a number with error in tex-ready format |
| `[ ]` | man/tikz.finalize.Rd | tikz.finalize | tikz.finalize | tikz.finalize | docs/api/tikz_finalize.md | unassigned | tikz.finalize |
| `[ ]` | man/tikz.init.Rd | tikz.init | tikz.init | tikz.init | docs/api/tikz_init.md | unassigned | tikz.init |
| `[ ]` | man/times-.raw_cf.Rd | *.raw_cf | *.raw_cf | *.raw_cf | docs/api/mul_raw_cf.md | unassigned | multiply two raw_cf objects |
| `[ ]` | man/unsymmetrise.cf.Rd | unsymmetrise.cf | unsymmetrise.cf | unsymmetrise.cf | docs/api/unsymmetrise_cf.md | unassigned | Unfold a correlation function which has been symmetrised |
| `[ ]` | man/uwerr.cf.Rd | uwerr.cf | uwerr.cf | uwerr.cf | docs/api/uwerr_cf.md | unassigned | uwerr.cf |
| `[ ]` | man/uwerr.raw_cf.Rd | uwerr.raw_cf | uwerr.raw_cf | uwerr.raw_cf | docs/api/uwerr_raw_cf.md | unassigned | Gamma method analysis on all time-slices in a 'raw_cf' object |
| `[ ]` | man/uwerr.Rd | uwerr | uwerr; uwerrprimary; uwerrderived | uwerr; uwerrprimary; uwerrderived | docs/api/uwerr.md | unassigned | Time Series Analysis With Gamma Method |
| `[ ]` | man/weight_shift_reweight.cf.Rd | weight_shift_reweight.cf | weight_shift_reweight.cf | weight_shift_reweight.cf | docs/api/weight_shift_reweight_cf.md | unassigned | Weight-shift-reweight a correlation function |
| `[ ]` | man/weight.cf.Rd | weight.cf | weight.cf | weight.cf | docs/api/weight_cf.md | unassigned | Weight a correlation function |
| `[ ]` | man/zetazp.Rd | zetazp | zetazp | zetazp | docs/api/zetazp.md | unassigned | Computes the running of Z_P from scale mu0 to scale mu2 |

## Native C/C++ files

| Status | Source path | Type | Proposed Python target | Policy | Test requirement | Owner | Notes |
|---|---|---|---|---|---|---|---|
| `[ ]` | src/alpha_s.c | c | src/pydron/native/alpha_s.py or src/pydron/_native/alpha_s.* | rewrite first, wrap only if parity/performance requires | R/native reference-output tests before `[E]` | unassigned | small numerical kernels with direct R wrappers |
| `[ ]` | src/alpha_s.h | h | src/pydron/native/alpha_s.py or src/pydron/_native/alpha_s.* | wrap/defer | R/native reference-output tests before `[E]` | unassigned | header for C implementation; keep with owning C file decision |
| `[?]` | src/cdh.c | c | src/pydron/native/finite_size.py or compiled extension | needs decision: wrap or rewrite | R/native reference-output tests before `[E]` | unassigned | GSL-heavy finite-size correction implementation |
| `[ ]` | src/cdh.h | h | src/pydron/native/finite_size.py or compiled extension | wrap/defer | R/native reference-output tests before `[E]` | unassigned | header for C implementation; keep with owning C file decision |
| `[ ]` | src/inv_cosh.c | c | src/pydron/native/inv_cosh.py or scipy-backed helper | rewrite first, wrap only if parity/performance requires | R/native reference-output tests before `[E]` | unassigned | small numerical kernels with direct R wrappers |
| `[ ]` | src/RcppExports.cpp | cpp | build-generated binding file only if wrapping | defer/regenerate | R/native reference-output tests before `[E]` | unassigned | generated Rcpp registration; regenerate only if native wrapping is chosen |
| `[ ]` | src/read_nissa_textcf_kernel.cpp | cpp | src/pydron/io/nissa.py or compiled parser | rewrite or wrap after I/O fixture review | R/native reference-output tests before `[E]` | unassigned | Rcpp text parser for NISSA correlator files |
| `[ ]` | src/tmcdh.c | c | src/pydron/native/finite_size.py or compiled extension | defer | R/native reference-output tests before `[E]` | unassigned | not registered in current R-to-native interface |

## Registered native callable routines

| Status | Source path | Registered name | C/C++ function | Args | Proposed Python target | Policy | Test requirement | Owner | Notes |
|---|---|---|---|---|---|---|---|---|---|
| `[?]` | src/RcppExports.cpp | _hadron_read_nissa_textcf_kernel | _hadron_read_nissa_textcf_kernel | 4 | build-generated binding file only if wrapping | wrap/rewrite/defer decision required | R-to-native parity fixture required | unassigned | registered in `src/RcppExports.cpp` |
| `[?]` | src/alpha_s.c | alphas | alphas | 5 | src/pydron/native/alpha_s.py or src/pydron/_native/alpha_s.* | wrap/rewrite/defer decision required | R-to-native parity fixture required | unassigned | registered in `src/RcppExports.cpp` |
| `[?]` | src/cdh.c | cdh_c | cdh_c | 13 | src/pydron/native/finite_size.py or compiled extension | wrap/rewrite/defer decision required | R-to-native parity fixture required | unassigned | registered in `src/RcppExports.cpp` |
| `[?]` | src/cdh.c | cdhnew_c | cdhnew_c | 11 | src/pydron/native/finite_size.py or compiled extension | wrap/rewrite/defer decision required | R-to-native parity fixture required | unassigned | registered in `src/RcppExports.cpp` |
| `[?]` | src/inv_cosh.c | invcosh | invcosh | 5 | src/pydron/native/inv_cosh.py or scipy-backed helper | wrap/rewrite/defer decision required | R-to-native parity fixture required | unassigned | registered in `src/RcppExports.cpp` |

## R-to-native interfaces

| R source path | Line | Interface | Native routine | Proposed port action |
|---|---|---|---|---|
| R/alpha_s.R | 27 | .Call | alphas | tie wrapper parity to registered native routine row |
| R/cdh.R | 61 | .Call | cdh_c | tie wrapper parity to registered native routine row |
| R/cdh.R | 132 | .Call | cdhnew_c | tie wrapper parity to registered native routine row |
| R/inv_cosh.R | 26 | .Call | invcosh | tie wrapper parity to registered native routine row |
| R/RcppExports.R | 5 | .Call | _hadron_read_nissa_textcf_kernel | tie wrapper parity to registered native routine row |

## Datasets

| Status | Source path | Dataset file name | R object | Class | Dimensions | Length | Approx. size bytes | Proposed Python handling | Test requirement | Owner |
|---|---|---|---|---|---|---|---|---|---|---|
| `[?]` | data/cA2.09.48_3pi_I3_0_A1u_1_pc.RData | cA2.09.48_3pi_I3_0_A1u_1_pc | cA2.09.48_3pi_I3_0_A1u_1_pc | list;cf;cf_meta;cf_boot;cf_principal_correlator | - | 19 | 794904 | load directly or convert to neutral fixture; decision pending | loader test plus checksum/object-shape parity | unassigned |
| `[?]` | data/correlatormatrix.RData | correlatormatrix | correlatormatrix | list;cf;cf_meta;cf_orig | - | 6 | 435160 | load directly or convert to neutral fixture; decision pending | loader test plus checksum/object-shape parity | unassigned |
| `[?]` | data/loopdata.RData | loopdata | loopdata | list;cf;cf_meta;cf_orig;cf_smeared | - | 12 | 64568 | load directly or convert to neutral fixture; decision pending | loader test plus checksum/object-shape parity | unassigned |
| `[?]` | data/plaq.sample.RData | plaq.sample | plaq.sample | numeric | - | 6352 | 50864 | load directly or convert to neutral fixture; decision pending | loader test plus checksum/object-shape parity | unassigned |
| `[?]` | data/pscor.sample.RData | pscor.sample | pscor.sample | data.frame | 15168x2 | 2 | 1153608 | load directly or convert to neutral fixture; decision pending | loader test plus checksum/object-shape parity | unassigned |
| `[?]` | data/samplecf.RData | samplecf | samplecf | list;cf;cf_meta;cf_orig | - | 6 | 205360 | load directly or convert to neutral fixture; decision pending | loader test plus checksum/object-shape parity | unassigned |

## Data metadata files

| Status | Source path | Role | Proposed Python handling | Owner | Notes |
|---|---|---|---|---|---|
| `[ ]` | data/datalist | R LazyData dataset listing | use to cross-check packaged dataset names | unassigned | lists 5 names; `cA2.09.48_3pi_I3_0_A1u_1_pc.RData` is present but not listed |

## Tests

| Status | Source path | R parse status | Test contexts/names | Detected tested R symbols | Proposed Python target | Owner | Notes |
|---|---|---|---|---|---|---|---|
| `[ ]` | tests/new_matrixfit_all_points.Rmd | ok | root documentation/workflow file | bootstrap.cf; matrixfit; new_matrixfit; residual_plot | tests/reference_workflows/test_new_matrixfit_all_points.py | unassigned | convert assertions and/or generate R reference fixtures |
| `[ ]` | tests/single_constant_model.Rmd | ok | root documentation/workflow file | bootstrap.cf; bootstrap.effectivemass; new_matrixfit; residual_plot | tests/reference_workflows/test_single_constant_model.py | unassigned | convert assertions and/or generate R reference fixtures |
| `[ ]` | tests/testthat.R | ok | root documentation/workflow file | - | tests/test_testthat.py | unassigned | convert assertions and/or generate R reference fixtures |
| `[ ]` | tests/testthat/test_bootstrapfit.R | ok | y errors; y errors cov; xy errors; xy errors cov; y errors with priors; y errors cov with priors; xy errors with priors; xy errors cov with priors | parametric.bootstrap; bootstrap.nlsfit | tests/test_bootstrapfit.py | unassigned | convert assertions and/or generate R reference fixtures |
| `[ ]` | tests/testthat/test_computeDisc.R | ok | cross_vs_diagonal | computeDisc | tests/test_compute_disc.py | unassigned | convert assertions and/or generate R reference fixtures |
| `[ ]` | tests/testthat/test_dummy.R | ok | root documentation/workflow file | - | tests/test_dummy.py | unassigned | convert assertions and/or generate R reference fixtures |
| `[ ]` | tests/testthat/test_extractSingleCor_cf.R | ok | extract_symmetrized | bootstrap.cf; mul.cf; extractSingleCor.cf | tests/test_extract_single_cor_cf.py | unassigned | convert assertions and/or generate R reference fixtures |
| `[ ]` | tests/testthat/test_new_matrixfit_higher_states.R | ok | two states with cov; two states | new_matrixfit | tests/test_new_matrixfit_higher_states.py | unassigned | convert assertions and/or generate R reference fixtures |
| `[ ]` | tests/testthat/test_new_matrixfit.R | ok | SingleModelJacobian; SingleModel; ShiftedModelPrediction; ShiftedModel; TwoStateModel | make_parlist; dmatrixChi; make_parind; make_sign_vec; make_ov_sign_vec; bootstrap.cf; matrixChi.shifted; dmatrixChi.shifted; takeTimeDiff.cf; bootstrap.gevp; gevp2cf | tests/test_new_matrixfit.py | unassigned | convert assertions and/or generate R reference fixtures |
| `[ ]` | tests/testthat/test_parlist.R | ok | parlist_1; parlist_4; parind_1; parind_4 | make_parlist; make_parind | tests/test_parlist.py | unassigned | convert assertions and/or generate R reference fixtures |
| `[ ]` | tests/testthat/test_removeTemporal.R | ok | equality | bootstrap.cf; matrixfit; old_removeTemporal.cf; removeTemporal.cf | tests/test_remove_temporal.py | unassigned | convert assertions and/or generate R reference fixtures |
| `[ ]` | tests/testthat/test_string2error.R | ok | 1 | string2error | tests/test_string2error.py | unassigned | convert assertions and/or generate R reference fixtures |
| `[ ]` | tests/testthat/test_tex_catwitherror.R | ok | small_error; borderline_error; another_borderline_error; even_nastier_borderline_error; scientific_notation; very_small_number; same_error; intermediate_error; large_error; similar_error; no_error; zero_error; zero_val_zero_err; NA; vector | tex.catwitherror | tests/test_tex_catwitherror.py | unassigned | convert assertions and/or generate R reference fixtures |

## Vignettes and workflows

| Status | Source path | Workflow topic | Proposed Python target | Owner | Notes |
|---|---|---|---|---|---|
| `[ ]` | vignettes/Two_Amplitudes_Model.Rmd | R vignette: Two Amplitudes Model | docs/examples/two_amplitudes_model.md or notebooks/two_amplitudes_model.ipynb | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | vignettes/gevp.Rmd | R vignette: gevp | docs/examples/gevp.md or notebooks/gevp.ipynb | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | vignettes/hankel.Rmd | R vignette: hankel | docs/examples/hankel.md or notebooks/hankel.ipynb | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | vignettes/hankel.bib | bibliography for Hankel/PGEVM docs | docs/examples/hankel.md or notebooks/hankel.ipynb | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | vignettes/jackknife_cov_and_missing_values.Rmd | R vignette: jackknife cov and missing values | docs/examples/jackknife_cov_and_missing_values.md or notebooks/jackknife_cov_and_missing_values.ipynb | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | vignettes/jackknife_error_normalization.Rmd | R vignette: jackknife error normalization | docs/examples/jackknife_error_normalization.md or notebooks/jackknife_error_normalization.ipynb | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | vignettes/multi_particle_fit_test.Rmd | R vignette: multi particle fit test | docs/examples/multi_particle_fit_test.md or notebooks/multi_particle_fit_test.ipynb | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | vignettes/pgevm.Rmd | R vignette: pgevm | docs/examples/pgevm.md or notebooks/pgevm.ipynb | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | vignettes/truncated_pgevm.Rmd | R vignette: truncated pgevm | docs/examples/truncated_pgevm.md or notebooks/truncated_pgevm.ipynb | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/3pt.R | executable/example analysis script | examples/_3pt.py or docs/migration/_3pt.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/analyse.eta_OS.4x4.R | executable/example analysis script | examples/analyse_eta_os_4x4.py or docs/migration/analyse_eta_os_4x4.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/analyse.eta_ss.R | executable/example analysis script | examples/analyse_eta_ss.py or docs/migration/analyse_eta_ss.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/analyse.nd-kaon.8x8.R | executable/example analysis script | examples/analyse_nd_kaon_8x8.py or docs/migration/analyse_nd_kaon_8x8.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/analyse.pion.2x2.R | executable/example analysis script | examples/analyse_pion_2x2.py or docs/migration/analyse_pion_2x2.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/analyse.pion0.2x2.R | executable/example analysis script | examples/analyse_pion0_2x2.py or docs/migration/analyse_pion0_2x2.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/analyse_pipi.R | executable/example analysis script | examples/analyse_pipi.py or docs/migration/analyse_pipi.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/averx/perform-fits.R | three-point weighted-average workflow | examples/averx_perform_fits.py or docs/migration/averx_perform_fits.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/averx/perform-weighted-median.R | three-point weighted-average workflow | examples/averx_perform_weighted_median.py or docs/migration/averx_perform_weighted_median.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/eta_extrapolate.R | executable/example analysis script | examples/eta_extrapolate.py or docs/migration/eta_extrapolate.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/feta.R | executable/example analysis script | examples/feta.py or docs/migration/feta.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/hdf5-example.R | executable/example analysis script | examples/hdf5_example.py or docs/migration/hdf5_example.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/mesons-cmi/a0.R | meson CMI analysis script or docs | examples/mesons_cmi_a0.py or docs/migration/mesons_cmi_a0.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/mesons-cmi/b1.R | meson CMI analysis script or docs | examples/mesons_cmi_b1.py or docs/migration/mesons_cmi_b1.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/mesons-cmi/man/b1.Rd | meson CMI analysis script or docs | examples/mesons_cmi_man_b1.py or docs/migration/mesons_cmi_man_b1.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/mesons-cmi/man/plot.cfit.Rd | meson CMI analysis script or docs | examples/mesons_cmi_man_plot_cfit.py or docs/migration/mesons_cmi_man_plot_cfit.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/mesons-cmi/man/rho.Rd | meson CMI analysis script or docs | examples/mesons_cmi_man_rho.py or docs/migration/mesons_cmi_man_rho.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/mesons-cmi/pion0.R | meson CMI analysis script or docs | examples/mesons_cmi_pion0.py or docs/migration/mesons_cmi_pion0.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/mesons-cmi/pionChiPTfit.R | meson CMI analysis script or docs | examples/mesons_cmi_pion_chi_ptfit.py or docs/migration/mesons_cmi_pion_chi_ptfit.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/mesons-cmi/plot.pionChiPTfit.R | meson CMI analysis script or docs | examples/mesons_cmi_plot_pion_chi_ptfit.py or docs/migration/mesons_cmi_plot_pion_chi_ptfit.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/mesons-cmi/summary.a0fit.R | meson CMI analysis script or docs | examples/mesons_cmi_summary_a0fit.py or docs/migration/mesons_cmi_summary_a0fit.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/mesons-cmi/summary.b1fit.R | meson CMI analysis script or docs | examples/mesons_cmi_summary_b1fit.py or docs/migration/mesons_cmi_summary_b1fit.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/mesons-cmi/summary.chiralfit.R | meson CMI analysis script or docs | examples/mesons_cmi_summary_chiralfit.py or docs/migration/mesons_cmi_summary_chiralfit.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/mesons-cmi/summary.pionChiPTfit.R | meson CMI analysis script or docs | examples/mesons_cmi_summary_pion_chi_ptfit.py or docs/migration/mesons_cmi_summary_pion_chi_ptfit.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/mesons-cmi/summary.rhofit.R | meson CMI analysis script or docs | examples/mesons_cmi_summary_rhofit.py or docs/migration/mesons_cmi_summary_rhofit.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/mesons-cmi/summarycfit.R | meson CMI analysis script or docs | examples/mesons_cmi_summarycfit.py or docs/migration/mesons_cmi_summarycfit.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/mesons-cmi/vector.R | meson CMI analysis script or docs | examples/mesons_cmi_vector.py or docs/migration/mesons_cmi_vector.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/old/OSChiPTFit.R | legacy executable analysis script | examples/old_oschi_ptfit.py or docs/migration/old_oschi_ptfit.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/old/ana.R | legacy executable analysis script | examples/old_ana.py or docs/migration/old_ana.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/old/analyseOS.R | legacy executable analysis script | examples/old_analyse_os.py or docs/migration/old_analyse_os.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/old/avercycle.R | legacy executable analysis script | examples/old_avercycle.py or docs/migration/old_avercycle.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/old/cfunction.R | legacy executable analysis script | examples/old_cfunction.py or docs/migration/old_cfunction.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/old/chiralfit.R | legacy executable analysis script | examples/old_chiralfit.py or docs/migration/old_chiralfit.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/old/chiralfit2.R | legacy executable analysis script | examples/old_chiralfit2.py or docs/migration/old_chiralfit2.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/old/cmfit.R | legacy executable analysis script | examples/old_cmfit.py or docs/migration/old_cmfit.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/old/cont_extr.R | legacy executable analysis script | examples/old_cont_extr.py or docs/migration/old_cont_extr.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/old/fit.R | legacy executable analysis script | examples/old_fit.py or docs/migration/old_fit.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/old/fit_fp.R | legacy executable analysis script | examples/old_fit_fp.py or docs/migration/old_fit_fp.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/old/gsl_fit.R | legacy executable analysis script | examples/old_gsl_fit.py or docs/migration/old_gsl_fit.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/old/kaon.R | legacy executable analysis script | examples/old_kaon.py or docs/migration/old_kaon.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/old/kaon.Rd | legacy executable analysis script | examples/old_kaon.py or docs/migration/old_kaon.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/old/man/avercycle.Rd | legacy executable analysis script | examples/old_man_avercycle.py or docs/migration/old_man_avercycle.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/old/man/cfunction.Rd | legacy executable analysis script | examples/old_man_cfunction.py or docs/migration/old_man_cfunction.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/old/nucleon.R | legacy executable analysis script | examples/old_nucleon.py or docs/migration/old_nucleon.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/old/oldplotutils.R | legacy executable analysis script | examples/old_oldplotutils.py or docs/migration/old_oldplotutils.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/old/pion.R | legacy executable analysis script | examples/old_pion.py or docs/migration/old_pion.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/old/pp2.R | legacy executable analysis script | examples/old_pp2.py or docs/migration/old_pp2.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/old/r0.R | legacy executable analysis script | examples/old_r0.py or docs/migration/old_r0.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/old/smearedpion.R | legacy executable analysis script | examples/old_smearedpion.py or docs/migration/old_smearedpion.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/old/strange.R | legacy executable analysis script | examples/old_strange.py or docs/migration/old_strange.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/old/summary.pionfit.R | legacy executable analysis script | examples/old_summary_pionfit.py or docs/migration/old_summary_pionfit.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/old/variational.R | legacy executable analysis script | examples/old_variational.py or docs/migration/old_variational.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/online_measurements_analysis_driver_template.R | executable/example analysis script | examples/online_measurements_analysis_driver_template.py or docs/migration/online_measurements_analysis_driver_template.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/online_measurements_status_template.Rmd | executable/example analysis script | examples/online_measurements_status_template.py or docs/migration/online_measurements_status_template.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/phaseshift/analyse.R | phase-shift workflow script | examples/phaseshift_analyse.py or docs/migration/phaseshift_analyse.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/phaseshift/finish.R | phase-shift workflow script | examples/phaseshift_finish.py or docs/migration/phaseshift_finish.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/phaseshift/fit-detEq-energylevels.R | phase-shift workflow script | examples/phaseshift_fit_det_eq_energylevels.py or docs/migration/phaseshift_fit_det_eq_energylevels.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/phaseshift/get-phaseshifts.R | phase-shift workflow script | examples/phaseshift_get_phaseshifts.py or docs/migration/phaseshift_get_phaseshifts.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/phaseshift/parameters.R | phase-shift workflow script | examples/phaseshift_parameters.py or docs/migration/phaseshift_parameters.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/phaseshift/phaseshift.R | phase-shift workflow script | examples/phaseshift_phaseshift.py or docs/migration/phaseshift_phaseshift.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/phaseshift/phaseshift.pipiswave.R | phase-shift workflow script | examples/phaseshift_phaseshift_pipiswave.py or docs/migration/phaseshift_phaseshift_pipiswave.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/phaseshift/pipi.R | phase-shift workflow script | examples/phaseshift_pipi.py or docs/migration/phaseshift_pipi.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/phaseshift/plot-Energies.R | phase-shift workflow script | examples/phaseshift_plot_energies.py or docs/migration/phaseshift_plot_energies.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/phaseshift/plot-deltaell-vs-qsq.R | phase-shift workflow script | examples/phaseshift_plot_deltaell_vs_qsq.py or docs/migration/phaseshift_plot_deltaell_vs_qsq.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/phaseshift/singlepi.R | phase-shift workflow script | examples/phaseshift_singlepi.py or docs/migration/phaseshift_singlepi.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/phaseshift/summary.R | phase-shift workflow script | examples/phaseshift_summary.py or docs/migration/phaseshift_summary.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/phaseshift/test-rootfinding.R | phase-shift workflow script | examples/phaseshift_test_rootfinding.py or docs/migration/phaseshift_test_rootfinding.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/putonlinetogether.sh | executable/example analysis script | examples/putonlinetogether.py or docs/migration/putonlinetogether.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/puttogether.sh | executable/example analysis script | examples/puttogether.py or docs/migration/puttogether.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/puttogether_reverse.sh | executable/example analysis script | examples/puttogether_reverse.py or docs/migration/puttogether_reverse.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/rho-phaseshift/analyse.R | phase-shift workflow script | examples/rho_phaseshift_analyse.py or docs/migration/rho_phaseshift_analyse.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/rho-phaseshift/average.data.R | phase-shift workflow script | examples/rho_phaseshift_average_data.py or docs/migration/rho_phaseshift_average_data.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/rho-phaseshift/detect_irrep_frame.R | phase-shift workflow script | examples/rho_phaseshift_detect_irrep_frame.py or docs/migration/rho_phaseshift_detect_irrep_frame.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/rho-phaseshift/fit.delta.R | phase-shift workflow script | examples/rho_phaseshift_fit_delta.py or docs/migration/rho_phaseshift_fit_delta.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/rho-phaseshift/infile-analyse.R | phase-shift workflow script | examples/rho_phaseshift_infile_analyse.py or docs/migration/rho_phaseshift_infile_analyse.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/rho-phaseshift/infile-fit.delta.R | phase-shift workflow script | examples/rho_phaseshift_infile_fit_delta.py or docs/migration/rho_phaseshift_infile_fit_delta.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/rho-phaseshift/phaseshift.rho.R | phase-shift workflow script | examples/rho_phaseshift_phaseshift_rho.py or docs/migration/rho_phaseshift_phaseshift_rho.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/rho-phaseshift/plot-Mrho.R | phase-shift workflow script | examples/rho_phaseshift_plot_mrho.py or docs/migration/rho_phaseshift_plot_mrho.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/rho-phaseshift/plot-delta.R | phase-shift workflow script | examples/rho_phaseshift_plot_delta.py or docs/migration/rho_phaseshift_plot_delta.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/rho-phaseshift/preliminary.fit.delta.R | phase-shift workflow script | examples/rho_phaseshift_preliminary_fit_delta.py or docs/migration/rho_phaseshift_preliminary_fit_delta.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/rho-phaseshift/summarise.R | phase-shift workflow script | examples/rho_phaseshift_summarise.py or docs/migration/rho_phaseshift_summarise.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/scatteringlength/analyse_pipi.R | scattering-length workflow script | examples/scatteringlength_analyse_pipi.py or docs/migration/scatteringlength_analyse_pipi.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/scatteringlength/analysis.R | scattering-length workflow script | examples/scatteringlength_analysis.py or docs/migration/scatteringlength_analysis.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/scatteringlength/chipt-fit.R | scattering-length workflow script | examples/scatteringlength_chipt_fit.py or docs/migration/scatteringlength_chipt_fit.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/scatteringlength/deltaE_from_ratios.R | scattering-length workflow script | examples/scatteringlength_delta_e_from_ratios.py or docs/migration/scatteringlength_delta_e_from_ratios.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/scatteringlength/fit_finite_range.R | scattering-length workflow script | examples/scatteringlength_fit_finite_range.py or docs/migration/scatteringlength_fit_finite_range.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/scatteringlength/fscor.R | scattering-length workflow script | examples/scatteringlength_fscor.py or docs/migration/scatteringlength_fscor.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/scatteringlength/gather-deltaE-values.R | scattering-length workflow script | examples/scatteringlength_gather_delta_e_values.py or docs/migration/scatteringlength_gather_delta_e_values.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/scatteringlength/get_finite_range_fits-ratio.R | scattering-length workflow script | examples/scatteringlength_get_finite_range_fits_ratio.py or docs/migration/scatteringlength_get_finite_range_fits_ratio.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/scatteringlength/get_summary-efm.R | scattering-length workflow script | examples/scatteringlength_get_summary_efm.py or docs/migration/scatteringlength_get_summary_efm.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/scatteringlength/get_summary-ratio.R | scattering-length workflow script | examples/scatteringlength_get_summary_ratio.py or docs/migration/scatteringlength_get_summary_ratio.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/scatteringlength/get_summary.R | scattering-length workflow script | examples/scatteringlength_get_summary.py or docs/migration/scatteringlength_get_summary.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/scatteringlength/plot-deltaE.R | scattering-length workflow script | examples/scatteringlength_plot_delta_e.py or docs/migration/scatteringlength_plot_delta_e.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/scatteringlength/plot-mpia0.R | scattering-length workflow script | examples/scatteringlength_plot_mpia0.py or docs/migration/scatteringlength_plot_mpia0.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/scatteringlength/plot-ratios.R | scattering-length workflow script | examples/scatteringlength_plot_ratios.py or docs/migration/scatteringlength_plot_ratios.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/scatteringlength/ratio.R | scattering-length workflow script | examples/scatteringlength_ratio.py or docs/migration/scatteringlength_ratio.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/scatteringlength/summary.R | scattering-length workflow script | examples/scatteringlength_summary.py or docs/migration/scatteringlength_summary.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/scatteringlength/test-data.R | scattering-length workflow script | examples/scatteringlength_test_data.py or docs/migration/scatteringlength_test_data.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | exec/scatteringlength/test-ratio.R | scattering-length workflow script | examples/scatteringlength_test_ratio.py or docs/migration/scatteringlength_test_ratio.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | inst/new_matrixfit.Rmd | installed R Markdown resource | docs/reference/new_matrixfit.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | inst/weighted_model.Rmd | installed R Markdown resource | docs/reference/weighted_model.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | notes/Matrixfit_Performance.Rmd | design/performance note | docs/design/matrixfit_performance.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | notes/gevp_review.md | design/performance note | docs/design/gevp_review.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | README.md | root documentation/workflow file | docs/migration/readme.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | NEWS.md | root documentation/workflow file | docs/migration/news.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | CONTRIBUTING.md | root documentation/workflow file | docs/migration/contributing.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | CODEX_PROMPT.md | root documentation/workflow file | docs/migration/codex_prompt.md | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | Weighted_Model.nb | Mathematica derivation notebook for weighted model | docs/design/weighted_model.md after manual review | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | check | repository workflow script | developer tooling decision; no algorithm port | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | test | repository workflow script | developer tooling decision; no algorithm port | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | document | repository workflow script | developer tooling decision; no algorithm port | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | install | repository workflow script | developer tooling decision; no algorithm port | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | cleanup | repository workflow script | developer tooling decision; no algorithm port | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | verify-exports | repository workflow script | developer tooling decision; no algorithm port | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |
| `[ ]` | configure.ac | native build configuration source | pyproject/scikit-build configuration if native wrapping is chosen | unassigned | convert to example/doc, reference fixture, or `[X]` with reason |

## Files still not parsed reliably

| Source path | Reason | Required follow-up | Owner |
|---|---|---|---|
| Weighted_Model.nb | Mathematica notebook inspected as text only, not as a Wolfram expression tree | Manually review or convert before using as documentation source | unassigned |
| inst/extdata/C2_bin.dat | Binary external data file; format not decoded in planning pass | Pair with reader fixture and reference checksum | unassigned |
| inst/extdata/C2_pi0.dat | Binary external data file; format not decoded in planning pass | Pair with reader fixture and reference checksum | unassigned |
| src/alpha_s.o | Compiled native artifact; not source-parsed | Do not port; regenerate via build if needed | unassigned |
| src/cdh.o | Compiled native artifact; not source-parsed | Do not port; regenerate via build if needed | unassigned |
| src/hadron.so | Compiled native artifact; not source-parsed | Do not port; regenerate via build if needed | unassigned |
| src/inv_cosh.o | Compiled native artifact; not source-parsed | Do not port; regenerate via build if needed | unassigned |
| src/RcppExports.o | Compiled native artifact; not source-parsed | Do not port; regenerate via build if needed | unassigned |
| src/read_nissa_textcf_kernel.o | Compiled native artifact; not source-parsed | Do not port; regenerate via build if needed | unassigned |
| src/tmcdh.o | Compiled native artifact; not source-parsed | Do not port; regenerate via build if needed | unassigned |

## First-pass issue breakdown

- [ ] **Planning verification:** Keep the R-backed inventory reproducible; if source changes, rerun the parser and update counts.
- [ ] **Package skeleton:** Add `pyproject.toml`, `src/pydron/`, import smoke tests, and no-op compatibility namespace.
- [ ] **R reference oracle:** Add deterministic fixture generation for high-value public functions now that Rscript is available.
- [ ] **Dataset strategy:** Choose direct `.RData` loading, conversion, or dual path; define checksums and object-shape tests from the metadata recorded here.
- [ ] **Native strategy:** Decide rewrite/wrap/defer for `alpha_s`, `invcosh`, `cdh`, and NISSA reader routines.
- [ ] **Core correlator model:** Port `raw_cf` and `cf` containers first, including arithmetic and summary compatibility tests.
- [ ] **Resampling/statistics:** Port bootstrap, jackknife, covariance, `UWerr`, and autocorrelation primitives behind shared sample classes.
- [ ] **Fit machinery:** Port nonlinear fit results and matrix-fit/R6 model hierarchy with R reference predictions.
- [ ] **Spectroscopy methods:** Port GEVP, Hankel/PGEVM, truncated PGEVM, and Lanczos workflows after core/fits are stable.
- [ ] **Documentation migration:** Convert vignettes, `inst/*.Rmd`, notes, root notebooks, and `man/*.Rd` examples into Python docs.
