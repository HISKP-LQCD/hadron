# pydron porting plan

This planning document fixes the inspected baseline for the Python port of the R package `hadron` to `pydron`. It is planning-only: no Python algorithms, package skeleton, or numerical behavior changes are part of this PR.

## Project goal

Port `hadron` to an object-oriented Python package while preserving the original scientific behavior and logical module structure. Public Python usage should be idiomatic, but every public R function, S3 method, R6 class, dataset, native routine, test, and workflow must remain traceable through `PORTING_INDEX.md` and the folder summaries.

## Non-goals

- No Python algorithm implementations in this planning PR.
- No Python package skeleton in this PR.
- No changes to numerical behavior, defaults, conventions, or reference data.
- No deletion, move, or rewrite of existing R/C/C++ source files.
- No broad redesign of the physics, statistics, fitting, plotting, or I/O workflow.
- No collapsing of unrelated R modules unless a later PR documents the exception and preserves traceability.
- No acceptance of a ported numerical item without R reference-output equivalence coverage.

## Baseline

| Item | Value | Notes |
|---|---|---|
| Inspection date | 2026-05-13 | Local repository inspection |
| Baseline commit | `42f8841c88e8ecb94a01a4abc3efc7a0a569fdc1` | `git rev-parse HEAD` |
| Branch at inspection | `port-plan` | Planning branch |
| R/Rscript version | `R version 4.3.3 (2024-02-29)` | `Rscript --version` |
| R package name/version | `hadron` / `3.4.1` | From `DESCRIPTION` |
| Compilation status | `NeedsCompilation: yes` | Native C/C++ and GSL are present |
| System requirements | Gnu Scientific Library version >= 1.8 | From `DESCRIPTION` |
| License | `GPL-3` | README says GPL 3 or later |

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

## R-backed inspection method

- Confirmed R with `Rscript --version` and recorded the exact version above.
- Parsed every `R/*.R` file with R `parse()`. The symbol inventory now tracks top-level function assignments and R6 declarations from the R parser rather than text-regex matches.
- Cross-checked `NAMESPACE` exports and `S3method()` registrations against parsed top-level definitions; all matched.
- Parsed `man/*.Rd` with `tools::parse_Rd()` and mapped names/aliases to source symbols where possible.
- Loaded every `data/*.RData` file in a clean R environment and recorded object names, classes, dimensions, lengths, and approximate object sizes without dumping object contents.
- Parsed R tests directly; R Markdown test files were inspected by parsing their R code chunks.
- Rechecked `.Call`, `.C`, `.Fortran`, Rcpp exports, and `src/RcppExports.cpp` native registration.

## Porting principles

1. Preserve the R package layout as the reference until parity is proven.
2. Use `R/foo.R` -> `src/pydron/foo.py` as the default mapping rule. If an R filename contains dots or other non-importable characters, normalize to snake_case and record the exception in the index.
3. Build Python internals around explicit objects for data, samples, fits, and analysis results; avoid unstructured dictionaries as the long-term model.
4. Keep the Python public API idiomatic first, with R-style compatibility wrappers where they improve traceability or migration.
5. Treat R reference-output equivalence as the gate for numerical behavior.
6. Keep native-code decisions explicit: rewrite, wrap, defer, or intentionally not port.
7. Prefer small, reviewable PRs with a fixed owner and status label per indexed item.

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

## Proposed Python package skeleton

```text
src/pydron/
├── __init__.py
├── alpha_s.py
├── cdh.py
├── cf.py
├── raw_cf.py
├── bootstrap_nlsfit.py
├── boot_ts_array.py
├── jackknife.py
├── timeseries.py
├── autocorrelation.py
├── covariance.py
├── effective_mass.py
├── fits.py
├── matrixfit.py
├── new_matrixfit.py
├── gevp.py
├── hankel.py
├── hankel_truncated.py
├── lanczos.py
├── luscher_method.py
├── gamma.py
├── momentum_utils.py
├── looptools.py
├── pcac.py
├── datasets.py
├── io/
│   ├── __init__.py
│   ├── binary.py
│   ├── cmi.py
│   ├── cvc.py
│   ├── cyprus.py
│   ├── gradient_flow.py
│   ├── hdf5.py
│   └── nissa.py
├── native/
│   ├── __init__.py
│   ├── alpha_s.py
│   ├── finite_size.py
│   └── inv_cosh.py
├── plotting/
│   ├── __init__.py
│   └── utils.py
└── compatibility/
    ├── __init__.py
    └── r_names.py
```

The skeleton is a target proposal only. This PR does not create it. `PORTING_INDEX.md` contains one row per R source file and parsed symbol with proposed normalized Python paths.

## Object-oriented design

| Candidate Python class | Primary R sources | Role |
|---|---|---|
| CorrelationFunction | `R/cf.R` | Processed correlator container, arithmetic, resampling, summaries |
| RawCorrelationFunction | `R/raw_cf.R` | Raw correlator data, metadata, conversion to `CorrelationFunction` |
| BootstrapSample | `R/cf.R`, `R/bootstrap.nlsfit.R`, `R/boot_ts_array.R` | Bootstrap samples, block bootstrap policy, error functions |
| JackknifeSample | `R/cf.R`, `R/jackknifeafterboot.R` | Jackknife samples and after-bootstrap jackknife diagnostics |
| TimeSeries | `R/timeseries.R`, `R/computeacf.R`, `R/UWerr.R` | Time-series storage, autocorrelation, UWerr integration |
| FitResult | `R/bootstrap.nlsfit.R`, `R/cosh_nlsfit.R`, `R/fit.plateau2cf.R` | Generic fit result, predictions, residuals, summaries |
| MatrixFitResult | `R/matrixfit.R`, `R/new_matrixfit.R` | Matrix-fit outputs and model metadata |
| MatrixModel and subclasses | `R/new_matrixfit.R` | OO fit model hierarchy: single, shifted, weighted, two-state, N-particle, constant |
| GEVPResult | `R/gevp.R` | Generalized eigenvalue outputs and amplitude conversion |
| HankelResult | `R/hankel.R`, `R/hankel.truncated.R` | Hankel/PGEVM spectra, coefficients, resampling |
| AutocorrelationResult | `R/computeacf.R`, `R/UWerr.R` | Autocorrelation windows, integrated tau, diagnostic plots |
| HadronDataset | `data/`, `inst/extdata/` | Packaged datasets and sample external resources |

## API naming policy

- Normal Python usage should expose classes, methods, and snake_case functions: for example `CorrelationFunction.bootstrap()` and `bootstrap_nlsfit()`.
- R symbols containing dots get snake_case wrappers: `bootstrap.cf` -> `bootstrap_cf`, `fit.plateau2cf` -> `fit_plateau2cf`, `removeTemporal.cf` -> `remove_temporal_cf`.
- S3 operator methods should become Python dunder methods where they belong to core classes, with compatibility wrappers such as `add_cf()` and `mul_raw_cf()` where useful.
- Public wrappers should carry docstrings naming the original R symbol and source path.
- Internal helper names may be reorganized only when the index explains the target owner and tests still cover the original behavior.

## Milestones

| Milestone | Scope | Exit criteria |
|---|---|---|
| 0. Inventory and planning | Replace placeholders with R-backed plan/index/summaries | This PR: exact R-backed counts, mappings, limitations, and issue breakdown recorded |
| 1. Python package skeleton | Add packaging and empty modules under `src/pydron/` | Imports work, no algorithms, metadata and tooling configured |
| 2. R reference-output infrastructure | Install/validate R oracle, freeze deterministic outputs | Selected R calls emit fixtures with seeds/tolerances |
| 3. Core data model | Port containers for raw/processed correlators and datasets | Construction, validation, serialization, and R object parity tested |
| 4. Statistical primitives | Bootstrap, jackknife, UWerr/autocorrelation, covariance | Numerical equivalence tests cover representative scalar/vector/matrix cases |
| 5. Correlation-function analysis | `cf`/`raw_cf` methods, effective masses, temporal handling | Public methods and compatibility wrappers match R outputs |
| 6. Fit machinery | NLS, cosh, plateau, matrixfit, new R6 model family | Fit results, predictions, residuals, and summaries parity-tested |
| 7. Spectroscopy methods | GEVP, Hankel, truncated PGEVM, Lanczos | Reference workflows reproduce R values on fixed fixtures |
| 8. I/O and physics utilities | Readers, CVC/Cyprus/NISSA, alpha_s, CDH, gamma/momentum utilities | Fixture I/O and numerical helpers pass equivalence tests |
| 9. Documentation and migration guide | Vignettes, examples, API mapping, migration notes | Docs link Python API to R symbols and examples run |
| 10. Release readiness | CI, packaging, dependency policy, coverage, release notes | Public API frozen for first release candidate |

## Definition of done for a ported item

- The `PORTING_INDEX.md` row has owner, status, Python target, compatibility wrapper decision, and test requirement filled in.
- The Python implementation has an idiomatic API and, where useful, an R-style compatibility wrapper.
- Numerical outputs match frozen R reference outputs within documented tolerances.
- Edge cases from the original R tests and vignettes are covered or explicitly marked `[X]` with a reason.
- Native-code policy is settled for the item if it depends on C/C++ or GSL.
- Documentation names the original R source path and symbol.

## Contributor workflow rules

- Start each implementation PR from the fixed baseline or a documented successor baseline.
- Keep one logical subsystem per PR and avoid broad API decisions without updating the index.
- Do not edit original R/C/C++ sources unless the PR is explicitly about upstream reference maintenance.
- Add or update R reference-output fixtures before marking numerical items `[E]`.
- Keep compatibility wrappers thin and traceable; do not let wrappers become divergent implementations.
- Mark uncertain items `[?]` and ask for a decision rather than silently choosing API, tolerance, or native-build policy.

## Open technical decisions

| Decision | Current status | Why it matters |
|---|---|---|
| Native code rewrite vs wrap vs defer | `[?]` | `cdh.c` uses GSL; `read_nissa_textcf_kernel.cpp` may benefit from a compiled parser; small kernels may be easier to rewrite |
| `.RData` handling | `[?]` | R metadata extraction now works, but loading policy still needs direct-load vs converted-fixture decision |
| Pythonic-first vs R-compatible-first API | `[?]` | Plan favors Pythonic-first plus wrappers, but wrapper breadth must be decided before public API freeze |
| Plotting backend | `[?]` | R plotting helpers need matplotlib or a data-first plotting abstraction |
| Typing strictness | `[?]` | Need policy for NumPy typing, optional runtime validation, and public protocol classes |
| Dependency policy | `[?]` | NumPy/SciPy/h5py/matplotlib and optional native/GSL dependencies affect installability |
| Reference fixture format | `[?]` | R can now emit metadata; next PR must choose fixture layout and tolerances |

## Final inventory summary

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

## Files still not parsed reliably

| Path | Reason | Required follow-up |
|---|---|---|
| Weighted_Model.nb | Mathematica notebook inspected as text only, not as a Wolfram expression tree | Manually review or convert before using as documentation source |
| inst/extdata/C2_bin.dat | Binary external data file; format not decoded in planning pass | Pair with reader fixture and reference checksum |
| inst/extdata/C2_pi0.dat | Binary external data file; format not decoded in planning pass | Pair with reader fixture and reference checksum |
| src/alpha_s.o | Compiled native artifact; not source-parsed | Do not port; regenerate via build if needed |
| src/cdh.o | Compiled native artifact; not source-parsed | Do not port; regenerate via build if needed |
| src/hadron.so | Compiled native artifact; not source-parsed | Do not port; regenerate via build if needed |
| src/inv_cosh.o | Compiled native artifact; not source-parsed | Do not port; regenerate via build if needed |
| src/RcppExports.o | Compiled native artifact; not source-parsed | Do not port; regenerate via build if needed |
| src/read_nissa_textcf_kernel.o | Compiled native artifact; not source-parsed | Do not port; regenerate via build if needed |
| src/tmcdh.o | Compiled native artifact; not source-parsed | Do not port; regenerate via build if needed |
