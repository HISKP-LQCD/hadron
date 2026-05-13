# pydron porting plan

> Annotation: This is a planning scaffold for the Python port of `hadron` to `pydron`.
> It is intentionally detailed, but it is not yet an inventory of the repository.
> Codex should replace placeholder sections after inspecting the files locally.

## 1. Project goal

Port the R package `hadron` to a Python package named `pydron`, preserving the logical structure, scientific behaviour, and reproducibility of the original package.

The port must be object oriented internally, but it should retain traceability to the original R API through compatibility wrappers and explicit source-to-target mappings.

## 2. Non-goals for the first port

- Do not redesign the physics or statistical methodology.
- Do not silently change default numerical conventions.
- Do not delete or move original R sources before parity is documented.
- Do not optimize algorithms before equivalence tests exist.
- Do not collapse unrelated R modules into a single Python file without documenting the reason.
- Do not start broad implementation work before the repository inventory is complete.

## 3. Source of truth

| Item | Value |
|---|---|
| Python port repository | `alessionegro99/pydron` |
| Upstream R repository | `HISKP-LQCD/hadron` |
| Baseline branch | `master`, unless changed deliberately |
| Baseline commit | `TBD: fill with git rev-parse HEAD` |
| Reference R version | `TBD` |
| Reference Python version | `TBD` |
| License policy | Preserve upstream GPL compatibility |

> Annotation: Fill the exact baseline commit before porting begins. Numerical parity is only meaningful relative to a fixed source commit.

## 4. Porting principles

1. Preserve logical module structure.
   - Default rule: `R/foo.R` maps to `src/pydron/foo.py`.
   - Exceptions must be documented in `PORTING_INDEX.md` and in the relevant folder `SUMMARY.md`.

2. Use object-oriented Python for core scientific concepts.
   - Public classes should represent correlation functions, raw correlators, bootstrap samples, jackknife samples, fit results, matrix fits, GEVP results, Hankel results, time series, and I/O containers where appropriate.

3. Keep compatibility wrappers.
   - R functions with names like `bootstrap.cf` should map to Pythonic methods such as `CorrelationFunction.bootstrap()`.
   - Where useful, add snake_case wrapper functions such as `bootstrap_cf()` to preserve discoverability.

4. Numerical routines require reference tests.
   - Each numerical port must compare against frozen R output, unless a written exception is added.

5. Every public item must be tracked.
   - No function, method, class, dataset, native routine, test, or vignette workflow should be ported without a corresponding checklist row.

## 5. Status labels

Use these exact labels in all checklist files.

| Label | Meaning |
|---|---|
| `[ ]` | Not started |
| `[~]` | In progress |
| `[T]` | Translated, not fully tested |
| `[E]` | Equivalence-tested against R output |
| `[D]` | Documented |
| `[X]` | Intentionally not ported, with reason |
| `[?]` | Needs decision |

## 6. Proposed Python package layout

```text
pydron/
├── R/                         # Original R sources retained as reference
├── src/
│   └── pydron/
│       ├── __init__.py
│       ├── cf.py              # from R/cf.R
│       ├── raw_cf.py          # from R/raw_cf.R
│       ├── bootstrap.py       # bootstrap and resampling utilities
│       ├── jackknife.py
│       ├── timeseries.py
│       ├── autocorrelation.py
│       ├── covariance.py
│       ├── effective_mass.py
│       ├── fits.py
│       ├── matrixfit.py
│       ├── gevp.py
│       ├── hankel.py
│       ├── luscher.py
│       ├── gamma.py
│       ├── momentum.py
│       ├── io/
│       │   ├── __init__.py
│       │   ├── hdf5.py
│       │   ├── cvc.py
│       │   ├── cyprus.py
│       │   └── nissa.py
│       ├── plotting/
│       │   ├── __init__.py
│       │   └── limits.py
│       └── utils/
│           ├── __init__.py
│           └── naming.py
├── tests/
│   ├── fixtures/
│   ├── reference_r_outputs/
│   └── test_*.py
├── docs/
├── PLAN.md
├── PORTING_INDEX.md
└── pyproject.toml
```

> Annotation: The layout above is a target proposal. Codex should adjust it only after inventorying real file groups and documenting each deviation.

## 7. Core object model

| Python class | Purpose | Likely R source |
|---|---|---|
| `RawCorrelationFunction` | Raw correlator data before derived analysis | `R/raw_cf.R` |
| `CorrelationFunction` | Correlator object and operations | `R/cf.R` |
| `BootstrapSample` | Bootstrap samples and derived statistics | bootstrap-related R files |
| `JackknifeSample` | Jackknife samples and derived statistics | jackknife-related R files |
| `TimeSeries` | Time-series statistics and autocorrelation | time-series/autocorrelation R files |
| `FitResult` | Generic fitted result container | fit-related R files |
| `NLSFitResult` | Nonlinear least-squares fit result | nlsfit-related R files |
| `MatrixFitResult` | Matrix-fit result container | matrix-fit R files |
| `GEVPResult` | Generalized eigenvalue problem output | `R/gevp.R` |
| `HankelResult` | Hankel-analysis output | Hankel-related R files |
| `AutocorrelationResult` | Autocorrelation analysis output | autocorrelation R files |
| `HadronDataset` | Loaded example/reference data | `data/` and I/O utilities |

## 8. API naming policy

For each R symbol, record both names when applicable.

| R style | Python method/function | Compatibility wrapper |
|---|---|---|
| `bootstrap.cf` | `CorrelationFunction.bootstrap()` | `bootstrap_cf()` |
| `concat.cf` | `CorrelationFunction.concat()` | `concat_cf()` |
| `effective.mass` or similar | `CorrelationFunction.effective_mass()` | `effective_mass()` |
| `computeacf` | `TimeSeries.autocorrelation()` | `compute_acf()` |

> Annotation: This table is illustrative. Codex must generate the actual mapping from parsed R symbols.

## 9. Folder strategy

### `R/`

Original R implementation. Treat this as the algorithmic reference. Do not edit except for comments needed by planning or inventory scripts.

### `src/`

Currently contains native source code in the R package. During porting, decide per file whether to rewrite in Python/NumPy/SciPy, wrap as a compiled extension, or defer.

### `tests/`

Convert R tests into `pytest` tests. Add frozen R reference outputs for numerical parity.

### `data/`

Keep original datasets. Decide whether Python should load `.RData` directly, convert to portable fixtures, or support both.

### `inst/`

Inventory installed package resources and decide which are needed in Python packaging.

### `vignettes/`

Convert R Markdown workflows into Python examples, notebooks, or documentation pages.

### `exec/`, `hooks/`, `notes/`

Inventory and classify as scripts, developer tooling, obsolete material, or documentation input.

## 10. Milestones

### Milestone 0: Inventory and planning

- [ ] Extract top-level repository tree.
- [ ] Extract all files under `R/`.
- [ ] Extract all function definitions from `R/*.R`.
- [ ] Extract all S3 methods.
- [ ] Extract all R6 classes.
- [ ] Extract all native C/C++ source files and callable routines from `src/`.
- [ ] Extract all datasets from `data/`.
- [ ] Extract all tests from `tests/`.
- [ ] Extract all vignettes and R Markdown workflows.
- [ ] Generate `PORTING_INDEX.md`.
- [ ] Generate one `SUMMARY.md` per top-level folder.

### Milestone 1: Python package skeleton

- [ ] Add `pyproject.toml`.
- [ ] Add `src/pydron/__init__.py`.
- [ ] Add package metadata.
- [ ] Add `pytest` configuration.
- [ ] Add linting with `ruff`.
- [ ] Add formatting rules.
- [ ] Add static typing policy with `mypy` or `pyright`.
- [ ] Add CI.
- [ ] Add documentation skeleton.

### Milestone 2: Reference-output infrastructure

- [ ] Add script to run selected R functions and export deterministic reference outputs.
- [ ] Add `tests/reference_r_outputs/`.
- [ ] Add tolerances for floating-point comparisons.
- [ ] Document random seed policy.
- [ ] Document platform-dependent numerical tolerances.

### Milestone 3: Core data model

- [ ] Port raw correlation-function data containers.
- [ ] Port processed correlation-function containers.
- [ ] Port bootstrap sample containers.
- [ ] Port jackknife sample containers.
- [ ] Port fit-result containers.
- [ ] Port matrix-result containers.

### Milestone 4: Statistical primitives

- [ ] Blocking.
- [ ] Bootstrap.
- [ ] Jackknife.
- [ ] Autocorrelation.
- [ ] Covariance and correlation matrices.
- [ ] Error propagation.

### Milestone 5: Correlation-function analysis

- [ ] Effective masses.
- [ ] Correlator concatenation.
- [ ] Temporal shifts and transformations.
- [ ] Symmetrization and averaging utilities.
- [ ] Plotting helpers.

### Milestone 6: Fit machinery

- [ ] Nonlinear least-squares bootstrap fits.
- [ ] Plateau fits.
- [ ] Cosh/sinh/exponential fits.
- [ ] Matrix fits.
- [ ] Fit summaries.
- [ ] Fit plotting.

### Milestone 7: Spectroscopy methods

- [ ] GEVP.
- [ ] Hankel method.
- [ ] Truncated Hankel method.
- [ ] Lüscher-related utilities.

### Milestone 8: I/O and physics utilities

- [ ] HDF5 utilities.
- [ ] CVC read utilities.
- [ ] Cyprus read utilities.
- [ ] NISSA read utilities.
- [ ] Momentum utilities.
- [ ] Gamma matrices.
- [ ] Loop tools.
- [ ] Form-factor utilities.

### Milestone 9: Documentation and migration guide

- [ ] Convert core vignettes.
- [ ] Add API documentation.
- [ ] Add examples using sample data.
- [ ] Add R-to-Python migration guide.
- [ ] Add contributor guide.

### Milestone 10: Release readiness

- [ ] All public symbols inventoried.
- [ ] High-priority numerical routines equivalence-tested.
- [ ] CI green on supported Python versions.
- [ ] Documentation builds.
- [ ] License files reviewed.
- [ ] Initial version tagged.

## 11. Definition of done for a ported item

A function, class, method, native routine, dataset, test, or vignette workflow is considered ported only when:

- It has a row in `PORTING_INDEX.md`.
- Its source path and Python target path are recorded.
- Its Python implementation exists, unless explicitly marked `[X]`.
- It has at least one test.
- Numerical routines compare against R reference output or document why this is impossible.
- Public API behaviour is documented.
- Known deviations from R are documented.
- The relevant folder `SUMMARY.md` is updated.

## 12. Contributor workflow

Each contributor should take one small checklist item or one tightly related group of functions.

Every implementation pull request must update:

- The Python implementation file.
- The relevant Python tests.
- `PORTING_INDEX.md`.
- The relevant folder `SUMMARY.md`.
- Documentation if the item is public API.

Pull requests should avoid mixing unrelated areas. For example, do not combine `R/cf.R`, native C bindings, and vignette conversion in one PR.

## 13. Open decisions

| Decision | Options | Status | Owner | Notes |
|---|---|---|---|---|
| Native code policy | rewrite / wrap / defer | `[?]` | unassigned | Decide per source file |
| `.RData` loading | direct / converted fixtures / both | `[?]` | unassigned | Needed for sample data |
| Public API style | Pythonic-first / R-compatible-first / hybrid | `[?]` | unassigned | Current recommendation: hybrid |
| Plotting backend | matplotlib only / optional backends | `[?]` | unassigned | Keep minimal initially |
| Type strictness | gradual / strict | `[?]` | unassigned | Start gradual |
| SciPy dependency policy | required / optional | `[?]` | unassigned | Likely required for fits |

## 14. Codex inventory requirements

Codex must update this file after local inspection with:

- Exact count of R source files.
- Exact count of parsed R symbols.
- Exact count of S3 methods.
- Exact count of R6 classes.
- Exact count of native source files and detected routines.
- Exact count of tests.
- Exact count of datasets.
- List of files that could not be parsed reliably.
- Proposed first ten implementation issues.
