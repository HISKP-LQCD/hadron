# Exact prompt to give Codex

```text
You are working inside the repository `pydron`, a Python port of the R package `hadron`.

Task: create a planning-only pull request. Do not implement Python algorithms yet.

Important context:
- The repository currently follows the upstream R package layout.
- The goal is to port `hadron` from R to Python while preserving the logical structure.
- The Python implementation should be object oriented internally.
- The public API should be Pythonic, but every public R function/class/method must remain traceable through explicit mapping rows and, where useful, compatibility wrappers.
- The repository already contains preliminary scaffold files: `PLAN.md`, `PORTING_INDEX.md`, and folder-level `SUMMARY.md` files. Replace placeholder content with inventory-based content after inspecting the repository. Do not preserve scaffold text that is incomplete, vague, or incorrect.

Hard constraints:
1. Do not delete, move, or rewrite existing source files.
2. Do not implement algorithms in this PR.
3. Do not add large generated artifacts.
4. Do not change numerical behaviour.
5. Do not collapse unrelated modules unless the reason is documented.
6. Every generated checklist must be based on repository inspection, not guesses.
7. If a file cannot be parsed reliably, list it explicitly with the reason.

Required outputs:
1. Update root `PLAN.md`.
2. Update root `PORTING_INDEX.md`.
3. Update or create these folder summaries:
   - `R/SUMMARY.md`
   - `src/SUMMARY.md`
   - `tests/SUMMARY.md`
   - `data/SUMMARY.md`
   - `inst/SUMMARY.md`
   - `vignettes/SUMMARY.md`
   - `exec/SUMMARY.md`
   - `hooks/SUMMARY.md`
   - `notes/SUMMARY.md`
4. If additional top-level folders exist, create a `SUMMARY.md` in each only if it is useful for port coordination. Do not create summaries for hidden build folders or dependency caches.

Repository inspection requirements:
1. List the top-level tree.
2. Inspect every file in `R/`.
3. Extract every function definition from `R/*.R`.
4. Extract every S3 method from `R/*.R`.
5. Extract every R6 class from `R/*.R`.
6. Inspect NAMESPACE, DESCRIPTION, and roxygen comments if present to identify public exports.
7. Inspect every native C/C++ file in `src/`.
8. Identify native callable routines, registration files, and R-to-native interfaces.
9. Inspect every dataset in `data/`.
10. Inspect every test file in `tests/`.
11. Inspect every vignette, R Markdown document, notebook, or workflow document in `vignettes/` and the repository root.
12. Inspect `inst/`, `exec/`, `hooks/`, and `notes/` enough to classify their role in the Python port.

`PLAN.md` requirements:
1. State the project goal and non-goals.
2. Record the exact baseline commit using `git rev-parse HEAD`.
3. Define porting principles.
4. Define status labels exactly:
   - `[ ]` not started
   - `[~]` in progress
   - `[T]` translated, not fully tested
   - `[E]` equivalence-tested against R output
   - `[D]` documented
   - `[X]` intentionally not ported, with reason
   - `[?]` needs decision
5. Define the Python package skeleton under `src/pydron/`.
6. Keep the default mapping rule: `R/foo.R` -> `src/pydron/foo.py`, unless a documented exception is necessary.
7. Define the object-oriented design, including candidate classes for correlation functions, raw correlators, bootstrap samples, jackknife samples, time series, fit results, matrix fits, GEVP, Hankel, autocorrelation, and datasets.
8. Define the API naming policy:
   - Pythonic class methods or snake_case functions for normal usage.
   - Compatibility wrappers for R-style names where useful.
   - For R names containing dots, propose snake_case wrappers.
9. Define milestones:
   - inventory and planning
   - Python package skeleton
   - R reference-output infrastructure
   - core data model
   - statistical primitives
   - correlation-function analysis
   - fit machinery
   - spectroscopy methods such as GEVP and Hankel
   - I/O and physics utilities
   - documentation and migration guide
   - release readiness
10. Define the definition of done for a ported item.
11. Define contributor workflow rules.
12. List open technical decisions, including:
   - native code rewrite vs wrap vs defer
   - `.RData` handling
   - Pythonic-first vs R-compatible-first API
   - plotting backend
   - typing strictness
   - dependency policy for NumPy/SciPy/h5py/matplotlib
13. Include a final inventory summary with exact counts.

`PORTING_INDEX.md` requirements:
1. It must be exhaustive enough to coordinate multiple contributors.
2. Include an inventory summary table with counts.
3. Include one row per R source file.
4. Include one row per parsed R symbol:
   - source path
   - R symbol
   - category: exported function, internal helper, S3 method, R6 class, generic, method, constant, dataset helper, plotting helper, I/O helper, test helper, unknown
   - proposed Python path
   - proposed Python API name
   - compatibility wrapper name if applicable
   - status label
   - owner
   - test requirement
   - notes
5. Include one row per native C/C++ file and one row per callable native routine where detectable.
6. Include one row per dataset.
7. Include one row per test file.
8. Include one row per vignette/workflow.
9. Include a section for files that could not be parsed reliably.
10. Include a first-pass issue breakdown suitable for GitHub Issues.

Folder `SUMMARY.md` requirements:
Each folder summary must contain:
1. Folder purpose.
2. Porting relevance.
3. Inventory table.
4. Mapping table where applicable.
5. Checklist of items to port or classify.
6. Testing requirements.
7. Known issues and ambiguities.
8. Owner/status fields.

Special handling by folder:
- `R/SUMMARY.md`: map R files and symbols to Python targets.
- `src/SUMMARY.md`: classify native C/C++ files as rewrite, wrap, defer, or not needed.
- `tests/SUMMARY.md`: map R tests to pytest tests and reference-output fixtures.
- `data/SUMMARY.md`: classify datasets and loading/conversion strategy.
- `inst/SUMMARY.md`: classify installed package resources.
- `vignettes/SUMMARY.md`: map R workflows to Python examples or documentation.
- `exec/SUMMARY.md`: classify executable scripts.
- `hooks/SUMMARY.md`: classify repository hooks and whether they matter for Python development.
- `notes/SUMMARY.md`: classify notes as documentation inputs, design references, or obsolete material.

Parsing guidance:
- Use scripts if helpful, but do not commit throwaway scripts unless they are clean and useful.
- Grep for common R patterns such as `<- function`, `= function`, `setGeneric`, `setMethod`, `UseMethod`, `R6Class`, and roxygen `@export`.
- Cross-check exports against `NAMESPACE` if present.
- Cross-check native routines against registration files, Rcpp exports, `.Call`, `.C`, `.Fortran`, and source files.
- Prefer conservative classification. If unsure, mark `[?]` and explain.

After editing, report in your final message:
1. Files changed.
2. Exact baseline commit.
3. Number of R files inventoried.
4. Number of R symbols found.
5. Number of S3 methods found.
6. Number of R6 classes found.
7. Number of native files and routines found.
8. Number of datasets found.
9. Number of tests found.
10. Files that could not be parsed reliably.
11. Suggested next PR after this planning PR.
```
