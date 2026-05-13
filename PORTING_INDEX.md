# pydron porting index

> Annotation: This is a scaffold. Codex must replace placeholder rows with an exhaustive repository-derived inventory.
> The purpose of this file is coordination: every source item must have exactly one tracking row or a clearly linked group row.

## Status labels

| Label | Meaning |
|---|---|
| `[ ]` | Not started |
| `[~]` | In progress |
| `[T]` | Translated, not fully tested |
| `[E]` | Equivalence-tested against R output |
| `[D]` | Documented |
| `[X]` | Intentionally not ported, with reason |
| `[?]` | Needs decision |

## Inventory summary

| Category | Count | Notes |
|---|---:|---|
| R source files | TBD | Fill by inspection |
| R functions | TBD | Fill by parsing `R/*.R` |
| S3 methods | TBD | Fill by parsing `R/*.R` |
| R6 classes | TBD | Fill by parsing `R/*.R` |
| Native C/C++ files | TBD | Fill by inspecting `src/` |
| Native callable routines | TBD | Fill by inspecting registration and source files |
| Datasets | TBD | Fill by inspecting `data/` |
| Test files | TBD | Fill by inspecting `tests/` |
| Vignettes/workflows | TBD | Fill by inspecting `vignettes/` and root notebooks |
| Unparsed or ambiguous files | TBD | List below |

## R source files

> Annotation: Codex should generate one row per R source file. The examples below are seed mappings only.

| Status | Source path | Proposed Python target | Main responsibility | Owner | Notes |
|---|---|---|---|---|---|
| `[ ]` | `R/cf.R` | `src/pydron/cf.py` | Correlation-function object and methods | unassigned | Seed row; verify by inspection |
| `[ ]` | `R/raw_cf.R` | `src/pydron/raw_cf.py` | Raw correlation-function object and methods | unassigned | Seed row; verify by inspection |
| `[ ]` | `R/gevp.R` | `src/pydron/gevp.py` | Generalized eigenvalue problem analysis | unassigned | Seed row; verify by inspection |
| `[ ]` | `R/hankel.R` | `src/pydron/hankel.py` | Hankel analysis | unassigned | Seed row; verify by inspection |
| `[ ]` | `R/timeseries.R` | `src/pydron/timeseries.py` | Time-series analysis | unassigned | Seed row; verify by inspection |

## R symbols

> Annotation: Codex must parse functions, S3 methods, and R6 classes into this table. Do not rely on these example rows.

| Status | Source path | R symbol | Category | Proposed Python API | Compatibility wrapper | Test requirement | Owner | Notes |
|---|---|---|---|---|---|---|---|---|
| `[ ]` | `R/cf.R` | `TBD` | function/method/class | `TBD` | `TBD` | R reference output required if numerical | unassigned | Replace after parsing |
| `[ ]` | `R/raw_cf.R` | `TBD` | function/method/class | `TBD` | `TBD` | R reference output required if numerical | unassigned | Replace after parsing |

## Native source files and routines

> Annotation: Each native source file should be classified as `rewrite`, `wrap`, `defer`, or `not needed`.

| Status | Source path | Routine or file | Language | Proposed Python target | Policy | Test requirement | Owner | Notes |
|---|---|---|---|---|---|---|---|---|
| `[?]` | `src/TBD` | `TBD` | C/C++ | `TBD` | rewrite/wrap/defer | R or native reference output required | unassigned | Replace after inspection |

## Datasets

| Status | Source path | Dataset name | Proposed Python handling | Test requirement | Owner | Notes |
|---|---|---|---|---|---|---|
| `[ ]` | `data/TBD` | `TBD` | direct load / converted fixture / both | Loader test required | unassigned | Replace after inspection |

## Tests

| Status | Source path | Test scope | Proposed Python target | Owner | Notes |
|---|---|---|---|---|---|
| `[ ]` | `tests/TBD` | `TBD` | `tests/test_TBD.py` | unassigned | Replace after inspection |

## Vignettes and workflows

| Status | Source path | Workflow topic | Proposed Python target | Owner | Notes |
|---|---|---|---|---|---|
| `[ ]` | `vignettes/TBD` | `TBD` | `docs/examples/TBD.md` or notebook | unassigned | Replace after inspection |

## Unparsed or ambiguous files

| Source path | Reason | Required follow-up | Owner |
|---|---|---|---|
| `TBD` | `TBD` | `TBD` | unassigned |

## First implementation candidates

> Annotation: Codex should replace this section with concrete issues after inventory. Prefer core objects and high-value tests first.

- [ ] `TBD`: first core data model issue.
- [ ] `TBD`: first bootstrap/statistics issue.
- [ ] `TBD`: first reference-output issue.
