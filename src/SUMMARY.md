# Folder summary: `src/`

## Folder purpose

Native source and R build metadata for compiled parts of the original package.

## Porting relevance

Classify each native file and registered callable routine as rewrite, wrap, defer, or not needed before implementing Python algorithms.

## Inventory table

| Category | Count | Notes |
|---|---|---|
| Files excluding summary | 19 | Includes Makevars, source/header, compiled artifacts |
| C/C++ source/header files | 8 | Inspected |
| Registered native routines | 5 | From `RcppExports.cpp` |
| R-to-native calls | 5 | All `.Call` |
| Compiled artifacts | 7 | Not source-port targets |

## Mapping table

| Status | Source path | Proposed Python target | Porting action | Owner | Notes |
|---|---|---|---|---|---|
| `[ ]` | src/alpha_s.c | src/pydron/native/alpha_s.py or src/pydron/_native/alpha_s.* | rewrite first, wrap only if parity/performance requires | unassigned | small numerical kernels with direct R wrappers |
| `[ ]` | src/alpha_s.h | src/pydron/native/alpha_s.py or src/pydron/_native/alpha_s.* | wrap/defer | unassigned | header for C implementation; keep with owning C file decision |
| `[?]` | src/cdh.c | src/pydron/native/finite_size.py or compiled extension | needs decision: wrap or rewrite | unassigned | GSL-heavy finite-size correction implementation |
| `[ ]` | src/cdh.h | src/pydron/native/finite_size.py or compiled extension | wrap/defer | unassigned | header for C implementation; keep with owning C file decision |
| `[ ]` | src/inv_cosh.c | src/pydron/native/inv_cosh.py or scipy-backed helper | rewrite first, wrap only if parity/performance requires | unassigned | small numerical kernels with direct R wrappers |
| `[ ]` | src/RcppExports.cpp | build-generated binding file only if wrapping | defer/regenerate | unassigned | generated Rcpp registration; regenerate only if native wrapping is chosen |
| `[ ]` | src/read_nissa_textcf_kernel.cpp | src/pydron/io/nissa.py or compiled parser | rewrite or wrap after I/O fixture review | unassigned | Rcpp text parser for NISSA correlator files |
| `[ ]` | src/tmcdh.c | src/pydron/native/finite_size.py or compiled extension | defer | unassigned | not registered in current R-to-native interface |
| `[X]` | src/alpha_s.o | not ported | build artifact | unassigned | regenerate from source if needed |
| `[X]` | src/cdh.o | not ported | build artifact | unassigned | regenerate from source if needed |
| `[X]` | src/hadron.so | not ported | build artifact | unassigned | regenerate from source if needed |
| `[X]` | src/inv_cosh.o | not ported | build artifact | unassigned | regenerate from source if needed |
| `[X]` | src/RcppExports.o | not ported | build artifact | unassigned | regenerate from source if needed |
| `[X]` | src/read_nissa_textcf_kernel.o | not ported | build artifact | unassigned | regenerate from source if needed |
| `[X]` | src/tmcdh.o | not ported | build artifact | unassigned | regenerate from source if needed |
| `[?]` | src/Makevars | pyproject/scikit-build config if native wrapping is chosen | defer | unassigned | R build metadata with GSL flags |
| `[?]` | src/Makevars.in | pyproject/scikit-build config if native wrapping is chosen | defer | unassigned | R build metadata with GSL flags |
| `[?]` | src/Makevars.win | pyproject/scikit-build config if native wrapping is chosen | defer | unassigned | R build metadata with GSL flags |

## Checklist

- [ ] Decide native policy for each registered `.Call` routine before porting dependent R wrappers.
- [ ] Keep compiled artifacts out of Python source mapping.
- [ ] Use R/native fixtures for `alphas`, `invcosh`, `cdh_c`, `cdhnew_c`, and NISSA reader behavior.
- [ ] Document GSL dependency impact if wrapping `cdh.c` or `tmcdh.c`.

## Testing requirements

| Source item | Required Python test or validation | Reference data needed | Owner/status | Notes |
|---|---|---|---|---|
| Registered `.Call` routines | R-to-native parity fixture and Python wrapper/rewrite comparison | yes | `[?]` / unassigned | 5 registered routines |
| Internal C/C++ helpers | Covered through registered-routine tests if retained | yes if wrapped | `[?]` / unassigned | Source-level declarations/helpers listed below |
| Build metadata | Build/install smoke tests if native extensions are added | no | `[?]` / unassigned | Makevars/configure.ac inform dependency policy |

## Known issues and ambiguities

| Source path | Issue | Proposed resolution | Owner/status |
|---|---|---|---|
| src/alpha_s.o | Compiled artifact is present in checkout | Do not port; regenerate or ignore in Python build | `[X]` / unassigned |
| src/cdh.o | Compiled artifact is present in checkout | Do not port; regenerate or ignore in Python build | `[X]` / unassigned |
| src/hadron.so | Compiled artifact is present in checkout | Do not port; regenerate or ignore in Python build | `[X]` / unassigned |
| src/inv_cosh.o | Compiled artifact is present in checkout | Do not port; regenerate or ignore in Python build | `[X]` / unassigned |
| src/RcppExports.o | Compiled artifact is present in checkout | Do not port; regenerate or ignore in Python build | `[X]` / unassigned |
| src/read_nissa_textcf_kernel.o | Compiled artifact is present in checkout | Do not port; regenerate or ignore in Python build | `[X]` / unassigned |
| src/tmcdh.o | Compiled artifact is present in checkout | Do not port; regenerate or ignore in Python build | `[X]` / unassigned |
| src/tmcdh.c | C/GSL helper is not registered in current R interface | Defer until a caller is found or mark intentionally not ported | `[?]` / unassigned |
| src/cdh.c | GSL-heavy implementation may be expensive to rewrite exactly | Choose wrap vs rewrite after fixture benchmarks | `[?]` / unassigned |

## Source-level native declarations and helpers

| Status | Source path | Routine/helper | Return type | Kind | Proposed Python target | Owner | Notes |
|---|---|---|---|---|---|---|---|
| `[?]` | src/alpha_s.c | alphas | SEXP | definition | src/pydron/native/alpha_s.py or src/pydron/_native/alpha_s.* | unassigned | registered callable |
| `[?]` | src/alpha_s.c | alphas_c | double | definition | src/pydron/native/alpha_s.py or src/pydron/_native/alpha_s.* | unassigned | source-level helper/declaration |
| `[?]` | src/alpha_s.h | alphas_c | double | declaration | src/pydron/native/alpha_s.py or src/pydron/_native/alpha_s.* | unassigned | source-level helper/declaration |
| `[?]` | src/cdh.c | my_errhandler | void | definition | src/pydron/native/finite_size.py or compiled extension | unassigned | source-level helper/declaration |
| `[?]` | src/cdh.c | g1 | double | definition | src/pydron/native/finite_size.py or compiled extension | unassigned | source-level helper/declaration |
| `[?]` | src/cdh.c | g1array | int | definition | src/pydron/native/finite_size.py or compiled extension | unassigned | source-level helper/declaration |
| `[?]` | src/cdh.c | cdh_c | SEXP | declaration | src/pydron/native/finite_size.py or compiled extension | unassigned | registered callable |
| `[?]` | src/cdh.c | cdhnew_c | SEXP | declaration | src/pydron/native/finite_size.py or compiled extension | unassigned | registered callable |
| `[?]` | src/cdh.h | g1 | double | declaration | src/pydron/native/finite_size.py or compiled extension | unassigned | source-level helper/declaration |
| `[?]` | src/cdh.h | g1array | int | declaration | src/pydron/native/finite_size.py or compiled extension | unassigned | source-level helper/declaration |
| `[?]` | src/inv_cosh.c | invcosh | SEXP | definition | src/pydron/native/inv_cosh.py or scipy-backed helper | unassigned | registered callable |
| `[?]` | src/RcppExports.cpp | read_nissa_textcf_kernel | NumericMatrix | definition | build-generated binding file only if wrapping | unassigned | source-level helper/declaration |
| `[?]` | src/RcppExports.cpp | _hadron_read_nissa_textcf_kernel | SEXP | definition | build-generated binding file only if wrapping | unassigned | registered callable |
| `[?]` | src/RcppExports.cpp | alphas | SEXP | declaration | build-generated binding file only if wrapping | unassigned | registered callable |
| `[?]` | src/RcppExports.cpp | cdh_c | SEXP | declaration | build-generated binding file only if wrapping | unassigned | registered callable |
| `[?]` | src/RcppExports.cpp | cdhnew_c | SEXP | definition | build-generated binding file only if wrapping | unassigned | registered callable |
| `[?]` | src/RcppExports.cpp | invcosh | SEXP | definition | build-generated binding file only if wrapping | unassigned | registered callable |
| `[?]` | src/RcppExports.cpp | R_init_hadron | void | definition | build-generated binding file only if wrapping | unassigned | source-level helper/declaration |
| `[?]` | src/read_nissa_textcf_kernel.cpp | read_nissa_textcf_kernel | NumericMatrix | declaration | src/pydron/io/nissa.py or compiled parser | unassigned | source-level helper/declaration |
| `[?]` | src/tmcdh.c | w | double | definition | src/pydron/native/finite_size.py or compiled extension | unassigned | source-level helper/declaration |
| `[?]` | src/tmcdh.c | Jb | double complex | definition | src/pydron/native/finite_size.py or compiled extension | unassigned | source-level helper/declaration |
| `[?]` | src/tmcdh.c | Jb1 | double complex | definition | src/pydron/native/finite_size.py or compiled extension | unassigned | source-level helper/declaration |
| `[?]` | src/tmcdh.c | x1 | double | definition | src/pydron/native/finite_size.py or compiled extension | unassigned | source-level helper/declaration |
| `[?]` | src/tmcdh.c | calc_R_int | void | definition | src/pydron/native/finite_size.py or compiled extension | unassigned | source-level helper/declaration |
