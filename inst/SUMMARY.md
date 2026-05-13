# Folder summary: `inst/`

## Folder purpose

Resources installed by the R package: example data, figures, and R Markdown documents.

## Porting relevance

Classify which resources become package data, test fixtures, or documentation inputs for Python.

## Inventory table

| Category | Count | Notes |
|---|---|---|
| Files excluding summary | 13 | Installed resources |
| R Markdown resources | 2 | Documentation/workflow inputs |
| External data files | 10 | Reader fixtures |
| Figures | 1 | Package image |

## Mapping table

| Status | Source path | Proposed Python target/handling | Porting action | Owner | Notes |
|---|---|---|---|---|---|
| `[ ]` | inst/extdata/C2_bin.dat | tests/fixtures/C2_bin.dat | retain as test fixture or package data if needed | unassigned | root documentation/workflow file |
| `[ ]` | inst/extdata/C2_pi0.dat | tests/fixtures/C2_pi0.dat | retain as test fixture or package data if needed | unassigned | root documentation/workflow file |
| `[ ]` | inst/extdata/gradflow.000124 | tests/fixtures/gradflow.000124 | retain as test fixture or package data if needed | unassigned | root documentation/workflow file |
| `[ ]` | inst/extdata/gradflow.000125 | tests/fixtures/gradflow.000125 | retain as test fixture or package data if needed | unassigned | root documentation/workflow file |
| `[ ]` | inst/extdata/newdisc.0.1373.0.006.k0v4.10 | tests/fixtures/newdisc.0.1373.0.006.k0v4.10 | retain as test fixture or package data if needed | unassigned | root documentation/workflow file |
| `[ ]` | inst/extdata/outprcvn.dddd.00.0000 | tests/fixtures/outprcvn.dddd.00.0000 | retain as test fixture or package data if needed | unassigned | root documentation/workflow file |
| `[ ]` | inst/extdata/output.data | tests/fixtures/output.data | retain as test fixture or package data if needed | unassigned | root documentation/workflow file |
| `[ ]` | inst/extdata/testfile000.dat | tests/fixtures/testfile000.dat | retain as test fixture or package data if needed | unassigned | root documentation/workflow file |
| `[ ]` | inst/extdata/testfile001.dat | tests/fixtures/testfile001.dat | retain as test fixture or package data if needed | unassigned | root documentation/workflow file |
| `[ ]` | inst/extdata/testfile002.dat | tests/fixtures/testfile002.dat | retain as test fixture or package data if needed | unassigned | root documentation/workflow file |
| `[ ]` | inst/figures/hadron.png | docs/assets/hadron.png | documentation asset | unassigned | root documentation/workflow file |
| `[ ]` | inst/new_matrixfit.Rmd | docs/reference/new_matrixfit.md | convert to docs/example | unassigned | installed R Markdown resource |
| `[ ]` | inst/weighted_model.Rmd | docs/reference/weighted_model.md | convert to docs/example | unassigned | installed R Markdown resource |

## Checklist

- [ ] Use `inst/extdata` files as fixtures for reader ports where possible.
- [ ] Convert `inst/*.Rmd` after the corresponding matrix-fit APIs exist.
- [ ] Keep image assets only if needed by migrated documentation.
- [ ] Document binary fixture formats before relying on them in tests.

## Testing requirements

| Source item | Required Python test or validation | Reference data needed | Owner/status | Notes |
|---|---|---|---|---|
| inst/extdata | reader fixture tests with checksums and parsed shape comparisons | yes | `[ ]` / unassigned | 10 external data files |
| inst/*.Rmd | example/doc execution once APIs exist | yes where numerical | `[ ]` / unassigned | 2 documents |
| inst/figures | documentation asset check | no | `[ ]` / unassigned | 1 PNG |

## Known issues and ambiguities

| Source path | Issue | Proposed resolution | Owner/status |
|---|---|---|---|
| inst/extdata/C2_bin.dat | Binary or domain-specific external format not decoded in planning pass | Pair with reader implementation and expected metadata | `[?]` / unassigned |
| inst/extdata/C2_pi0.dat | Binary or domain-specific external format not decoded in planning pass | Pair with reader implementation and expected metadata | `[?]` / unassigned |
