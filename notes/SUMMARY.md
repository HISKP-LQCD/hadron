# Folder summary: `notes/`

## Folder purpose

Project notes, performance notes, and design references.

## Porting relevance

Use notes to inform documentation and design choices, but do not treat them as tested algorithm implementations.

## Inventory table

| Category | Count | Notes |
|---|---|---|
| Files excluding summary | 2 | Design/performance notes |
| R Markdown notes | 1 | Performance note |
| Markdown notes | 1 | Design/review note |

## Mapping table

| Status | Source path | Role | Proposed Python handling | Owner | Notes |
|---|---|---|---|---|---|
| `[ ]` | notes/Matrixfit_Performance.Rmd | design/performance note | docs/design/matrixfit_performance.md | unassigned | use as design/reference input; not algorithm source by itself |
| `[ ]` | notes/gevp_review.md | design/performance note | docs/design/gevp_review.md | unassigned | use as design/reference input; not algorithm source by itself |

## Checklist

- [ ] Review notes before matrix-fit and GEVP/Hankel design PRs.
- [ ] Extract durable design decisions into docs when corresponding APIs are implemented.
- [ ] Mark stale or superseded notes `[X]` only after owner review.

## Testing requirements

| Source item | Required Python test or validation | Reference data needed | Owner/status | Notes |
|---|---|---|---|---|
| Matrixfit_Performance.Rmd | benchmark/reference note review | maybe | `[ ]` / unassigned | use during fit machinery planning |
| gevp_review.md | design consistency review | no | `[ ]` / unassigned | use during spectroscopy planning |

## Known issues and ambiguities

| Source path | Issue | Proposed resolution | Owner/status |
|---|---|---|---|
