# Folder summary: `data/`

## Folder purpose

Packaged R datasets and the LazyData listing.

## Porting relevance

Use R-backed object metadata to decide whether Python loads `.RData` directly, converts to neutral fixtures, or supports both.

## Inventory table

| Category | Count | Notes |
|---|---|---|
| Files excluding summary | 7 | 6 `.RData` plus `datalist` |
| Dataset files loaded | 6 | All succeeded in clean R environments |
| Objects inspected | 6 | Metadata only, no content dumps |
| Names in datalist | 5 | `cA2...` file is present but not listed |

## Mapping table

| Status | Source path | R object/item | Class | Dimensions | Length | Approx. size bytes | Proposed Python target/handling | Owner | Notes |
|---|---|---|---|---|---|---|---|---|---|
| `[?]` | data/cA2.09.48_3pi_I3_0_A1u_1_pc.RData | cA2.09.48_3pi_I3_0_A1u_1_pc | list;cf;cf_meta;cf_boot;cf_principal_correlator | - | 19 | 794904 | HadronDataset loader or converted fixture | unassigned | listed in datalist: no |
| `[?]` | data/correlatormatrix.RData | correlatormatrix | list;cf;cf_meta;cf_orig | - | 6 | 435160 | HadronDataset loader or converted fixture | unassigned | listed in datalist: yes |
| `[?]` | data/loopdata.RData | loopdata | list;cf;cf_meta;cf_orig;cf_smeared | - | 12 | 64568 | HadronDataset loader or converted fixture | unassigned | listed in datalist: yes |
| `[?]` | data/plaq.sample.RData | plaq.sample | numeric | - | 6352 | 50864 | HadronDataset loader or converted fixture | unassigned | listed in datalist: yes |
| `[?]` | data/pscor.sample.RData | pscor.sample | data.frame | 15168x2 | 2 | 1153608 | HadronDataset loader or converted fixture | unassigned | listed in datalist: yes |
| `[?]` | data/samplecf.RData | samplecf | list;cf;cf_meta;cf_orig | - | 6 | 205360 | HadronDataset loader or converted fixture | unassigned | listed in datalist: yes |
| `[ ]` | data/datalist | dataset listing | text | - | 5 | - | package-data manifest cross-check | unassigned | metadata file |

## Checklist

- [ ] Choose `.RData` direct loading vs conversion policy.
- [ ] Add loader tests for every packaged dataset retained in Python.
- [ ] Record checksums/object-shape expectations in the reference fixture PR.
- [ ] Reconcile `data/datalist` with the extra `cA2...RData` file.

## Testing requirements

| Source item | Required Python test or validation | Reference data needed | Owner/status | Notes |
|---|---|---|---|---|
| `.RData` datasets | loader test, object-shape check, checksum/reference metadata | yes | `[?]` / unassigned | all metadata extraction succeeded |
| `data/datalist` | manifest consistency test | no | `[ ]` / unassigned | lists 5 datasets |

## Known issues and ambiguities

| Source path | Issue | Proposed resolution | Owner/status |
|---|---|---|---|
| data/datalist | Does not list `cA2.09.48_3pi_I3_0_A1u_1_pc` | Confirm whether it should be packaged in Python | `[?]` / unassigned |
