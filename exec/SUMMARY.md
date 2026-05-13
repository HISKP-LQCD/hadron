# Folder summary: `exec/`

## Folder purpose

Executable R analysis scripts, shell helpers, and legacy workflow material.

## Porting relevance

Classify scripts as Python examples, migration references, developer tools, obsolete material, or intentionally not ported.

## Inventory table

| Category | Count | Notes |
|---|---|---|
| Files excluding summary | 99 | Executable examples/scripts/docs |
| R scripts | 89 | Analysis workflows |
| Shell scripts | 3 | File-combining helpers |
| R Markdown | 1 | Status/report template |
| Rd docs | 6 | Legacy script documentation |

## Mapping table

| Status | Source path | Current role | Proposed Python target | Owner | Notes |
|---|---|---|---|---|---|
| `[ ]` | exec/3pt.R | executable/example analysis script | examples/_3pt.py or docs/migration/_3pt.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/analyse.eta_OS.4x4.R | executable/example analysis script | examples/analyse_eta_os_4x4.py or docs/migration/analyse_eta_os_4x4.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/analyse.eta_ss.R | executable/example analysis script | examples/analyse_eta_ss.py or docs/migration/analyse_eta_ss.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/analyse.nd-kaon.8x8.R | executable/example analysis script | examples/analyse_nd_kaon_8x8.py or docs/migration/analyse_nd_kaon_8x8.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/analyse.pion.2x2.R | executable/example analysis script | examples/analyse_pion_2x2.py or docs/migration/analyse_pion_2x2.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/analyse.pion0.2x2.R | executable/example analysis script | examples/analyse_pion0_2x2.py or docs/migration/analyse_pion0_2x2.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/analyse_pipi.R | executable/example analysis script | examples/analyse_pipi.py or docs/migration/analyse_pipi.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/averx/perform-fits.R | three-point weighted-average workflow | examples/averx_perform_fits.py or docs/migration/averx_perform_fits.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/averx/perform-weighted-median.R | three-point weighted-average workflow | examples/averx_perform_weighted_median.py or docs/migration/averx_perform_weighted_median.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/eta_extrapolate.R | executable/example analysis script | examples/eta_extrapolate.py or docs/migration/eta_extrapolate.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/feta.R | executable/example analysis script | examples/feta.py or docs/migration/feta.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/hdf5-example.R | executable/example analysis script | examples/hdf5_example.py or docs/migration/hdf5_example.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/mesons-cmi/a0.R | meson CMI analysis script or docs | examples/mesons_cmi_a0.py or docs/migration/mesons_cmi_a0.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/mesons-cmi/b1.R | meson CMI analysis script or docs | examples/mesons_cmi_b1.py or docs/migration/mesons_cmi_b1.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[?]` | exec/mesons-cmi/man/b1.Rd | meson CMI analysis script or docs | examples/mesons_cmi_man_b1.py or docs/migration/mesons_cmi_man_b1.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[?]` | exec/mesons-cmi/man/plot.cfit.Rd | meson CMI analysis script or docs | examples/mesons_cmi_man_plot_cfit.py or docs/migration/mesons_cmi_man_plot_cfit.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[?]` | exec/mesons-cmi/man/rho.Rd | meson CMI analysis script or docs | examples/mesons_cmi_man_rho.py or docs/migration/mesons_cmi_man_rho.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/mesons-cmi/pion0.R | meson CMI analysis script or docs | examples/mesons_cmi_pion0.py or docs/migration/mesons_cmi_pion0.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/mesons-cmi/pionChiPTfit.R | meson CMI analysis script or docs | examples/mesons_cmi_pion_chi_ptfit.py or docs/migration/mesons_cmi_pion_chi_ptfit.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/mesons-cmi/plot.pionChiPTfit.R | meson CMI analysis script or docs | examples/mesons_cmi_plot_pion_chi_ptfit.py or docs/migration/mesons_cmi_plot_pion_chi_ptfit.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/mesons-cmi/summary.a0fit.R | meson CMI analysis script or docs | examples/mesons_cmi_summary_a0fit.py or docs/migration/mesons_cmi_summary_a0fit.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/mesons-cmi/summary.b1fit.R | meson CMI analysis script or docs | examples/mesons_cmi_summary_b1fit.py or docs/migration/mesons_cmi_summary_b1fit.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/mesons-cmi/summary.chiralfit.R | meson CMI analysis script or docs | examples/mesons_cmi_summary_chiralfit.py or docs/migration/mesons_cmi_summary_chiralfit.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/mesons-cmi/summary.pionChiPTfit.R | meson CMI analysis script or docs | examples/mesons_cmi_summary_pion_chi_ptfit.py or docs/migration/mesons_cmi_summary_pion_chi_ptfit.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/mesons-cmi/summary.rhofit.R | meson CMI analysis script or docs | examples/mesons_cmi_summary_rhofit.py or docs/migration/mesons_cmi_summary_rhofit.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/mesons-cmi/summarycfit.R | meson CMI analysis script or docs | examples/mesons_cmi_summarycfit.py or docs/migration/mesons_cmi_summarycfit.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/mesons-cmi/vector.R | meson CMI analysis script or docs | examples/mesons_cmi_vector.py or docs/migration/mesons_cmi_vector.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[?]` | exec/old/OSChiPTFit.R | legacy executable analysis script | examples/old_oschi_ptfit.py or docs/migration/old_oschi_ptfit.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[?]` | exec/old/ana.R | legacy executable analysis script | examples/old_ana.py or docs/migration/old_ana.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[?]` | exec/old/analyseOS.R | legacy executable analysis script | examples/old_analyse_os.py or docs/migration/old_analyse_os.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[?]` | exec/old/avercycle.R | legacy executable analysis script | examples/old_avercycle.py or docs/migration/old_avercycle.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[?]` | exec/old/cfunction.R | legacy executable analysis script | examples/old_cfunction.py or docs/migration/old_cfunction.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[?]` | exec/old/chiralfit.R | legacy executable analysis script | examples/old_chiralfit.py or docs/migration/old_chiralfit.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[?]` | exec/old/chiralfit2.R | legacy executable analysis script | examples/old_chiralfit2.py or docs/migration/old_chiralfit2.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[?]` | exec/old/cmfit.R | legacy executable analysis script | examples/old_cmfit.py or docs/migration/old_cmfit.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[?]` | exec/old/cont_extr.R | legacy executable analysis script | examples/old_cont_extr.py or docs/migration/old_cont_extr.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[?]` | exec/old/fit.R | legacy executable analysis script | examples/old_fit.py or docs/migration/old_fit.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[?]` | exec/old/fit_fp.R | legacy executable analysis script | examples/old_fit_fp.py or docs/migration/old_fit_fp.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[?]` | exec/old/gsl_fit.R | legacy executable analysis script | examples/old_gsl_fit.py or docs/migration/old_gsl_fit.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[?]` | exec/old/kaon.R | legacy executable analysis script | examples/old_kaon.py or docs/migration/old_kaon.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[?]` | exec/old/kaon.Rd | legacy executable analysis script | examples/old_kaon.py or docs/migration/old_kaon.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[?]` | exec/old/man/avercycle.Rd | legacy executable analysis script | examples/old_man_avercycle.py or docs/migration/old_man_avercycle.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[?]` | exec/old/man/cfunction.Rd | legacy executable analysis script | examples/old_man_cfunction.py or docs/migration/old_man_cfunction.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[?]` | exec/old/nucleon.R | legacy executable analysis script | examples/old_nucleon.py or docs/migration/old_nucleon.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[?]` | exec/old/oldplotutils.R | legacy executable analysis script | examples/old_oldplotutils.py or docs/migration/old_oldplotutils.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[?]` | exec/old/pion.R | legacy executable analysis script | examples/old_pion.py or docs/migration/old_pion.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[?]` | exec/old/pp2.R | legacy executable analysis script | examples/old_pp2.py or docs/migration/old_pp2.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[?]` | exec/old/r0.R | legacy executable analysis script | examples/old_r0.py or docs/migration/old_r0.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[?]` | exec/old/smearedpion.R | legacy executable analysis script | examples/old_smearedpion.py or docs/migration/old_smearedpion.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[?]` | exec/old/strange.R | legacy executable analysis script | examples/old_strange.py or docs/migration/old_strange.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[?]` | exec/old/summary.pionfit.R | legacy executable analysis script | examples/old_summary_pionfit.py or docs/migration/old_summary_pionfit.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[?]` | exec/old/variational.R | legacy executable analysis script | examples/old_variational.py or docs/migration/old_variational.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/online_measurements_analysis_driver_template.R | executable/example analysis script | examples/online_measurements_analysis_driver_template.py or docs/migration/online_measurements_analysis_driver_template.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/online_measurements_status_template.Rmd | executable/example analysis script | examples/online_measurements_status_template.py or docs/migration/online_measurements_status_template.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/phaseshift/analyse.R | phase-shift workflow script | examples/phaseshift_analyse.py or docs/migration/phaseshift_analyse.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/phaseshift/finish.R | phase-shift workflow script | examples/phaseshift_finish.py or docs/migration/phaseshift_finish.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/phaseshift/fit-detEq-energylevels.R | phase-shift workflow script | examples/phaseshift_fit_det_eq_energylevels.py or docs/migration/phaseshift_fit_det_eq_energylevels.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/phaseshift/get-phaseshifts.R | phase-shift workflow script | examples/phaseshift_get_phaseshifts.py or docs/migration/phaseshift_get_phaseshifts.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/phaseshift/parameters.R | phase-shift workflow script | examples/phaseshift_parameters.py or docs/migration/phaseshift_parameters.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/phaseshift/phaseshift.R | phase-shift workflow script | examples/phaseshift_phaseshift.py or docs/migration/phaseshift_phaseshift.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/phaseshift/phaseshift.pipiswave.R | phase-shift workflow script | examples/phaseshift_phaseshift_pipiswave.py or docs/migration/phaseshift_phaseshift_pipiswave.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/phaseshift/pipi.R | phase-shift workflow script | examples/phaseshift_pipi.py or docs/migration/phaseshift_pipi.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/phaseshift/plot-Energies.R | phase-shift workflow script | examples/phaseshift_plot_energies.py or docs/migration/phaseshift_plot_energies.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/phaseshift/plot-deltaell-vs-qsq.R | phase-shift workflow script | examples/phaseshift_plot_deltaell_vs_qsq.py or docs/migration/phaseshift_plot_deltaell_vs_qsq.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/phaseshift/singlepi.R | phase-shift workflow script | examples/phaseshift_singlepi.py or docs/migration/phaseshift_singlepi.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/phaseshift/summary.R | phase-shift workflow script | examples/phaseshift_summary.py or docs/migration/phaseshift_summary.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/phaseshift/test-rootfinding.R | phase-shift workflow script | examples/phaseshift_test_rootfinding.py or docs/migration/phaseshift_test_rootfinding.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/putonlinetogether.sh | executable/example analysis script | examples/putonlinetogether.py or docs/migration/putonlinetogether.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/puttogether.sh | executable/example analysis script | examples/puttogether.py or docs/migration/puttogether.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/puttogether_reverse.sh | executable/example analysis script | examples/puttogether_reverse.py or docs/migration/puttogether_reverse.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/rho-phaseshift/analyse.R | phase-shift workflow script | examples/rho_phaseshift_analyse.py or docs/migration/rho_phaseshift_analyse.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/rho-phaseshift/average.data.R | phase-shift workflow script | examples/rho_phaseshift_average_data.py or docs/migration/rho_phaseshift_average_data.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/rho-phaseshift/detect_irrep_frame.R | phase-shift workflow script | examples/rho_phaseshift_detect_irrep_frame.py or docs/migration/rho_phaseshift_detect_irrep_frame.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/rho-phaseshift/fit.delta.R | phase-shift workflow script | examples/rho_phaseshift_fit_delta.py or docs/migration/rho_phaseshift_fit_delta.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/rho-phaseshift/infile-analyse.R | phase-shift workflow script | examples/rho_phaseshift_infile_analyse.py or docs/migration/rho_phaseshift_infile_analyse.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/rho-phaseshift/infile-fit.delta.R | phase-shift workflow script | examples/rho_phaseshift_infile_fit_delta.py or docs/migration/rho_phaseshift_infile_fit_delta.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/rho-phaseshift/phaseshift.rho.R | phase-shift workflow script | examples/rho_phaseshift_phaseshift_rho.py or docs/migration/rho_phaseshift_phaseshift_rho.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/rho-phaseshift/plot-Mrho.R | phase-shift workflow script | examples/rho_phaseshift_plot_mrho.py or docs/migration/rho_phaseshift_plot_mrho.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/rho-phaseshift/plot-delta.R | phase-shift workflow script | examples/rho_phaseshift_plot_delta.py or docs/migration/rho_phaseshift_plot_delta.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/rho-phaseshift/preliminary.fit.delta.R | phase-shift workflow script | examples/rho_phaseshift_preliminary_fit_delta.py or docs/migration/rho_phaseshift_preliminary_fit_delta.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/rho-phaseshift/summarise.R | phase-shift workflow script | examples/rho_phaseshift_summarise.py or docs/migration/rho_phaseshift_summarise.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/scatteringlength/analyse_pipi.R | scattering-length workflow script | examples/scatteringlength_analyse_pipi.py or docs/migration/scatteringlength_analyse_pipi.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/scatteringlength/analysis.R | scattering-length workflow script | examples/scatteringlength_analysis.py or docs/migration/scatteringlength_analysis.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/scatteringlength/chipt-fit.R | scattering-length workflow script | examples/scatteringlength_chipt_fit.py or docs/migration/scatteringlength_chipt_fit.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/scatteringlength/deltaE_from_ratios.R | scattering-length workflow script | examples/scatteringlength_delta_e_from_ratios.py or docs/migration/scatteringlength_delta_e_from_ratios.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/scatteringlength/fit_finite_range.R | scattering-length workflow script | examples/scatteringlength_fit_finite_range.py or docs/migration/scatteringlength_fit_finite_range.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/scatteringlength/fscor.R | scattering-length workflow script | examples/scatteringlength_fscor.py or docs/migration/scatteringlength_fscor.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/scatteringlength/gather-deltaE-values.R | scattering-length workflow script | examples/scatteringlength_gather_delta_e_values.py or docs/migration/scatteringlength_gather_delta_e_values.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/scatteringlength/get_finite_range_fits-ratio.R | scattering-length workflow script | examples/scatteringlength_get_finite_range_fits_ratio.py or docs/migration/scatteringlength_get_finite_range_fits_ratio.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/scatteringlength/get_summary-efm.R | scattering-length workflow script | examples/scatteringlength_get_summary_efm.py or docs/migration/scatteringlength_get_summary_efm.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/scatteringlength/get_summary-ratio.R | scattering-length workflow script | examples/scatteringlength_get_summary_ratio.py or docs/migration/scatteringlength_get_summary_ratio.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/scatteringlength/get_summary.R | scattering-length workflow script | examples/scatteringlength_get_summary.py or docs/migration/scatteringlength_get_summary.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/scatteringlength/plot-deltaE.R | scattering-length workflow script | examples/scatteringlength_plot_delta_e.py or docs/migration/scatteringlength_plot_delta_e.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/scatteringlength/plot-mpia0.R | scattering-length workflow script | examples/scatteringlength_plot_mpia0.py or docs/migration/scatteringlength_plot_mpia0.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/scatteringlength/plot-ratios.R | scattering-length workflow script | examples/scatteringlength_plot_ratios.py or docs/migration/scatteringlength_plot_ratios.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/scatteringlength/ratio.R | scattering-length workflow script | examples/scatteringlength_ratio.py or docs/migration/scatteringlength_ratio.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/scatteringlength/summary.R | scattering-length workflow script | examples/scatteringlength_summary.py or docs/migration/scatteringlength_summary.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/scatteringlength/test-data.R | scattering-length workflow script | examples/scatteringlength_test_data.py or docs/migration/scatteringlength_test_data.md | unassigned | classify as example, migration input, developer script, or `[X]` |
| `[ ]` | exec/scatteringlength/test-ratio.R | scattering-length workflow script | examples/scatteringlength_test_ratio.py or docs/migration/scatteringlength_test_ratio.md | unassigned | classify as example, migration input, developer script, or `[X]` |

## Checklist

- [ ] Group scripts by workflow before porting; do not port one-off legacy scripts blindly.
- [ ] Promote only maintained workflows to Python examples.
- [ ] Mark obsolete or domain-specific scripts `[X]` with reason after owner review.
- [ ] Use shell helpers only as migration references unless still needed.

## Testing requirements

| Source item | Required Python test or validation | Reference data needed | Owner/status | Notes |
|---|---|---|---|---|
| Maintained R scripts | converted example smoke tests when dependencies exist | yes where numerical | `[ ]` / unassigned | current top-level and workflow directories |
| Legacy `exec/old` scripts | classification review; tests only if revived | maybe | `[?]` / unassigned | legacy material |
| Shell scripts | developer-tool smoke tests if retained | no | `[?]` / unassigned | 3 shell helpers |

## Known issues and ambiguities

| Source path | Issue | Proposed resolution | Owner/status |
|---|---|---|---|
| exec/old/** | Legacy status unclear from file names alone | Owner decision before porting or marking `[X]` | `[?]` / unassigned |
| exec/**/*.Rd | Documentation for executable scripts, not package API docs | Use only if the corresponding workflow is ported | `[?]` / unassigned |
