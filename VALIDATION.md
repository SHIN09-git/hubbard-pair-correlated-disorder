# Public-package verification — 22 October 2026

The public package was checked in a new directory before publication. The manuscript and scientific inputs were not revised.

| Check | Result |
| --- | --- |
| Original reproducibility manifest | All 901 listed file hashes match; the archive has 902 files including its manifest. |
| Paper and reproduction payload | All 910 entries in `RELEASE_MANIFEST.json` match. |
| Saved-data reconstruction | Tables, post hoc analysis, structural checks and figures completed successfully; the saved primary and post hoc numerical comparisons match at the runner's stated tolerances. |
| Unit checks | 24 tests passed. |
| 25-point scattering-state/resolvent comparison | Maximum relative discrepancy 6.89561 × 10⁻⁵, below the recorded 7 × 10⁻⁵ threshold. |
| Additional checks | Three-point resolvent refinement, infinite-chain finite-time diagnostic and higher-order momentum-support check completed successfully. Their scope is defined in the manuscript and supplied code. |
| Manuscript build | Two pdfLaTeX passes succeeded using the supplied bibliography; 18 pages, with identical page-by-page extracted text to the preserved PDF. No undefined references or overfull-box warnings in the final log. |
| New navigation | All relative file links in the new README and paper instructions resolve. |
| Citation metadata | `CITATION.cff` passes the official CFF 1.2.0 schema. |
| Public-file selection | Submission correspondence and internal review packages were excluded; no matches in the selected files for the checked credential formats. |

The saved-data runner took 42.31 seconds of wall time on the release workstation. This is not a runtime guarantee. Its original run manifest records process timings and hashes. The public package retains the earlier validation records, and the commands in the root README generate fresh records in any new output directory.

These checks verify packaging, reproducibility of the stated calculations and consistency with the supplied manuscript. They do not constitute a new independent experiment or an external peer review. No original dynamics were rerun and no additional independent samples were produced.

Preserved manuscript SHA-256:

```text
623729ee6168b27ef23e3739935d6abf6c39866f16e82638bb123ee170246be5
```

Original data/code ZIP SHA-256:

```text
87cfc5759d7b84ee8e2ae664caad18118fe4ae35abe2d32b75ba73e30214bc17
```
