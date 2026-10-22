# PHY-PAIR-02: reproducible research attachment

This archive accompanies the article *Second-order dissociation of a Hubbard bound pair in a common spatially correlated lattice potential* by JIN HANZHE. Start with `PAPER_TO_FILES.md` for the exact paper-item / input-file / command mapping. The article itself contains the proofs. This attachment contains self-produced observations, configurations, seeds, predictions, scientific code and figure data. It excludes third-party full texts, previous manuscripts, review correspondence and before/ backups.

## Quick start

Use Python 3.12 and install `requirements.txt` in a virtual environment. The original Windows confirmation manifests record Python 3.12.14, NumPy 2.3.5 and SciPy 1.18.1. They are the source for these versions; the independent review's different environment is not substituted. Figure rendering used Matplotlib 3.11.2. Times New Roman gives the supplied figure appearance; if absent, Matplotlib can substitute a font without changing numerical values. TeX is not needed to reproduce the research data.

From this directory:

```text
python -m pip install -r requirements.txt
python reproduce.py --target verify --output reproductions/integrity
python reproduce.py --target all --output reproductions/saved_data
```

Every output directory must be new. `all` reconstructs saved-data tables, post hoc statistics, the original structural/kernel checks and figures; it does not create new dynamics or change any frozen input. Outputs and logs go under the selected directory. Thread counts are set to one by the runner. Its manifest records input hashes, environment, commands, exit status, wall time, parent CPU time and available direct-child CPU time; a wrapper child's CPU figure is not a total for its descendants. Individual scientific scripts retain their own calculation timings.

For a narrower task, select `tables`, `posthoc`, `checks` or `figures`. See the mapping table before running individual scripts: some older entry points have fixed output filenames and should be run in the isolated snapshot provided by the runner.

## Contents and interpretation

- `results/confirmation/`: 752 saved trajectories from two frozen runs, including configurations, seeds and run manifests. Each trajectory includes `times`, `bound_populations`, input parameters, allowed-mode count, noise hash and numerical invariants. It is an observable time series, not the full wave function.
- `results/leakage/`: 48 auxiliary leakage trajectories, retained as auxiliary evidence, not counted in the main table.
- `results/sample_diagnostic/`, `window_pilot/`: saved finite-ring second-order predictions and the arrays required to trace them.
- `results/analysis/`, `reanalysis/`: frozen primary estimates and later diagnostics. The post hoc estimates do not replace the frozen slope estimator or establish an exponential lifetime.
- `configs/`, `src/`, `tests/`, root scientific scripts: model, masks, seeds, original estimator and optional full replay entry point.
- `figures/`: current manuscript figures, renderer and four byte-identical frozen data files; `data/figure_data/` also retains the original exported tables and arrays.
- `checks/`, `tools/`: original kernel, normalization and support checks plus the existing three-point refinement and infinite-chain window diagnostic.
- `PROVENANCE.md`, `MANIFEST.json`: line-ending explanation and per-file SHA-256 hashes. The normalized LF protocol exactly matches both original confirmation manifests.
- `literature/`: current bibliographic metadata, DOI/arXiv access routes and recorded reading versions/hashes, without third-party article files.

Rerunning all original dynamics is optional: `python REPRODUCE_RESEARCH.py` creates a separate source snapshot and results. It reuses the original seeds, so it adds no independent statistical samples. That expensive path is not invoked by the default runner and is not required for these definition/wording revisions.

## AI workflow

`AI_USE_AND_VALIDATION.md` gives the reported model identifiers, task allocation, human verification scope and links to the executable validation materials. The main manuscript contains both a research-workflow paragraph and an AI disclosure. Generative-model access is not needed to rerun the numerical procedures in this archive.

## Release status

Public release v1.0 on Zenodo, 2026-10-15. The DOI is shown on the Zenodo record and cited in the data availability statement of the article. The author is JIN HANZHE (family name JIN, given name HANZHE), Independent Researcher, Chongqing, China. Author contact, contribution, funding, competing-interest and AI-use statements appear in the article.

Licence: Creative Commons Attribution 4.0 International (CC BY 4.0); see `LICENSE.md`. Third-party software dependencies are not included and keep their own licences.

All scientific inputs, code and run records are unchanged from the preceding reproducibility release; only this README, `LICENSE.md` and the notice in `MANIFEST.json` were updated for the public release. `REPRODUCTION_VALIDATION.json` retains the earlier saved-data reconstruction results. Run records keep the local working paths of the original computations. Internal historical drafts and review records remain outside this archive.
