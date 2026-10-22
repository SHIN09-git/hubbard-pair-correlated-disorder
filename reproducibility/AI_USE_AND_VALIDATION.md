# AI use and validation record

This record accompanies *Second-order dissociation of a Hubbard bound pair in a common spatially correlated lattice potential*. It describes the workflow reported for this project and the validation materials supplied with the article. It is not a transcript of all model interactions.

## Human author and model roles

The sole human author is JIN HANZHE (family name: JIN; given name: HANZHE), Independent Researcher, Chongqing, China. The following allocation was confirmed by the author on 2026-10-14:

| Contributor or model identifier | Reported role |
| --- | --- |
| JIN HANZHE | Substantial participation in result interpretation and verification of the principal conclusions; critical manuscript review; final approval and responsibility for the work. |
| `gpt-6-pro` | Initial research idea, principal scientific review and final pre-submission review. |
| `gpt-6-astra-max` | Main research work. |
| `opus5.5` | Auxiliary scientific review and manuscript drafting. |

The final pre-submission review involved `gpt-6-pro` and JIN HANZHE. This describes an internal review before submission; it is not a claim of journal peer review or independent external human review. The model identifiers are preserved as reported by the author. Immutable model build identifiers and a complete prompt history are not supplied in this archive.

## Research workflow

Generative AI was used for analytic derivations, numerical implementation and execution, analysis, literature checking, scientific review and manuscript preparation. The resulting research uses the stated Hubbard Hamiltonian, explicit analytic formulas and executable numerical procedures. It does not use an inference-time generative model as a predictor of the physical observables. Reproduction of the saved-data results therefore does not require access to the named models or their original sessions.

The author-directed scope was to audit the core contribution and prepare a minimal research paper without expanding the model. Checks addressed state normalization, exchange symmetry, momentum folding, continuum endpoints and special limits. Literature comparisons were made proposition by proposition, with novelty claims narrowed when a result or continuum factor was already known. The distinction between second-order closure and all-order protection was retained. Fixed-variance and fixed-height spectra were both reported. The principal numerical estimator, windows and seeds were frozen before the confirmatory calculations; later diagnostics remain identified as post hoc.

## Validation materials and limits

The manuscript contains the explicit proofs and describes its numerical checks in appendix C. The companion reproducibility archive provides:

- `PAPER_TO_FILES.md`, which maps claims, tables, figures and diagnostics to their source files and commands;
- the simulation, analysis and figure code, frozen configurations and seeds;
- 752 confirmatory and 48 auxiliary observable trajectories, together with the saved finite-ring second-order predictions;
- structural and kernel checks, including normalization, symmetry, folding and scattering-state/resolvent comparisons;
- the infinite-chain finite-time diagnostic, retained as post hoc evidence;
- run records, environment information and per-file hashes.

Start with the archive's `README.md`. `python reproduce.py --target verify --output reproductions/integrity` checks its recorded file hashes. `python reproduce.py --target all --output reproductions/saved_data` reconstructs saved-data tables, diagnostics, checks and figures; every output directory must be new. This does not rerun all original dynamics or generate independent statistical samples.

The author's stated verification scope is substantial participation in interpreting the results, checking the principal conclusions and critically reviewing the manuscript. This record does not assert that the author personally rederived every formula, manually inspected every numerical sample or independently repeated all computations. Agreement between AI-assisted routes alone is not treated as proof of correctness; the explicit arguments, inputs, executable checks and their limitations are provided for scrutiny. Responsibility for the manuscript rests with the human author.

This author/declaration revision changes no scientific formulas, figures, numerical tables, frozen observations or scientific software. The existing source-comparison and reproduction records remain applicable to those unchanged materials.
