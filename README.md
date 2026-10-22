# Hubbard-pair dissociation in a correlated lattice potential

**Second-order dissociation of a Hubbard bound pair in a common spatially correlated lattice potential**  
**JIN HANZHE** — Independent Researcher, Chongqing, China  
Research manuscript / preprint, 2026. This is not an accepted journal publication.

[Read the paper (PDF)](paper/manuscript.pdf) · [LaTeX source](paper/source/main.tex) · [Evidence and reproduction map](reproducibility/PAPER_TO_FILES.md) · [Archived data and code](https://doi.org/10.5281/zenodo.23390681)

## What the paper establishes

For two repulsively interacting bosons in a one-dimensional Bose–Hubbard chain, both particles experience the same weak, static, stationary Gaussian potential. The paper derives a closed-form second-order dissociation kernel and uses its positivity and folded momentum support to classify when high-pass disorder closes dissociation or elastic backscattering. It also derives the endpoint powers and the fixed-variance plateau, with explicit proofs and finite-ring numerical checks.

The results concern the coefficient of the disorder strength squared. A zero second-order coefficient is not a claim of all-order protection or an infinite lifetime. Fixed on-site variance and fixed pass-band height are treated separately. Post hoc diagnostics are identified separately from the frozen primary estimator. The paper describes the connection to existing dissociation formulas and the limits of its literature comparison.

## Read and inspect

| Item | Location |
| --- | --- |
| Complete paper, including proofs and checks (18 pages) | [`paper/manuscript.pdf`](paper/manuscript.pdf) |
| Editable manuscript, bibliography and figures | [`paper/source/`](paper/source/) |
| Plain-text abstract | [`paper/abstract.txt`](paper/abstract.txt) |
| Exact mapping from paper items to evidence and commands | [`reproducibility/PAPER_TO_FILES.md`](reproducibility/PAPER_TO_FILES.md) |
| Frozen observations, configurations, seeds and scientific programs | [`reproducibility/`](reproducibility/) |
| AI roles, human review scope and validation limits | [`reproducibility/AI_USE_AND_VALIDATION.md`](reproducibility/AI_USE_AND_VALIDATION.md) |
| Origin and file integrity of this public edition | [`RELEASE_MANIFEST.json`](RELEASE_MANIFEST.json) |
| Checks performed for this public release | [`VALIDATION.md`](VALIDATION.md) |

The reproducibility directory is the byte-for-byte content of the existing [Zenodo data/code release](https://doi.org/10.5281/zenodo.23390681). It contains 752 confirmatory and 48 auxiliary observable trajectories, not full wave-function histories. Historical run records retain their original working paths; the documented reproduction entry points resolve inputs relative to the archive and work after relocation. The canonical paper and bibliography for this edition are under `paper/`.

## Reproduce the saved-data results

Use Python 3.12 and a fresh virtual environment. From the repository root:

```text
cd reproducibility
python -m pip install -r requirements.txt
python reproduce.py --target verify --output reproductions/integrity
python reproduce.py --target all --output reproductions/saved_data
```

Each output directory must be new. The first command checks the frozen hashes. The `all` target reconstructs tables, post hoc statistics, structural/kernel checks and figures from the saved inputs. It does not rerun all original dynamics or add independent samples. Narrower targets and the optional expensive dynamics replay are documented in [`reproducibility/README.md`](reproducibility/README.md).

The original records specify Python 3.12.14, NumPy 2.3.5, SciPy 1.18.1 and Matplotlib 3.11.2. Missing Times New Roman may change the font used in regenerated figures without changing the numbers. No generative-model access is required to run the numerical procedures. For the PDF build, see [`paper/README.md`](paper/README.md).

## Citation and licence

For the research manuscript, cite JIN HANZHE (2026), *Second-order dissociation of a Hubbard bound pair in a common spatially correlated lattice potential*, preprint, this repository. [`CITATION.cff`](CITATION.cff) supplies machine-readable author and version information. The data/code archive has its own identifier and should be cited separately: [10.5281/zenodo.23390681](https://doi.org/10.5281/zenodo.23390681).

The manuscript and original repository materials are available under [CC BY 4.0](LICENSE.md). Third-party software dependencies retain their own licences.

## Authorship and AI disclosure

The sole author is **JIN HANZHE** (family name **JIN**, given name **HANZHE**), contact: jinhanzhe09@gmail.com. Generative AI was used extensively in the research and writing. As reported by the author, `gpt-6-pro` supplied the initial idea and principal review; `gpt-6-astra-max` carried out the main research; `opus5.5` assisted with review and drafted the manuscript. The human author participated substantially in interpreting the results, checking the principal conclusions and critically reviewing the manuscript, approved it, and accepts responsibility. The manuscript and the [workflow record](reproducibility/AI_USE_AND_VALIDATION.md) state the detailed scope and limitations. No specific funding or competing interests were declared.
