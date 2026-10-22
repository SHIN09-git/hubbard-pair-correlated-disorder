# Manuscript and source

`manuscript.pdf` is the complete 18-page research manuscript, including the analytic proofs and numerical-check appendices. It is preserved byte for byte from the author's final 2026-10-15 manuscript. Public release does not imply journal acceptance. The manuscript's scientific text, figures and declared author information have not been altered for this repository.

## Build the PDF

The `source/` directory contains `main.tex`, `main.bbl`, `references.bib` and the two original figure PDFs. Use a current TeX Live or MiKTeX installation with the standard packages listed at the top of `main.tex`.

From the repository root:

```text
cd paper/source
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

The generated `source/main.pdf` is separate from the preserved `paper/manuscript.pdf`. The supplied `main.bbl` lets these commands build the complete bibliography without running BibTeX. If editing and regenerating the bibliography with BibTeX, obtain the `iopart-num` bibliography style from its official distribution; it is not bundled here. A rebuilt PDF may differ in timestamps and embedded font metadata even when its text and pagination agree.

The paper cites the independently archived data and code at https://doi.org/10.5281/zenodo.23390681 . All editable proof content is included in `source/main.tex`; there are no local-drive links needed to read it.
