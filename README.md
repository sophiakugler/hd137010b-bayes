# Bayesian-Stats-Project

Bayesian analysis of the K2 transit photometry of HD 137010, a tenth-magnitude
K-dwarf hosting an Earth-sized planet candidate.

Nicolas Volken, Ian Röthlisberger, Denis Titov, Sophia Kugler

Overleaf Project: https://www.overleaf.com/project/6acf34e287704edad40bbd7f

## Reproducible analysis (4 November 2026)

Use the existing environment from the repository root:

```powershell
.venv\Scripts\python.exe tests\test_analysis.py
.venv\Scripts\python.exe src\analysis.py --steps 20000 --output results\new_run
.venv\Scripts\python.exe src\refine.py --source results\refined_20261104 --output results\extended_run --max-steps 120000
```

The executed refined run is in `results/refined_20261104`. Notebooks 04 and 05 inspect these outputs. Original notebooks 01–03 and 06 and legacy arrays are preserved. Each new run must use a new output directory. Complete chains, seeds, input hashes, diagnostics, plots, derived samples and summary tables are saved. Large generated arrays stay local through `.gitignore`.

See [the audit](docs/PROJECT_STATUS.md) for priors, parameter ordering, reproducibility details and scientific caveats. `src/refine.py` extends the dated pilot run, while `src/independent_circular.py` can run a separate validation ensemble. The extra independent run started during this session was interrupted at the user's request to shorten the work; its checkpoint is not part of the reported results.

Check convergence and period bounds before interpretation. Circular reproduction of the supplied Venner samples is approximate. No TESS non-detection likelihood is used. Supplied `data/posteriors.csv` was already tracked; confirm permission before distributing it. No commit or push was made. Overleaf was not edited.
