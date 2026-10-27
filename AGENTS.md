# AGENTS.md

Notes for coding agents working in this repository.  Read `README.md` for what
the toolbox does and `CONTRIBUTING.md` for the channel and noise file format.

## Layout

- `src/` holds the three public functions: `replay.m`, `noisegen.m`, `unpack.m`.
- `examples/` holds runnable scripts that use them.
- `tests/` holds `matlab.unittest` classes, one per function.
- `packaging/buildToolbox.m` builds the `.mltbx`; `VERSION` holds its version.

## Tests

```bash
matlab -batch "addpath('src'); assertSuccess(runtests('tests'))"
```

CI (`.github/workflows/ci.yml`) runs the same suite on R2021a and on the latest
release, with the Signal Processing Toolbox, and builds and installs the
`.mltbx` on the latest release; it also runs weekly.  Keep the code working on R2021a:
do not use functions or syntax introduced after it.  The tests use
`matlab.unittest`.  The toolbox supports MATLAB only; `src/` relies on
MATLAB's `resample` with the `'Dimension'` option, which Octave lacks.

## Parity with Python

`uwa-channels/python` implements the same functions, and
`uwa-channels/calibrations` checks every week that both give the same output on
the released Zenodo files.  A change to the behavior of a function here needs
the same change in the Python package.  Do not change argument order, index
conventions (`array_index`, `start`), or output length in one implementation
only.

## File format

The `.mat` files on Zenodo are fixed.  Code must read them as they are;
`CONTRIBUTING.md` is the specification.  Do not change the format to suit the
code.

## Conventions

- Each function starts with a help block: an H1 line, the call forms, Inputs,
  Outputs, Example, required files and toolboxes, See also, Author, Email,
  License, and a dated Revision history.  Add a line to the Revision history
  when you change the behavior of a function.
- Two-space indentation.
- Commit messages are one imperative sentence that says what changed and why,
  for example "Start the f_resamp phase ramp at zero, and test that it does".

## Releases

A tag `vX.Y.Z` that matches `VERSION` runs `.github/workflows/release.yml`,
which builds the `.mltbx` on the latest release, checks that it installs on
R2021a too, and attaches it to a GitHub release.  When you change
`VERSION`, change `version` in `CITATION.cff` to the same value.

## Citation

`CITATION.cff` cites the arXiv paper (arXiv:2609.03207) and the Zenodo data
(doi:10.5281/zenodo.21287414).  Use these when a README or a help block needs
a citation.
