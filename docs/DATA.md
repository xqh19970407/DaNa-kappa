# Data, structures and provenance

All conductivity labels refer to 300 K and units W m^-1 K^-1. The target is
`trace(kappa_P + kappa_C)/3` for Wigner/RTA lattice transport; no electronic
conductivity is included. CSVs are convenient for analysis; matching JSONL
files preserve numeric and Boolean types. All paths in these tables are
relative to `data/`.

| File | Meaning |
|---|---|
| `grace_results.csv`, `.jsonl` | All 4,724 GRACE numerical results, including 13 zero-valued records |
| `grace_low2.csv` | 1,712 records satisfying **0 < κ < 2** |
| `dft_results.csv`, `.jsonl` | 60 completed DFT results: 59 independent-static-v3, one older-v2 SrHCl result |
| `dft_low2.csv` | The paper's 58 DFT results below 2 |
| `dft_parameters_58.csv` | Material-specific cutoffs, PAW labels, k meshes, force-task counts and observed costs |
| `frozen_cost_candidates_6170.csv` | Original 6,170 candidates with frozen historical selection flags |
| `registered_dft_ids.csv` | 123 historical registered IDs for optional repeated-submission exclusion |

The 60 completed structures include two above 2 (SrHCl and YTeI). Registered
IDs include unsuccessful and unfinished candidates and are not completed
conductivity labels. `convergence_validated=false` identifies uniform
screening results whose individual convergence studies remain outstanding.
GRACE values are potential-derived labels, not DFT ground truth. Read the
quality fields before using records for training; zero values cannot be
used in a logarithmic target.

## CIF roles

* `source_cif_path`: original source input used by the direct DANA predictor.
* `cif_path`: corresponding relaxed primitive geometry for the GRACE/DFT label.
* `calculation_cell_cif_path`: conventional GRACE cell used to build force
  supercells, when included.

The archive contains 1,712 GRACE-relaxed primitive CIFs, 1,712 corresponding
conventional CIFs, 1,713 source CIFs and 60 DFT-relaxed CIFs (5,197 files).
Source CIFs overlap across the two result cohorts; many structures also
share a composition. Neither file count nor composition is a count of
globally unique materials. Use source aliases and structural comparison.

`cif_sha256` identifies the relaxed CIF. `DANA_source_kappa_W_mK` is the frozen
random-split model prediction on the **source** structure. `GRACE_kappa_W_mK`
in the DFT table is the corresponding potential result, not an evaluation
on the DFT-relaxed geometry. Population/coherence columns sum to the total.

RbAu retains conditional acceptance of a small acoustic instability. Denser
paths also reveal small negative frequencies for RbTlS and BaTlP. These flags
are retained; passage of the original mesh screen does not prove stability
at every wave vector or finite temperature. A field named `strict_pass`
means passage of the archived numerical tolerance, not an exact positivity
test. See the DFT workflow for the actual tolerances.

## Historical funnel versus retrospective filtering

The 1,712-record set includes all positive GRACE results below 2 from the
original screening. A retrospective DANA < 8 filter keeps all 58 accepted
low-κ DFT structures but intersects that GRACE set in 1,693 records.
These scopes must not be interchanged. DANA < 2 keeps only 50 of the 58
DFT low-κ structures. The figure's stricter cutoff is not a full-recall
claim. The DANA analysis was performed after the original GRACE/DFT selection.

For the frozen 1,953,259-source snapshot: DANA < 8 retains 1,627,135 records;
source physical rules retain 47,588; intersection with the original cost
selection retains 4,520, including 3,410 GRACE numerical results. The portable
screening script recomputes structure-based descriptors; exact historical
counts require the original snapshot and frozen flags.

## Obtain other source structures

* [Materials Project](https://materialsproject.org/),
  [API documentation](https://docs.materialsproject.org/downloading-data/using-the-api):
  use numeric `mp-...` aliases from `aliases` or `numeric_alias`.
* [Crystoria data record](https://doi.org/10.5281/zenodo.23572554),
  [dataset](https://huggingface.co/datasets/kexinpang/Crystoria),
  [paper](https://arxiv.org/abs/2610.38195): use the canonical structure IDs.
  The current public dataset and the project's archived 1,798,882-record
  source snapshot can differ; preserve version and CIF hashes when comparing.

The full-pool predictions archive contains 154,377 MP and 1,798,882 Crystoria
records. Bulk source structures are not mirrored here. Source-specific
terms apply; see `NOTICE.md`. The project archives are published in
[Zenodo version v1](https://zenodo.org/records/24520420), dated November 25,
2026, DOI [10.5281/zenodo.24520420](https://doi.org/10.5281/zenodo.24520420).
Direct file URLs, sizes, SHA256 and published MD5 checksums are listed
in `configs/assets.json`; use `scripts/get_assets.py` to download them.
The dataset citation is available in `docs/DATA_CITATION.bib`.
