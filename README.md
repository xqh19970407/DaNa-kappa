# DANA-κ

**Distribution-Aware Neighborhood Attention for lattice thermal conductivity prediction and materials screening.**

This release provides CIF-to-conductivity inference, physical and
computational-budget filters, GRACE phonon transport, portable VASP scripts,
and crystal–conductivity pairs at 300 K.

**Data and pretrained models:** [Zenodo, version v1](https://zenodo.org/records/24520420),
published November 25, 2026
([DOI: 10.5281/zenodo.24520420](https://doi.org/10.5281/zenodo.24520420)).
The archive filenames retain the software release
identifier `v0.1.0`.

| Released cohort | Records | Included structures |
|---|---:|---|
| GRACE conductivity results | 4,724 | CIFs for 1,712 records with **0 < κ < 2** |
| Completed DFT results | 60 | All relaxed structures, including **58 with κ < 2** |
| Full-source DANA predictions | 1,953,259 | Separate archive; bulk source structures via database links |
| Frozen DANA checkpoints | 3 | Random, space-group and low-κ OOD splits |

Units are W m⁻¹ K⁻¹. The scalar is the diagonal average of the Wigner
population-plus-coherence tensor, excluding electronic conductivity.
Records are source entries, not globally deduplicated compositions.
Thirteen zero-valued GRACE results are excluded from the positive low-κ set.
See [data and provenance](docs/DATA.md).

## Installation and inference

Use the pinned environment in [INSTALL.md](docs/INSTALL.md), then:

```bash
python -m pip install --no-deps -e .
python scripts/get_aliegnn.py
# Download and extract the published checkpoints.
python scripts/get_assets.py models
python -m dana_kappa.predict --manifest examples/manifest.csv \
  --checkpoint checkpoints/dana_random_s123.pt --device cpu \
  --batch-size 3 --output predictions.csv
```

ALiEGNN source is downloaded separately from a pinned author commit; no
editing of installed packages is required. Supply `--upstream PATH` for an
existing author source checkout. Use `--device cuda:0` with a compatible
GPU installation. A manifest contains `material_id,cif_path`, with CIF paths
relative to the manifest directory.

## Screening and transport

```bash
python -m dana_kappa.screen predictions.csv --metadata examples/metadata.csv \
  --dana-threshold 8 --output selected.csv
python -m dana_kappa.grace examples/BaPd.source.cif \
  --output runs/BaPd_grace --device cpu --threads 4
python -m dana_kappa.vasp --phase prepare --cif examples/CsKAu2.source.cif \
  --workdir runs/CsKAu2_dft --potcar-root /path/to/your/PAW_PBE \
  --paw-map configs/paw_example_CsKAu2.json
```

GRACE and VASP use separate environments. VASP requires the user's licensed
executable and PAW files. [WORKFLOWS.md](docs/WORKFLOWS.md) explains subsequent
phases, archived-force-constant transport, duplicate exclusion and budget
formulas. Scripts do not manage remote accounts or submit jobs automatically.

## Model and validation

DANA adds two-pass neighborhood slot attention and sum/mean/max/standard-
deviation aggregation to an ALiEGNN atom/line-graph backbone, with LayerNorm
in embedding and line-graph updates and a smooth pseudo-Huber loss.

| Frozen checkpoint | Test records | MAE of ln κ |
|---|---:|---:|
| Random | 1,403 | 0.369672 |
| Space group | 1,393 | 0.479707 |
| Low-κ OOD | 1,914 | 1.181701 |
| Arithmetic mean | — | **0.677027** |

The public adapter strictly loads the original weights. Public test sets
were used during architecture exploration: these are development results,
not a new blind evaluation or a claim of consistent improvement across
retrained seeds. See [MODEL_CARD.md](docs/MODEL_CARD.md) and
[validation results](docs/validation.json).

The GRACE workflow was checked from a source BaPd CIF through relaxation,
force constants and conductivity. DFT input generation and postprocessing
were checked using archived force results; a new full VASP calculation
was not executed for this release.

## Downloads and reuse

Small CIFs and tables are included in GitHub. The complete release is available
from [Zenodo](https://zenodo.org/records/24520420):

| Archive | Contents | Size |
|---|---|---:|
| [CIFs and results](https://zenodo.org/records/24520420/files/dana-kappa-cifs-results-v0.1.0.zip?download=1) | Structures and GRACE/DFT label tables | 5.4 MB |
| [DANA models](https://zenodo.org/records/24520420/files/dana-kappa-models-v0.1.0.zip?download=1) | Three split-specific checkpoints | 260.0 MB |
| [Full-pool predictions](https://zenodo.org/records/24520420/files/dana-kappa-fullpool-predictions-v0.1.0.zip?download=1) | 1,953,259 source-record predictions | 208.4 MB |
| [DFT force constants](https://zenodo.org/records/24520420/files/dana-kappa-dft58-force-constants-v0.1.0.zip?download=1) | FC2/FC3 and reference outputs for 58 low-conductivity phases | 286.2 MB |

Download and extract individual archives from the repository root:

```bash
python scripts/get_assets.py models
python scripts/get_assets.py predictions --dest downloads/fullpool
python scripts/get_assets.py force_constants --dest downloads/dft58
# Optional: obtain the same CIFs/tables already included in this repository.
python scripts/get_assets.py data --dest downloads/cifs-results
```

The script verifies SHA256 before extraction. Published MD5 values and SHA256
checksums are listed in [configs/assets.json](configs/assets.json). A previously
downloaded ZIP can be supplied with `--archive /path/to/file.zip`.
See [data citation](docs/DATA_CITATION.bib) when citing the dataset.

Code: GPL-3.0-only. Project-generated labels and DANA weights: CC BY 4.0.
Source structures and external models retain their terms. Read
[NOTICE.md](NOTICE.md); VASP/PAW files and GRACE weights are not redistributed.

中文：本仓库提供DANA推理、规则与成本筛选、GRACE输运、VASP接口，以及1712条
低热导率GRACE记录和60项DFT完成结果的结构。四个数据与模型包已在
[Zenodo v1](https://zenodo.org/records/24520420)公开，可用上述命令直接下载。
源结构、GRACE优化结构和DFT优化结构分别保存，不能混用标签。
