# Environments

Run these commands from the repository root. Keep inference, GRACE and the
historical DFT postprocessor in separate environments.

## DANA inference

The checked configuration uses Python 3.10, PyTorch 2.2.1 and DGL 2.1.0:

```bash
conda create -n dana-inference python=3.10 -y
conda activate dana-inference
python -m pip install torch==2.2.1 --index-url https://download.pytorch.org/whl/cpu
python -m pip install torchdata==0.7.1
python -m pip install dgl==2.1.0 -f https://data.dgl.ai/wheels/repo.html
python -m pip install -e '.[inference]'
python scripts/get_aliegnn.py
```

For CUDA 12.1, replace the Torch index with
`https://download.pytorch.org/whl/cu121` and use the DGL wheel page
`https://data.dgl.ai/wheels/torch-2.2/cu121/repo.html`. Match Torch/DGL builds.
The adapter was checked with jarvis-tools 2025.5.30, e3nn 0.5.8 and Pydantic 1.x.

## GRACE and current transport API

```bash
conda create -n dana-grace python=3.11 -y
conda activate dana-grace
python -m pip install 'numpy<2' tensorflow==2.16.2 tensorpotential==0.4.5
python -m pip install -e '.[transport]'
```

This pins phono3py 3.24.0 and phonopy 2.47.1. The external model is
`GRACE_2L_OAM_28Jan25`, downloaded by TensorPotential on first use under its
own terms. Its weights are not included here. For GPU execution, provide a
compatible TensorFlow CUDA runtime and use
`CUDA_VISIBLE_DEVICES=0 python -m dana_kappa.grace ... --device gpu`.

## Historical DFT postprocessor

```bash
conda create -n dana-dft python=3.9 -y
conda activate dana-dft
python -m pip install 'numpy<2' ase==3.26.0 h5py spglib PyYAML \
  phonopy==2.15.1 phono3py==2.3.2
python -m pip install --no-deps -e .
```

The runner expects `phono3py` beside the active Python executable; override
with `--phono3py /path/to/phono3py` if necessary. Configure your licensed VASP
6.5.0, MPI environment and PAW_PBE directory separately. Newer backend versions
can change generated displacements and numerical results. Archived force
constants can also be read by the current transport API; bitwise identity
across versions is not promised.

## Lightweight tests

```bash
python -m pip install -e . pytest
python -m pytest -q
```

Tests cover tensor reduction, a cost example and dataset consistency.
Full checkpoint and transport checks are recorded in `docs/validation.json`.
