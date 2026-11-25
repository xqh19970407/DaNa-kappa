# DANA-kappa frozen inference model

**Name:** Distribution-Aware Neighborhood Attention (DANA-κ).
**Input:** periodic crystal structure with element identities.
**Output:** ln of the 300 K total lattice conductivity in W m^-1 K^-1;
exponentiation is performed by the inference interface.

## Architecture

The pinned ALiEGNN backbone alternates four atom/line-graph blocks, then
four atom-only updates (hidden width 256). DANA modifies incoming messages
before node reduction:

1. Four learned 32-channel slots per layer gather neighboring edge messages;
   each edge then reads the slots and receives a gated residual correction.
2. Corrected messages are reduced by sum, mean, maximum and smoothed
   population standard deviation, concatenated for the node MLP.
3. Embedding MLPs and line-graph gated updates use LayerNorm; atom-message
   and node MLPs do not acquire extra LayerNorm layers.
4. Mean pooling and a linear head predict ln κ. Training used
   `mean(sqrt(residual**2 + 0.2**2) - 0.2)` without target standardization.

The adapter keeps the author's registered module names, including dormant
baseline branches, so the existing state dictionaries load **strictly**.
These dormant branches are not extra active forward layers.

## Checkpoints and data

The [published Zenodo models archive](https://zenodo.org/records/24520420/files/dana-kappa-models-v0.1.0.zip?download=1)
contains `dana_random_s123.pt`,
`dana_space_group_s123.pt` and `dana_ood_s123.pt` under `checkpoints/`.
They were trained independently with seed 123 on the respective splits of
the Phonix-derived ALiEGNN benchmark (6,966 records overall). They are not
an ensemble. The full-source screening uses only the random-split weight.
Per-weight SHA256 values are inside the archive's `MANIFEST_MODELS.json`.

The author baseline/data source is
[Wang et al., Fast and Accurate Prediction of Lattice Thermal Conductivity
via Machine Learning Surrogates](https://arxiv.org/abs/2605.11610),
[source repository](https://github.com/dannyzekunren/Dataset_thermoconductivity_pred)
at commit `4b41cb2986402288c94d9b42a4565e59b48edef4`.
The underlying database is
[Ohnishi et al., npj Computational Materials](https://doi.org/10.1038/s41524-026-02033-w).
The source benchmark, its target and split definition must accompany
comparisons; NC's 103-crystal potential benchmark is a different task.

## Evaluation and limits

The release adapter was checked against all 4,710 archived test predictions:
random 0.36967181, space-group 0.47970736, OOD 1.18170142 MAE_ln;
arithmetic mean 0.67702686. Maximum per-record deviation from the archived
predictions was below 1.7e-6 in ln κ on the checked CPU implementation.
See `validation.json` for counts and tolerances.

Public tests were used as feedback during architecture exploration, so
these are development results. Re-loading frozen checkpoints reproduces
their predictions; it does not establish retraining stability or isolate
the causal benefit of each module. Paired ablations and repeated training
are needed for such claims. Low-κ extrapolation remains the weakest split.

Predictions do not certify dynamical stability and do not calculate force
constants, band gaps, phonon lifetimes or electronic conductivity. Inputs
need not be relaxed by GRACE or DFT first. Source, potential-relaxed and
DFT-relaxed CIFs can give different predictions.

The simple CLI prepares graphs sequentially and batches neural-network
evaluation. It is not the original many-worker full-pool timing pipeline;
do not use its one-material startup time to claim the paper's bulk throughput.
