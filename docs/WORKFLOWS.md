# Reproducing the multilevel workflow

## 1. Source CIF to DANA prediction

Follow `INSTALL.md`; download the DANA models archive and obtain the pinned
ALiEGNN source. Run the example in the root README. CIFs are used as supplied:
the direct predictor does not call GRACE or perform structural relaxation.
Use the random-split checkpoint for the released full-source screening.
Example expected predictions are in `examples/expected_predictions.csv`.
These are predictions of the total lattice conductivity, not force constants.

## 2. Source-property rules, cost and duplicates

`dana_kappa.screen` joins source metadata by `material_id`, recomputes
primitive size/symmetry from each CIF, then applies a DANA threshold and the
frozen central T2 cost proxy. Supply these source-property fields:

| Source | Required metadata |
|---|---|
| MP | `source_collection=mp`, `band_gap_eV`, `energy_above_hull_eV_atom`, `is_metal`, `is_magnetic` |
| Crystoria | `source_collection=crystoria`, `source_nonmetal` |

The primitive must contain at most 20 atoms. The MP rule requires gap > 0.1 eV,
hull distance <= 0.1 eV/atom, nonmetal and nonmagnetic source flags. Crystoria
uses its source nonmetal label and the historical M1 element exclusion:
V, Cr, Mn, Fe, Co, Ni, Ce–Yb and Pa–No. This compositional magnetic-risk rule
is a screening approximation, not a DFT proof of a nonmagnetic ground state.
Missing source properties do not count as a physical-rule pass.

The T2 proxy uses primitive atom count `n`, point-group order `g` and
conventional-cell centering multiplicity `c` (P=1; A/B/C/I=2; R=3; F=4):

```text
N3 = 8 c n; N2 = 64 c n
d3 = ceil(6 n (1 + 6 N3) / g^0.9)
d2 = ceil(6 n / g^0.6)
full_FC3_memory_GiB = 216 N3^3 / 2^30
relative_force_work = (d3 N3 + d2 N2) / (146*64 + 2*512)
relative_transport_work = (n/2)^3 * 48 / g_reciprocal
```

Here `g_reciprocal=g` for centrosymmetric groups and `2g` otherwise.
Retain `d3 <= 2000`, FC3 memory <= 4 GiB, force work <= 100 and transport
work <= 1000. These empirical estimates refer to conventional GRACE cells;
they are not measured runtimes or exact phono3py displacement counts.
For CsKAu2 they give `d3=567`, force proxy `7.296296`.

```bash
python -m dana_kappa.screen predictions.csv --metadata examples/metadata.csv \
  --dana-threshold 8 --output selected.csv \
  --exclude-ids data/registered_dft_ids.csv
# Optional: compare periodic structures within candidates and to existing labels.
python -m pip install -e '.[structure]'
python -m dana_kappa.deduplicate selected.csv \
  --existing data/dft_results.csv --output unique.csv
```

Exact-ID exclusion also excludes historical failed/unfinished candidates to
avoid repeated submissions. Structural matching uses pymatgen's species-aware
`StructureMatcher` with `ltol=0.2`, `stol=0.3`, `angle_tol=5`, primitive-cell
reduction and volume scaling. These settings should be reported when used;
composition alone is not a structural duplicate criterion. Literature
novelty checks require a separate literature review and are not automated here.

Frozen original counts and later retrospective DANA intersections are
documented in `DATA.md`; do not claim that all original DFT candidates were
prospectively chosen by the later DANA predictor.

## 3. GRACE forces and Wigner conductivity

In the GRACE environment:

```bash
python -m dana_kappa.grace examples/BaPd.source.cif \
  --output runs/BaPd_grace --device cpu --threads 4
```

The external pretrained model is GRACE-2L-OAM (`GRACE_2L_OAM_28Jan25`),
accessed through TensorPotential 0.4.5. The source is standardized to a
conventional cell, then relaxed with two FIRE/FrechetCellFilter stages,
each at `fmax=1e-4 eV/A`, at most 300 steps. Stage 1 uses ASE FixSymmetry;
stage 2 removes that constraint, retaining stage 1 if stage 2 changes its
space group. The cell filter varies the three normal strains. FixSymmetry's
default tolerance and the 1e-5 symmetry monitor are distinct.

FC2/FC3 supercells are conventional-cell 4^3/2^3 multiples with 0.03 A
displacements. GRACE supplies the forces; phono3py symmetrizes FC2/FC3 in
real space. The harmonic screen permits Gamma's three lowest modes down to
-0.01 THz and otherwise rejects negative frequencies. An imaginary result
is recorded as `status=imaginary`, not a completed conductivity value.

Transport uses a 19^3 mesh, tetrahedron integration, Wigner/RTA, isotope
scattering and a 1 m boundary mean free path. No NAC, four-phonon scattering
or finite-temperature force-constant renormalization is included.
Outputs include relaxed geometry, displacement YAML, cached force arrays,
FC2/FC3 HDF5, conductivity HDF5, `kappa.csv` and `result.json`.
Use a separate output directory for each structure/model/settings combination;
existing cached forces are reused when rerunning the same calculation.

## 4. DFT validation

Use a GRACE-relaxed primitive structure as the starting geometry. In the
historical DFT environment, provide licensed PAW_PBE files and VASP 6.5.0:

```bash
python -m dana_kappa.vasp --phase prepare \
  --cif data/structures/grace_low2/canonical__2a0991efbd81fb0dcf09c32c9a4d3454c6a59e58783eceff9786961a8deea776.cif \
  --workdir runs/CsKAu2_dft --potcar-root /path/to/PAW_PBE \
  --paw-map configs/paw_example_CsKAu2.json
python -m dana_kappa.vasp --phase relax --workdir runs/CsKAu2_dft \
  --vasp-command 'mpiexec -np 64 vasp_std'
python -m dana_kappa.vasp --phase displacements --workdir runs/CsKAu2_dft
python -m dana_kappa.vasp --phase fc2 --workdir runs/CsKAu2_dft
python -m dana_kappa.vasp --phase fc3 --workdir runs/CsKAu2_dft
python -m dana_kappa.vasp --phase collect --workdir runs/CsKAu2_dft
```

`--phase all` executes these phases sequentially. Each force phase runs
displacements sequentially, each with the requested MPI command. Edit the
generic `configs/slurm_dft_example.sh` for your scheduler; no account or
partition is assumed. The 24 h request in that example is editable.
For input generation without licensed PAW files, use `--encut`; this does
not permit force execution without a matching POTCAR.

The independent-static-v3 defaults are:

| Part | Settings |
|---|---|
| Electronic | PAW-PBE, non-spin-polarized, no SOC/U/dispersion; ENCUT=1.3 max ENMAX, EDIFF=1e-8 eV, PREC=Accurate, ALGO=Normal |
| Occupations | ISMEAR=0, SIGMA=0.01 eV, NELM=200, NELMIN=5 |
| Projection | LREAL=False, ADDGRID=True, LASPH=True |
| Relaxation | IBRION=2, ISIF=3, PSTRESS=0, POTIM=0.2, NSW=10 per segment; maximum 200 ionic steps; EDIFFG=-1e-4 eV/A |
| Static acceptance | Independent initial/final static; final force <= 1e-3 eV/A, stress component <= 0.05 GPa, energy rise <= 0.01 eV/atom, retained space group |
| Electronic meshes | Gamma 6^3 primitive static/relaxation; 2^3 FC2; 3^3 FC3 |
| Force supercells | Primitive 4^3 FC2 and 2^3 FC3; displacement 0.03 A; IBRION=-1, NSW=0, ISYM=0 |
| Parallelism | Default 64 MPI ranks, one thread/rank, NCORE=8 |
| Transport | phono3py 2.3.2 / phonopy 2.15.1; 19^3, tetrahedra, Wigner/RTA, isotope, 1 m boundary |

Material-specific values for the paper's 58 cases are in
`data/dft_parameters_58.csv`. PAW choices affect ENCUT; the historical range
is 286.1456–525.1077 eV. The force optimization target and independent final
static acceptance threshold are different.

The archived DFT harmonic policy (`configs/harmonic_policy.json`) has a
baseline tolerance of -0.05 THz for the three closest-to-zero Gamma modes
and -0.02 THz otherwise. A separately recorded exception accepts only tiny
near-Gamma lowest acoustic instabilities that pass ratio, location,
translation-projection and sum-rule tests. Frequencies are never clipped.
Passing this policy is not proof of complete dynamical or thermal stability.

This portable runner preserves scientific stage defaults but omits the
original site's SSH, daemon, job-array and account-recovery machinery.
VASP XML parsing, displacement generation, harmonic screening and conductivity
postprocessing were tested with archived outputs. Fresh licensed VASP
relaxation/force execution remains a user-side integration step.

## 5. Transport from the 58-material archive

Extract the force-constants archive into `downloads/dft58/`:

```bash
python scripts/get_assets.py force_constants --dest downloads/dft58
```

With the current
transport environment, point to a single extracted material folder:

```bash
python -m dana_kappa.transport --input downloads/dft58/MATERIAL_FOLDER \
  --output runs/MATERIAL_transport --temperatures 100 150 200 250 300 350 400
```

The reference HDF5 and force constants are included. Temperature sweeps reuse
fixed force constants, not temperature-renormalized interactions. Group
velocities, linewidths and mode conductivities can be read from reference
HDF5. Population mode contributions are distinct from coherence transport.

## References

* [ALiEGNN and surrogate benchmark](https://arxiv.org/abs/2605.11610)
* [GRACE](https://doi.org/10.1038/s41524-026-01979-1)
* [Foundation-potential thermal conductivity benchmark](https://doi.org/10.1038/s41467-026-76391-w)
* [Phonopy/phono3py implementation](https://doi.org/10.1088/1361-648X/acd831)
* [Wigner transport](https://doi.org/10.1103/PhysRevX.12.041011)
* [Original relaxation workflow, k_SRME](https://github.com/MPA2suite/k_SRME)
