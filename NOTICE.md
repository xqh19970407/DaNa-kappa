# Attribution and terms

The workflow code in this repository is licensed under **GPL-3.0-only**; see
[LICENSE](LICENSE). The two-stage GRACE relaxation workflow was adapted from
[MPA2suite/k_SRME](https://github.com/MPA2suite/k_SRME), also GPL-3.0.
The DANA adapter and portable screening/DFT interfaces are maintained here.

DANA uses ALiEGNN from
[Dataset_thermoconductivity_pred](https://github.com/dannyzekunren/Dataset_thermoconductivity_pred),
pinned to commit `4b41cb2986402288c94d9b42a4565e59b48edef4`.
Its source is obtained separately with `scripts/get_aliegnn.py`; it is not
redistributed in this repository. The inspected upstream snapshot did not
contain a standalone license file. Consult the authors for uses requiring
additional permission; our license does not relicense their source.

Project-generated DANA checkpoint weights, conductivity labels and derived
tables in the accompanying release are provided under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
Source crystal structures retain their original provenance and applicable
source terms: [Materials Project](https://materialsproject.org/about/terms)
and [Crystoria](https://doi.org/10.5281/zenodo.23572554). Dataset records identify
the source collection and aliases. These terms do not grant rights over
third-party source databases, software or pretrained potentials.

GRACE/TensorPotential and the `GRACE_2L_OAM_28Jan25` pretrained potential are
external dependencies, obtained through their own distribution channels.
The inspected TensorPotential distribution uses the Academic Software License;
check the terms for the version and model you obtain. GRACE source and GRACE
weights are not bundled here. VASP and PAW datasets must be supplied by a
licensed user. No VASP executable, source, POTCAR, account or access key is
included.

Please cite the DANA work, the ALiEGNN benchmark, GRACE, phono3py/phonopy and
the appropriate source databases when using the corresponding components.
See the component references in [docs/WORKFLOWS.md](docs/WORKFLOWS.md).
