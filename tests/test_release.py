import csv
from pathlib import Path

import numpy as np

from dana_kappa.screen import cost_proxy
from dana_kappa.transport import conductivity_rows

ROOT=Path(__file__).resolve().parents[1]


def test_wigner_tensor_trace_and_temperature_axes():
    # Off-diagonal tensor entries must not enter the scalar conductivity.
    p=np.array([[[1.,2.,3.,100.,200.,300.],[4.,5.,6.,100.,200.,300.]]])
    c=p*.1
    rows=conductivity_rows([100,300],p,c,p+c)
    assert np.allclose([r['kappa_W_mK'] for r in rows],[2.2,5.5])
    assert [r['temperature_K'] for r in rows]==[100,300]


def test_published_conventional_cell_cost_example():
    r=cost_proxy(4,'m-3m','Fm-3m')
    assert r['estimated_FC3_displacements']==567
    assert np.isclose(r['relative_force_work'],7.296296296296297)
    assert r['T2_pass']


def test_released_label_counts_and_component_closure():
    with (ROOT/'data/grace_results.csv').open() as f:g=list(csv.DictReader(f))
    with (ROOT/'data/dft_results.csv').open() as f:d=list(csv.DictReader(f))
    assert len(g)==4724 and sum(0<float(r['kappa_W_mK'])<2 for r in g)==1712
    assert len(d)==60 and sum(0<float(r['kappa_W_mK'])<2 for r in d)==58
    for r in g+d:
        assert np.isclose(float(r['kappa_P_W_mK'])+float(r['kappa_C_W_mK']),float(r['kappa_W_mK']))
        if r['cif_path']:assert (ROOT/'data'/r['cif_path']).is_file()
