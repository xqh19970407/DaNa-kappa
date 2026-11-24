"""Wigner/RTA transport from existing phono3py YAML and force constants."""
import argparse
import csv
import json
import os
from pathlib import Path


def conductivity_rows(temperatures,population,coherence,total):
    import numpy as np
    def tensors(value):
        array=np.asarray(value)
        return array[0] if array.ndim==3 else array
    return [dict(temperature_K=float(t),kappa_P_W_mK=float(np.mean(p[:3])),
                 kappa_C_W_mK=float(np.mean(c[:3])),kappa_W_mK=float(np.mean(all_components[:3])))
            for t,p,c,all_components in zip(temperatures,tensors(population),tensors(coherence),tensors(total))]


def evaluate(ph3,output,temperatures=(300,),mesh=(19,19,19)):
    import numpy as np
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=True)
    ph3.mesh_numbers=list(mesh)
    old=Path.cwd()
    try:
        os.chdir(output)
        ph3.init_phph_interaction(symmetrize_fc3q=False)
        ph3.run_thermal_conductivity(temperatures=list(temperatures),is_isotope=True,
             conductivity_type='wigner',boundary_mfp=1e6,write_kappa=True)
    finally:os.chdir(old)
    k=ph3.thermal_conductivity
    # API arrays include an integration-sigma axis; average xx/yy/zz only.
    rows=conductivity_rows(temperatures,k.kappa_P_RTA,k.kappa_C,k.kappa_TOT_RTA)
    with (output/'kappa.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (output/'result.json').write_text(json.dumps(dict(target='Wigner total RTA trace/3',mesh=list(mesh),
        isotope=True,boundary_mfp_um=1000000,NAC=False,rows=rows),indent=2)+'\n')
    return rows


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',type=Path,required=True)
    ap.add_argument('--output',type=Path,default=Path('transport'))
    ap.add_argument('--temperatures',type=float,nargs='+',default=[300])
    ap.add_argument('--mesh',type=int,nargs=3,default=[19,19,19])
    args=ap.parse_args()
    from phono3py import load
    p=args.input.resolve()
    ph3=load(str(p/'phono3py_disp.yaml'),fc2_filename=str(p/'fc2.hdf5'),
             fc3_filename=str(p/'fc3.hdf5'),is_nac=False)
    print(json.dumps(evaluate(ph3,args.output,args.temperatures,args.mesh)))


if __name__=='__main__':main()
