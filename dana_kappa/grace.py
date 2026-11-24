"""CIF -> two-stage GRACE relaxation -> FC2/FC3 -> Wigner conductivity.

Relaxation workflow adapted from MPA2suite/k_SRME (GPL-3.0-only).
GRACE/TensorPotential and pretrained models remain external dependencies.
"""
import argparse
import json
import os
from pathlib import Path


def relaxed(atoms,calculator,fmax=1e-4,steps=300):
    from ase.constraints import FixSymmetry
    from ase.filters import FrechetCellFilter
    from ase.optimize import FIRE
    import spglib
    def sg(a):return int(spglib.get_symmetry_dataset((a.cell.array,a.get_scaled_positions(),a.numbers),symprec=1e-5).number)
    original=sg(atoms);atoms.calc=calculator;atoms.set_constraint(FixSymmetry(atoms))
    cell=FrechetCellFilter(atoms,mask=[True,True,True,False,False,False])
    first=FIRE(cell,logfile='relax_stage1.log');ok1=first.run(fmax=fmax,steps=steps)
    one=atoms.copy();one.calc=calculator
    atoms.constraints=None
    second=FIRE(cell,logfile='relax_stage2.log');ok2=second.run(fmax=fmax,steps=steps)
    chosen,ok=(atoms,ok2) if sg(atoms)==sg(one) else (one,ok1)
    if not ok or sg(chosen)!=original:raise RuntimeError('GRACE relaxation did not converge with the retained space group.')
    chosen.constraints=None
    return chosen


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('cif',type=Path);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--model',default='GRACE_2L_OAM_28Jan25')
    ap.add_argument('--threads',type=int,default=4);ap.add_argument('--device',choices=['cpu','gpu'],default='cpu')
    ap.add_argument('--fmax',type=float,default=1e-4);ap.add_argument('--steps',type=int,default=300)
    ap.add_argument('--mesh',type=int,nargs=3,default=[19,19,19])
    ap.add_argument('--temperature',type=float,default=300)
    args=ap.parse_args();source=args.cif.resolve();output=args.output.resolve();output.mkdir(parents=True,exist_ok=True)
    if args.device=='cpu':os.environ['CUDA_VISIBLE_DEVICES']='-1'
    os.environ['OMP_NUM_THREADS']=str(args.threads)
    import numpy as np
    from ase import Atoms
    from ase.io import read,write
    import spglib,tensorflow as tf
    tf.config.threading.set_intra_op_parallelism_threads(args.threads)
    tf.config.threading.set_inter_op_parallelism_threads(1)
    from tensorpotential.calculator import grace_fm
    from phono3py import Phono3py
    from phonopy.structure.atoms import PhonopyAtoms
    from phonopy.file_IO import write_force_constants_to_hdf5
    from phono3py.file_IO import write_fc3_to_hdf5
    from .transport import evaluate
    atoms=read(source)
    lattice,positions,numbers=spglib.standardize_cell((atoms.cell.array,atoms.get_scaled_positions(),atoms.numbers),
        to_primitive=False,no_idealize=True,symprec=1e-5)
    atoms=Atoms(numbers=numbers,cell=lattice,scaled_positions=positions,pbc=True)
    calc=grace_fm(args.model);old=Path.cwd()
    try:
        os.chdir(output)
        if (output/'relaxed.extxyz').exists():atoms=read(output/'relaxed.extxyz')
        else:
            atoms=relaxed(atoms,calc,args.fmax,args.steps);atoms.calc=None
            write(output/'relaxed.extxyz',atoms)
        p=Phono3py(PhonopyAtoms(numbers=atoms.numbers,cell=atoms.cell.array,scaled_positions=atoms.get_scaled_positions()),
            supercell_matrix=np.eye(3,dtype=int)*2,phonon_supercell_matrix=np.eye(3,dtype=int)*4,
            primitive_matrix='auto',symprec=1e-5,calculator='vasp',log_level=0)
        p.mesh_numbers=args.mesh;p.generate_displacements(distance=.03)
        p.save('phono3py_disp.yaml',settings={'force_sets':False,'force_constants':False,'displacements':True})
        def forces(cells,name):
            path=output/(name+'.npy')
            if path.exists():return np.load(path)
            result=[]
            for i,cell in enumerate(cells):
                a=Atoms(numbers=cell.numbers,cell=cell.cell,scaled_positions=cell.scaled_positions,pbc=True)
                a.calc=calc;result.append(a.get_forces())
                if (i+1)%20==0:print(name,i+1,len(cells),flush=True)
            result=np.asarray(result);np.save(path,result);return result
        p.phonon_forces=forces(p.phonon_supercells_with_displacements,'forces_fc2')
        p.produce_fc2(symmetrize_fc2=True)
        write_force_constants_to_hdf5(p.fc2[p.phonon_primitive.p2s_map],filename='fc2.hdf5',p2s_map=p.phonon_primitive.p2s_map)
        p.init_phph_interaction(symmetrize_fc3q=False);p.run_phonon_solver()
        frequency=p.get_phonon_data()[0]
        if np.any(frequency[0,:3]<-.01) or np.any(frequency[0,3:]<0) or np.any(frequency[1:]<0):
            (output/'result.json').write_text(json.dumps({'status':'imaginary','minimum_frequency_THz':float(frequency.min())})+'\n')
            print('Imaginary phonons: FC3 and conductivity were skipped.');return
        p.forces=forces(p.supercells_with_displacements,'forces_fc3');p.produce_fc3(symmetrize_fc3r=True)
        write_fc3_to_hdf5(p.fc3[p.primitive.p2s_map],filename='fc3.hdf5',p2s_map=p.primitive.p2s_map)
        print(json.dumps(evaluate(p,output,[args.temperature],args.mesh)))
    finally:os.chdir(old)


if __name__=='__main__':main()
