"""Portable, resumable VASP finite-displacement workflow.

Uses user-provided VASP and PAW files. No SSH, account, daemon or scheduler
assumptions. Defaults follow the independent-static DFT validation protocol.
"""
import argparse
import json
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys

BASE={'PREC':'Accurate','EDIFF':'1E-8','GGA':'PE','ISPIN':1,'ISMEAR':0,'SIGMA':.01,
      'ALGO':'Normal','NELM':200,'NELMIN':5,'LREAL':'.FALSE.','ADDGRID':'.TRUE.',
      'LASPH':'.TRUE.','LCHARG':'.FALSE.','LWAVE':'.FALSE.','NCORE':8}
STATIC={'IBRION':-1,'NSW':0,'ISIF':2,'ISTART':0,'ICHARG':2,'ISYM':2,'SYMPREC':'1E-5','PSTRESS':0}
RELAX={'IBRION':2,'NSW':10,'ISIF':3,'ISTART':0,'ICHARG':2,'ISYM':2,'SYMPREC':'1E-5',
       'PSTRESS':0,'POTIM':.2,'EDIFFG':-.0001}


def dump(path,value):path.write_text(json.dumps(value,indent=2)+'\n')


def incar(directory,encut,extra,mesh):
    directory.mkdir(parents=True,exist_ok=True)
    values=dict(BASE,ENCUT=f'{encut:.4f}',**extra)
    (directory/'INCAR').write_text(''.join(f'{k} = {v}\n' for k,v in values.items()))
    (directory/'KPOINTS').write_text('Gamma mesh\n0\nGamma\n%d %d %d\n0 0 0\n'%tuple(mesh))


def metrics(directory):
    """Read a finished VASP result and its actual SCF convergence."""
    import numpy as np
    from ase.io import read
    import xml.etree.ElementTree as ET
    xml=directory/'vasprun.xml';tree=ET.parse(xml).getroot();calcs=tree.findall('calculation')
    if not calcs:raise RuntimeError('No VASP calculation found in '+str(xml))
    count=len(calcs[-1].findall('scstep'))
    if count>=200:raise RuntimeError('Electronic convergence was not reached: '+str(directory))
    atoms=read(xml,index=-1)
    return atoms,dict(energy_eV=float(atoms.get_potential_energy()),
        maximum_force_eV_A=float(np.linalg.norm(atoms.get_forces(apply_constraint=False),axis=1).max()),
        maximum_stress_GPa=float(abs(atoms.get_stress()).max()*160.2176634),
        ionic_steps=len(calcs),electronic_iterations=count)


def execute(directory,command):
    if (directory/'vasprun.xml').exists():
        try:return metrics(directory)
        except Exception:pass
    if not (directory/'POTCAR').exists():raise FileNotFoundError('Provide licensed POTCAR files before running VASP.')
    with (directory/'vasp.stdout').open('w') as log:
        subprocess.run(shlex.split(command),cwd=directory,stdout=log,stderr=subprocess.STDOUT,check=True)
    return metrics(directory)


def prepare(args):
    root=args.workdir.resolve();root.mkdir(parents=True,exist_ok=True)
    if (root/'settings.json').exists():return json.loads((root/'settings.json').read_text())
    from ase import Atoms
    from ase.io import read,write
    import spglib
    if not args.cif:raise ValueError('--cif is required for a new work directory.')
    atoms=read(args.cif)
    lattice,positions,numbers=spglib.standardize_cell((atoms.cell.array,atoms.get_scaled_positions(),atoms.numbers),
        to_primitive=True,no_idealize=True,symprec=1e-5)
    atoms=Atoms(numbers=numbers,cell=lattice,scaled_positions=positions,pbc=True)
    order=list(dict.fromkeys(atoms.get_chemical_symbols()))
    symbols=atoms.get_chemical_symbols()
    atoms=atoms[[i for symbol in order for i,value in enumerate(symbols) if value==symbol]]
    paw=json.loads(args.paw_map.read_text()) if args.paw_map else {s:s for s in order}
    chunks=[];enmax=[]
    if args.potcar_root:
        for symbol in order:
            data=(args.potcar_root/paw[symbol]/'POTCAR').read_bytes();chunks.append(data)
            enmax.extend(float(x) for x in re.findall(r'ENMAX\s*=\s*([\d.]+)',data.decode(errors='replace')))
        (root/'POTCAR').write_bytes(b''.join(chunks))
    encut=args.encut or (round(1.3*max(enmax),4) if enmax else None)
    if encut is None:raise ValueError('Use --potcar-root, or --encut for input preparation without PAW files.')
    write(root/'POSCAR.initial',atoms,format='vasp',direct=True,vasp5=True)
    sg=int(spglib.get_symmetry_dataset((atoms.cell.array,atoms.get_scaled_positions(),atoms.numbers),symprec=1e-5).number)
    settings=dict(encut_eV=encut,n_primitive=len(atoms),spacegroup=sg,PAW_labels=[paw[s] for s in order],
                  fc2_dim=[4,4,4],fc3_dim=[2,2,2],displacement_A=.03)
    dump(root/'settings.json',settings)
    return settings


def copy_paw(root,dest):
    if (root/'POTCAR').exists():shutil.copyfile(root/'POTCAR',dest/'POTCAR')


def relax(root,settings,command):
    from ase.io import read,write
    import spglib
    initial=root/'initial_static';incar(initial,settings['encut_eV'],STATIC,[6,6,6])
    shutil.copyfile(root/'POSCAR.initial',initial/'POSCAR');copy_paw(root,initial)
    _,reference=execute(initial,command)
    source=initial/'POSCAR';steps=0
    for segment in range(20):
        work=root/'relax'/('segment_%03d'%segment);incar(work,settings['encut_eV'],RELAX,[6,6,6])
        shutil.copyfile(source,work/'POSCAR');copy_paw(root,work)
        atoms,info=execute(work,command);steps+=info['ionic_steps']
        symmetry=spglib.get_symmetry_dataset((atoms.cell.array,atoms.get_scaled_positions(),atoms.numbers),symprec=1e-5)
        if int(symmetry.number)!=settings['spacegroup']:raise RuntimeError('Space group changed during DFT relaxation.')
        source=work/'CONTCAR'
        converged='reached required accuracy' in (work/'OUTCAR').read_text(errors='replace')
        if converged:
            final=root/'final_static';incar(final,settings['encut_eV'],STATIC,[6,6,6])
            shutil.copyfile(source,final/'POSCAR');copy_paw(root,final)
            atoms,check=execute(final,command)
            if check['maximum_force_eV_A']>.001:raise RuntimeError('Fresh static force exceeds 0.001 eV/A.')
            if (check['energy_eV']-reference['energy_eV'])/settings['n_primitive']>.01:
                raise RuntimeError('Fresh static energy increased by more than 0.01 eV/atom.')
            if check['maximum_stress_GPa']<=.05:
                dump(root/'geometry.json',check);write(root/'structure_final.cif',atoms)
                shutil.copyfile(final/'POSCAR',root/'POSCAR');return
            # A failed stress check needs another fresh calculation next segment.
            (final/'vasprun.xml').unlink()
        if steps>=200:break
    raise RuntimeError('No accepted zero-pressure geometry within 200 ionic steps.')


def displacements(root,settings,phono):
    subprocess.run([phono,'-d','-c','POSCAR','--dim=2 2 2','--dim-fc2=4 4 4','--amplitude=0.03'],cwd=root,check=True)
    for group,pattern,mesh in [('fc2','POSCAR_FC2-[0-9]*',[2,2,2]),('fc3','POSCAR-[0-9]*',[3,3,3])]:
        for poscar in sorted(root.glob(pattern)):
            dest=root/'tasks'/group/('disp-'+poscar.name.rsplit('-',1)[1])
            incar(dest,settings['encut_eV'],{'IBRION':-1,'NSW':0,'ISYM':0},mesh)
            shutil.copyfile(poscar,dest/'POSCAR');copy_paw(root,dest)


def force_group(root,group,settings,command,phono,policy):
    for task in sorted((root/'tasks'/group).glob('disp-*')):execute(task,command)
    xmls=[str(p.resolve()) for p in sorted((root/'tasks'/group).glob('disp-*/vasprun.xml'))]
    if not xmls:raise RuntimeError('No displaced force calculations for '+group)
    subprocess.run([phono,'--cf'+group[-1]]+xmls,cwd=root,check=True)
    if group=='fc2':
        from .harmonic_policy import harmonic_check
        result=harmonic_check(root,settings,json.loads(policy.read_text()))
        if not result['pass']:raise RuntimeError('Imaginary phonons: '+result['reason'])


def collect(root,phono):
    command=[phono,'-c','phono3py_disp.yaml','--sym-fc2','--sym-fc3r','--mesh','19','19','19',
             '--br','--wigner','--isotope','--boundary-mfp','1000000','--ts','300','--loglevel','1']
    poscar=root/'POSCAR';hidden=root/'POSCAR.transport_input'
    if poscar.exists():poscar.rename(hidden)
    try:subprocess.run(command,cwd=root,check=True)
    finally:
        if hidden.exists():hidden.rename(poscar)
    import h5py
    with h5py.File(root/'kappa-m191919.hdf5') as f:
        result={name:float(f[key][0,:3].mean()) for name,key in [('kappa_W_mK','kappa_TOT_RTA'),
                  ('kappa_P_W_mK','kappa_P_RTA'),('kappa_C_W_mK','kappa_C')]}
    result.update(temperature_K=300,target='Wigner total RTA trace/3',convergence_validated=False)
    dump(root/'kappa_result.json',result);print(json.dumps(result))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--phase',choices=['prepare','relax','displacements','fc2','fc3','collect','all'],default='prepare')
    ap.add_argument('--cif',type=Path);ap.add_argument('--workdir',type=Path,required=True)
    ap.add_argument('--potcar-root',type=Path);ap.add_argument('--paw-map',type=Path);ap.add_argument('--encut',type=float)
    ap.add_argument('--vasp-command',default='mpiexec -np 64 vasp_std')
    ap.add_argument('--phono3py',default=str(Path(sys.executable).with_name('phono3py')))
    ap.add_argument('--harmonic-policy',type=Path,default=Path(__file__).resolve().parent/'defaults/harmonic_policy.json')
    args=ap.parse_args();root=args.workdir.resolve();settings=prepare(args)
    for phase in ['relax','displacements','fc2','fc3','collect']:
        if args.phase not in [phase,'all']:continue
        if phase=='relax':relax(root,settings,args.vasp_command)
        elif phase=='displacements':displacements(root,settings,args.phono3py)
        elif phase in ['fc2','fc3']:force_group(root,phase,settings,args.vasp_command,args.phono3py,args.harmonic_policy)
        else:collect(root,args.phono3py)
    print('Finished',args.phase,root)


if __name__=='__main__':main()
