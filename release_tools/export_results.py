#!/usr/bin/env python3
"""Export portable CIF/label pairs from the project's archived result tables."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil


def read_csv(p):
    with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def write_table(p,rows):
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--dataset',type=Path,required=True,help='Root of the archived crystal-kappa pair dataset')
    ap.add_argument('--dft-registry',type=Path,required=True)
    ap.add_argument('--dft-parameters',type=Path,required=True)
    ap.add_argument('--dana-predictions',type=Path,help='Optional CSV with material_id,DANA_kappa_W_mK')
    ap.add_argument('--output',type=Path,default=Path('data'))
    args=ap.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=True)
    def copy(src,rel):
        dst=out/rel;dst.parent.mkdir(parents=True,exist_ok=True)
        if not dst.exists():shutil.copyfile(src,dst)
        return str(rel),sha(dst)
    grace=[];by_id={}
    for row in read_csv(args.dataset/'02_grace/labels.csv'):
        mid=row['material_id'];value=float(row['kappa_W_mK']);included=0<value<2
        cif=calc=source='';digest=''
        if included:
            cif,digest=copy(args.dataset/row['cif_path'],Path('structures/grace_low2')/(mid+'.cif'))
            calc,_=copy(args.dataset/row['calculation_cell_cif_path'],Path('calculation_cells/grace_low2')/(mid+'.cif'))
            source,_=copy(args.dataset/row['database_input_cif_path'],Path('source_cifs')/(mid+'.cif'))
        item=dict(material_id=mid,formula=row['formula'],spacegroup=int(row['spacegroup_number']),
            source_collection=row['source_collection'],aliases=row['aliases'],temperature_K=300,
            kappa_W_mK=value,kappa_P_W_mK=float(row['kappa_particle_W_mK']),
            kappa_C_W_mK=float(row['kappa_coherence_W_mK']),
            cif_path=cif,calculation_cell_cif_path=calc,source_cif_path=source,cif_sha256=digest,
            n_primitive=int(row['n_atoms_model_input']),n_calculation_cell=int(row['n_atoms_calculation_cell']),
            fc2_supercell_atoms=int(row['fc2_supercell_atoms']),fc3_supercell_atoms=int(row['fc3_supercell_atoms']),
            model='GRACE_2L_OAM_28Jan25',target='Wigner_total_RTA_trace_over_3',NAC=False,SOC=False,
            isotope=True,boundary_mfp_um=1000000,train_eligible=row['train_eligible']=='True',
            quality_flags=row['quality_notes'],convergence_validated=False,
            low2=included,cif_included=included)
        grace.append(item);by_id[mid]=row
        for alias in json.loads(row['aliases']):by_id[alias]=row
    write_table(out/'grace_results.csv',grace)
    low=[r for r in grace if r['low2']];write_table(out/'grace_low2.csv',low)
    params={r['formula']:r for r in read_csv(args.dft_parameters)}
    dft=[]
    for row in json.loads(args.dft_registry.read_text()):
        ids=row['material_ids'];mid=ids[0]
        public_id='dft_'+row['formula']+'_sg'+str(row['reported_spacegroup'])+'_'+row['geometry_sha256'][:8]
        cif,digest=copy(row['structure_path'],Path('structures/dft')/(public_id+'.cif'))
        src=''
        g=by_id.get(mid)
        if g:
            src,_=copy(args.dataset/g['database_input_cif_path'],Path('source_cifs')/(mid+'.cif'))
        p=params.get(row['formula'],{})
        item=dict(structure_id=public_id,material_id=mid,aliases=json.dumps(ids),formula=row['formula'],
            spacegroup=int(row['reported_spacegroup']),n_primitive=int(row['n_primitive']),temperature_K=300,
            kappa_W_mK=float(row['kappa_W_mK']),kappa_P_W_mK='',kappa_C_W_mK='',
            GRACE_kappa_W_mK=row['grace_kappa_W_mK'],
            cif_path=cif,source_cif_path=src,cif_sha256=digest,
            pressure_GPa=row.get('pressure_GPa'),
            phonon_mesh_strict_pass=row['strict_harmonic_pass'],conditional_acoustic_acceptance=row['acoustic_tolerance_applied'],
            dense_path_min_frequency_THz=p.get('dense_path_min_THz',''),
            workflow='independent_static_v3' if row['workflow_policy']=='zero_pressure_symmetry_static_v3' else 'ionic_stress_v2',
            ENCUT_eV=p.get('ENCUT_eV',''),PAW_labels=p.get('PAW_labels_in_POSCAR_order',''),
            active_allocation_hours=p.get('active_allocation_hours',''),
            target='Wigner_total_RTA_trace_over_3',NAC=False,SOC=False,isotope=True,boundary_mfp_um=1000000,
            convergence_validated=False,low2=float(row['kappa_W_mK'])<2)
        import h5py
        with h5py.File(row['hdf5_path']) as result:
            item['kappa_P_W_mK']=float(result['kappa_P_RTA'][0,:3].mean())
            item['kappa_C_W_mK']=float(result['kappa_C'][0,:3].mean())
        dft.append(item)
    write_table(out/'dft_results.csv',dft);write_table(out/'dft_low2.csv',[r for r in dft if r['low2']])
    if args.dana_predictions:
        index={r['material_id']:float(r['DANA_kappa_W_mK']) for r in read_csv(args.dana_predictions)}
        for group in [grace,dft]:
            for row in group:row['DANA_source_kappa_W_mK']=index.get(row['material_id'])
        for name,group in [('grace_results',grace),('grace_low2',low),('dft_results',dft),('dft_low2',[r for r in dft if r['low2']])]:
            write_table(out/(name+'.csv'),group)
    for name,rows in [('grace_results',grace),('dft_results',dft)]:
        with (out/(name+'.jsonl')).open('w') as f:
            for row in rows:
                item=dict(row)
                for key in ['aliases','quality_flags']:
                    if key in item:item[key]=json.loads(item[key])
                f.write(json.dumps(item,ensure_ascii=False,allow_nan=False)+'\n')
    summary=dict(temperature_K=300,unit='W m^-1 K^-1',cif_paths_relative_to='data/',
        grace_successful_records=len(grace),grace_positive_below2=len(low),
        grace_zero_records=sum(r['kappa_W_mK']==0 for r in grace),
        dft_completed_records=len(dft),dft_positive_below2=sum(r['low2'] for r in dft),
        dft_primary_protocol=sum(r['workflow']=='independent_static_v3' for r in dft),
        note='Source records are not a globally deduplicated set of compounds. Low2 uses strict 0 < kappa < 2.')
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary))


if __name__=='__main__':main()
