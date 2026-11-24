"""Physical-property and empirical-cost filters after DANA prediction."""
import argparse
import csv
import json
import math
from pathlib import Path

POINT_ORDERS={'1':1,'-1':2,'2':2,'m':2,'2/m':4,'222':4,'mm2':4,'mmm':8,
 '4':4,'-4':4,'4/m':8,'422':8,'4mm':8,'-42m':8,'4/mmm':16,
 '3':3,'-3':6,'32':6,'3m':6,'-3m':12,'6':6,'-6':6,'6/m':12,
 '622':12,'6mm':12,'-6m2':12,'6/mmm':24,'23':12,'m-3':24,'432':24,'-43m':24,'m-3m':48}
CENTROSYMMETRIC={'-1','2/m','mmm','4/m','4/mmm','-3','-3m','6/m','6/mmm','m-3','m-3m'}
CENTRING={'P':1,'A':2,'B':2,'C':2,'I':2,'F':4,'R':3}


def flag(value):
    if value is True or str(value).lower() in ['true','1','yes']:return True
    if value is False or str(value).lower() in ['false','0','no']:return False
    return None


def number(value):
    try:return float(value)
    except (TypeError,ValueError):return None


def cost_proxy(n_primitive,pointgroup,spacegroup_symbol):
    """Frozen central T2 proxy, based on conventional-cell force supercells."""
    n=int(n_primitive);g=POINT_ORDERS[pointgroup];c=CENTRING[spacegroup_symbol[0]]
    n3=8*c*n;n2=64*c*n
    d3=math.ceil(6*n*(1+6*n3)/(g**.9));d2=math.ceil(6*n/(g**.6))
    memory=216*n3**3/2**30
    force=(d3*n3+d2*n2)/(146*64+2*512)
    reciprocal=g if pointgroup in CENTROSYMMETRIC else 2*g
    transport=(n/2)**3*48/reciprocal
    return dict(estimated_FC2_displacements=d2,estimated_FC3_displacements=d3,
        full_FC3_array_GiB=memory,relative_force_work=force,relative_transport_work=transport,
        T2_pass=d3<=2000 and memory<=4 and force<=100 and transport<=1000)


def structure_metadata(path):
    from ase.io import read
    import spglib
    atoms=read(path);cell=(atoms.cell.array,atoms.get_scaled_positions(),atoms.numbers)
    sym=spglib.get_symmetry_dataset(cell,symprec=1e-5)
    primitive=spglib.standardize_cell(cell,to_primitive=True,no_idealize=True,symprec=1e-5)
    return dict(n_primitive=len(primitive[1]),spacegroup=int(sym.number),
                spacegroup_symbol=sym.international,pointgroup=sym.pointgroup,elements=atoms.get_chemical_symbols())


def physical_pass(row,elements):
    from ase.data import chemical_symbols
    source=str(row.get('source_collection',row.get('dataset',''))).lower().replace(' ','').replace('_','')
    n=number(row.get('n_primitive'))
    if n is None or n>20:return False,'primitive_size'
    if source in ['mp','materialsproject']:
        gap=number(row.get('band_gap_eV'));hull=number(row.get('energy_above_hull_eV_atom'))
        ok=(gap is not None and gap>.1 and hull is not None and hull<=.1 and
            flag(row.get('is_metal')) is False and flag(row.get('is_magnetic')) is False)
        return ok,'MP_source_properties' if not ok else ''
    if source=='crystoria':
        risk={'V','Cr','Mn','Fe','Co','Ni'}|set(chemical_symbols[58:71])|set(chemical_symbols[91:103])
        ok=flag(row.get('source_nonmetal')) is True and not(set(elements)&risk)
        return ok,'Crystoria_source_label_or_M1_rule' if not ok else ''
    return False,'unknown_source_properties'


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('predictions',type=Path)
    ap.add_argument('--metadata',type=Path,required=True,help='CSV with material_id, source properties, and optional cached symmetry')
    ap.add_argument('--dana-threshold',type=float,default=8.)
    ap.add_argument('--exclude-ids',type=Path,help='Optional CSV of already registered material_id values')
    ap.add_argument('--output',type=Path,default=Path('selected.csv'))
    ap.add_argument('--all-results',type=Path,default=Path('screening_all.csv'))
    args=ap.parse_args()
    def rows(path):
        with path.open(encoding='utf-8-sig') as f:return list(csv.DictReader(f))
    metadata={r['material_id']:r for r in rows(args.metadata)}
    excluded={r['material_id'] for r in rows(args.exclude_ids)} if args.exclude_ids else set()
    records=[]
    for prediction in rows(args.predictions):
        mid=prediction['material_id'];r=dict(metadata.get(mid,{}),**prediction)
        value=number(prediction.get('kappa_W_mK'));r['dana_pass']=value is not None and 0<value<args.dana_threshold
        if not r['dana_pass']:continue
        path=Path(r['cif_path'])
        if not path.is_absolute():path=args.predictions.parent/path
        meta=structure_metadata(path)
        for k,v in meta.items():
            if k!='elements':r[k]=v
        r['physical_pass'],r['exclusion_reason']=physical_pass(r,meta['elements'])
        r.update(cost_proxy(r['n_primitive'],r['pointgroup'],r['spacegroup_symbol']))
        r['already_registered']=mid in excluded
        r['selected']=r['physical_pass'] and r['T2_pass'] and not r['already_registered']
        r['cif_path']=str(path.resolve());records.append(r)
    def write(path,values):
        fields=list(dict.fromkeys(k for r in records for k in r))
        path.parent.mkdir(parents=True,exist_ok=True)
        with path.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(values)
    write(args.all_results,records);write(args.output,[r for r in records if r['selected']])
    print(json.dumps(dict(dana_retained=len(records),physical_retained=sum(r['physical_pass'] for r in records),
        selected=sum(r['selected'] for r in records),threshold=args.dana_threshold)))


if __name__=='__main__':main()
