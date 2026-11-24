"""Remove matching periodic structures within a candidate list or reference set."""
import argparse
import csv
from pathlib import Path


def unique_records(candidates, references=()):
    from pymatgen.analysis.structure_matcher import StructureMatcher
    from pymatgen.core import Structure
    matcher=StructureMatcher(ltol=.2,stol=.3,angle_tol=5,primitive_cell=True,scale=True)
    buckets={};kept=[];duplicates=[]
    for row,is_reference in [(r,True) for r in references]+[(r,False) for r in candidates]:
        structure=Structure.from_file(row['cif_path'])
        key=structure.composition.reduced_formula
        bucket=buckets.setdefault(key,[])
        match=next((mid for mid,old in bucket if matcher.fit(structure,old)),None)
        if match is None:
            bucket.append((row['material_id'],structure))
            if not is_reference:kept.append(row)
        elif not is_reference:
            duplicates.append(dict(material_id=row['material_id'],matches_material_id=match))
    return kept,duplicates


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('candidates',type=Path)
    ap.add_argument('--existing',type=Path,help='Reference CSV with material_id,cif_path')
    ap.add_argument('--output',type=Path,default=Path('unique.csv'))
    ap.add_argument('--duplicates',type=Path,default=Path('duplicates.csv'))
    args=ap.parse_args()
    def load(path):
        if path is None:return []
        with path.open() as f:records=list(csv.DictReader(f))
        for row in records:
            p=Path(row['cif_path']);row['cif_path']=str((path.parent/p).resolve())
        return records
    records=load(args.candidates);kept,duplicates=unique_records(records,load(args.existing))
    fields=list(records[0]) if records else ['material_id','cif_path']
    for path,rows,columns in [(args.output,kept,fields),
            (args.duplicates,duplicates,['material_id','matches_material_id'])]:
        path.parent.mkdir(parents=True,exist_ok=True)
        with path.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=columns);w.writeheader();w.writerows(rows)
    print(f'{len(kept)} retained; {len(duplicates)} structure matches removed')


if __name__=='__main__':main()
