"""CIF -> atom/line graphs -> frozen DANA -> 300 K lattice conductivity."""
import argparse
import csv
import json
import os
from pathlib import Path

from .upstream import configure


def make_graph(path):
    import torch
    from ase.io import read
    from jarvis.core.atoms import Atoms
    from jarvis.core.graphs import Graph
    atoms=read(path)
    crystal=Atoms(lattice_mat=atoms.cell.array.tolist(),coords=atoms.get_scaled_positions().tolist(),
                  elements=atoms.get_chemical_symbols(),cartesian=False)
    graph,line=Graph.atom_dgl_multigraph(crystal,neighbor_strategy='k-nearest',cutoff=8.,max_neighbors=12,
                                        atom_features='cgcnn',compute_line_graph=True,use_canonize=True)
    # Support stock jarvis-tools without editing its installed source.
    graph.ndata['coords']=torch.tensor(crystal.cart_coords,dtype=torch.get_default_dtype())
    graph.ndata['frac_coords']=torch.tensor(crystal.frac_coords,dtype=torch.get_default_dtype())
    return graph,line


def load_model(checkpoint,config,upstream=None,device='cpu'):
    configure(upstream)
    import torch
    from .model import create_model
    model=create_model(json.loads(Path(config).read_text()))
    state=torch.load(checkpoint,map_location='cpu',weights_only=True)
    model.load_state_dict(state.get('model',state),strict=True)
    return model.to(device).eval()


def predict_records(model,records,device='cpu',batch_size=32):
    import dgl
    import torch
    result=[]
    for start in range(0,len(records),batch_size):
        group=records[start:start+batch_size]
        graphs=[make_graph(row['cif_path']) for row in group]
        g=dgl.batch([x[0] for x in graphs]).to(device)
        line=dgl.batch([x[1] for x in graphs]).to(device)
        with torch.inference_mode():pred=model([g,line]).reshape(-1).cpu()
        for row,value in zip(group,pred.tolist()):
            result.append(dict(material_id=row['material_id'],cif_path=str(row['cif_path']),temperature_K=300,
                               ln_kappa=value,kappa_W_mK=float(torch.exp(torch.tensor(value,dtype=torch.float64)))))
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('structures',nargs='*',type=Path)
    ap.add_argument('--manifest',type=Path,help='CSV with material_id,cif_path; paths relative to CSV directory')
    ap.add_argument('--checkpoint',required=True,type=Path)
    ap.add_argument('--config',type=Path,default=Path(__file__).resolve().parent/'defaults/dana.json')
    ap.add_argument('--upstream',type=Path)
    ap.add_argument('--device',default='cpu');ap.add_argument('--batch-size',type=int,default=32)
    ap.add_argument('--threads',type=int,default=4);ap.add_argument('--output',type=Path,default=Path('predictions.csv'))
    args=ap.parse_args()
    import torch
    torch.set_num_threads(args.threads)
    if args.manifest:
        with args.manifest.open(encoding='utf-8-sig') as f:records=list(csv.DictReader(f))
        records=[dict(material_id=r['material_id'],cif_path=args.manifest.parent/r['cif_path']) for r in records]
    else:records=[dict(material_id=p.stem,cif_path=p) for p in args.structures]
    if not records:ap.error('Provide CIF paths or --manifest.')
    model=load_model(args.checkpoint,args.config,args.upstream,args.device)
    values=predict_records(model,records,args.device,args.batch_size)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    for row in values:
        row['cif_path']=os.path.relpath(Path(row['cif_path']).resolve(),args.output.resolve().parent)
    with args.output.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(values[0]));w.writeheader();w.writerows(values)
    print(json.dumps(dict(records=len(values),output=str(args.output),device=args.device)))


if __name__=='__main__':main()
