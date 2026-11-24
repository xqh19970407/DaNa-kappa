"""Auditable tolerance for small near-Gamma lowest acoustic imaginary modes.

Python 3.7 compatible. Original frequencies/force constants are never clipped.
"""
import hashlib
import itertools
import json
import math
from pathlib import Path
import uuid
import numpy as np


def policy_fingerprint(policy):
    return hashlib.sha256(json.dumps(policy,sort_keys=True,separators=(',',':'),ensure_ascii=True,allow_nan=False).encode()).hexdigest()


def validate_policy(policy):
    if policy.get('schema')!=1 or policy.get('name')!='gamma_acoustic_ratio_v1' or policy.get('lowest_branch_only') is not True:
        raise ValueError('Unknown harmonic screening policy')
    for n in ['strict_gamma_tolerance_THz','strict_other_tolerance_THz','maximum_imaginary_fraction','maximum_imaginary_THz',
              'maximum_q_over_shortest_reciprocal','minimum_translation_weight','maximum_gamma_acoustic_abs_THz',
              'maximum_asr_residual_eV_A2','maximum_asr_relative_residual']:
        v=policy.get(n)
        if isinstance(v,bool) or not isinstance(v,(float,int)) or not math.isfinite(v) or v<=0:
            raise ValueError('Invalid harmonic policy value: '+n)
    if (policy['maximum_imaginary_fraction']>=1 or policy['maximum_q_over_shortest_reciprocal']>=.5
            or policy['minimum_translation_weight']>1
            or policy['strict_gamma_tolerance_THz']!=.05 or policy['strict_other_tolerance_THz']!=.02):
        raise ValueError('Invalid ratio or changed legacy strict gate')


def relative_q_distances(qpoints,reciprocal_lattice,return_images=False):
    """Cartesian nearest reciprocal Gamma distance / shortest reciprocal vector."""
    import spglib
    basis=np.asarray(reciprocal_lattice,dtype=float)
    q=np.asarray(qpoints,dtype=float)
    if basis.shape!=(3,3) or not np.isfinite(basis).all() or abs(np.linalg.det(basis))<1e-20 or not np.isfinite(q).all():
        raise ValueError('Invalid reciprocal geometry')
    reduced=spglib.niggli_reduce(basis,eps=1e-8)
    if reduced is None:raise ValueError('Reciprocal lattice reduction failed')
    shifts=np.array(list(itertools.product(range(-2,3),repeat=3)),dtype=float)
    lengths=np.linalg.norm(shifts@reduced,axis=1)
    shortest=float(lengths[lengths>1e-12].min())
    fractional=q@basis@np.linalg.inv(reduced)
    folded=fractional-np.rint(fractional)
    cartesian=(folded[:,None,:]+shifts[None,:,:])@reduced
    lengths=np.linalg.norm(cartesian,axis=2)
    selected=np.argmin(lengths,axis=1)
    distances=lengths[np.arange(len(q)),selected]/shortest
    images=cartesian[np.arange(len(q)),selected]@np.linalg.inv(basis)
    return (distances,images) if return_images else distances


def _strict_result(f,q,policy):
    if f.ndim!=2 or f.shape[1]<3 or q.shape!=(len(f),3) or not np.isfinite(f).all() or not np.isfinite(q).all():
        return {'pass':False,'reason':'invalid_or_nonfinite_frequencies'},None,None
    if np.max(np.abs(f))<1e-8:return {'pass':False,'reason':'all_zero_frequencies'},None,None
    gamma=np.all(np.isclose(q-np.rint(q),0,atol=1e-8),axis=1)
    tolerance=np.full(f.shape,policy['strict_other_tolerance_THz'])
    for row in np.where(gamma)[0]:tolerance[row,np.argsort(np.abs(f[row]))[:3]]=policy['strict_gamma_tolerance_THz']
    bad=f < -tolerance
    return {'pass':not bool(bad.any()),'reason':'imaginary_modes' if bad.any() else 'passed',
        'minimum_THz':float(f.min()),'maximum_THz':float(f.max()),'negative_modes_beyond_tolerance':int(bad.sum()),
        'gamma_acoustic_tolerance_THz':.05,'other_tolerance_THz':.02},bad,gamma


def screen_frequencies(frequencies,qpoints,reciprocal_lattice,policy,translation_weights=None,asr=None):
    validate_policy(policy)
    f=np.asarray(frequencies,dtype=float);q=np.asarray(qpoints,dtype=float)
    strict,bad,gamma=_strict_result(f,q,policy)
    result=dict(strict,strict_pass=bool(strict['pass']),strict_reason=strict['reason'],acoustic_tolerance_applied=False,
        screening_policy_name=policy['name'],screening_policy_fingerprint=policy_fingerprint(policy),screening_policy=policy,
        interpretation='Conditional numerical tolerance is not proof of dynamical stability',tolerated_modes=[])
    if bad is None:return result
    maximum=float(f.max());ratio=abs(min(float(f.min()),0.))/maximum if maximum>0 else None
    result['imaginary_to_max_ratio']=ratio
    if strict['pass']:return result
    failures=[]
    if ratio is None or ratio>policy['maximum_imaginary_fraction']:failures.append('imaginary_fraction_exceeded')
    if abs(float(f.min()))>policy['maximum_imaginary_THz']:failures.append('absolute_imaginary_exceeded')
    indices=np.argwhere(bad)
    if np.any(indices[:,1]!=0):failures.append('not_only_lowest_branch')
    if not gamma.any() or any(np.max(np.sort(np.abs(f[i]))[:3])>policy['maximum_gamma_acoustic_abs_THz'] for i in np.where(gamma)[0]):
        failures.append('gamma_acoustic_zeros_not_verified')
    try:
        residual=float(asr['maximum_residual_eV_A2']);scale=float(asr['maximum_fc_eV_A2'])
        if (not math.isfinite(residual) or not math.isfinite(scale) or scale<=0 or residual<0
                or residual>policy['maximum_asr_residual_eV_A2'] or residual/scale>policy['maximum_asr_relative_residual']):
            failures.append('asr_residual_exceeded')
        result['asr']=asr
    except (TypeError,KeyError,ValueError):failures.append('asr_not_verified')
    weights=None if translation_weights is None else np.asarray(translation_weights,dtype=float)
    if weights is None or weights.shape!=f.shape:failures.append('acoustic_projection_not_verified')
    else:
        selected=weights[bad]
        if not np.isfinite(selected).all() or np.any(selected<policy['minimum_translation_weight']) or np.any(selected>1+1e-8):
            failures.append('not_rigid_translation_acoustic')
    try:
        relative=relative_q_distances(q[indices[:,0]],reciprocal_lattice)
        if np.any(relative>policy['maximum_q_over_shortest_reciprocal']+1e-12):failures.append('outside_gamma_neighborhood')
    except (ValueError,np.linalg.LinAlgError):relative=None;failures.append('reciprocal_geometry_not_verified')
    result['exception_rejection_reasons']=failures
    if failures:return result
    for k,(i,j) in enumerate(indices):
        result['tolerated_modes'].append({'q_index':int(i),'q_fractional':q[i].tolist(),'branch_zero_based':int(j),
            'frequency_THz':float(f[i,j]),'translation_weight':float(weights[i,j]),'q_over_shortest_reciprocal':float(relative[k])})
    result.update({'pass':True,'reason':'small_near_gamma_acoustic_tolerated','acoustic_tolerance_applied':True})
    return result


def atomic_json(path,value):
    path=Path(path);tmp=path.with_name(path.name+'.tmp.'+uuid.uuid4().hex)
    tmp.write_text(json.dumps(value,indent=2,sort_keys=True,ensure_ascii=True,allow_nan=False)+'\n');tmp.replace(path)


def harmonic_check(root,config,policy):
    """Fit the same FC2 as the original helper; evaluate new policy with eigenvectors."""
    import phono3py
    from phono3py import load
    from phono3py.file_IO import parse_FORCES_FC2
    from phonopy import Phonopy
    root=Path(root);p=load(str(root/'phono3py_disp.yaml'),produce_fc=False,is_nac=False)
    parse_FORCES_FC2(p.phonon_dataset,filename=str(root/'FORCES_FC2'))
    kwargs={'symmetrize_fc2':True,'is_compact_fc':False}
    if int(phono3py.__version__.split('.')[0])>=3:kwargs['fc_calculator']='traditional'
    p.produce_fc2(**kwargs)
    fc=np.asarray(p.fc2)
    if fc.ndim!=4 or fc.shape[0]!=fc.shape[1] or not np.isfinite(fc).all():raise RuntimeError('Invalid full FC2 array')
    asr={'maximum_residual_eV_A2':float(max(np.abs(fc.sum(axis=0)).max(),np.abs(fc.sum(axis=1)).max())),
         'maximum_fc_eV_A2':float(np.abs(fc).max()),'permutation_residual_eV_A2':float(np.abs(fc-fc.transpose(1,0,3,2)).max())}
    h=Phonopy(p.unitcell,p.phonon_supercell_matrix,primitive_matrix=p.primitive_matrix,symprec=1e-5);h.force_constants=fc
    h.run_mesh([19,19,19],is_gamma_center=True,is_mesh_symmetry=False);data=h.get_mesh_dict()
    f=np.asarray(data['frequencies']);q=np.asarray(data['qpoints'])
    if f.shape[1]!=3*config['n_primitive']:raise RuntimeError('Phonon primitive cell differs from submitted structure')
    reciprocal=np.linalg.inv(h.primitive.cell).T
    strict,bad,_=_strict_result(f,q,policy);weights=None
    if bad is not None and bad.any():
        masses=np.asarray(h.primitive.masses,dtype=float)
        if not np.isfinite(masses).all() or np.any(masses<=0):raise RuntimeError('Invalid primitive masses')
        translations=np.zeros((3*len(masses),3))
        for i,mass in enumerate(masses):translations[3*i:3*i+3,:]=np.eye(3)*np.sqrt(mass/masses.sum())
        weights=np.full(f.shape,np.nan)
        for i in np.unique(np.argwhere(bad)[:,0]):
            # Evaluate the nearest-Gamma image: Bloch eigenvector phases at q+G
            # differ by atomic position and would otherwise suppress translations.
            _,images=relative_q_distances(q[i:i+1],reciprocal,return_images=True)
            values,vectors=h.get_frequencies_with_eigenvectors(images[0])
            if np.max(np.abs(values-f[i]))>1e-6:raise RuntimeError('Eigenvector and mesh frequencies do not match')
            weights[i]=np.sum(np.abs(translations.conj().T@vectors)**2,axis=0)
    result=screen_frequencies(f,q,reciprocal,policy,translation_weights=weights,asr=asr)
    result['asr']=asr
    result['acoustic_projection_reference']='Eigenvectors at nearest reciprocal-Gamma image; stored mesh frequencies unchanged'
    result['force_constants_calculator']='traditional';result['negative_frequency_handling']='No clipping/absolute value; unchanged backend cutoff behavior'
    np.savez_compressed(str(root/'harmonic_mesh.npz'),frequencies=f,qpoints=q)
    atomic_json(root/'harmonic_gate.json',result)
    return result


def kappa_harmonic_metadata(root):
    gate=json.loads((Path(root)/'harmonic_gate.json').read_text())
    return {'screening_policy_name':gate.get('screening_policy_name','legacy_strict'),
            'screening_policy_fingerprint':gate.get('screening_policy_fingerprint'),
            'strict_pass':gate.get('strict_pass',gate.get('pass')),
            'acoustic_tolerance_applied':gate.get('acoustic_tolerance_applied',False),
            'minimum_THz':gate.get('minimum_THz'),'maximum_THz':gate.get('maximum_THz'),
            'imaginary_to_max_ratio':gate.get('imaginary_to_max_ratio'),
            'reason':gate.get('reason'),'interpretation':'Conditional tolerance, if applied, is not proof of dynamical stability'}
