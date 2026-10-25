"""Read prepared inputs and verify compatibility offline; no legacy atlas builder."""
from pathlib import Path
import gzip
import hashlib
import json
import numpy as np


def read_json(path):
    path=Path(path)
    return json.loads(gzip.decompress(path.read_bytes()) if path.suffix=='.gz' else path.read_bytes())


def load_inputs(folder):
    from xenium_runtime.io import fingerprint
    import anndata as ad
    from shapely.geometry import shape
    folder=Path(folder)
    manifest=read_json(folder/'BUNDLE_MANIFEST.json')
    for name,expected in manifest['files'].items():
        p=folder/name
        if not p.is_file() or p.stat().st_size!=expected['bytes']:
            raise ValueError(f'Missing or changed bundle file: {name}')
        with p.open('rb') as f:
            if hashlib.file_digest(f,'sha256').hexdigest()!=expected['sha256']:
                raise ValueError(f'Checksum mismatch: {name}')
    atlas=read_json(folder/'atlas.json.gz');contract=read_json(folder/'runtime_contract.json')
    if fingerprint({k:v for k,v in atlas.items() if k!='fingerprint'})!=atlas['fingerprint']:
        raise ValueError('Atlas fingerprint mismatch')
    if fingerprint(contract)!=atlas['provenance']['runtime_contract_fingerprint']:
        raise ValueError('Runtime contract mismatch')
    data=ad.read_h5ad(folder/'region.h5ad')
    cells=atlas['cells'];ids=[c['cell_id'] for c in cells]
    if len(cells)>1000 or len(set(ids))!=len(ids): raise ValueError('Invalid population or duplicate identities')
    if list(data.obs_names)!=[c['external_cell_id'] for c in cells]: raise ValueError('Expression/geometry ID order mismatch')
    if data.obs.internal_id.tolist()!=ids: raise ValueError('Internal ID mapping mismatch')
    if not np.allclose(data.obsm['spatial'],[[c['x'],c['y']] for c in cells]): raise ValueError('Spatial coordinate mismatch')
    known=set(ids)
    for c in cells:
        geom=shape(c['geometry'])
        if not geom.is_valid or geom.is_empty or geom.area<=0: raise ValueError('Invalid geometry')
        if not np.isclose(geom.area,c['area_um2']): raise ValueError('Geometry area mismatch')
        if any(e['cell_id'] not in known for e in c['contacts']): raise ValueError('Dangling contact')
        if c['role']=='roi' and c['identity'] is None: raise ValueError('Active cell has no supported identity')
        for s in c['signatures'].values():
            if s['method']!='AUCell' or not 0<=s['raw']<=1 or not 0<=s['percentile_within_type_depth']<=1:
                raise ValueError('Invalid AUCell evidence')
    return atlas,contract,data


def initialize(folder,scenario='static',seed=21):
    from xenium_runtime.v2.autonomous import AutonomousSimulation
    from xenium_runtime.v2.policy import ScoreAwareLocalPolicy
    atlas,contract,_=load_inputs(folder)
    return AutonomousSimulation(atlas,contract,scenario=scenario,seed=seed,policy=ScoreAwareLocalPolicy())


def validate_bundle(folder,smoke_steps=1):
    """A static ten-minute tick verifies numerical compatibility, not biology."""
    from xenium_runtime.io import write_json
    sim=initialize(folder,scenario='static')
    t=sim.tissue
    before=t.xy.copy()
    for _ in range(smoke_steps): sim.step()
    t.check_invariants();sim.check_biology()
    if sim.decisions: raise AssertionError('Static validation must make no agent decisions')
    if not np.allclose(before,t.xy,rtol=0,atol=1e-8): raise AssertionError('Static tissue moved')
    if not np.isfinite(t.xy).all() or not t.present.all(): raise AssertionError('Invalid initial state')
    if (t.acts != np.array([c['role']=='roi' for c in sim.atlas['cells']])).any():
        raise AssertionError('Context became an actor')
    result=dict(passed=True,total_cells=len(t.ids),acting_cells=int(t.acts.sum()),context_cells=int((~t.acts).sum()),
                static_minutes=sim.time,decisions=len(sim.decisions),http_calls=0,
                fidelity=sim.fidelity,scope='Format, ID joins, geometry, evidence ranges, initial mechanics and one static tick; not biological validation.')
    write_json(Path(folder)/'compatibility_report.json',result)
    return result
