"""Create an explicitly synthetic heterogeneous v3 experiment; no API calls."""
import argparse
import random
from pathlib import Path
from virtualtissue_xenium.v2.inputs import prepare, write_json
from virtualtissue_xenium.v2.controls import apply_control
from virtualtissue_xenium.v3.ontology import load_signatures


def create_demo(output):
    out=Path(output)
    if out.exists():
        raise ValueError('Choose a new demo directory; existing work is never overwritten.')
    genes,_=load_signatures('human')
    cycle=genes['cell_cycle'][:8]
    housekeeping=[f'SYNTHETIC_HK_{i:02}' for i in range(80)]
    panel=sorted({g for values in genes.values() for g in values}|set(housekeeping))
    rng=random.Random(20261031)
    cells=[]
    for i in range(90):
        for j,lineage in enumerate(('malignant','T_CD8','myeloid')):
            depth=(1,3,10)[i%3]
            expression={g:float(rng.randint(1,4)*depth) for g in rng.sample(housekeeping,25)}
            if lineage=='malignant' and i%2==0:
                expression.update({g:float(rng.randint(5,12)*depth) for g in cycle})
            cells.append(dict(cell_id=f'{lineage}_{i:03}',cell_type=lineage,lineage=lineage,sample_id='synthetic',
                x=float((i%15)*20+j*4),y=float((i//15)*20),expression=expression))
    context=dict(schema_version='3.0',ontology_version='3.0',platform='xenium',assay='xenium',
        tissue='synthetic melanoma-like geometry',organism='human',disease='synthetic engineering example',
        sampling=dict(design='subsampled',description='270 generated cells; artificial counts/XY, not a biological sample.',missing_populations=[]),
        perturbations=[],recognition=False,time_step_hours=.5,
        scoring=dict(method='depth_conditioned_auc',depth_bins=2,min_reference_cells=30,min_transcripts=20),
        dynamics=dict(crowding_enabled=False))
    bundle=prepare(cells,context,assayed_genes=panel)
    out.mkdir(parents=True)
    write_json(out/'input.json',dict(cells=cells,assayed_genes=panel,description='Synthetic engineering fixture.'))
    write_json(out/'context.json',context)
    write_json(out/'control.json',bundle)
    for kind in ('beta_zero','shuffle_expression_within_type'):
        write_json(out/(kind+'.json'),apply_control(bundle,kind,seed=11))
    write_json(out/'arms.json',{'arms':{'control':'control.json','beta_zero':'beta_zero.json',
                                    'within_type_shuffle':'shuffle_expression_within_type.json'}})
    write_json(out/'bounds.json',{'action_priors.malignant.PROLIFERATE':[.003,.009],
                                'scoring.beta.PROLIFERATE':[0.,1.5]})
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',required=True)
    print(create_demo(p.parse_args().output))
