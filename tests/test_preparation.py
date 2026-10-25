"""Scientific/regression checks independent of Colab UI and the full source data."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import unittest
import numpy as np
import pandas as pd
import scipy.sparse as sp
import shapely
import anndata as ad
import decoupler as dc
import xenium_prepare as p


class PreparationTests(unittest.TestCase):
    def test_aucell_known_enrichment_and_rank_invariance(self):
        # Strict descending, untied counts: first five genes are enriched maximally;
        # last five are below the cutoff. This checks direction and normalization.
        genes=np.array([f'g{i:02d}' for i in range(20)])
        counts=np.arange(20,0,-1,dtype=float)[None,:]
        net=pd.DataFrame({'source':['top']*5+['bottom']*5,'target':genes[:5].tolist()+genes[-5:].tolist()})
        results=[]
        for x in [counts, counts*731, np.log1p(counts)]:
            a=ad.AnnData(sp.csr_matrix(x),obs=pd.DataFrame(index=['cell']),var=pd.DataFrame(index=genes))
            dc.mt.aucell(a,net,tmin=5,n_up=5,empty=False)
            results.append(a.obsm['score_aucell'].iloc[0])
        for result in results:
            self.assertAlmostEqual(result['top'],1.0)
            self.assertAlmostEqual(result['bottom'],0.0)

    def test_aucell_ties_do_not_depend_on_batch_boundaries(self):
        genes=np.array([f'g{i:02d}' for i in range(20)])
        x=np.array([[5]*10+[1]*10,[2]*5+[4]*15],dtype=float)
        net=pd.DataFrame({'source':['set']*5,'target':genes[::4]})
        def score(mat,names):
            a=ad.AnnData(sp.csr_matrix(mat),obs=pd.DataFrame(index=names),var=pd.DataFrame(index=genes))
            dc.mt.aucell(a,net,tmin=5,n_up=5,empty=False)
            return a.obsm['score_aucell'].to_numpy()
        whole=score(x,['a','b'])
        separate=np.vstack([score(x[i:i+1],[name]) for i,name in enumerate(['a','b'])])
        np.testing.assert_array_equal(whole,separate)

    def test_region_contains_all_nearby_objects_and_respects_total_cap(self):
        xx,yy=np.meshgrid(np.arange(50)*12,np.arange(50)*12)
        xy=np.c_[xx.ravel(),yy.ravel()]
        shapes=shapely.buffer(shapely.points(xy),5)
        sample={'obs':pd.DataFrame(xy,columns=['x_centroid','y_centroid']),
                'polygons':shapes,'geometry_usable':np.ones(len(xy),bool),'tree':shapely.STRtree(shapes)}
        region=p.make_region(sample,[300,300],core_target=750,max_total=1000,halo_um=30)
        self.assertIsNotNone(region);self.assertLessEqual(len(region['selected']),1000)
        distance=np.linalg.norm(xy-[300,300],axis=1)
        self.assertEqual(set(region['core']),set(np.flatnonzero(distance<=region['radius_um']+1e-9)))
        # Independent distance check, including cells whose centroid is beyond the halo.
        union=shapely.union_all(shapes[region['core']])
        required=set(np.flatnonzero(shapely.distance(shapes,union)<29.9))
        self.assertTrue(required<=set(region['selected']))
        with self.assertRaises(ValueError): p.make_region(sample,[300,300],max_total=1001)
        with self.assertRaises(ValueError): p.make_region(sample,[300,300],halo_um=float('nan'))

    def test_missing_cycle_score_cannot_license_cycle_in_local_score_policy(self):
        from xenium_runtime.v2.policy import ScoreAwareLocalPolicy
        options={'ENTER_CYCLE':{},'REMAIN_QUIESCENT':{}}
        obs=dict(primary_label='Tumor',sim_state={'exhaustion':0,'pulses':0,'pdl1':0},fields_sim={'oxygen':1,'ifng':0},time_min=0)
        request=dict(cell_id=1,options=options,observation=obs)
        result=ScoreAwareLocalPolicy().choose([request])
        self.assertEqual(result[1]['choice'],'REMAIN_QUIESCENT')
        obs['t0_aucell_percentile']={'MITOTIC_CELL_CYCLE':.9}
        obs['t0_aucell_raw']={'MITOTIC_CELL_CYCLE':0}
        self.assertEqual(ScoreAwareLocalPolicy().choose([request])[1]['choice'],'REMAIN_QUIESCENT')

    def test_stale_annotations_refuse_export(self):
        with self.assertRaisesRegex(ValueError,'Annotations changed'):
            p.export_region({'annotation_fingerprint':'new','score_annotation_fingerprint':'old'},'R01')

if __name__=='__main__': unittest.main()
