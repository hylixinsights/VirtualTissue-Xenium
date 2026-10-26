"""Behavioral regressions for the independent 0.1.1 review findings."""
import copy
import gzip
import json
import os
import tempfile
import unittest
import warnings
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

import numpy as np
import pandas as pd
import virtualtissue_xenium as vtx
from virtualtissue_xenium import results, own_data
from jev_payload import build_batches, payload
from player import FIELDS, refresh_payload, write_payload
from recorded_simulation import RecordedSimulation
from recording_analysis import event_totals, population_limit
from xenium_runtime.v2.autonomous import AutonomousSimulation


def fixture(end=40, scenario="declared_recognition", seed=21):
    fields = {key: [0.0]*9 for key in FIELDS}
    frames = [dict(time_min=t, ids=[2] if t >= 20 else [1, 2],
                   x=[15] if t >= 20 else [10, 15], y=[10] if t >= 20 else [10, 10],
                   rho=[2] if t >= 20 else [2, 2],
                   states=[dict(viability="alive", cycle=None)]*(1 if t >= 20 else 2), fields=fields)
              for t in range(0, end+1, 20)]
    contract = dict(autonomous_scenarios={scenario: dict(title=scenario, recognition=scenario=="declared_recognition",
                    effector_competence=True, clearance_competence=True, fields=True, nerve="none")},
                    biology=dict(max_population=3), division=dict(division_blocked_local_pressure=.4),
                    run=dict(horizon_min=60), provider=dict(pinned_model="jev-1.14.0"))
    episode = dict(contract=contract, scenario=scenario, seed=seed, provider="local_score_rules",
                   atlas_fingerprint="atlas", engine_fingerprint="kernel", prompt_version="vtx-choice-2",
                   reporting_version="vtx-recording-2", patch_id="sample_R01", sample_name="Example sample",
                   bodies=[dict(id=1,label="Tumor",origin="roi",born_at_min=0),
                           dict(id=2,label="T_CD8_cytotoxic",origin="roi",born_at_min=0)],
                   frames=frames, field_grid=dict(origin_um=[0,0],spacing_um=10,shape=[3,3]),
                   events=[dict(time_min=20,kind="hypoxic_necrosis",cell_id=1,position_um=[10,10])], decisions=[],
                   summary=dict(time_min=end,events={"hypoxic_necrosis":1},decisions=0,final={}))
    original=dict(origin=[0,0],frames=[dict(t=f["time_min"]) for f in frames],cells={})
    return dict(folder=Path(f"run_{scenario}"),episode=episode,config=dict(zip_sha256="bundle",sources={},environment={}),
                report=dict(mode="local_preview",completed=False,usage={})),original


class RecordingTests(unittest.TestCase):
    def test_necrosis_survives_removed_frame_and_counts_once(self):
        result, original = fixture()
        p = refresh_payload(original,result["episode"])
        e = p["events"][0]
        self.assertEqual(e[4:7],[10,10,False])
        self.assertNotIn(1,result["episode"]["frames"][-1]["ids"])
        self.assertEqual(results.summary(result)["hypoxic_necrosis"],1)
        self.assertEqual(results.summary(result)["unique_deaths"],1)
        result["episode"]["events"].append(dict(time_min=30,kind="corpse_degraded",cell_id=1))
        self.assertEqual(event_totals(result["episode"])["unique_deaths"],1)

    def test_legacy_position_is_marked_approximate_and_missing_fields_are_not_zero(self):
        result, original = fixture()
        result["episode"]["events"][0].pop("position_um")
        for frame in result["episode"]["frames"]:
            frame["fields"].pop("ne",None)
        p = refresh_payload(original,result["episode"])
        self.assertEqual(p["events"][0][4:7],[10,10,True])
        self.assertEqual(p["field_grid"]["statuses"]["ne"],"unavailable")
        self.assertNotIn("ne",p["fields"][0])
        self.assertEqual(p["auto_field"],"oxygen")
        self.assertEqual(p["fields"][0]["stats"]["oxygen"],dict(min=0,max=0))
        self.assertEqual(len(p["field_grid"]["names"]),10)

    def test_zero_headroom_and_removed_bodies_still_reach_legacy_cap(self):
        result, _ = fixture()
        e = result["episode"]
        e["contract"]["biology"]["max_population"]=2
        self.assertEqual(population_limit(e)["reached_min"],0)
        self.assertEqual(population_limit(e)["allocated_bodies"],2)

    def test_block_causes_record_actual_gate_without_changing_state(self):
        sim=object.__new__(RecordedSimulation)
        sim.tissue=SimpleNamespace(ids=np.array([1,2]),xy=np.array([[1.,2.],[3.,4.]]))
        sim.p={"max_population":2};sim.div={"division_blocked_local_pressure":.4}
        sim.lp=np.array([.1,.5]);sim.time=1440;sim.events=[]
        sim.event("division_blocked_crowding",0,local_pressure_sim=.1)
        self.assertEqual(sim.events[-1]["kind"],"division_blocked_population_ceiling")
        self.assertEqual(sim.events[-1]["position_um"],[1.,2.])
        sim.p["max_population"]=3
        sim.event("division_blocked_crowding",1,local_pressure_sim=.5)
        self.assertEqual(sim.events[-1]["kind"],"division_blocked_local_pressure")
        sim.event("division_blocked_compression",0)
        self.assertEqual(sim.events[-1]["block_reasons"],["compression"])
        sim.event("growth_blocked_compression",0)
        counts=event_totals(dict(summary=dict(time_min=1440),events=sim.events,contract={}))
        self.assertEqual(counts["division_blocks_compression"],1)
        self.assertEqual(counts["growth_blocks_compression"],1)


class ComparisonTests(unittest.TestCase):
    def test_common_saved_frame_not_independent_endpoints(self):
        a,_=fixture(40);b,_=fixture(20,"recognition_off")
        a["episode"]["events"].append(dict(time_min=40,kind="division",cell_id=2))
        with patch.object(results,"load_results",side_effect=[a,b]), warnings.catch_warnings():
            warnings.simplefilter("ignore")
            table=results.compare_runs(["a","b"])
        self.assertEqual(table.attrs["common_time_min"],20)
        self.assertEqual(table.divisions.tolist(),[0,0])
        self.assertEqual(table.hypoxic_necrosis.tolist(),[1,1])
        self.assertNotEqual(table.reached_hours.iloc[0],table.reached_hours.iloc[1])

    def test_mismatched_seed_bundle_model_and_ceiling_are_rejected(self):
        for field in ("seed","bundle","model","ceiling"):
            a,_=fixture();b,_=fixture(scenario="recognition_off")
            if field=="seed":b["episode"]["seed"]=22
            if field=="bundle":b["config"]["zip_sha256"]="different"
            if field=="model":b["episode"]["contract"]["provider"]["pinned_model"]="jev-1.15.0"
            if field=="ceiling":b["episode"]["contract"]["biology"]["max_population"]=4
            with self.subTest(field=field),patch.object(results,"load_results",side_effect=[a,b]):
                with self.assertRaisesRegex(ValueError,"not matched"):
                    results.compare_runs(["a","b"])

    def test_sham_matches_and_missing_provenance_fails(self):
        a,_=fixture();b=copy.deepcopy(a)
        with patch.object(results,"load_results",side_effect=[a,b]):
            self.assertEqual(results.compare_runs(["a","b"]).divisions.tolist(),[0,0])
        b["config"]={}
        with patch.object(results,"load_results",side_effect=[a,b]),self.assertRaisesRegex(ValueError,"provenance"):
            results.compare_runs(["a","b"])

    def test_shared_replays_preserve_original_and_use_common_scales(self):
        a,oa=fixture();b,ob=fixture(20,"recognition_off")
        with tempfile.TemporaryDirectory() as td:
            for i,(result,old) in enumerate([(a,oa),(b,ob)]):
                old.update(patch="sample",scenario=result["episode"]["scenario"])
                result["replay"]=Path(td)/f"original{i}.html"
                write_payload(old,result["replay"])
            before=a["replay"].read_bytes()
            with patch.object(results,"load_results",side_effect=[a,b,a,b]),warnings.catch_warnings():
                warnings.simplefilter("ignore")
                exported=results.export_comparison_replays(["a","b"],Path(td)/"compared")
            self.assertEqual(len(exported),2)
            self.assertTrue(all(p.is_file() for p in exported))
            self.assertEqual(a["replay"].read_bytes(),before)


class PromptTests(unittest.TestCase):
    def test_non_cd8_payload_is_identical_and_cd8_cannot_change_its_batch(self):
        from test_provider import rows
        a=rows(3);a[0]["observation"].update(recognition_sim=True,excluded={"x":["recognition_declared:true"]})
        b=copy.deepcopy(a);b[0]["observation"].update(recognition_sim=False,excluded={"x":["recognition_declared:false"]})
        self.assertEqual(payload(a,"jev-1.14.0"),payload(b,"jev-1.14.0"))
        a[1]["observation"]["primary_label"]="T_CD8_cytotoxic"
        batches=build_batches(a,"jev-1.14.0",dict(max_questions=20,max_payload_bytes=100000))
        self.assertEqual([[r["cell_id"] for r in batch] for batch in batches],[[0,2],[1]])
        self.assertIn("frozen",payload(a,"jev-1.14.0")["state"]["semantics"])


class PreparationTests(unittest.TestCase):
    def test_unusable_manual_center_stops_before_scoring(self):
        prep=Mock();prep.make_region.return_value=None
        sample=dict(obs=pd.DataFrame(dict(compartment=["tumor_candidate"],qc_pass=[True])),output_dir="audit")
        with self.assertRaisesRegex(own_data.PreparationScopeError,"before AUCell"):
            own_data.preflight_regions(prep,sample,center=[1,2])
        prep.score_aucell.assert_not_called()
        prep.propose_regions.assert_not_called()

    def test_non_tumour_annotation_scope_warning(self):
        self.assertTrue(own_data.annotation_scope_warnings([dict(compartment="immune")]))
        self.assertFalse(own_data.annotation_scope_warnings([dict(compartment="tumor_candidate")]))


@unittest.skipUnless(os.environ.get("VTX_TEST_BUNDLE"),"set VTX_TEST_BUNDLE")
class RealBundleRegression(unittest.TestCase):
    def test_reporting_preserves_one_hour_biological_state_and_legacy_bundle(self):
        with tempfile.TemporaryDirectory() as td:
            wf=vtx.workflow
            bundle=wf.extract_bundle(os.environ["VTX_TEST_BUNDLE"],td)
            config=dict(scenario="declared_recognition",seed=21,horizon_hours=1,
                        max_population=1000,max_decisions=10000,model="jev-1.14.0")
            new=wf.create_sim(bundle,config)
            from xenium_runtime.v2.policy import ScoreAwareLocalPolicy
            old=AutonomousSimulation(bundle["atlas"],wf.effective_contract(bundle,config),
                scenario=config["scenario"],seed=21,policy=ScoreAwareLocalPolicy(),
                max_decisions=10000,stop_on_provider_failure=True)
            while not new.stop_reason:new.step()
            while not old.stop_reason:old.step()
            a,b=new.episode(),old.episode()
            self.assertEqual(a["final_state_fingerprint"],b["final_state_fingerprint"])
            self.assertEqual(a["decisions"],b["decisions"])
            self.assertEqual(a["metrics"],b["metrics"])
            from xenium_runtime.v2.episode import validate_episode
            self.assertTrue(validate_episode(a,bundle["atlas"])["verified"])


if __name__=="__main__":unittest.main()
