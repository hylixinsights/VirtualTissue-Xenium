"""Exercise installed command wiring in addition to module-level tests."""
from virtualtissue_xenium.v2.cli import main
from virtualtissue_xenium.v2.inputs import write_json, read_json
from test_v3_contracts import cells, context


def test_prepare_control_ensemble_and_sensitivity_commands(tmp_path):
    write_json(tmp_path/'cells.json', {'cells':cells(4), 'assayed_genes':list(cells(1)[0]['expression'])})
    write_json(tmp_path/'context.json', context())
    assert main(['prepare',str(tmp_path/'cells.json'),str(tmp_path/'context.json'),str(tmp_path/'base.json')])==0
    assert main(['control',str(tmp_path/'base.json'),str(tmp_path/'beta.json'),'--kind','beta_zero'])==0
    write_json(tmp_path/'arms.json',{'arms':{'control':'base.json','beta':'beta.json'}})
    assert main(['ensemble',str(tmp_path/'arms.json'),str(tmp_path/'ensemble'),
                 '--horizon-hours','1','--seeds','0','1','--bootstrap-replicates','100'])==0
    assert read_json(tmp_path/'ensemble/report.json')['n_completed_runs']==4
    write_json(tmp_path/'bounds.json',{'action_priors.malignant.PROLIFERATE':[.003,.009]})
    assert main(['sensitivity',str(tmp_path/'base.json'),str(tmp_path/'bounds.json'),str(tmp_path/'sensitivity'),
                 '--horizon-hours','1','--seeds','0','1','--n','2','--method','lhs',
                 '--bootstrap-replicates','100'])==0


def test_annotation_agreement_command(tmp_path):
    text='cell_id,cell_type,lineage\na,tumour,malignant\nb,T,T_CD8\n'
    (tmp_path/'a.csv').write_text(text)
    (tmp_path/'b.csv').write_text(text)
    assert main(['annotation-agreement',str(tmp_path/'a.csv'),str(tmp_path/'b.csv'),str(tmp_path/'agreement')])==0
