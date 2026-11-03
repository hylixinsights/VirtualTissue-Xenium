import math
from virtualtissue_xenium.v2.provider import JevProvider
from virtualtissue_xenium.v2.provider_audit import audit_journals
from test_v2_provider import rows, context


def test_independent_journals_measure_variance_and_batch_order(tmp_path):
    paths = [tmp_path/'a',tmp_path/'b']
    for i,path in enumerate(paths):
        p = [.6,.8][i]
        def transport(body, p=p):
            return dict(model=body['model'],answers={qid:dict(type='choice',choice='REST',
                probabilities={'REST':p,'ACTIVATE':1-p},confidence=.8) for qid in body['questions']},
                usage={'input_tokens':10,'output_tokens':3})
        provider=JevProvider('',path,{'replicate':i},transport=transport)
        provider.choose(rows(2)[::1 if i==0 else -1],context())
    report=audit_journals(paths)
    assert report['n_matched_questions']==2
    assert all(q['question_order_varied'] for q in report['questions'])
    assert all(math.isclose(q['probability_sample_variance']['REST'],.02) for q in report['questions'])
    assert report['usage']['input_tokens']==20
