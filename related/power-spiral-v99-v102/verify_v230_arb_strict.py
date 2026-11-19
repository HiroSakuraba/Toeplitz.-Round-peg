"""Fail-closed entry point for recovered v101 Arb eight-face certificates.

The historical verify_v101_arb.py prints FAILED on unresolved cells but exits
with status 0. This adapter turns *any* unresolved cell or nonpositive sign
margin into an ordinary Python assertion failure, suitable for GitHub CI.
It requires python-flint; this is genuinely independent of iv101.py.
"""
from verify_v101_arb import TESTS, adaptive_face, PREC


def certify_all():
    results=[]
    for box,name,kind,val,target,sign,alpha in TESTS:
        result=adaptive_face(box,kind,val,target,sign,alpha)
        assert result['bad']==0, (name,'unresolved',result['bad_examples'])
        assert result['margin'] is not None and result['margin']>0, (name,result)
        results.append((name,result['nodes'],result['leaves'],result['bad']))
    assert len(results)==8
    return tuple(results)


if __name__=='__main__':
    for name,nodes,leaves,bad in certify_all():
        print(f'{name:9s} nodes={nodes:6d} leaves={leaves:6d} unresolved={bad}')
    print(f'v230: all eight recovered v101 Arb face certificates certified at {PREC} bits')
