"""Exact rational cap-residue counterexample. No Floer solver."""
from fractions import Fraction as Q
from itertools import combinations

E = ((1, 0), (0, 1), (-1, 0), (0, -1))


def add(a, b):
    return tuple(x+y for x,y in zip(a,b))


def sub(a, b):
    return tuple(x-y for x,y in zip(a,b))


def scale(a, s):
    return tuple(s*x for x in a)


def det(a, b):
    return a[0]*b[1]-a[1]*b[0]


def integral(path):
    return sum((det(a,b)/2 for a,b in zip(path,path[1:])), Q(0))


def gadget(n):
    if type(n) is not int or n < 1:
        raise ValueError("positive integer number of turns required")
    m,d=4*n,Q(1,16*n)
    u=[scale(E[j%4],Q(3,4)-j*d) for j in range(m+1)]
    v=[scale(E[j%4],Q(3,4)-(j+2)*d) for j in range(-1,m)]
    p=[(Q(1),Q(0)),*u,(Q(0),Q(0)),*reversed(v),
       (Q(0),Q(-1)),(Q(-1),Q(-1)),(Q(-1),Q(1)),
       (Q(0),Q(1)),(Q(1),Q(1))]
    return p


def fixture(k):
    if type(k) is not int or k < 2:
        raise ValueError("scale denominator must be an integer >= 2")
    e=Q(1,k)
    p=[scale(sub(a,(1,0)),e) for a in gadget(k*k)]
    p.extend([(Q(0),Q(4)),(Q(4),Q(4)),(Q(4),Q(0))])
    return p


def joint_formula(n):
    if type(n) is not int or n < 1:
        raise ValueError("positive integer number of turns required")
    return Q(19*n,24)-Q(1,384*n)-Q(1,2)


def scaled_payment(k):
    if type(k) is not int or k < 2:
        raise ValueError("scale denominator must be an integer >= 2")
    return joint_formula(k*k)/k**2


def on(a,b,p):
    return (det(sub(b,a),sub(p,a)) == 0 and
            all(min(x,y)<=z<=max(x,y) for x,y,z in zip(a,b,p)))


def intersects(a,b,c,d):
    x,y=det(sub(b,a),sub(c,a)),det(sub(b,a),sub(d,a))
    z,w=det(sub(d,c),sub(a,c)),det(sub(d,c),sub(b,c))
    return ((x*y<0 and z*w<0) or on(a,b,c) or on(a,b,d)
            or on(c,d,a) or on(c,d,b))


def simple(p):
    n=len(p)
    if n<3 or len(set(p))!=n:
        return False
    for i in range(n):
        a,b,c=p[i-1],p[i],p[(i+1)%n]
        if det(sub(b,a),sub(c,b)) == 0 and (on(a,b,c) or on(b,c,a)):
            return False
    edges=list(zip(p,p[1:]+p[:1]))
    for i,j in combinations(range(n),2):
        if j==i+1 or (i==0 and j==n-1):
            continue
        a,b=edges[i];c,d=edges[j]
        if intersects(a,b,c,d):
            return False
    return True


def winding(p, x):
    """Closed polyline passed with its endpoints equal."""
    if len(p)<2 or p[0]!=p[-1]:
        raise ValueError("expected closed path")
    count=0
    for a,b in zip(p,p[1:]):
        if on(a,b,x):
            raise ValueError("winding is undefined on trace")
        z=det(sub(b,a),sub(x,a))
        if a[1]<=x[1]<b[1] and z>0:
            count+=1
        if b[1]<=x[1]<a[1] and z<0:
            count-=1
    return count


def payment(k, check_simple=True):
    p=fixture(k)
    if check_simple and not simple(p):
        raise AssertionError("non-simple fixture")
    e=Q(1,k)
    a,b,c,d=(Q(0),Q(0)),(-e,Q(0)),(-e,e),(Q(0),e)
    ia,ib,ic,id_=[p.index(v) for v in (a,b,c,d)]
    assert ia==0<ib<ic<id_
    # Start diagonal (b,d) rotates to (a,c). Both paths follow the loop.
    assert sub(a,((b[0]+d[0])/2,(b[1]+d[1])/2)) == (
        -(b[1]-d[1])/2,(b[0]-d[0])/2)
    joint=integral(p[ia:ib+1])+integral(p[ic:id_+1])
    assert joint==scaled_payment(k)
    area=abs(integral(p+[p[0]]))
    assert abs(area-16)<=4*e*e and area>=15
    assert Q(1,2)<joint<1
    return {"k":k,"vertices":len(p),"side":e,"area":area,
            "action_residue":(-joint)%area,"action_gap":joint}


def endpoint_model(b,x,mixing=0):
    b,x,mixing=map(Q,(b,x,mixing))
    if not 0<b<=Q(1,2) or not 0<=x<=1 or not 0<=mixing<=1:
        raise ValueError("outside endpoint-model domain")
    f=lambda t:(1-mixing)*min(t,b)+mixing*b*t
    return f(x),b-f(1-x)


if __name__=="__main__":
    for k in (2,3,4):
        print(payment(k))
    for k in (8,32,1024):
        print({"k":k,"closed_form_gap":scaled_payment(k)})
