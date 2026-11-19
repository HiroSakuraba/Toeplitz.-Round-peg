"""Numerical-only v247 grid optimization. Never a mathematical certificate.

Reads rational v245/v247 JSON, explores small single-coordinate grid moves,
and ranks candidates with a floating-point 10^4 edge-chart square search.
Every accepted theorem must be replayed by the independent Fraction-based
v243 square engine and v246 global inradius verifier.

Requires numpy. A singular or ambiguous numerical input is rejected from
ranking; singular charts MUST be handled by the exact proof engine.
"""
import argparse
import json
from fractions import Fraction as F
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SOURCES={"v245":"research/v245-improved-ten-gon.json",
         "v247":"research/v247-refined-ten-gon.json"}
N=10
IDS=np.indices((N,N,N,N)).reshape(4,-1).T
IDS=IDS[np.any(IDS!=IDS[:,0:1],axis=1)]
CONTACTS=(1,6,8)
REFERENCE_C=np.array([0.000302622629,0.0000002143102])


def read(name):
    data=json.loads((ROOT/SOURCES[name]).read_text(encoding="utf-8"))
    return np.array([[float(F(x)),float(F(y))] for x,y in data["vertices"]])


def distance_squared_to_edges(P,c):
    e=np.roll(P,-1,axis=0)-P
    t=np.clip(np.sum((c-P)*e,axis=1)/np.sum(e*e,axis=1),0,1)
    return np.sum((P+t[:,None]*e-c)**2,axis=1)


def candidate_centre(P):
    """Three-contact numerical heuristic; exact disk is checked elsewhere."""
    e=np.roll(P,-1,axis=0)-P
    A=[];b=[]
    for i in CONTACTS:
        n=np.array([-e[i,1],e[i,0]])/np.linalg.norm(e[i])
        if np.dot(REFERENCE_C-P[i],n)<0:n=-n
        A.append([n[0],n[1],-1])
        b.append(np.dot(n,P[i]))
    try:x,y,r=np.linalg.solve(np.array(A),np.array(b))
    except np.linalg.LinAlgError:return None
    c=np.array([x,y])
    if not (np.isfinite(c).all() and .95<r<1.2):return None
    d2=distance_squared_to_edges(P,c)
    if min(d2)<r*r-1e-7:return None
    return c,float(min(d2))


def numerical_square_max(P):
    e=np.roll(P,-1,axis=0)-P
    n=np.stack([-e[:,1],e[:,0]],axis=1)[IDS]
    nx,ny=n[:,:,0],n[:,:,1]
    A=np.empty((len(IDS),4,4))
    A[:,:,0]=nx;A[:,:,1]=ny
    A[:,:,2]=np.stack([nx[:,0],ny[:,1],-nx[:,2],-ny[:,3]],axis=1)
    A[:,:,3]=np.stack([ny[:,0],-nx[:,1],-ny[:,2],nx[:,3]],axis=1)
    b=np.sum(n*P[IDS],axis=2)
    use=np.abs(np.linalg.det(A))>1e-13
    if np.count_nonzero(~use)>200:return None
    try:solutions=np.linalg.solve(A[use],b[use,:,None])[:,:,0]
    except np.linalg.LinAlgError:return None
    c=solutions[:,:2];u=solutions[:,2:]
    corners=np.empty((len(c),4,2))
    corners[:,0]=c+u
    corners[:,1]=c+np.column_stack([-u[:,1],u[:,0]])
    corners[:,2]=c-u
    corners[:,3]=c+np.column_stack([u[:,1],-u[:,0]])
    assigned=IDS[use]
    edge=e[assigned];origins=P[assigned]
    t=np.sum((corners-origins)*edge,axis=2)/np.sum(edge*edge,axis=2)
    nr=np.stack([-e[:,1],e[:,0]],axis=1)[assigned]
    residual=np.abs(np.sum((corners-origins)*nr,axis=2))
    ok=np.all((t>=-1e-9)&(t<=1+1e-9)&(residual<1e-8),axis=1)
    sides2=2*np.sum(u*u,axis=1)
    sides2[~ok]=0
    k=int(np.argmax(sides2))
    if sides2[k]<=0:return None
    return float(sides2[k]),int(np.count_nonzero(sides2>1e-12)),assigned[k].tolist()


def score(P):
    disk=candidate_centre(P)
    sq=numerical_square_max(P)
    if disk is None or sq is None:return None
    return {"approximate_ratio":float(np.sqrt(sq[0]/disk[1])),
            "square_side_squared":sq[0],
            "disk_radius_squared":disk[1],
            "positive_oriented_assignments":sq[1],
            "max_assignment":sq[2]}


def greedy(start,grid,steps):
    P=start.copy()
    moves=[]
    for _ in range(steps):
        current=score(P)
        if current is None:raise RuntimeError("invalid numeric baseline")
        best=current["approximate_ratio"];chosen=None
        for i in range(N):
            for axis in (0,1):
                for sign in (-1,1):
                    trial=P.copy();trial[i,axis]+=sign/grid
                    s=score(trial)
                    if s is not None and s["approximate_ratio"]<best-1e-10:
                        best=s["approximate_ratio"];chosen=(i,axis,sign)
        if chosen is None:break
        i,axis,sign=chosen;P[i,axis]+=sign/grid
        moves.append({"vertex":i,"axis":axis,"delta_grid_units":sign,
                      "approximate_ratio":best})
    return P,moves


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start",choices=tuple(SOURCES),default="v247")
    parser.add_argument("--steps",type=int,default=0)
    parser.add_argument("--grid",type=int,default=1200000)
    args=parser.parse_args()
    if args.steps<0 or args.steps>1000 or args.grid<=0:
        parser.error("steps must be 0..1000 and grid >0")
    scores={name:score(read(name)) for name in SOURCES}
    if any(s is None for s in scores.values()):
        raise RuntimeError("numerical audit could not rank known candidates")
    result={"status":"NUMERICAL HEURISTIC ONLY; not a certificate",
            "known_candidate_scores":scores}
    if args.steps:
        P,moves=greedy(read(args.start),args.grid,args.steps)
        result["grid_search_moves"]=moves
        result["unverified_candidate_vertices"]=P.tolist()
        result["unverified_candidate_score"]=score(P)
        result["warning"]="Rationalize and run v243/v246 exact proofs before claiming improvement."
    print(json.dumps(result,indent=2,sort_keys=True))
