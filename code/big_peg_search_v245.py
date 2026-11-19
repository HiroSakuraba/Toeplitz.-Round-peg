"""v245 numerical-only candidate generator for exact Big Peg certification.

All discovered shapes are UNVERIFIED until the Fraction-based
code/big_peg_exact_v243.py maximum-square and disk checker accepts them.
Requires numpy. The initial v245 candidate was found using coarse and fine
integer-grid neighborhood search around Pochuev's rational decagon.

Run:
  python code/big_peg_search_v245.py                    # replay approximate scores
  python code/big_peg_search_v245.py --steps 40 --seed 245  # explore nearby grids

This is a heuristic, NOT a complete square census and NOT a theorem.
"""
import argparse
import json
from pathlib import Path
from fractions import Fraction as Q

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"research"/"v244-pochuev-original-ten-gon.json"
CANDIDATE=ROOT/"research"/"v245-improved-ten-gon.json"
DENOM=12000
N=10

# All ordered four-edge assignments other than the obvious all-one-edge
# singular charts. Other singular systems reject a search candidate; they
# must never be treated as evidence for nonexistence of squares.
ASSIGNMENTS=np.indices((N,N,N,N)).reshape(4,-1).T
ASSIGNMENTS=ASSIGNMENTS[np.any(ASSIGNMENTS!=ASSIGNMENTS[:,0:1],axis=1)]


def _grid(path):
    raw=json.loads(Path(path).read_text(encoding="utf-8"))
    return np.array([[int(Q(v)*DENOM) for v in xy] for xy in raw["vertices"]],dtype=np.int64)


def _center(path):
    raw=json.loads(Path(path).read_text(encoding="utf-8"))
    return np.array([float(Q(v)) for v in raw["disk"]["center"]])


def numerical_max_side_squared(points):
    """Approximate candidate ranking only; skips known singular charts."""
    P=np.asarray(points,float)/DENOM
    E=np.roll(P,-1,axis=0)-P
    normals=np.column_stack((-E[:,1],E[:,0]))
    B=P[ASSIGNMENTS]
    edges=E[ASSIGNMENTS]
    nor=normals[ASSIGNMENTS]
    nx,ny=nor[:,:,0],nor[:,:,1]
    M=np.empty((len(ASSIGNMENTS),4,4))
    M[:,:,:2]=nor
    M[:,:,2]=np.stack((nx[:,0],ny[:,1],-nx[:,2],-ny[:,3]),axis=1)
    M[:,:,3]=np.stack((ny[:,0],-nx[:,1],-ny[:,2],nx[:,3]),axis=1)
    rhs=np.sum(nor*B,axis=2)
    try:
        sol=np.linalg.solve(M,rhs[...,None])[...,0]
    except np.linalg.LinAlgError:
        return None
    c,u=sol[:,:2],sol[:,2:]
    corners=np.empty((len(sol),4,2))
    corners[:,0]=c+u
    corners[:,1]=c+np.stack((-u[:,1],u[:,0]),axis=1)
    corners[:,2]=c-u
    corners[:,3]=c+np.stack((u[:,1],-u[:,0]),axis=1)
    projection=np.sum((corners-B)*edges,axis=2)/np.sum(edges**2,axis=2)
    residual=np.abs(np.sum(nor*(corners-B),axis=2))
    finite=np.isfinite(sol).all(axis=1)
    valid=(finite & np.all((projection>=-1e-10)&(projection<=1+1e-10)
                           &(residual<1e-8),axis=1))
    side2=2*np.sum(u*u,axis=1)
    side2[~valid]=0
    return (float(np.max(side2)),int(np.count_nonzero(side2>1e-10)))


def approximate_disk_radius_squared(points,centre):
    """Radius of largest disk with this FIXED centre; no global inradius claim."""
    P=np.asarray(points,float)/DENOM
    A=P; E=np.roll(P,-1,axis=0)-A
    t=np.clip(np.sum((centre-A)*E,axis=1)/np.sum(E*E,axis=1),0,1)
    closest=A+t[:,None]*E
    return float(np.min(np.sum((closest-centre)**2,axis=1)))


def numerical_score(points,centre):
    side=numerical_max_side_squared(points)
    if side is None or side[0]<=0:
        return None
    rad2=approximate_disk_radius_squared(points,centre)
    return (side[0]/rad2,side[0],rad2,side[1])


def greedy_grid_neighborhood(points,centre,max_steps=10):
    """Coordinate descent on the rational 1/12000 grid; unverified outputs."""
    points=points.copy()
    history=[]
    for iteration in range(max_steps):
        current=numerical_score(points,centre)
        if current is None:
            raise RuntimeError("numerical baseline is unusable")
        best=(current[0],None,None,None)
        for i in range(N):
            for j in range(2):
                for move in (-1,1):
                    trial=points.copy()
                    trial[i,j]+=move
                    result=numerical_score(trial,centre)
                    if result is not None and result[0]<best[0]-1e-12:
                        best=(result[0],i,j,move)
        if best[1] is None:
            break
        _,i,j,move=best
        points[i,j]+=move
        history.append({"vertex":i,"axis":j,"delta_numerator":move,
                        "numerical_ratio_squared":best[0]})
    return points,history


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--steps",type=int,default=0,help="numerical greedy trial steps")
    args=parser.parse_args()
    if args.steps<0:
        parser.error("steps must be nonnegative")
    baseline=_grid(BASE)
    candidate=_grid(CANDIDATE)
    out={
        "status":"NUMERICAL DIAGNOSTICS ONLY, not a rigorous certificate",
        "original_source_origin_disk":numerical_score(baseline,np.zeros(2)),
        "v245_rational_candidate":numerical_score(candidate,_center(CANDIDATE)),
        "grid_denominator":DENOM,
    }
    if args.steps:
        trial,history=greedy_grid_neighborhood(candidate,_center(CANDIDATE),args.steps)
        out["experimental_steps"]=history
        out["new_candidate_integer_grid"]=trial.tolist()
        out["new_candidate_numerical_score"]=numerical_score(trial,_center(CANDIDATE))
        out["warning"]="Must run v243 exact census and disk containment before any new claim."
    print(json.dumps(out,indent=2))


if __name__=="__main__":
    main()
