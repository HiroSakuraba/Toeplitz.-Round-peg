"""v231 experimental numerical discovery of square charts on a sampled power spiral.

Numerical candidates are NOT proofs of a square on the smooth spiral.
Only rationalized polygon chart roots can be checked exactly; a future
interval/degree bridge may certify a square on the continuous spiral.
"""
from math import sin,cos,hypot,floor,ceil,factorial
from fractions import Fraction as Q
from dataclasses import dataclass
from collections import defaultdict

@dataclass(frozen=True)
class Sample:
    t: float
    x: float
    y: float


def spiral_samples(p, inner_radius, step=0.12):
    """Rational t=u^2 makes theta=u^(-2p) rational for p=2,5/2,3."""
    p=Q(str(p))
    if 2*p not in (Q(4),Q(5),Q(6)):
        raise ValueError('certified sampler currently supports p=2,5/2,3')
    if not 0<inner_radius<1 or step<=0:raise ValueError('invalid sampling')
    theta_max=inner_radius**(-float(p))
    n=max(2,int((theta_max-1)/step)+1)
    samples=[]
    for i in range(n+1):
        th=1+(theta_max-1)*i/n
        uq=Q(round(th**(-1/(2*float(p)))*2**22),2**22)
        t=uq*uq
        theta=1/(uq**int(2*p))
        samples.append(Sample(t,float(t)*cos(float(theta)),float(t)*sin(float(theta))))
    return tuple(samples)


def dist_segment_sq(p,a,b):
    vx=b.x-a.x;vy=b.y-a.y; w=(p[0]-a.x)*vx+(p[1]-a.y)*vy
    u=max(0.,min(1.,w/(vx*vx+vy*vy)))
    return (p[0]-a.x-u*vx)**2+(p[1]-a.y-u*vy)**2

class SpatialEdges:
    def __init__(self,pts,cell=None):
        self.pts=pts
        lengths=[hypot(pts[i+1].x-pts[i].x,pts[i+1].y-pts[i].y) for i in range(len(pts)-1)]
        self.max_edge=max(lengths)
        self.cell=cell or max(0.017,3*self.max_edge)
        self.grid=defaultdict(list)
        for i in range(len(pts)-1):
            a=pts[i];b=pts[i+1]
            x0=floor(min(a.x,b.x)/self.cell);x1=floor(max(a.x,b.x)/self.cell)
            y0=floor(min(a.y,b.y)/self.cell);y1=floor(max(a.y,b.y)/self.cell)
            for x in range(x0,x1+1):
                for y in range(y0,y1+1):self.grid[x,y].append(i)
    def nearest_edges(self,p,limit=4,max_distance=None):
        x=floor(p[0]/self.cell);y=floor(p[1]/self.cell)
        inds=set()
        for dx in (-1,0,1):
            for dy in (-1,0,1):inds.update(self.grid.get((x+dx,y+dy),()))
        ds=sorted((dist_segment_sq(p,self.pts[i],self.pts[i+1]),i) for i in inds)
        d2=(max_distance or self.max_edge*2.0)**2
        return tuple((d,i) for d,i in ds if d<=d2)[:limit]


def square_residual(pts,sg):
    (x0,y0),(x1,y1),(x2,y2),(x3,y3)=pts
    ex=x1-x0;ey=y1-y0
    return (x2-x1+sg*ey,y2-y1-sg*ex,x3-x0+sg*ey,y3-y0-sg*ex)


def solve_chart(samples,ids,sg,eps=1e-9):
    a=[(samples[i].x,samples[i].y) for i in ids]
    d=[(samples[i+1].x-samples[i].x,samples[i+1].y-samples[i].y) for i in ids]
    b=[-z for z in square_residual(a,sg)]
    cols=[]
    for k in range(4):
        q=list(a);q[k]=(q[k][0]+d[k][0],q[k][1]+d[k][1]);f=square_residual(q,sg)
        cols.append(tuple(f[j]+b[j] for j in range(4)))
    M=[[cols[j][i] for j in range(4)]+[b[i]] for i in range(4)]
    for k in range(4):
        piv=max(range(k,4),key=lambda i:abs(M[i][k]))
        if abs(M[piv][k])<1e-14:return None
        M[k],M[piv]=M[piv],M[k]
        z=M[k][k]
        for j in range(k,5):M[k][j]/=z
        for i in range(4):
            if i!=k:
                fac=M[i][k]
                for j in range(k,5):M[i][j]-=fac*M[k][j]
    u=[M[i][4] for i in range(4)]
    if any(z<-eps or z>1+eps for z in u):return None
    result=[(a[i][0]+u[i]*d[i][0],a[i][1]+u[i]*d[i][1]) for i in range(4)]
    if max(map(abs,square_residual(result,sg)))>1e-7:return None
    return tuple(u),tuple(result)


def discover(p,inner=0.28,step=0.12,max_pairs=1_000_000, topk=5, min_margin=0.07):
    samples=spiral_samples(p,inner,step)
    spatial=SpatialEdges(samples)
    n=len(samples)-1
    attempts=0;near=0;best=None
    # Scan edges from outer to inner. No v100 Miranda box or center is provided.
    for i in range(n):
        a=samples[i];A=((a.x+samples[i+1].x)/2,(a.y+samples[i+1].y)/2)
        for j in range(i+2,n):
            attempts+=1
            if attempts>max_pairs:return best or {'status':'budget','pairs':attempts,'near':near,'samples':n}
            c=samples[j];C=((c.x+samples[j+1].x)/2,(c.y+samples[j+1].y)/2)
            vx=(C[0]-A[0])/2;vy=(C[1]-A[1])/2
            center=((A[0]+C[0])/2,(A[1]+C[1])/2)
            B=(center[0]-vy,center[1]+vx);D=(center[0]+vy,center[1]-vx)
            bs=spatial.nearest_edges(B,limit=topk)
            ds=spatial.nearest_edges(D,limit=topk)
            if not bs or not ds:continue
            near+=1
            for bdist,k in bs:
                for ddist,l in ds:
                    for ids,sg in (((i,k,j,l),-1),((i,l,j,k),1)):
                        if len(set(ids))<4:continue
                        x=solve_chart(samples,ids,sg)
                        if x is not None:
                            u,pts=x
                            ss=hypot(pts[1][0]-pts[0][0],pts[1][1]-pts[0][1])
                            if ss>0.18 and min(hypot(*pt) for pt in pts)>0.20:
                                row={'status':'found','p':p,'inner':inner,'step':step,'samples':n,'pairs':attempts,'near':near,'ids':ids,'sigma':sg,'parameters':u,'points':pts,'side':ss,'samples_data':samples}
                                margin=min(min(v,1-v) for v in u)
                                if best is None or margin>min(min(v,1-v) for v in best['parameters']):best=row
                                if margin>=min_margin:return row
    return best or {'status':'not_found','pairs':attempts,'near':near,'samples':n}



# The certification path uses rational interval trigonometry and a 4-dimensional
# homotopy degree, not a near-zero numerical residual or a priori Miranda box.

def _atan_interval(z,terms=34):
    z=Q(z)
    acc=Q(0)
    for k in range(terms):acc+=(-1)**k*z**(2*k+1)/Q(2*k+1)
    nxt=z**(2*terms+1)/Q(2*terms+1)
    return (acc,acc+nxt) if terms%2==0 else (acc-nxt,acc)


def _pi_interval():
    a,b=_atan_interval(Q(1,5));c,d=_atan_interval(Q(1,239))
    return 16*a-4*d,16*b-4*c


def _round_fraction(v,scale=2**40):
    x=v*scale
    return Q((2*x.numerator+x.denominator)//(2*x.denominator),scale)


def _trig_rational_enclosure(theta,terms=20):
    """Rational approximants plus proven per-coordinate absolute error.

    Machin pi enclosure; modular reduction with exact rational midpoint;
    Taylor's theorem with |derivative|<=1, plus argument uncertainty.
    """
    theta=Q(theta)
    pl,pu=_pi_interval();pm=(pl+pu)/2
    k=round(theta/(2*pm))
    raw=theta-2*k*pm
    x=_round_fraction(raw,2**34)
    arg_err=abs(raw-x)+abs(k)*(pu-pl)
    term_s=x;term_c=Q(1);s=x;c=Q(1)
    for j in range(1,terms+1):
        term_s*=-x*x/Q((2*j)*(2*j+1))
        term_c*=-x*x/Q((2*j-1)*(2*j))
        s+=term_s;c+=term_c
    absx=abs(x)
    err_s=absx**(2*terms+2)/Q(factorial(2*terms+2))+arg_err
    err_c=absx**(2*terms+1)/Q(factorial(2*terms+1))+arg_err
    return (s,c,max(err_s,err_c))


def _certified_point(sample,p):
    p=Q(str(p));t=Q(sample.t)
    # t = u^2; therefore theta=t^(-p) is rational for requested half-integer p.
    num=t.numerator; den=t.denominator
    from math import isqrt
    un=isqrt(num);ud=isqrt(den)
    if un*un!=num or ud*ud!=den:
        raise ValueError('sampler t must be a rational square')
    u=Q(un,ud);theta=u**(-int(2*p))
    ss,cc,trig_e=_trig_rational_enclosure(theta)
    x=_round_fraction(t*cc);y=_round_fraction(t*ss)
    error=max(abs(x-t*cc),abs(y-t*ss))+t*trig_e
    return (x,y),error


def _R_residual(v,sg):
    (x0,y0),(x1,y1),(x2,y2),(x3,y3)=v
    return (x2-x1+sg*(y1-y0),y2-y1-sg*(x1-x0),
            x3-x0+sg*(y1-y0),y3-y0-sg*(x1-x0))


def _R_matrix(pts_a,pts_b,sg):
    base=_R_residual(pts_a,sg)
    cols=[]
    for j in range(4):
        z=list(pts_a);z[j]=pts_b[j]
        r=_R_residual(z,sg)
        cols.append(tuple(r[k]-base[k] for k in range(4)))
    return tuple(tuple(cols[j][i] for j in range(4)) for i in range(4)),tuple(-b for b in base)


def _inverse_and_root(A,b):
    M=[list(A[i])+[Q(int(i==j)) for j in range(4)]+[b[i]] for i in range(4)]
    for k in range(4):
        pivot=next((i for i in range(k,4) if M[i][k]),None)
        if pivot is None:raise ValueError('singular chart')
        M[k],M[pivot]=M[pivot],M[k]
        fac=M[k][k]
        M[k]=[v/fac for v in M[k]]
        for i in range(4):
            if i==k:continue
            fac=M[i][k]
            M[i]=[v-fac*w for v,w in zip(M[i],M[k])]
    return tuple(row[-1] for row in M),max(sum(abs(v) for v in row[4:8]) for row in M)


def _side_floor_holds(starts, ends, root, error, floor):
    """Bound the side uniformly over the first two source intervals."""
    vertices = []
    allowance = 4*error
    for a, b, v in zip(starts[:2], ends[:2], root[:2]):
        vertices.append(tuple(x+v*(y-x) for x, y in zip(a, b)))
        allowance += max(v, 1-v) * sum(abs(y-x) for x, y in zip(a, b))
    (x0, y0), (x1, y1) = vertices
    return (x1-x0)**2 + (y1-y0)**2 > (Q(floor)+allowance)**2


def _chord_error(t0,t1,p):
    """Sup coordinate error for linear interpolation in rational radius t.

    gamma''(t)=e^(i t^-p)[-p^2 t^(-2p-1)+i p(p-1)t^(-p-1)]
    and interpolation error <= |dt|^2 sup |gamma''|/8.
    The ceil exponents are a safe conservative rational bound for 0<t<=1.
    """
    p=Q(str(p));a=min(t0,t1);dt=abs(t1-t0)
    if not (0<a<=1):raise ValueError('radii outside (0,1]')
    m1=ceil(2*p+1);m2=ceil(p+1)
    return dt*dt/Q(8)*(p*p*a**(-m1)+abs(p*(p-1))*a**(-m2))


def certify_on_smooth_spiral(candidate, radius_floor=Q(1,5),side_floor=Q(1,8)):
    """A posteriori Brouwer-degree certificate, entirely independent of Box A.

    The numerical search supplies only four edge indices. The certificate
    uses exact rational vertices enclosing actual transcendental spiral
    endpoints, a deterministic second-derivative interpolation bound,
    and an exact affine boundary margin. If successful, every homotopy
    retains nonzero degree on [0,1]^4; hence the TRUE smooth spiral
    has a square on these four source arcs with the asserted size/radii.
    """
    if candidate['status']!='found':raise ValueError('numerical candidate required')
    P=Q(str(candidate['p']));ids=candidate['ids'];sg=candidate['sigma']
    samples=candidate['samples_data']
    ends=[(_certified_point(samples[j],P),_certified_point(samples[j+1],P)) for j in ids]
    pa=tuple(pair[0][0] for pair in ends);pb=tuple(pair[1][0] for pair in ends)
    A,b=_R_matrix(pa,pb,sg)
    try:u,inv_norm=_inverse_and_root(A,b)
    except ValueError:return {'certified':False,'reason':'singular-rational-chart'}
    margin=min(min(z,1-z) for z in u)
    if margin<=0:return {'certified':False,'reason':'rational-root-exits-box'}
    ep=[]
    for idx,((_,err0),(_,err1)) in zip(ids,ends):
        err=max(err0,err1)+_chord_error(samples[idx].t,samples[idx+1].t,P)
        ep.append(err)
    E=max(ep);boundary=margin/inv_norm
    # Each component of the square residual involves at most four coordinates.
    if 4*E>=boundary:
        return {'certified':False,'reason':'degree-margin-insufficient',
                'ratio':float(4*E/boundary),'margin':float(margin),
                'bound':float(E),'edge_indices':ids}
    square=tuple((pa[i][0]+u[i]*(pb[i][0]-pa[i][0]),
                  pa[i][1]+u[i]*(pb[i][1]-pa[i][1])) for i in range(4))
    assert _R_residual(square,sg)==(0,0,0,0)
    side2=(square[1][0]-square[0][0])**2+(square[1][1]-square[0][1])**2
    # The affine and true zeros need not have the same source parameters.
    if not _side_floor_holds(pa,pb,u,E,side_floor):
        return {'certified':False,'reason':'side-floor-insufficient','edge_indices':ids}
    if any(min(samples[j].t,samples[j+1].t)<=Q(radius_floor) for j in ids):
        return {'certified':False,'reason':'radius-floor-insufficient','edge_indices':ids}
    return {'certified':True,'p':str(P),'inner_search_radius':candidate['inner'],
            'edge_indices':ids,'orientation':sg,'degree_boundary_margin':float(boundary),
            'max_residual_perturbation':float(4*E),'max_coordinate_error':float(E),
            'minimum_rational_source_margin':float(margin),
            'polygon_side_length':float(side2**Q(1,2)) if False else float(side2)**0.5,
            'certified_side_strictly_greater_than':str(Q(side_floor)),
            'certified_min_radius_strictly_greater_than':str(Q(radius_floor)),
            'verified_rational_root':True,
            'exact_separator_4E_lt_1over20000_lt_boundary': 4*E < Q(1,20000) < boundary,
            'scope':'one real smooth power-spiral square by a posteriori degree; not arbitrary Jordan'}


def _solve4_float(A,b):
    m=[list(row)+[b[i]] for i,row in enumerate(A)]
    for j in range(4):
        k=max(range(j,4),key=lambda i:abs(m[i][j]))
        if abs(m[k][j])<1e-12:return None
        m[j],m[k]=m[k],m[j]
        s=m[j][j]
        for l in range(j,5):m[j][l]/=s
        for i in range(4):
            if i==j:continue
            s=m[i][j]
            for l in range(j,5):m[i][l]-=s*m[j][l]
    return tuple(m[i][4] for i in range(4))


def smooth_spiral(t,p):
    theta=t**(-p)
    return (t*cos(theta),t*sin(theta))


def smooth_spiral_derivative(t,p):
    theta=t**(-p)
    return (cos(theta)+p*theta*sin(theta),sin(theta)-p*theta*cos(theta))


def newton_square(candidate,max_iterations=18):
    """Numerical locator ONLY; all proof later uses rational degree certificate."""
    p=float(candidate['p'])
    return refine_square(candidate, lambda t:smooth_spiral(t,p),
                         lambda t:smooth_spiral_derivative(t,p), max_iterations)


def refine_square(candidate, point, tangent, max_iterations=18):
    """Refine a chord-chart square on a differentiable parameterized arm."""
    samples=candidate['samples_data'];ids=candidate['ids'];sigma=candidate['sigma']
    ts=[float(samples[j].t)*(1-u)+float(samples[j+1].t)*u for j,u in zip(ids,candidate['parameters'])]
    for k in range(max_iterations):
        pts=[point(t) for t in ts];F=square_residual(pts,sigma)
        err=max(map(abs,F))
        if err<5e-14:return tuple(ts),err,k
        dirs=[tangent(t) for t in ts]
        zero=(0.,0.);cols=[]
        for j in range(4):
            z=[zero]*4;z[j]=dirs[j]
            cols.append(square_residual(z,sigma))
        A=tuple(tuple(cols[j][i] for j in range(4)) for i in range(4))
        dt=_solve4_float(A,tuple(-z for z in F))
        if dt is None:raise ValueError('singular Newton step')
        done=False
        for m in range(16):
            scale=2.**(-m)
            nxt=[x+scale*y for x,y in zip(ts,dt)]
            if any(not 0.05<x<1 for x in nxt):continue
            e=max(map(abs,square_residual([point(x) for x in nxt],sigma)))
            if e<err:
                ts=nxt;done=True;break
        if not done:raise ValueError('Newton line search failed')
    raise ValueError('Newton iteration cap exceeded')


def local_arcs_at_newton(candidate,half_angle=0.02):
    p=float(candidate['p']);ts,err,it=newton_square(candidate)
    if min(ts)<0.18 or any(abs(a-b)<1e-8 for i,a in enumerate(ts) for b in ts[i+1:]):
        raise ValueError('candidate has unsuitable radial geometry')
    samples=[]
    for t in ts:
        dt=half_angle/(p*t**(-p-1))
        ul=Q(round(((t-dt)**.5)*2**30),2**30)
        ur=Q(round(((t+dt)**.5)*2**30),2**30)
        for uq in (ul,ur):
            r=uq*uq;th=uq**(-int(2*Q(str(candidate['p']))))
            samples.append(Sample(r,float(r)*cos(float(th)),float(r)*sin(float(th))))
    local={'status':'found','p':candidate['p'],'inner':candidate['inner'],
           'samples_data':tuple(samples),'ids':(0,2,4,6),
           'sigma':candidate['sigma'],'newton_residual':err,
           'half_angle':half_angle,'newton_iterations':it}
    return local


def certify_discovered_square(candidate,half_angle=0.02):
    local=local_arcs_at_newton(candidate,half_angle)
    out=certify_on_smooth_spiral(local)
    out['newton_residual']=local['newton_residual']
    out['source_global_edges']=candidate['ids']
    return out

# Stored proof witnesses were discovered by the entirely separate, generic
# numerical search and Newton routine. Replaying them needs NO floating point
# search, and makes all proof inequalities depend only on Fraction arithmetic.
# u integers encode t=(u/2^30)^2 and theta=(u/2^30)^(-2p).
STORED_WITNESSES = {
    Q(2): (-1,(818567386,819260064,621854659,622029681,
               535397360,535480134,786041032,786606421)),
    Q(5,2): (-1,(855020987,855569551,706339249,706513403,
                 598989141,599053883,802365853,802740285)),
    Q(3): (1,(855337675,855702480,662377952,662438811,
              764627192,764793505,909516016,910077141)),
}


STORED_STRONG_FLOORS = {
    Q(2):(Q(1,4),Q(6,25)),
    Q(5,2):(Q(3,5),Q(3,10)),
    Q(3):(Q(7,10),Q(3,8)),
}

def replay_stored_witness(p):
    p=Q(str(p))
    sigma,us=STORED_WITNESSES[p]
    src=tuple(Sample(Q(u,2**30)**2,0.,0.) for u in us)
    row={'status':'found','p':p,'inner':None,'samples_data':src,
         'ids':(0,2,4,6),'sigma':sigma}
    side_floor,radius_floor=STORED_STRONG_FLOORS[p]
    return certify_on_smooth_spiral(row,side_floor=side_floor,radius_floor=radius_floor)


def blind_discover_and_certify(p, radius_ladder=(.55,.45,.36,.30,.22),step=.22):
    """Discovery needs no prior square position, v100 phase box, or stored witness.

    Search original power spiral's polygonal approximation using spatial hash,
    refine numerically, then prove on four actual smooth spiral arcs using degree.
    If no degree certificate is found within budget: INCONCLUSIVE (not no square).
    """
    for radius in radius_ladder:
        candidate=discover(p,radius,step,max_pairs=30000,topk=4,min_margin=.14)
        if candidate['status']!='found':continue
        try:
            for half in (.01,.005,.002):
                result=certify_discovered_square(candidate,half)
                if result['certified']:
                    result['search_pairs']=candidate['pairs']
                    result['search_segments']=candidate['samples']
                    result['searched_from_radius']=radius
                    result['local_half_angle']=half
                    return result
        except (ValueError,OverflowError,ZeroDivisionError):
            pass
    return {'certified':False,'reason':'search-budget-inconclusive','p':str(p)}


def report():
    return {'stored_exact_replay':[replay_stored_witness(p) for p in (2,Q(5,2),3)],
            'blind_discovery':[blind_discover_and_certify(p) for p in (2,Q(5,2),3)],
            'scope':'three particular non-Liouville spiral exponents; not a general Toeplitz theorem'}

if __name__=='__main__':
    import json
    out=report()
    print(json.dumps(out,indent=2))
    if not all(x['certified'] for x in out['stored_exact_replay']+out['blind_discovery']):
        raise SystemExit('FAIL: no certified square in one or more positive controls')
