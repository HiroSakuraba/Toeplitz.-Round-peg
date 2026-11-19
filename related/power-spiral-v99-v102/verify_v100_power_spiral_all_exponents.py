from decimal import Decimal as D, getcontext, localcontext, ROUND_FLOOR, ROUND_CEILING
import math, time
getcontext().prec=60

# Conservative rational decimal constants
PI_LO=D(333)/D(106)
PI_HI=D(355)/D(113)

class IV:
    __slots__=('lo','hi')
    def __init__(self,lo,hi=None):
        self.lo=D(lo); self.hi=D(lo if hi is None else hi)
        if self.lo>self.hi: raise ValueError((self.lo,self.hi))
    def __add__(self,o):
        o=o if isinstance(o,IV) else IV(o)
        with localcontext() as c:
            c.rounding=ROUND_FLOOR; lo=self.lo+o.lo
            c.rounding=ROUND_CEILING; hi=self.hi+o.hi
        return IV(lo.next_minus(),hi.next_plus())
    __radd__=__add__
    def __neg__(self): return IV(-self.hi,-self.lo)
    def __sub__(self,o): return self+(-(o if isinstance(o,IV) else IV(o)))
    def __rsub__(self,o): return IV(o)-self
    def __mul__(self,o):
        o=o if isinstance(o,IV) else IV(o)
        with localcontext() as c:
            c.rounding=ROUND_FLOOR
            ls=[self.lo*o.lo,self.lo*o.hi,self.hi*o.lo,self.hi*o.hi]
            lo=min(ls)
            c.rounding=ROUND_CEILING
            hs=[self.lo*o.lo,self.lo*o.hi,self.hi*o.lo,self.hi*o.hi]
            hi=max(hs)
        return IV(lo.next_minus(),hi.next_plus())
    __rmul__=__mul__
    def recip(self):
        if self.lo<=0<=self.hi: raise ZeroDivisionError((self.lo,self.hi))
        with localcontext() as c:
            c.rounding=ROUND_FLOOR; lo=D(1)/self.hi
            c.rounding=ROUND_CEILING; hi=D(1)/self.lo
        if self.hi<0: lo,hi=hi,lo
        return IV(min(lo,hi).next_minus(),max(lo,hi).next_plus())
    def __truediv__(self,o): return self*(o if isinstance(o,IV) else IV(o)).recip()
    def __rtruediv__(self,o): return IV(o)/self
    def __repr__(self): return f'[{self.lo},{self.hi}]'

def dec_exp_bounds(x):
    # Explicit Taylor enclosure for exp(x), 0<=x<=1.6 in this verifier.
    if x < 0:
        raise ValueError('exp positive evaluator got negative')
    term=D(1); sm=D(1); N=42
    for n in range(1,N+1):
        term=term*x/D(n); sm+=term
    first=term*x/D(N+1)
    ratio=x/D(N+2)
    upper=sm+first/(D(1)-ratio)
    return sm.next_minus(),upper.next_plus()

def exp_iv(x):
    lo,_=dec_exp_bounds(x.lo); _,hi=dec_exp_bounds(x.hi)
    return IV(lo,hi)

def sin_scalar(x):
    # |x| <= 0.8 here; use odd alternating series through x^17
    neg=x<0
    a=-x if neg else x
    s=D(0); term=a
    vals=[]
    for k in range(10):
        if k==0: term=a
        else: term = term*a*a/D((2*k)*(2*k+1))
        s = s + term if k%2==0 else s-term
        vals.append(s)
    # alternating: even k upper, odd lower for a>=0
    lo=max(vals[k] for k in range(len(vals)) if k%2==1)
    hi=min(vals[k] for k in range(len(vals)) if k%2==0)
    # widen
    lo=lo.next_minus(); hi=hi.next_plus()
    return (-hi,-lo) if neg else (lo,hi)

def cos_scalar(x):
    a=-x if x<0 else x
    s=D(0); term=D(1); vals=[]
    for k in range(10):
        if k==0: term=D(1)
        else: term=term*a*a/D((2*k-1)*(2*k))
        s=s+term if k%2==0 else s-term
        vals.append(s)
    lo=max(vals[k] for k in range(len(vals)) if k%2==1)
    hi=min(vals[k] for k in range(len(vals)) if k%2==0)
    return lo.next_minus(),hi.next_plus()

def sin_iv(x):
    # domain [0,0.8], monotone increasing
    lo,_=sin_scalar(x.lo); _,hi=sin_scalar(x.hi)
    return IV(lo,hi)

def cos_iv(x):
    # [0,0.8], decreasing
    lo,_=cos_scalar(x.hi); _,hi=cos_scalar(x.lo)
    return IV(lo,hi)

def sinc_scalar(x):
    if x==0:return D(1),D(1)
    lo,hi=sin_scalar(x)
    return (lo/x).next_minus(),(hi/x).next_plus()

def sinc_iv(x):
    # sinc decreasing on [0,0.8]
    lo,_=sinc_scalar(x.hi); _,hi=sinc_scalar(x.lo)
    return IV(lo,hi)

def coshC_scalar(x):
    # (cosh x-1)/x^2 = sum x^(2k)/(2k+2)!, positive, increasing
    a=x
    term=D(1)/D(2); s=term
    # term k=0. recurrence to k+1 multiply x^2/((2k+3)(2k+4))
    for k in range(10):
        term=term*a*a/D((2*k+3)*(2*k+4))
        s+=term
    # tail tiny; geometric upper
    r=a*a/D((23)*(24))
    upper=s + term*(a*a/D(23*24))/(D(1)-r)
    return s.next_minus(),upper.next_plus()

def coshC_iv(x):
    lo,_=coshC_scalar(x.lo); _,hi=coshC_scalar(x.hi)
    return IV(lo,hi)

def L1_scalar(z):
    # log(1+z)/z, |z|<0.7
    if z==0:return D(1),D(1)
    if z>0:
        if z < D('0.5'):
            s=D(0); vals=[]; power=D(1)
            for n in range(80):
                term=power/D(n+1)
                s=s+term if n%2==0 else s-term
                vals.append(s); power*=z
            lo=max(vals[n] for n in range(len(vals)) if n%2==1)
            hi=min(vals[n] for n in range(len(vals)) if n%2==0)
            return lo.next_minus(),hi.next_plus()
        # log(1+z)=2*atanh(z/(2+z)); positive series with small argument.
        t=z/(D(2)+z)
        sm=D(0); power=t; N=45
        for n in range(N):
            sm += power/D(2*n+1); power*=t*t
        tail=power/D(2*N+1)/(D(1)-t*t)
        lo=(D(2)*sm/z).next_minus()
        hi=(D(2)*(sm+tail)/z).next_plus()
        return lo,hi
    a=-z
    s=D(0); power=D(1)
    N=180
    for n in range(N):
        s+=power/D(n+1); power*=a
    # next n=N, tail <= a^N/(N+1)/(1-a)
    tail=power/D(N+1)/(D(1)-a)
    return s.next_minus(),(s+tail).next_plus()

def L1_iv(z):
    # decreasing for z>-1
    lo,_=L1_scalar(z.hi); _,hi=L1_scalar(z.lo)
    return IV(lo,hi)

def Efun_scalar(x):
    # (1-e^-x)/x = sum (-1)^n x^n/(n+1)!, x in [0,2]
    if x==0:return D(1),D(1)
    s=D(0); vals=[]; power=D(1)
    for n in range(36):
        term=power/D(math.factorial(n+1))
        s=s+term if n%2==0 else s-term
        vals.append(s)
        power*=x
    # terms decreasing for x<2 after n=0; here max <1.9
    lo=max(vals[n] for n in range(len(vals)) if n%2==1)
    hi=min(vals[n] for n in range(len(vals)) if n%2==0)
    return lo.next_minus(),hi.next_plus()

def Efun_iv(x):
    # decreasing in x>=0
    lo,_=Efun_scalar(x.hi); _,hi=Efun_scalar(x.lo)
    return IV(lo,hi)

def atan_pos_scalar(x):
    # x>=0 arbitrary, use reciprocal if >1
    if x==0:return D(0),D(0)
    if x>D(1):
        lo,hi=atan_pos_scalar(D(1)/x)
        return (PI_LO/D(2)-hi).next_minus(),(PI_HI/D(2)-lo).next_plus()
    s=D(0); vals=[]; power=x
    for n in range(80):
        term=power/D(2*n+1)
        s=s+term if n%2==0 else s-term
        vals.append(s)
        power*=x*x
    lo=max(vals[n] for n in range(len(vals)) if n%2==1)
    hi=min(vals[n] for n in range(len(vals)) if n%2==0)
    return lo.next_minus(),hi.next_plus()

def atan_scalar(x):
    if x>=0:return atan_pos_scalar(x)
    lo,hi=atan_pos_scalar(-x)
    return (-hi).next_minus(),(-lo).next_plus()

def atan_iv(x):
    lo,_=atan_scalar(x.lo); _,hi=atan_scalar(x.hi)
    return IV(lo,hi)

PI=IV(PI_LO,PI_HI)

def evalE(cI,dI,qI):
    x=cI*qI; y=dI*qI
    # x,y are mathematically nonnegative; outward arithmetic may create a tiny negative zero-neighborhood.
    x=IV(max(D(0),x.lo),x.hi); y=IV(max(D(0),y.lo),y.hi)
    S=sinc_iv(y); C=coshC_iv(x)
    wp=cI*cI*qI*C + dI*S
    wm=cI*cI*qI*C - dI*S
    zp=qI*wp; zm=qI*wm
    if zm.lo<=D('-0.995'):
        raise ValueError(('L1 domain',zp,zm))
    Lp=L1_iv(zp); Lm=L1_iv(zm)
    sp=(IV(1)+wp*Lp/cI)/2
    sm=(IV(1)+wm*Lm/cI)/2
    a=cI*(IV(1)-qI)
    gp=sp*Efun_iv(a*sp)/Efun_iv(a)
    gm=sm*Efun_iv(a*sm)/Efun_iv(a)
    rho=exp_iv(x)
    sy=sin_iv(y); cy=cos_iv(y)
    rap=(rho*(cy-sy)-IV(1))/(IV(1)+rho*(cy+sy))
    ram=(rho*(cy+sy)-IV(1))/(IV(1)+rho*(cy-sy))
    ap=atan_iv(rap); am=atan_iv(ram)
    Delta=3*PI-y
    Ep=D(5)/D(2)*PI + ap - Delta*gp
    Em=D(3)/D(2)*PI - am - Delta*gm
    return Ep,Em

def parts(lo,hi,n):
    lo=D(str(lo));hi=D(str(hi)); step=(hi-lo)/D(n)
    return [(lo+step*i,lo+step*(i+1)) for i in range(n)]

def cert_face(box,kind,val,target,sign,n1,n2,alpha=D(0)):
    cL,cU,dL,dU,qL,qU=map(lambda z:D(str(z)),box)
    if kind=='c':
        A=parts(dL,dU,n1); B=parts(qL,qU,n2)
    else:
        A=parts(cL,cU,n1); B=parts(qL,qU,n2)
    global_lo=D('1e99'); global_hi=D('-1e99'); bad=0; worst=None
    for aa in A:
        for bb in B:
            if kind=='c': cI=IV(val); dI=IV(*aa); qI=IV(*bb)
            else: cI=IV(*aa); dI=IV(val); qI=IV(*bb)
            Ep,Em=evalE(cI,dI,qI)
            if target=='Ep': T=Ep
            elif target=='Em': T=Em
            else: T=Ep+alpha*Em
            if T.lo<global_lo: global_lo=T.lo; worst=(aa,bb,T)
            if T.hi>global_hi: global_hi=T.hi
            if sign>0 and T.lo<=0: bad+=1
            if sign<0 and T.hi>=0: bad+=1
    return global_lo,global_hi,bad,worst



# ------------------------------------------------------------------
# v100 adaptive face certificates.
# ------------------------------------------------------------------
def adaptive_face(box,kind,val,target,sign,alpha=D(0),maxdepth=18):
    cL,cU,dL,dU,qL,qU=[D(str(x)) for x in box]
    if kind=='c': stack=[(dL,dU,qL,qU,0)]
    else: stack=[(cL,cU,qL,qU,0)]
    leaves=0; bad=[]; margin=D('1e99'); nodes=0; maxdep=0
    while stack:
        alo,ahi,qlo,qhi,dep=stack.pop(); nodes+=1; maxdep=max(maxdep,dep)
        if kind=='c': cI=IV(val); dI=IV(alo,ahi)
        else: cI=IV(alo,ahi); dI=IV(val)
        Ep,Em=evalE(cI,dI,IV(qlo,qhi))
        T=Ep if target=='Ep' else Em if target=='Em' else Ep+alpha*Em
        ok=T.lo>0 if sign>0 else T.hi<0
        if ok:
            leaves+=1
            mm=T.lo if sign>0 else -T.hi
            if mm<margin: margin=mm
            continue
        if dep>=maxdepth:
            bad.append((alo,ahi,qlo,qhi,T.lo,T.hi,dep)); continue
        aw=ahi-alo; qw=qhi-qlo
        fulla=(dU-dL) if kind=='c' else (cU-cL); fullq=qU-qL
        if aw/fulla >= qw/fullq:
            mid=(alo+ahi)/2
            stack.append((alo,mid,qlo,qhi,dep+1)); stack.append((mid,ahi,qlo,qhi,dep+1))
        else:
            mid=(qlo+qhi)/2
            stack.append((alo,ahi,qlo,mid,dep+1)); stack.append((alo,ahi,mid,qhi,dep+1))
    return {'nodes':nodes,'leaves':leaves,'bad':len(bad),'margin':margin,'maxdepth':maxdep,'bad_examples':bad[:3]}

# Import and execute v99 exact small-p verifier as a regression/theorem input.
from pathlib import Path as _Path
import runpy as _runpy
HERE=_Path(__file__).resolve().parent
v99=_runpy.run_path(str(HERE/'verify_v99_small_exponent_log_spiral.py'))

# Box A over-covers q in [0,1/3].
BOX_A=('1','1.9','0.45','0.72','0','0.333333333333333334')
TESTS_A=[
    ('A_c_low','c','1','Em',+1,D(0)),
    ('A_c_high','c','1.9','Em',-1,D(0)),
    ('A_d_low','d','0.45','Ep',+1,D(0)),
    ('A_d_high','d','0.72','Ep',-1,D(0)),
]

# Box B over-covers q in [1/3,200/201].
BOX_B=('0.75','1.5','0.4','0.72','0.333333333333333333','0.995024875621890548')
ALPHA=D(57)/D(400)
TESTS_B=[
    ('B_c_low','c','0.75','Em',+1,D(0)),
    ('B_c_high','c','1.5','Em',-1,D(0)),
    ('B_d_low','d','0.4','G',+1,ALPHA),
    ('B_d_high','d','0.72','G',-1,ALPHA),
]

summaries=[]
for name,kind,val,target,sign,alpha in TESTS_A:
    r=adaptive_face(BOX_A,kind,val,target,sign,alpha,maxdepth=18)
    assert r['bad']==0, (name,r)
    summaries.append((name,r))
for name,kind,val,target,sign,alpha in TESTS_B:
    r=adaptive_face(BOX_B,kind,val,target,sign,alpha,maxdepth=18)
    assert r['bad']==0, (name,r)
    summaries.append((name,r))

print('v100 adaptive enclosure summaries:')
for name,r in summaries:
    print(name,'nodes=',r['nodes'],'leaves=',r['leaves'],'maxdepth=',r['maxdepth'],'certified_margin=',r['margin'])
print('v100 Box A covers 0 <= q <= 1/3 and proves a square for every p >= 2.')
print('v100 Box B covers 1/3 <= q <= 200/201 and proves a square for 1/200 <= p <= 2.')
print('v100 v99 covers 0 < p <= 1/200.')
print('v100 conclusion: every power spiral t*exp(i/t^p), p>0, contains a non-graceful 0132 square.')
print('v100 verifier: PASS')
