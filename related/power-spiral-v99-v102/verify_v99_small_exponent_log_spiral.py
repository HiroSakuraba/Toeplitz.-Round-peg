"""
v99 exact-rational verifier.

Checks:
  A. p=0/logarithmic-spiral corner sign bounds on the explicit box;
  B. exact rational margins for the four monotonicity lemmas;
  C. p-derivative bound and persistence through 0 <= p <= 1/200;
  D. numerical locator for the limiting logarithmic square (not theorem input).

No interval library is used.  The proof bounds use exact Fraction arithmetic,
alternating Taylor inequalities, and rational bounds on pi.
"""
from fractions import Fraction as F
from pathlib import Path
import math
import mpmath as mp

PI_LO=F(333,106)
PI_HI=F(355,113)

rlo=F(119,50)
rhi=F(62,25)
elo=F(27,50)
ehi=F(61,100)

# ------------------------------------------------------------
# Exact elementary bound helpers.
# ------------------------------------------------------------
def sin_bounds(x,n=5):
    vals=[]
    s=F(0)
    for k in range(n+1):
        term=x**(2*k+1)/F(math.factorial(2*k+1))
        s += term if k%2==0 else -term
        vals.append(s)
    lo=max(vals[k] for k in range(len(vals)) if k%2==1)
    hi=min(vals[k] for k in range(len(vals)) if k%2==0)
    return lo,hi

def cos_bounds(x,n=5):
    vals=[]
    s=F(0)
    for k in range(n+1):
        term=x**(2*k)/F(math.factorial(2*k))
        s += term if k%2==0 else -term
        vals.append(s)
    lo=max(vals[k] for k in range(len(vals)) if k%2==1)
    hi=min(vals[k] for k in range(len(vals)) if k%2==0)
    return lo,hi

def atanh_bounds(x,n=14):
    s=F(0)
    for k in range(n+1):
        s += x**(2*k+1)/F(2*k+1)
    tail=x**(2*n+3)/(F(2*n+3)*(1-x*x))
    return s,s+tail

def log_bounds(x,n=18):
    z=(x-1)/(x+1)
    lo,hi=atanh_bounds(z,n)
    return 2*lo,2*hi

def atan01_bounds(x,n=80):
    vals=[]
    s=F(0)
    for k in range(n+1):
        term=x**(2*k+1)/F(2*k+1)
        s += term if k%2==0 else -term
        vals.append(s)
    lo=max(vals[k] for k in range(len(vals)) if k%2==1)
    hi=min(vals[k] for k in range(len(vals)) if k%2==0)
    return lo,hi

def atan_bounds(x,n=80):
    if x<=1:
        return atan01_bounds(x,n)
    loz,hiz=atan01_bounds(1/x,n)
    return PI_LO/2-hiz,PI_HI/2-loz

# ------------------------------------------------------------
# p=0 corner bounds.
# ------------------------------------------------------------
def corner_bounds(r,e):
    sinlo,sinhi=sin_bounds(e)
    coslo,coshi=cos_bounds(e)

    u_lo=2*r*sinlo/(r*r+1)
    u_hi=2*r*sinhi/(r*r+1)
    atlo=atanh_bounds(u_lo)[0]
    athi=atanh_bounds(u_hi)[1]

    logrl,logrh=log_bounds(r)

    barg_lo=(r*r-1)/(2*r*coshi)
    barg_hi=(r*r-1)/(2*r*coslo)
    blo=atan_bounds(barg_lo)[0]
    bhi=atan_bounds(barg_hi)[1]

    Dlo=3*PI_LO-e
    Dhi=3*PI_HI-e

    H1lo=PI_LO+blo-Dhi*athi/logrl
    H1hi=PI_HI+bhi-Dlo*atlo/logrh

    A=(1+r*r)/2
    Blo=A*A-r*r*sinhi*sinhi
    Bhi=A*A-r*r*sinlo*sinlo

    logBlo=log_bounds(Blo)[0]
    logBhi=log_bounds(Bhi)[1]

    s2lo,s2hi=sin_bounds(2*e,7)
    c2lo,c2hi=cos_bounds(2*e,7)

    ylo=r*r*s2lo/(1+r*r*c2hi)
    yhi=r*r*s2hi/(1+r*r*c2lo)

    philo=atan_bounds(ylo)[0]
    phihi=atan_bounds(yhi)[1]

    H2lo=2*PI_LO-Dhi*logBhi/(4*logrl)-phihi/2
    H2hi=2*PI_HI-Dlo*logBlo/(4*logrh)-philo/2

    return H1lo,H1hi,H2lo,H2hi

LL=corner_bounds(rlo,elo)
UR=corner_bounds(rhi,ehi)

assert LL[0] > F(46,1000)
assert LL[2] > F(18,1000)
assert UR[1] < -F(37,1000)
assert UR[3] < -F(19,1000)

# ------------------------------------------------------------
# Monotonicity quantitative margins.
#
# H1_r > 0 is structural:
# beta_r>0 and atanh(u)/log r decreases for r>1.
#
# H1_e < 0:
# Delta F_e - F - beta_e > 0.
# ------------------------------------------------------------
sElo,sEhi=sin_bounds(ehi)
cElo,cEhi=cos_bounds(ehi)

cpmin=2*rhi*cElo/(rhi*rhi+1)
logrhi_lo,logrhi_hi=log_bounds(rhi)
Fe_low=cpmin/logrhi_hi
Dmin=3*PI_LO-ehi

u_max=2*rlo*sEhi/(rlo*rlo+1)
F_upper=atanh_bounds(u_max)[1]/log_bounds(rlo)[0]
betae_upper=(sEhi/cElo)/2

H1e_margin=Dmin*Fe_low-F_upper-betae_upper
assert H1e_margin > F(4)

# ------------------------------------------------------------
# H2_r < 0:
# G_r>0 and phi_r>0.
# Prove numerator of G_r is positive.
# ------------------------------------------------------------
Brmin=rlo*(rlo*rlo+1-2*sEhi*sEhi)

slo_lo,slo_hi=sin_bounds(elo)
Ahi=(1+rhi*rhi)/2
Bmax=Ahi*Ahi-rhi*rhi*slo_lo*slo_lo

Lmin=log_bounds(rlo)[0]
lhs=Lmin*Brmin/Bmax
rhs=log_bounds(Bmax)[1]/rlo
Gr_num_margin=lhs-rhs
assert Gr_num_margin > F(9,100)

# ------------------------------------------------------------
# H2_e > 0:
# G - Delta G_e - phi_e/2 > 0.
# Here G_e<0 and phi_e/2<1.
# ------------------------------------------------------------
Alo=(1+rlo*rlo)/2
Bmin=Alo*Alo-rlo*rlo*sEhi*sEhi
G_lower=log_bounds(Bmin)[0]/(4*logrhi_hi)

celo_lo,celo_hi=cos_bounds(elo)
Beabs_min=2*rlo*rlo*slo_lo*celo_lo
Geabs_low=Beabs_min/(4*Bmax*logrhi_hi)

H2e_margin=G_lower+Dmin*Geabs_low-1
assert H2e_margin > F(69,100)

# ------------------------------------------------------------
# Positive-p persistence.
# Universal bound from v98:
# |dH_i/dp| <= Delta log(rho)/4.
# Use Delta<3*pi and log rho < rho-1.
# ------------------------------------------------------------
Cbound=3*PI_HI*(rhi-1)/4
assert Cbound==F(7881,2260)

pmax=F(1,200)
drift=Cbound*pmax
assert drift==F(7881,452000)

margins=[
    F(46,1000)-drift,
    F(37,1000)-drift,
    F(18,1000)-drift,
    F(19,1000)-drift,
]
assert margins[2]==F(51,90400)
assert min(margins)>0

print("v99 LL H1 lower:",float(LL[0]))
print("v99 LL H2 lower:",float(LL[2]))
print("v99 UR H1 upper:",float(UR[1]))
print("v99 UR H2 upper:",float(UR[3]))
print("v99 H1_e negativity margin:",float(H1e_margin))
print("v99 G_r numerator positivity margin:",float(Gr_num_margin))
print("v99 H2_e positivity margin:",float(H2e_margin))
print("v99 exponent derivative bound:",Cbound)
print("v99 retained margins at p=1/200:",[float(x) for x in margins])
print("v99 exact p=0 Miranda and small-p persistence audit: PASS")

# ------------------------------------------------------------
# Numerical locator only.
# ------------------------------------------------------------
mp.mp.dps=70

def H0(rho,eps):
    D=3*mp.pi-eps
    A=(1+rho*rho)/2
    S=rho*mp.sin(eps)
    Rp=mp.sqrt(A+S)
    Rm=mp.sqrt(A-S)
    z2=-rho*mp.e**(1j*eps)
    pp=(1+1j)/2
    qq=(1-1j)/2
    Qp=qq+pp*z2
    Qm=pp+qq*z2
    k=D/mp.log(rho)
    fp=-k*mp.log(Rp)-mp.arg(Qp)
    fm=-k*mp.log(Rm)-mp.arg(Qm)
    return fp-fm,(fp+fm)/2+2*mp.pi

rho0,eps0=mp.findroot(H0,(mp.mpf("2.429"),mp.mpf("0.57425")))
k0=(3*mp.pi-eps0)/mp.log(rho0)

print("v99 numerical log-spiral rho0:",mp.nstr(rho0,45))
print("v99 numerical log-spiral eps0:",mp.nstr(eps0,45))
print("v99 numerical log-spiral pitch k0:",mp.nstr(k0,45))
