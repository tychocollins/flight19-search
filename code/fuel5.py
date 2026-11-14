"""Five-plane minimum-order fuel-out time. Inputs as ranges (unverified handbook figures treated as priors)."""
import numpy as np, sys
rng=np.random.default_rng(3); N=200000
burn_scale=float(sys.argv[1]) if len(sys.argv)>1 else 1.0
def plane(lead, cap_lo, cap_hi, eng):
    F=rng.uniform(cap_lo,cap_hi,N)
    F-=rng.uniform(12,16,N)*(1.1 if eng else 1)                  # start/taxi/climb 14:10-14:25
    bomb=rng.uniform(85,100,N)*0.5; F-=bomb                         # bombing ~0.5 h
    cruise=rng.uniform(52,60,N)*(1.04 if eng else 1)*burn_scale     # gph auto-lean
    F-=cruise*2.55                                                  # cruise to 17:28
    night=cruise*(1.0 if lead else rng.uniform(1.08,1.15,N))        # formation penalty after sunset
    unusable=rng.uniform(3,5,N); gauge=rng.normal(0,3,N)
    usable_to_10=F-unusable-(10+gauge)                              # gallons until gauge reads 10
    return 17+28/60+usable_to_10/night                              # clock hours
T=np.vstack([plane(True,320,335,True)]+[plane(False,330,340,False) for _ in range(4)])
t=T.min(0); first=T.argmin(0)
q=np.percentile(t,[5,25,50,75,95])
f=lambda h:f"{int(h)}:{int(round((h%1)*60)):02d}"
print(f"burn x{burn_scale}: ditch time 5/25/50/75/95% = "+" ".join(f(x) for x in q), f"| lead plane first in {np.mean(first==0):.0%}")
print(f"  P(19:05-19:35) {np.mean((t>=19+5/60)&(t<=19+35/60)):.0%}  P(<19:05) {np.mean(t<19+5/60):.0%}  P(>19:35) {np.mean(t>19+35/60):.0%}")
