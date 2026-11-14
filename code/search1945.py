"""Negative-information update from the documented 6-10 Dec 1945 search boxes (NARA M1657 CG 7ND log).

For each simulated ditch point:
  1. drift a floating-evidence marker (raft/debris/slick) from ditch time to each search day
     (Gulf Stream inside its axis band, weak NE ambient current elsewhere, NW-wind leeway after the front);
  2. for every box searched that day that contains the drifted marker, multiply the weight by
     (1 - q_day * (1 - exp(-C))):  q_day = chance something detectable was still afloat, C = coverage.
Run with SCEN=low|mid|high to scale q.
"""
import sys, numpy as np
from shapely.geometry import shape, box as sbox
from shapely import contains_xy
from shapely.ops import unary_union
import json

rng = np.random.default_rng(7)
SCEN = sys.argv[1] if len(sys.argv) > 1 else "mid"
QSCALE = {"low": 0.5, "mid": 1.0, "high": 1.8}[SCEN]
Q = {6: 0.40, 7: 0.30, 8: 0.20, 9: 0.15, 10: 0.10}          # P(detectable evidence afloat)
SEA6 = 0.6                                                    # heavy seas 6 Dec reduce coverage
C0 = -82.0                                                    # "coast" (land cut separately)

# (day, lat_s, lat_n, lon_w, lon_e, coverage C, source)
BOX = [
 (6, 29.0, 30.5, C0, -78.5, .8*SEA6, "Jax NAOTC coast-78-30W"), (6, 28, 29, C0, -79, .8*SEA6, "Daytona"),
 (6, 27, 28, C0, -79, .8*SEA6, "Vero"), (6, 27, 29, -79, -77, .8*SEA6, "Banana River"),
 (6, 25, 27, C0, -78, .8*SEA6, "Ft Lauderdale (38-40 planes)"), (6, 25, 27, -78, -77, .5*SEA6, "CGAS Miami"),
 (6, 25, 28, -77, -75, .5*SEA6, "Master Field"), (6, 24, 25, -76, -73, .5*SEA6, "RAF Nassau"),
 (6, 29.65, 29.85, -78, -77, .3, "Solomons surface track (no flying wx)"),
 (7, 29.5, 31.5, -78, -76, .6, "Jax PBYs"), (7, 28.5, 29.5, -78, -77, .6, "Daytona"),
 (7, 30.0, 31.5, -79, -77, .8, "USS Solomons air group - 'completed, negative'"),
 (7, 29.5, 31.5, -79, -77, .6, "Banana River"), (7, 26.5, 28, C0, -77.5, .8, "Ft Lauderdale"),
 (7, 26, 27, -77, -76, .3, "Miami ATC (1 PBY, fuel-limited)"), (7, 25, 26, -78, -76, .5, "RAF Nassau"),
 (7, 28, 30, -79, -77, .6, "NAS Opa Locka"), (7, 28, 29.5, -77, -76, .5, "Miami PBM"),
 (7, 27, 28, -77, -76, .5, "CGAS Miami"), (7, 27, 28.5, -78, -77, .5, "Vero Dumbo"),
 (7, 30, 31.5, C0, -79, .8, "Mayport/Jax/Cecil"),
 (8, 28.5, 30, -76, -75, .6, "Jax PBYs"), (8, 28.5, 30, -74, -73, .6, "Banana River PBYs"),
 (8, 28.5, 31, C0, -79, .8, "Jax land planes"), (8, 25.5, 26, -76, -75.5, .5, "RAF Nassau"),
 (8, 28.3, 28.7, -75.7, -75.3, .5, "CGC Pandora square search"),
 (9, 27.8, 28.6, -76.8, -75.6, .8, "raft-report concentration (army planes, ships)"),
 (10, 27.5, 29.5, -78, -74, .6, "'concentrating search' 74-78W 27-30/29-30N"),
 (10, 28.92, 31.08, -73.58, -72.5, .6, "Solomons east of 73-35W"),
]

F = [("out_ships_isl_rec_bax_rad_fw.npz", .5*.65), ("out_ships_isl_rec_bax_rad_zz_fw.npz", .5*.35),
     ("out_ships_isl_rec_bax_rad_fe.npz", .5*.65), ("out_ships_isl_rec_bax_rad_zz_fe.npz", .5*.35)]
import os
if os.environ.get("MODE")=="final":
    F = [("out_ships_isl_rec_bax_rad_fe.npz", .5), ("out_ships_isl_rec_bax_rad_zz_fe.npz", .5)]
L = [np.load(f) for f, _ in F] if os.environ.get("MODE")!="v3" else []
CUT = 1160 if os.environ.get("MODE")=="final" else 1e9
if os.environ.get("MODE")=="v3":
    _v=np.load("work/final_v3_prior.npz"); lat=_v["lat"]; lon=_v["lon"]
else:
    lat = np.concatenate([x["lat"] for x in L]); lon = np.concatenate([x["lon"] for x in L])
if os.environ.get("MODE")=="v3":
    tf = np.full(len(lat), 1152.0); W = _v["w"]
else:
  tf = np.concatenate([x["t_fuel"] for x in L])
  W = np.concatenate([(x["w"]*(x["t_fuel"]<=CUT)) / (x["w"]*(x["t_fuel"]<=CUT)).sum() * p for x, (_, p) in zip(L, F)])
n = len(lat)

# drift: hourly steps from ditch to noon of each day
la, lo = lat.copy(), lon.copy()
t = tf / 60.0  # minutes->hours                       # hours after midnight 5 Dec
lee_dir = np.radians(140 + rng.normal(0, 30, n)); lee_k = rng.uniform(0.5, 1.5, n) * 0.6   # kt, wind from NW
pos = {}
for day in range(6, 11):
    tgt = (day - 5) * 24 + 12
    while True:
        go = t < tgt
        if not go.any(): break
        dt = np.minimum(1.0, tgt - t) * go
        axis = -79.6 + np.clip(la - 27, 0, None) * 0.45      # Gulf Stream axis veers NE north of 27N
        ingulf = np.abs(lo - axis) * 60 * np.cos(np.radians(la)) < 25
        cu = np.where(ingulf, 2.5, 0.3); cd = np.radians(np.where(ingulf, 15 + np.clip(la - 29, 0, None) * 12, 30))
        vn = cu * np.cos(cd) + lee_k * np.cos(lee_dir); ve = cu * np.sin(cd) + lee_k * np.sin(lee_dir)
        la = la + vn * dt / 60; lo = lo + ve * dt / 60 / np.cos(np.radians(la)); t = t + dt
    pos[day] = (la.copy(), lo.copy())

lik = np.ones(n); hit = {}
for day, s, nn, w_, e, C, src in BOX:
    a, o = pos[day]; m = (a >= s) & (a <= nn) & (o >= w_) & (o <= e)
    pod = min(0.95, QSCALE * Q[day]) * (1 - np.exp(-C))
    lik[m] *= 1 - pod; hit[src + f" ({day} Dec)"] = (W[m].sum(), pod)
W2 = W * lik; W2 /= W2.sum()
print(f"[{SCEN}] prob. mass removed by 1945 search: {1 - (W*lik).sum()/W.sum():.3f}")
for k, (pm, pod) in sorted(hit.items(), key=lambda x: -x[1][0] * x[1][1])[:8]:
    print(f"  {k:60s} prior-in-box {pm:.3f}  POD {pod:.2f}")

# land mask for areas, HPD stats
land = unary_union([shape(f["geometry"]) for f in json.load(open("inputs/ne_50m_land.geojson"))["features"]]).intersection(sbox(-84, 24, -69, 34))
def hpd(w):
    xe = np.linspace(-84, -69, 301); ye = np.linspace(24.5, 33, 171)
    H, _, _ = np.histogram2d(lon, lat, bins=[xe, ye], weights=w); H = H.T
    yc = (ye[:-1] + ye[1:]) / 2
    cell = (0.05 * 60 * np.cos(np.radians(yc)))[:, None] * (0.05 * 60) * np.ones_like(H)
    o = np.argsort(H.ravel())[::-1]; c = np.cumsum(H.ravel()[o]); ca = np.cumsum(cell.ravel()[o])
    a50 = ca[np.searchsorted(c, .5)]; a90 = ca[np.searchsorted(c, .9)]
    best1k = c[np.searchsorted(ca, 1000)]
    return a50, a90, best1k, H, xe, ye
for tag, w in [("before", W), ("after", W2)]:
    a50, a90, b, *_ = hpd(w)
    print(f"  {tag:6s}: centre {np.average(lat,weights=w):.2f}N {np.average(lon,weights=w):.2f}W | 50% {a50:,.0f} sq nm | 90% {a90:,.0f} | best 1,000 sq nm holds {b:.1%} | west of 76.8W {w[lon<-76.8].sum():.1%}")
np.savez(f"work/search_{SCEN}{'_'+os.environ.get('MODE','') if os.environ.get('MODE') else ''}.npz", lat=lat, lon=lon, w=W2, w0=W)
