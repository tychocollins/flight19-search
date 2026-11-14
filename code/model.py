"""Flight 19 Monte Carlo dead-reckoning reconstruction.

Inputs (sources in comments):
- Departure 14:10 from 26°03'N 80°07'W; planned legs 091°T 56 nm to Hens & Chickens Shoals,
  091°T 67 nm, 346°T 73 nm, 241°T 120 nm (Board of Investigation, Finding 2; NHHC H-Gram 057).
- Bombing ~30 min at the shoals; bomb-drop permission ~15:00.
- ~15:40 first trouble ("lost after that last turn"); Cox advises fly north.
- ~16:45 Taylor: 030° for 45 min; ~7 min later changes to 090° ("for 10 minutes", ~16:56+).
- ~17:24-17:36 Taylor: fly 270° until landfall or out of fuel.
- 18:04 Taylor: turns back to 090°.
- 17:50 HF/DF fix: within 100-mile radius of 29°15'N 79°00'W (Finding 28).
- Fuel: full tanks, ~5.5 h endurance (NHHC) -> ~19:40; one account estimates fuel to ~20:00.
- TBM-3 cruise ~145-153 mph (~126-133 kt); route timing implies ~150 kt for the exercise.
- Wind: surface 22 kt at 16:00 (NHHC); an account cites 30-35 kt aloft from the west; a front
  arriving in the evening (high seas reported ~18:20).
All times local (R zone, UTC-5), in minutes after midnight.
"""
import sys, json, numpy as np

R_NM = 3440.065
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
SCEN = sys.argv[1] if len(sys.argv) > 1 else "base"
TAG = SCEN + ("_isl" if "--island" in sys.argv else "") + ("_e270" if "--early270" in sys.argv else "") + ("_rec" if "--record" in sys.argv else "") + ("_bax" if "--baxter" in sys.argv else "") + ("_rad" if "--radio" in sys.argv else "") + ("_zz" if "--zigzag" in sys.argv else "") + ("_fw" if "--fuelwide" in sys.argv else "") + ("_fe" if "--fuelearly" in sys.argv else "") + __import__("os").environ.get("TAGX", "")
N = 200_000

def T(h, m): return h * 60 + m

import os as _osm
if _osm.environ.get("MAGVAR"):
    # 1945.9 magnetic declination (IGRF via ppigrf), east positive; headings flown are magnetic
    import ppigrf as _pg, datetime as _dt
    _glat = np.arange(24, 34.01, 1.0); _glon = np.arange(-84, -68.99, 1.0)
    _G = np.zeros((len(_glat), len(_glon)))
    for _i, _a in enumerate(_glat):
        for _j, _o in enumerate(_glon):
            _e, _n, _u = _pg.igrf(_o, _a, 0, _dt.datetime(1945, 12, 5)); _G[_i, _j] = np.degrees(np.arctan2(np.ravel(_e)[0], np.ravel(_n)[0]))
    from scipy.interpolate import RegularGridInterpolator as _RGI
    _decl = _RGI((_glat, _glon), _G, bounds_error=False, fill_value=None)
else:
    _decl = None
def step(lat, lon, hdg_deg, tas, wdir_from, wspd, dt_min):
    """Advance positions (arrays) by dt minutes. Heading magnetic if MAGVAR set (else true), wind FROM direction (true)."""
    if _decl is not None: hdg_deg = hdg_deg + _decl(np.c_[lat, lon])
    h = np.radians(hdg_deg); w = np.radians(wdir_from + 180.0)  # wind TO direction
    vn = tas * np.cos(h) + wspd * np.cos(w)
    ve = tas * np.sin(h) + wspd * np.sin(w)
    dn = vn * dt_min / 60.0; de = ve * dt_min / 60.0
    lat2 = lat + np.degrees(dn / R_NM)
    lon2 = lon + np.degrees(de / (R_NM * np.cos(np.radians(lat))))
    return lat2, lon2

def hav(lat1, lon1, lat2, lon2):
    p = np.pi / 180
    a = np.sin((lat2 - lat1) * p / 2) ** 2 + np.cos(lat1 * p) * np.cos(lat2 * p) * np.sin((lon2 - lon1) * p / 2) ** 2
    return 2 * R_NM * np.arcsin(np.sqrt(a))

U = lambda a, b: rng.uniform(a, b, N)
Nrm = lambda m, s: rng.normal(m, s, N)

# --- per-run parameters ---
tas = U(125, 160)
wdir = U(230, 300); wspd = U(15, 35)           # aloft wind FROM W-SW
if SCEN == "calm": wspd = U(0, 10)
if SCEN in ("exhibit", "exhibit_decay"):
    # NAS Fort Lauderdale winds aloft, 1300 EWT (1200 EST) 5 Dec 1945 (Board exhibit):
    # 2000 ft 242/32, 3000 240/36, 4000 246/33, 5000 245/34 kt. Measurement + spatial/temporal uncertainty.
    wdir = Nrm(243, 8); wspd = Nrm(34, 4)
    if SCEN == "exhibit_decay":
        # surface winds at Ft Lauderdale fell from ~20 to 6 after 1650 ("DRK NW"); allow aloft wind to drop/veer later
        wdir2 = U(230, 280); wspd2 = U(15, 35); t_wx = U(T(16, 30), T(17, 30))
if SCEN == "ships":
    # measured aloft wind (NAS Ft Lauderdale 1200 EST) until the evening, then ICOADS ship-obs wind field
    wdir = Nrm(243, 8); wspd = Nrm(34, 4)
    import pandas as _pd
    _o = _pd.read_csv("icoads/near.csv")
    _o["t"] = (_o.day - 5) * 1440 + _o.hr_utc * 60 - 300        # minutes after 00:00 EST 5 Dec
    _o = _o[(_o.t >= T(14, 0)) & (_o.t <= T(23, 0)) & _o.wdir_deg.notna() & _o.w_ms.notna()]
    _o = _o[_o.w_ms * 1.944 > 2]
    OB = _o[["lat", "lon", "t", "wdir_deg"]].to_numpy(float); OBs = (_o.w_ms * 1.944).to_numpy(float)
    # extra land/coast obs: NAS Ft Lauderdale surface (Board exhibit) 1650 W/15, 1750 W/6, 1850 W/6 (arrow ~west)
    OB = np.vstack([OB, [[26.07, -80.15, T(16, 50), 255], [26.07, -80.15, T(17, 50), 255], [26.07, -80.15, T(18, 50), 270]]])
    OBs = np.concatenate([OBs, [15, 6, 6]])
    t_wx = U(T(16, 15), T(17, 15))                    # hand-over from measured aloft wind to ship field
    import os as _os
    if _os.environ.get("HOLD"): t_wx = np.full(N, 1e9)   # sensitivity: measured SW aloft wind held all evening
    veer = Nrm(15, 8); fac = Nrm(float(_os.environ.get("FAC", 1.3)), 0.15)            # surface -> 2-5k ft: veer and strengthen
    dnoise = Nrm(0, 15)
    def ship_wind(la, lo, t):
        dlat = la[:, None] - OB[None, :, 0]; dlon = (lo[:, None] - OB[None, :, 1]) * np.cos(np.radians(29))
        dist = np.sqrt(dlat ** 2 + dlon ** 2) * 60 + np.abs(t - OB[None, :, 2]) / 60 * 40 + 20
        wgt = 1 / dist ** 2
        th = np.radians(OB[:, 3]); uo = -OBs * np.sin(th); vo = -OBs * np.cos(th)
        u = (wgt * uo).sum(1) / wgt.sum(1); v = (wgt * vo).sum(1) / wgt.sum(1)
        spd = np.hypot(u, v); dfrom = (np.degrees(np.arctan2(-u, -v)) + 360) % 360
        return (dfrom + veer + dnoise) % 360, spd * fac
if SCEN == "era5":
    # ERA5 (ECMWF reanalysis, Google ARCO copy) 900 hPa u/v, 1945-12-05 19Z..12-06 01Z, 1-degree grid 25-33N 82-73W (era5/era5_900.npz)
    _E = np.load("era5/era5_900.npz")
    from scipy.interpolate import RegularGridInterpolator as _RGI5
    _tt = (_E["hour_utc"] - 5) * 60.0; _la = _E["lat"][::-1]; _lo = _E["lon"]
    _Ui = _RGI5((_tt, _la, _lo), _E["u"][:, ::-1, :], bounds_error=False, fill_value=None)
    _Vi = _RGI5((_tt, _la, _lo), _E["v"][:, ::-1, :], bounds_error=False, fill_value=None)
    efac = Nrm(1.0, 0.15); edir = Nrm(0, 12)          # 1945 reanalysis uncertainty (per simulated flight)
    def era5_wind(la, lo, t):
        tc = np.clip(t, _tt[0], _tt[-1]); lac = np.clip(la, 25, 33); loc = np.clip(lo, -82, -73)
        pts = np.c_[np.full(len(la), tc), lac, loc]
        u = _Ui(pts); v = _Vi(pts)
        return ((np.degrees(np.arctan2(-u, -v)) + edir) % 360), np.hypot(u, v) * efac
if SCEN == "wx1945":
    # 1945-12-05 01:30 EST Weather Bureau map: post-frontal W-NW winds over Florida (front ~78.5W);
    # 1945-12-06 01:30 map: 999 mb low near 30.5N 78.5W -> SW flow on its SE side in the evening.
    wdir = U(255, 315); wspd = U(15, 30)
    wdir2 = U(200, 250); wspd2 = U(25, 40); t_wx = U(T(16, 45), T(17, 30))
hbias = Nrm(0, 4)                                # systematic compass/variation error
seg_noise = lambda: Nrm(0, 6)
dt = 1.0

t_bomb_end = U(T(14, 55), T(15, 10))
t_030 = Nrm(T(16, 45), 6)
t_090a = t_030 + 7 + Nrm(0, 2)
t_270 = U(T(17, 5), T(17, 36))                    # west leg begins somewhere in this window
if "--record" in sys.argv:
    # NARA M1657: Board finding 33 (270 assumed prior to 1715R); Pensacola HF/DF log 2219Z "flying 270";
    # finding 35 (1806R suggestion to turn east); finding 37 (landings after 1904R); finding 39 (fuel to ~2000R)
    t_270 = np.maximum(t_090a + 3, U(T(16, 45), T(17, 15)))
    t_east = Nrm(T(18, 6), 2)
    t_fuel = np.clip(Nrm(T(19, 55), 12), T(19, 5), T(20, 15))
    if "--fuelearly" in sys.argv:
        # Houma log (Exh. 23): 1802R "we may have to ditch any minute"; still heard 1832 (Houma), 1844 (PE), 1904 (Board)
        t_fuel = U(T(19, 5), T(19, 35))
    elif "--fuelwide" in sys.argv:
        t_fuel = U(T(19, 5), T(20, 5))   # last heard 19:04 (finding 37) .. Board ~20:00; ASRTU-4 estimate 19:30
elif "--early270" in sys.argv:
    t_270 = np.maximum(t_090a + 3, U(T(16, 40), T(17, 36)))
if "--record" not in sys.argv:
    t_east = Nrm(T(18, 4), 3)
    t_fuel = np.clip(Nrm(T(19, 40), 15), T(19, 10), T(20, 15))

# planned-route pieces, in time
v_out = tas  # leg flown at exercise speed
t_shoal = T(14, 10) + 56 / v_out * 60
t_leg1b_end = t_bomb_end + 67 / v_out * 60
t_leg2_end = t_leg1b_end + 73 / v_out * 60

# confusion phase headings (from end of planned leg 2 or first trouble until t_030):
# Cox advised "sun on port wing, fly north"; model two sub-segments around north.
conf1 = Nrm(0, 35) % 360; conf2 = Nrm(10, 35) % 360

lat = np.full(N, 26 + 3 / 60); lon = np.full(N, -(80 + 7 / 60))
fix_lat, fix_lon = 29.25, -79.0
lat1750 = np.full(N, np.nan); lon1750 = np.full(N, np.nan)
track_idx = rng.choice(N, 60, replace=False)
_bd = np.load("bahamas_dist.npz"); _bl, _bo, _BD = _bd["lats"], _bd["lons"], _bd["D"]
def dist_bahamas(la, lo):
    i = np.clip(np.round((la - _bl[0]) / 0.05).astype(int), 0, len(_bl) - 1)
    j = np.clip(np.round((lo - _bo[0]) / 0.05).astype(int), 0, len(_bo) - 1)
    return _BD[i, j]
dmin_island = np.full(N, 1e9)
dmin_walker = np.full(N, 1e9); d_canav_1800 = np.full(N, np.nan)
d_pe_1733 = np.full(N, np.nan); d_pe_1850 = np.full(N, np.nan)
USE_ISLAND = "--island" in sys.argv
tracks = []

b1, b2 = conf1 + hbias, conf2 + hbias
n1, n2, n3, n4, n5, n6 = (seg_noise() for _ in range(6))
t = T(14, 10)
traj = [[] for _ in track_idx]
while True:
    hd = np.full(N, np.nan)
    hd = np.where(t < t_shoal, 91 + hbias, hd)
    bomb = (t >= t_shoal) & (t < t_bomb_end)            # orbiting the shoals
    hd = np.where((t >= t_bomb_end) & (t < t_leg1b_end), 91 + hbias + n1, hd)
    hd = np.where((t >= t_leg1b_end) & (t < t_leg2_end), 346 + hbias + n2, hd)
    mid = (t_leg2_end + t_030) / 2
    hd = np.where((t >= t_leg2_end) & (t < mid), b1, hd)
    hd = np.where((t >= mid) & (t < t_030), b2, hd)
    hd = np.where((t >= t_030) & (t < t_090a), 30 + hbias + n3, hd)
    hd = np.where((t >= t_090a) & (t < t_270), 90 + hbias + n4, hd)
    hd = np.where((t >= t_270) & (t < t_east), 270 + hbias + n5, hd)
    if "--zigzag" in sys.argv:
        ph = ((t - t_east) % 28) < 20
        hd = np.where(t >= t_east, np.where(ph, 90, 270) + hbias + n6, hd)
    elif SCEN == "west_after_1804":
        hd = np.where(t >= t_east, 270 + hbias + n6, hd)
    elif SCEN == "zigzag":
        # "back and forth, but mostly east": alternate 20 min east / 8 min west
        ph = ((t - t_east) % 28) < 20
        hd = np.where(t >= t_east, np.where(ph, 90, 270) + hbias + n6, hd)
    else:
        hd = np.where(t >= t_east, 90 + hbias + n6, hd)
    hd = np.nan_to_num(hd, nan=0.0)
    alive = t < t_fuel
    spd = np.where(bomb, 0.0, tas)
    # front: after 18:00 wind backs toward NW and strengthens a little in 'front' scenario
    wd, ws = wdir, wspd
    if SCEN == "ships":
        sd, ss = ship_wind(lat, lon, t)
        late = t >= t_wx
        wd = np.where(late, sd, wdir); ws = np.where(late, ss, wspd)
    if SCEN == "era5":
        wd, ws = era5_wind(lat, lon, t)
    if SCEN in ("wx1945", "exhibit_decay"):
        late = t >= t_wx
        wd = np.where(late, wdir2, wdir); ws = np.where(late, wspd2, wspd)
    if SCEN == "front" and t >= T(18, 0):
        wd, ws = U(300, 340), wspd + 10
    nlat, nlon = step(lat, lon, hd, spd, wd, ws, dt)
    lat = np.where(alive, nlat, lat); lon = np.where(alive, nlon, lon)
    if T(15, 35) <= t <= T(16, 15):
        dmin_island = np.minimum(dmin_island, dist_bahamas(lat, lon))
    if T(15, 35) <= t <= T(16, 5):
        dmin_walker = np.minimum(dmin_walker, hav(lat, lon, 27.27, -78.40))
    if t == T(17, 33):
        d_pe_1733 = hav(lat, lon, 26.09, -80.12)
    if t == T(18, 50):
        d_pe_1850 = hav(lat, lon, 26.09, -80.12)
    if t == T(18, 0):
        d_canav_1800 = hav(lat, lon, 28.46, -80.54)
    if t == T(17, 50):
        lat1750, lon1750 = lat.copy(), lon.copy()
    if t % 5 == 0:
        for k, i in enumerate(track_idx):
            traj[k].append((float(lat[i]), float(lon[i])))
    t += dt
    if t > T(20, 20): break

d_fix = hav(lat1750, lon1750, fix_lat, fix_lon)
w = np.exp(-0.5 * (d_fix / 50.0) ** 2)            # fix radius 100 ~ 2 sigma
if "--radio" in sys.argv:
    # Port Everglades (ASRTU-4, Exhibit 8): strength 3 at 1733; only garbled at 1844 and nothing after.
    # Assumption-light constraint: farther from Port Everglades at 18:50 than at 17:33.
    _f = 1 / (1 + np.exp(-(d_pe_1850 - d_pe_1733) / 20.0))
    RADIO_EVID = float((w * _f).sum() / w.sum())
    _g = 1 / (1 + np.exp(-(d_pe_1850 - d_pe_1733 - 40) / 20.0))   # stricter: >=40 nm farther
    RADIO_EVID_STRICT = float((w * _g).sum() / w.sum())
    w = w * _f
if "--baxter" in sys.argv:
    # Baxter testimony (record p.119-120): near Walker Cay when Taylor thought he was over the Keys (~1545R);
    # ~100-150 miles off Cape Canaveral at ~1800R. Treated as soft (they are his DR estimates).
    w = w * np.exp(-0.5 * (dmin_walker / 25.0) ** 2)
    mid = 125 / 1.15  # statute -> nautical
    w = w * np.exp(-0.5 * ((d_canav_1800 - mid) / 35.0) ** 2)
if USE_ISLAND:
    w = w * np.exp(-0.5 * (dmin_island / 20.0) ** 2)  # over/near Bahamas islands ~15:35-16:15
w_in = (d_fix <= 100)
W = w / w.sum()

def wq(x, q):
    o = np.argsort(x); c = np.cumsum(W[o]); return float(x[o][np.searchsorted(c, q)])

mlat, mlon = float((W * lat).sum()), float((W * lon).sum())
dist_mean = hav(lat, lon, mlat, mlon)
r50, r90 = wq(dist_mean, 0.5), wq(dist_mean, 0.9)
ess = 1 / (W ** 2).sum()
out = dict(scenario=TAG, evidence=float(w.mean()), island_med_nm=round(float(np.median(dmin_island)),1), n=N, ess=round(float(ess)), share_within_fix=round(float(w_in.mean()), 3),
           mean_lat=round(mlat, 3), mean_lon=round(mlon, 3),
           lat_q=[round(wq(lat, q), 2) for q in (0.05, 0.5, 0.95)],
           lon_q=[round(wq(lon, q), 2) for q in (0.05, 0.5, 0.95)],
           r50_nm=round(r50), r90_nm=round(r90),
           dist_from_fix_center_at_1750_med=round(wq(d_fix, 0.5)))
if "--radio" in sys.argv:
    out["radio_evidence"] = round(RADIO_EVID, 3); out["radio_evidence_strict"] = round(RADIO_EVID_STRICT, 3)
print(json.dumps(out))
np.savez_compressed(f"out_{TAG}.npz", lat=lat, lon=lon, w=W, lat1750=lat1750, lon1750=lon1750,
                    tas=tas, wdir=wdir, wspd=wspd, t_270=t_270, t_fuel=t_fuel)
json.dump([traj[k] for k in range(len(track_idx))], open(f"tracks_{TAG}.json", "w"))
json.dump([float(W[i]) for i in track_idx], open(f"trackw_{TAG}.json", "w"))
