"""Map of reconstructed Flight 19 ditch area."""
import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from shapely.geometry import shape, box

land = json.load(open("inputs/ne_50m_land.geojson"))
view = box(-84, 24.5, -69, 33)
fig, ax = plt.subplots(figsize=(11, 7.2), dpi=150)
ax.set_facecolor("#eef3f7")
for f in land["features"]:
    g = shape(f["geometry"]).intersection(view)
    if g.is_empty: continue
    for p in getattr(g, "geoms", [g]):
        if p.geom_type == "Polygon":
            x, y = p.exterior.xy; ax.fill(x, y, color="#d9d4c7", ec="#9a9384", lw=0.5, zorder=1)

_S=np.load("work/search_mid_v3.npz"); W=_S["w"]; lat=_S["lat"]; lon=_S["lon"]  # posterior after 1945 search negative info
class _D(dict): pass
d = _D(lat=lat, lon=lon, w=W)
H, xe, ye = np.histogram2d(lon, lat, bins=[np.linspace(-84, -69, 151), np.linspace(24.5, 33, 86)], weights=W)
Hs = H.T
# HPD thresholds
flat = np.sort(Hs.ravel())[::-1]; c = np.cumsum(flat)
t50 = flat[np.searchsorted(c, 0.5)]; t90 = flat[np.searchsorted(c, 0.9)]
xc = (xe[:-1] + xe[1:]) / 2; yc = (ye[:-1] + ye[1:]) / 2
ax.contourf(xc, yc, Hs, levels=[t90, t50, Hs.max() + 1], colors=["#f2b8a2", "#d9603b"], alpha=0.85, zorder=2)
ax.contour(xc, yc, Hs, levels=[t90, t50], colors=["#b5452a"], linewidths=0.8, zorder=3)

_B=np.load("inputs/bathy/bands.npy"); _Z=np.select([_B=="L",_B=="a",_B=="b",_B=="c",_B=="d"],[-10,400,800,2000,4000])
_cs=ax.contour(-81+np.arange(241)/30, 27+np.arange(150)/30, _Z, levels=[1000, 3000], colors=["#5a7d9a"], linewidths=[0.7, 0.7], linestyles=["--", ":"], zorder=2.5)
ax.clabel(_cs, fmt={1000: "1,000 m", 3000: "3,000 m"}, fontsize=6.5)
# sample tracks (weighted top ones)
tr = json.load(open("work/tracks_era5_isl_rec_bax_rad_fe_m.json")); tw = np.array(json.load(open("work/trackw_era5_isl_rec_bax_rad_fe_m.json")))
for k in np.argsort(tw)[::-1][:25]:
    p = np.array(tr[k]); ax.plot(p[:, 1], p[:, 0], color="#3b5b7a", lw=0.6, alpha=0.45, zorder=4)

# planned route
R = 3440.065
def leg(lat, lon, hdg, nm):
    h = np.radians(hdg); return lat + np.degrees(nm * np.cos(h) / R), lon + np.degrees(nm * np.sin(h) / (R * np.cos(np.radians(lat))))
pts = [(26.05, -80.1167)]
for h, n in [(91, 56), (91, 67), (346, 73), (241, 120)]:
    pts.append(leg(*pts[-1], h, n))
pts = np.array(pts); ax.plot(pts[:, 1], pts[:, 0], color="k", lw=1.6, ls="--", zorder=5, label="Planned route (Nav. Problem No. 1)")

# 17:50 fix circle (100 nm radius)
th = np.linspace(0, 2 * np.pi, 200); r = 100
ax.plot(-79 + np.degrees(r * np.sin(th) / (R * np.cos(np.radians(29.25)))), 29.25 + np.degrees(r * np.cos(th) / R),
        color="#2a7a3b", lw=1.4, zorder=5, label="17:50 radio fix (100-mile radius)")
ax.plot(-79, 29.25, "+", color="#2a7a3b", ms=10, zorder=6)

ax.plot(-(80 + 25 / 60), 28 + 59 / 60, "*", color="#7a2a6b", ms=12, zorder=6, label="Gaines Mills flames, 19:50 (PBM Mariner)")
ax.plot(-78.5, 30.5, "o", color="#1f4e9e", ms=9, mfc="none", mew=2, zorder=6, label="999 mb storm low, 01:30 Dec 6 (Weather Bureau map)")
ax.annotate("cold front 01:30 Dec 5", xy=(-78.6, 27.5), xytext=(-78.4, 27.4), fontsize=7, color="#1f4e9e")
ax.plot([-78.4, -78.6, -78.9], [32.5, 27.5, 24.6], color="#1f4e9e", lw=1, ls=":", zorder=4)
for s, mk, lab in [("era5_isl_rec_bax_rad_fe_m", "D", "Centre: held east, early fuel-out"), ("era5_isl_rec_bax_rad_zz_fe_m", "s", "Centre: zig-zag, early fuel-out")]:
    z = np.load(f"work/out_{s}.npz"); m = (z["w"] * z["lon"]).sum(), (z["w"] * z["lat"]).sum()
    ax.plot(*m, mk, color="#222", ms=7, mfc="white", zorder=6, label=lab)
import pandas as _pd
_o = _pd.read_csv("inputs/icoads/icoads_1945-12-05_06_near.csv"); _o = _o[(_o.day == 5) & (_o.hr_utc >= 19)]
_u = _o[["lat", "lon"]].drop_duplicates()
ax.plot(_u.lon, _u.lat, "^", color="#2a6f97", ms=7, zorder=6, label="Ship weather reports, evening 5 Dec (ICOADS)")
ax.plot(-78.40, 27.27, "x", color="#333", ms=7, mew=2, zorder=6, label="Walker Cay (Baxter: position ~15:45)")
ax.plot([], [], color="#5a7d9a", lw=0.8, ls="--", label="1,000 m depth (Blake Plateau edge)")
ax.plot([], [], color="#5a7d9a", lw=0.8, ls=":", label="3,000 m depth (basin)")
ax.plot([], [], color="#d9603b", lw=8, label="Most likely ditch area (50%)")
ax.plot([], [], color="#f2b8a2", lw=8, label="Likely ditch area (90%)")
ax.set_xlim(-84, -69); ax.set_ylim(24.5, 33); ax.set_aspect(1 / np.cos(np.radians(29)))
ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")
ax.set_title("Flight 19: where the Avengers most likely ran out of fuel (Monte Carlo reconstruction)", fontsize=11)
ax.legend(loc="upper left", fontsize=7.5, framealpha=0.92)
ax.text(-83.8, 24.7, "NARA M1657 + ERA5 1945 winds at 900 hPa (checked vs Miami/Jacksonville soundings) + magnetic variation (IGRF 1945).\nHeld east / zig-zag 50/50; five-plane fuel-out 19:05-19:20; HF/DF ellipse; 6-10 Dec 1945 search applied. Depths: ETOPO1 2-min.", fontsize=6.8, color="#444")
ax.grid(color="white", lw=0.6)
plt.tight_layout(); plt.savefig("maps/flight19_ditch_map.png")
print("saved")
