"""Combine the two ERA5 runs (held east / zigzag) into the final prior.
- keeps simulated flights whose fuel ran out by 19:20 EST (five-plane fuel model + last transmission at 19:04)
- weights held east / zigzag 50/50, scaled by each run's mean likelihood ("evidence") printed by model.py
- replaces the 17:50 radio-fix circle (sigma 50 nm) with an east-west ellipse (sigma 35 nm N-S, 65 nm E-W)
Writes work/final_v3_prior.npz, which search1945.py (MODE=v3) reads."""
import json, sys, numpy as np
runs = [("work/out_era5_isl_rec_bax_rad_fe_m.npz", "work/evidence_east.json"),
        ("work/out_era5_isl_rec_bax_rad_zz_fe_m.npz", "work/evidence_zigzag.json")]
lat, lon, w = [], [], []
for npz, evf in runs:
    d = np.load(npz); ev = json.load(open(evf))["evidence"]
    keep = d["t_fuel"] <= 19 * 60 + 20
    ww = d["w"] * keep; ww = ww / ww.sum() * 0.5 * ev
    dn = (d["lat1750"] - 29.25) * 60; de = (d["lon1750"] + 79) * 60 * np.cos(np.radians(29.25))
    ww *= np.exp(-0.5 * ((dn / 35) ** 2 + (de / 65) ** 2)) / np.exp(-0.5 * (dn ** 2 + de ** 2) / 50 ** 2)
    lat.append(d["lat"]); lon.append(d["lon"]); w.append(ww)
lat, lon, w = map(np.concatenate, (lat, lon, w)); w /= w.sum()
np.savez("work/final_v3_prior.npz", lat=lat, lon=lon, w=w)
print(f"prior centre {np.average(lat, weights=w):.2f}N {np.average(lon, weights=w):.2f}W")
