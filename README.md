# Flight 19 search package (14 November 2026)

Probability map for where the five TBM Avengers of Flight 19 (5 December 1945) most likely ditched, with the inputs, code and a draft search order. Prepared by Tycho, Vestigia Labs. The full written case is in the case file (Claude Docs): https://claude.ai/artifact/FcxLtnS4DTtdjNepUXTHFK

## Headline
- Probability-weighted centre: **30.13°N, 76.53°W** (about 140 nm east-northeast of the 1945 radio fix at 29°15'N 79°00'W), straddling the Blake Escarpment.
- 50% area: about 6,400 sq nm (6,000–6,800 depending on grid smoothing). 90% area: about 20,500–21,200 sq nm.
- Depth: 26% on the Blake Plateau (600–1,000 m), 21% on the escarpment (1,000–3,000 m), 52% in the Blake-Bahama Basin (3,000–5,000 m). Essentially none shallower than 600 m.
- Draft search order (assumed fleet sweep 330 sq nm/day; detection 90% flat ground, 60% escarpment): ~18% chance of finding the wreck in 7 days, ~30% in 14, ~50% in 30, ~72% in 60.

## Contents
| Path | What it is |
| --- | --- |
| `maps/flight19_map_simple.png` | Simple map for general readers: take-off, 1945 radio fix, most likely area |
| `maps/flight19_ditch_map.png` | Detailed map: 50%/90% areas, simulated tracks, planned route, 17:50 fix, weather, depth contours |
| `data/flight19_probability_grid.csv` | 2-arcminute cells: probability, 50%/90% flags, depth band, search rank, cumulative days and chance of success |
| `data/flight19_zones.geojson`, `.kml` | 50%/90% polygons, week-1 and weeks-1–2 search boxes, centre point, 17:50 fix (open the KML in Google Earth) |
| `data/flight19_summary.json` | Headline numbers and assumptions |
| `code/model.py` | Monte Carlo flight reconstruction (200,000 flights per run) |
| `code/fuel5.py` | Five-plane minimum-order fuel-out model |
| `code/search1945.py` | Negative information from the 6–10 Dec 1945 search boxes (debris drift + detection odds) |
| `code/plot.py`, `code/plot_simple.py`, `code/export_pkg.py` | Maps and export |
| `code/icoads_parse.py` | ICOADS IMMA1 parser |
| `inputs/era5/` | ERA5 900 hPa winds (knots), 1945-12-05 19Z to 12-06 01Z, 25–33°N, 82–73°W |
| `inputs/igra/` | Miami and Jacksonville soundings, 4–6 Dec 1945 |
| `inputs/icoads/` | ICOADS ship observations near the flight area, 5–6 Dec 1945 |
| `inputs/bathy/` | ETOPO1 depth bands, 2-arcminute, 27–32°N, 81–73°W |

The Navy Board record itself (NARA microfilm M1657, 549 pages) is not included: https://catalog.archives.gov/id/146850136

## Reproducing the final result
Python 3 with numpy, scipy, pandas, matplotlib, shapely, ppigrf. Model runs need `data/ne_50m_land.geojson` (Natural Earth) in the working directory (not included) and `inputs/bahamas_dist.npz`.
```
MAGVAR=1 TAGX=_m python3 model.py era5 1 --island --record --baxter --radio --fuelearly
MAGVAR=1 TAGX=_m python3 model.py era5 1 --island --record --baxter --radio --fuelearly --zigzag
# combine 50/50 with fuel-out <= 19:20, apply 17:50 ellipse -> final_v3_prior.npz (see case file, Method)
MODE=v3 python3 search1945.py mid
python3 export_pkg.py && python3 plot.py
```

## Evidence tiers
1. Primary records: NARA M1657 radio logs, findings, Coast Guard search plans; ICOADS 1945 ship reports; NAS Fort Lauderdale winds aloft; IGRA soundings; Weather Bureau maps.
2. Physics reconstructions: ERA5 reanalysis (ECMWF, via Google's public ARCO-ERA5 store); IGRF 1945 magnetic variation.
3. Bounded assumptions: ERA5 error ±15% speed / ±12°; fuel 320–340 gal per aircraft; burn ranges; radio-fix ellipse 35 × 65 nm; debris detection odds; survey sweep rate.
4. Neutral: 50/50 held east vs zigzag after 18:06.
Excluded as unverified: "tanks checked full"; Cuban interference quoted from the Port Everglades log; 2016 blog report on the Solomons deck log.

## Limits
A probability map, not a location. Not peer-reviewed. The spread is set by unrecorded crew behaviour after 18:07 and five hours of dead reckoning; only a new 1945 measurement (a radar plot or raw direction-finding bearings) could shrink it substantially. Any find must be reported to the Naval History and Heritage Command (Sunken Military Craft Act).
