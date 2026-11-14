"""Parse ICOADS IMMA1 core fields and extract ship obs near Flight 19, 5-6 Dec 1945."""
import sys, csv

def fld(s, a, b):  # 1-indexed inclusive columns
    v = s[a - 1:b].strip()
    return v

def num(v, scale=1.0):
    try: return int(v) * scale
    except ValueError: return None

rows = []
for line in open(sys.argv[1], errors="ignore"):
    yr, mo, dy, hr = num(fld(line, 1, 4)), num(fld(line, 5, 6)), num(fld(line, 7, 8)), num(fld(line, 9, 12), 0.01)
    lat, lon = num(fld(line, 13, 17), 0.01), num(fld(line, 18, 23), 0.01)
    if None in (yr, mo, dy, lat, lon): continue
    if lon > 180: lon -= 360
    if not (dy in (5, 6) and 24 <= lat <= 34 and -82 <= lon <= -68): continue
    rows.append(dict(
        day=dy, hr_utc=hr, lat=lat, lon=lon, id=fld(line, 35, 43), ii=fld(line, 33, 34),
        wdir_deg=num(fld(line, 47, 49)), wi=fld(line, 50, 50), w_ms=num(fld(line, 51, 53), 0.1),
        ww=fld(line, 57, 58), slp_hpa=num(fld(line, 60, 64), 0.1), at_c=num(fld(line, 70, 73), 0.1),
        deck=fld(line, 119, 121) if len(line) > 121 else ""))
w = csv.DictWriter(open("near.csv", "w"), fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print(len(rows), "obs in box 24-34N, 68-82W on 5-6 Dec")
