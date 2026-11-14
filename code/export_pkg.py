"""Export handover data: probability grid CSV, zones GeoJSON/KML, summary JSON (light 2 nm smoothing of the Monte Carlo grid)."""
import numpy as np, json, matplotlib
from scipy.ndimage import gaussian_filter
matplotlib.use("Agg"); import matplotlib.pyplot as plt
B=np.load("inputs/bathy/bands.npy"); d=np.load("work/search_mid_v3.npz"); lat,lon,w=d["lat"],d["lon"],d["w"]
i=np.round((lat-27)*30).astype(int); j=np.round((lon+81)*30).astype(int); ok=(i>=0)&(i<150)&(j>=0)&(j<241)
P=np.zeros(B.shape); np.add.at(P,(i[ok],j[ok]),w[ok]); P=gaussian_filter(P,1.0); P[B=="L"]=0; P/=P.sum()
yc=27+np.arange(150)/30; xc=-81+np.arange(241)/30
A=((2*np.cos(np.radians(yc)))[:,None]*2)*np.ones(B.shape)
rate={"a":1.0,"b":1.0,"c":0.6,"d":1.0,"L":1e-9}; pod={"a":0.9,"b":0.9,"c":0.6,"d":0.9,"L":0}
R=np.vectorize(rate.get)(B); D=np.vectorize(pod.get)(B); days=A/(330*R)
val=P*D/days; o=np.argsort(val.ravel())[::-1]
rank=np.empty(B.size,int); rank[o]=np.arange(1,B.size+1)
cumd=np.empty(B.size); cumd[o]=np.cumsum(days.ravel()[o]); cump=np.empty(B.size); cump[o]=np.cumsum((P*D).ravel()[o])
oh=np.argsort(P.ravel())[::-1]; ch=np.cumsum(P.ravel()[oh]); hpd=np.empty(B.size); hpd[oh]=ch
names={"a":"<600m","b":"plateau_600-1000m","c":"escarpment_1000-3000m","d":"basin_>3000m","L":"land"}
Y,X=np.meshgrid(yc,xc,indexing="ij")
with open("data/flight19_probability_grid.csv","w") as f:
    f.write("lat,lon,cell_area_sqnm,probability,cumulative_hpd,in_50pct_zone,in_90pct_zone,depth_band,search_rank,cum_search_days,cum_prob_find\n")
    for k in np.argsort(rank):
        if P.ravel()[k] < 1e-9: continue
        f.write(f"{Y.ravel()[k]:.4f},{X.ravel()[k]:.4f},{A.ravel()[k]:.2f},{P.ravel()[k]:.6e},{hpd[k]:.4f},{int(hpd[k]<=0.5)},{int(hpd[k]<=0.9)},{names[B.ravel()[k]]},{rank[k]},{cumd[k]:.2f},{cump[k]:.4f}\n")
HPD=hpd.reshape(B.shape); feats=[]
cs=plt.contour(xc,yc,HPD,levels=[0.5,0.9])
for lev,segs in zip([0.5,0.9],cs.allsegs):
    polys=[s for s in segs if len(s)>10]
    feats.append({"type":"Feature","properties":{"name":f"{int(lev*100)}% most likely ditch area","probability":lev},
      "geometry":{"type":"MultiPolygon","coordinates":[[[[round(x,4),round(y,4)] for x,y in np.vstack([s,s[:1]])]] for s in polys]}})
for T,name in [(7,"Week 1 search cells (bounding box)"),(14,"Weeks 1-2 search cells (bounding box)")]:
    sel=(cumd<=T).reshape(B.shape); ii,jj=np.where(sel)
    a0,a1,o0,o1=[round(float(v),3) for v in (yc[ii].min(),yc[ii].max(),xc[jj].min(),xc[jj].max())]
    feats.append({"type":"Feature","properties":{"name":name,"days":T,"prob_find":round(float(cump[(cumd<=T)].max()),3)},
      "geometry":{"type":"Polygon","coordinates":[[[o0,a0],[o1,a0],[o1,a1],[o0,a1],[o0,a0]]]}})
cen=(float(np.average(lat,weights=w)),float(np.average(lon,weights=w)))
feats.append({"type":"Feature","properties":{"name":"Probability-weighted centre"},"geometry":{"type":"Point","coordinates":[round(cen[1],3),round(cen[0],3)]}})
feats.append({"type":"Feature","properties":{"name":"17:50 HF/DF fix centre (Board finding 28), radius 100 statute miles"},"geometry":{"type":"Point","coordinates":[-79.0,29.25]}})
json.dump({"type":"FeatureCollection","features":feats},open("data/flight19_zones.geojson","w"))
def kml_poly(name,rings,color):
    pm=f'<Placemark><name>{name}</name><Style><LineStyle><color>{color}</color><width>2</width></LineStyle><PolyStyle><color>40{color[2:]}</color></PolyStyle></Style><MultiGeometry>'
    for r in rings: pm+='<Polygon><outerBoundaryIs><LinearRing><coordinates>'+" ".join(f"{x},{y},0" for x,y in r)+'</coordinates></LinearRing></outerBoundaryIs></Polygon>'
    return pm+'</MultiGeometry></Placemark>'
k='<?xml version="1.0" encoding="UTF-8"?><kml xmlns="http://www.opengis.net/kml/2.2"><Document><name>Flight 19 search zones</name>'
cols={"50% most likely ditch area":"ff3b60d9","90% most likely ditch area":"ffa2b8f2","Week 1 search cells (bounding box)":"ff00aa00","Weeks 1-2 search cells (bounding box)":"ff008800"}
for f_ in feats:
    g=f_["geometry"]; n=f_["properties"]["name"]
    if g["type"]=="MultiPolygon": k+=kml_poly(n,[p[0] for p in g["coordinates"]],cols[n])
    elif g["type"]=="Polygon": k+=kml_poly(n,g["coordinates"],cols[n])
    else: k+=f'<Placemark><name>{n}</name><Point><coordinates>{g["coordinates"][0]},{g["coordinates"][1]},0</coordinates></Point></Placemark>'
open("data/flight19_zones.kml","w").write(k+'</Document></kml>')
bands={names[b]:round(float(P[B==b].sum()),3) for b in "abcd"}
wk={f["properties"]["name"]:f["geometry"]["coordinates"][0][:3] for f in feats if f["geometry"]["type"]=="Polygon"}
summ={"centre_lat":round(cen[0],2),"centre_lon":round(cen[1],2),
 "area_50pct_sqnm":int(A.ravel()[hpd<=0.5].sum()),"area_90pct_sqnm":int(A.ravel()[hpd<=0.9].sum()),
 "area_note":"50% area is ~6,000-6,800 sq nm depending on grid resolution/smoothing of the Monte Carlo output; 6,400 quoted in the case file",
 "depth_band_probability":bands,"search_curve_prob_find_by_day":{str(T):round(float(cump[cumd<=T].max()),3) for T in [3,7,14,21,30,45,60]},
 "search_boxes":wk,
 "search_assumptions":{"fleet_sweep_sqnm_per_day":330,"detection_flat":0.9,"detection_escarpment":0.6,"escarpment_speed_factor":0.6},
 "model":"ERA5 900 hPa winds; magnetic headings (IGRF 1945); held east/zigzag 50/50; five-plane fuel-out 19:05-19:20; 17:50 fix as 35x65 nm ellipse; 1945 search negative information; ETOPO1 2-min depths"}
json.dump(summ,open("data/flight19_summary.json","w"),indent=2); print(json.dumps(summ,indent=1))
