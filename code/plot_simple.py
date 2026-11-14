"""Plain-language map: departure, 1945 radio fix, most likely area."""
import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from shapely.geometry import shape, box
from scipy.ndimage import gaussian_filter
land=json.load(open("data/ne_50m_land.geojson")); view=box(-83,24.3,-72.5,32.8)
fig,ax=plt.subplots(figsize=(9,7.4),dpi=170); ax.set_facecolor("#e8f0f6")
for f in land["features"]:
    g=shape(f["geometry"]).intersection(view)
    for p in getattr(g,"geoms",[g]):
        if p.geom_type=="Polygon": x,y=p.exterior.xy; ax.fill(x,y,color="#e4dfd2",ec="#a69f8f",lw=0.6,zorder=1)
d=np.load("search_mid_v3.npz"); lat,lon,w=d["lat"],d["lon"],d["w"]
xe=np.linspace(-83,-72.5,211); ye=np.linspace(24.3,32.8,171)
H,_,_=np.histogram2d(lon,lat,bins=[xe,ye],weights=w); H=gaussian_filter(H.T,1.2)
fl=np.sort(H.ravel())[::-1]; c=np.cumsum(fl)/fl.sum(); t50=fl[np.searchsorted(c,.5)]; t90=fl[np.searchsorted(c,.9)]
xc=(xe[:-1]+xe[1:])/2; yc=(ye[:-1]+ye[1:])/2
ax.contourf(xc,yc,H,levels=[t90,t50,H.max()+1],colors=["#f4c3ae","#d9603b"],zorder=2)
ax.contour(xc,yc,H,levels=[t90,t50],colors=["#c4613f","#9c3d1e"],linewidths=[0.8,1.2],zorder=3)
# departure
ax.plot(-80.12,26.07,"o",color="#222",ms=7,zorder=5)
ax.annotate("Fort Lauderdale\n(take-off, 2:10 pm)",(-80.12,26.07),xytext=(-82.8,26.6),fontsize=10,color="#222",bbox=dict(boxstyle="round,pad=0.3",fc="white",ec="#bbb",lw=0.6),arrowprops=dict(arrowstyle="-",color="#555",lw=0.8),zorder=6)
# 1945 fix
ax.plot(-79.0,29.25,"X",color="#2a6f3b",ms=10,mew=1.5,zorder=5)
ax.annotate("Where most accounts place them\n(1945 radio fix, 5:50 pm)\nabout 140 nautical miles from\nthe most likely area",(-79.0,29.25),xytext=(-79.75,27.35),fontsize=9.5,color="#2a6f3b",bbox=dict(boxstyle="round,pad=0.3",fc="white",ec="#9cc3a6",lw=0.6),arrowprops=dict(arrowstyle="-",color="#2a6f3b",lw=0.8),zorder=6)
# arrow fix -> centre
ax.annotate("",xy=(-76.75,30.05),xytext=(-78.8,29.35),arrowprops=dict(arrowstyle="-|>",color="#444",lw=1.4,mutation_scale=14),zorder=6)
# zone labels
ax.text(-76.3,31.35,"Most likely area\n(even chance the planes are here)",fontsize=10.5,fontweight="bold",color="#7a2c12",ha="center",zorder=6)
ax.annotate("Wider area\n(9-in-10 chance)",(-75.3,29.2),xytext=(-74.3,28.0),fontsize=9.5,color="#a2512f",ha="center",arrowprops=dict(arrowstyle="-",color="#a2512f",lw=0.8),zorder=6)
ax.text(-81.5,28.2,"FLORIDA",fontsize=11,color="#8a8372",fontweight="bold",rotation=-62,zorder=4)
ax.text(-76.6,25.6,"BAHAMAS",fontsize=10,color="#8a8372",fontweight="bold",zorder=4)
ax.text(-75.4,26.9,"ATLANTIC OCEAN",fontsize=11,color="#7d97ad",fontweight="bold",zorder=4)
# scale bar 100 nm at 25N
x0=-75.6; y0=24.75; dl=100/(60*np.cos(np.radians(25)))
ax.plot([x0,x0+dl],[y0,y0],color="#333",lw=2.5,zorder=6); ax.text(x0+dl/2,y0+0.12,"100 nautical miles",ha="center",fontsize=8.5,color="#333",zorder=6)
ax.set_xlim(-83,-72.5); ax.set_ylim(24.3,32.8); ax.set_aspect(1/np.cos(np.radians(28.5)))
ax.set_xticks([]); ax.set_yticks([])
for s in ax.spines.values(): s.set_color("#999")
ax.set_title("Where Flight 19 most likely came down, 5 December 1945",fontsize=13,pad=10)
fig.text(0.5,0.015,"Centre about 30.1°N, 76.5°W. Based on the 1945 Navy Board record, 1945 weather data and the 1945 search.",fontsize=8,color="#555",ha="center")
plt.tight_layout(rect=(0,0.03,1,1)); plt.savefig("flight19_map_simple.png"); print("ok")
