import cartopy
import cartopy.crs as ccrs 
import cartopy.feature as cfeature 
from matplotlib.pyplot import *
import numpy as np 
import pandas as pd 
from netCDF4 import Dataset
from wrf import getvar, interplevel

fileDir = "/home/twsand/fskao/wrfOUT43v1/NC20210415/wrfout_d04_2021-04-15_12:00:00"
wrfdata = Dataset(fileDir)
lon = np.array(wrfdata["XLONG"])
lat = np.array(wrfdata["XLAT"])
t2 = np.array(wrfdata["T2"]) - 273

def findArgmin(data, value):
    IdxMinflat = np.argmin(abs(data - value))
    idxMin = np.unravel_index(IdxMinflat, data.shape)
    return idxMin

xitouLon, xitouLat = 120.7838246, 23.6759616
xitouGridargX = findArgmin(lon, xitouLon)
xitouGridargY = findArgmin(lat, xitouLat)

fig = figure(dpi=300)
ax = axes(projection=ccrs.PlateCarree())
# draw T2
t2 = ax.contourf(lon[0], lat[0], t2[0], vmin=0, vmax=30, cmap='coolwarm')
colorbar(t2, pad=0.1)
# draw xitou
scatter(x=xitouLon, y = xitouLat, s=5, c='red')
scatter(x=lon[xitouGridargX], y = lat[xitouGridargY], s=5, c='black')
ax.set_extent([119, 123, 21.5, 26], crs=ccrs.PlateCarree())
ax.gridlines(draw_labels=True, alpha=0.5, linestyle='--')
# draw geography
ax.coastlines()
ax.set_title("Taiwan [T2] (K)")
savefig("Taiwan.jpg")
