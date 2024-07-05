# -*- coding: utf-8 -*-
"""
Created on Thu May  9 19:50:07 2024

@author: kevinhayley
"""
%matplotlib auto
import os
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt

import flopy 
from flopy.utils.triangle import Triangle
import shapely

def shapely_linestring2_coords_list(linestring,decimals=2):
    # decimals=number of decimals in output coordinates
    xy=linestring.xy
    out=np.zeros((len(xy[0]),2))
    out[:,0]=np.around(xy[0],decimals)
    out[:,1]=np.around(xy[1],decimals)
    return out

def ele_gdf(tri):
    polys=[]
    for e in tri.ele:
        verts=[]
        for i in range(1,4):
            verts.append((tri.node[e[i]]['x'],tri.node[e[i]]['y']))
        verts.append((tri.node[e[1]]['x'],tri.node[e[1]]['y']))
        polys.append(shapely.Polygon(verts))
    gdf=gpd.GeoDataFrame({'ele':tri.ele['icell'],'geometry':polys})
    return gdf


model_poly=gpd.read_file('Polygon.shp')
# shp=model_poly.loc[0,'geometry']
# new_coords=[]
# for c in shp.boundary.coords:
#     new_coords.append((round(c[0]),round(c[1])))
# model_poly.loc[0,'geometry']=shapely.geometry.Polygon(new_coords)
# model_poly.to_file('Polygon.shp')    

lines=gpd.read_file('lines1.shp')
# for i in lines.index:
#     shp=lines.loc[i,'geometry']
#     new_coords=[]
#     for c in shp.coords:
#         new_coords.append((round(c[0]),round(c[1])))
#     lines.loc[i,'geometry']=shapely.geometry.LineString(new_coords)    
# lines.to_file('Lines1.shp')    





points=gpd.read_file('points.shp')
refinement_zone=gpd.read_file('Polygon2.shp')

wellpts=np.zeros((len(points.index),2))
for i in points.index:
    wellpts[i,0]=points.loc[i,'geometry'].x
    wellpts[i,1]=points.loc[i,'geometry'].y
    

workspace = os.getcwd()

tri = Triangle(maximum_area=100000, angle=25, nodes=wellpts, model_ws=workspace)
tri.add_polygon(model_poly.iloc[0]['geometry'])


for line in lines['geometry']:
#    tri.add_LineString(line)
    tri.add_linestring(shapely_linestring2_coords_list(line,decimals=12))
#        lines.append(shapely_linestring2_coords_list(line,decimals=8))


tri.build(verbose=True)
fig = plt.figure(figsize=(10, 10))
ax = plt.subplot(1, 1, 1, aspect="equal")
pc = tri.plot(ax=ax)
plt.title('inital_mesh')




# first refinement iteration Refine at zone defined by polygon 
ele=[]
gdf=ele_gdf(tri)
for i in gdf.index:
    if (gdf.loc[i,'geometry'].touches(refinement_zone.loc[0,'geometry']) or
        gdf.loc[i,'geometry'].within(refinement_zone.loc[0,'geometry']) or
        gdf.loc[i,'geometry'].overlaps(refinement_zone.loc[0,'geometry'])):
        ele.append(gdf.loc[i,'ele'])

tri.refine_ele(ele,20000,iteration=1,verbose=True)

fig = plt.figure(figsize=(10, 10))
ax = plt.subplot(1, 1, 1, aspect="equal")
pc = tri.plot(ax=ax)
plt.title('first refinement')

#second refinement at pumping wells


# second refinement around pumping bores
nds_2_refine=np.array(wellpts)
nd_numbers=tri.node[np.logical_and(np.isin(tri.node['x'],nds_2_refine[:,0]),np.isin(tri.node['y'],nds_2_refine[:,1]))]['ivert']
tri.refine_nds(nd_numbers,2000,iteration=2,verbose=True)
fig = plt.figure(figsize=(10, 10))
ax = plt.subplot(1, 1, 1, aspect="equal")
pc = tri.plot(ax=ax)
plt.title('second refinement')

