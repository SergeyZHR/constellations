#!/usr/bin/python3
# -*- coding: UTF-8 -*-
import numpy as np
#from astropy.io import fits
from astropy.table import Table
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time, TimeDelta
import pymesh
#from astropy.utils import iers
import warnings
from astropy.utils.exceptions import AstropyWarning
#iers.conf.auto_download = False
warnings.simplefilter('ignore', AstropyWarning)


def make(Cosntellation):
	Clines_f = map(lambda s: s[:-1].split(), Clines )				#read constellation lines information
	line_arr = pd.DataFrame(list(Clines_f))
	line_arr=np.array(line_arr[line_arr.loc[:,0] ==Cosntellation ].iloc[0].dropna()[2:]).astype(int)
	#print(line_arr)

	dat = Table.read('hip.fit', format='fits')		
	df = dat.to_pandas().set_index('HIP')				# i llove pandas :)
	df = df.loc[np.unique(line_arr)]
	#print(df)

	RA_center = ( np.max(df['_RA_icrs']) + np.min(df['_RA_icrs']) )/2.  			#constellation center
	if np.max(df['_RA_icrs']) - np.min(df['_RA_icrs'])>180:
		RA_center = ( np.max(df['_RA_icrs'])-360 + np.min(df['_RA_icrs']) )/2.  
		if RA_center<0:
			RA_center+=360

	DE_center = ( np.max(df['_DE_icrs']) + np.min(df['_DE_icrs']) )/2.			

	#print(RA_center,DE_center)

	loc = EarthLocation(lat=DE_center*u.deg, lon=0*u.deg)			#convert coordinates
	fake_time = Time('2019-09-21T00:01:40.2', format='isot', scale='utc',location = loc) + TimeDelta(RA_center/360*86164, format='sec')
	fakeFrame = AltAz(obstime=fake_time,location=loc)	
	cords = SkyCoord(list(df['_RA_icrs'])*u.deg, list(df['_DE_icrs'])*u.deg).transform_to(fakeFrame)
	df['phi'] = np.radians(cords.az.deg)+np.pi/2
	df['z'] = 90 - cords.alt.deg

	# df['phi'] = np.radians(df['_RA_icrs'])
	# df['z'] = 90 - df['_DE_icrs']

	df['rho'] = df['z']
	df['x'] = df['rho']*np.cos(df['phi'])
	df['y'] = df['rho']*np.sin(df['phi'])

	scale= np.max([df['x'],df['y']]) - np.min([df['x'],df['y']])
	scale = 45
	df['x']*=230/scale			#calculate scale
	df['y']*=230/scale
	# import matplotlib.pyplot as plt
	# plt.plot(df['x'],df['y'],'o')
	# plt.show()
	print('coverted'+Cosntellation)

	star_list = []
	for hip, star in df.iterrows():  #make stars
		r = 8*1.3**(-star['Vmag'])
		if r>15: r=15
		sphere = pymesh.generate_icosphere(r,(star['x'],star['y'],0),refinement_order=2)
		sphere = pymesh.form_mesh(sphere.vertices*np.array((1,1,0.75)), sphere.faces);
		star_list.append(sphere)

	mesh = pymesh.merge_meshes(star_list)

	print('stars'+Cosntellation)


	# print(df)
	cil_list=[]
	for i in range(0,len(line_arr),2):    #make lines
		cilinder = pymesh.generate_cylinder((df.loc[line_arr[i],'x'],df.loc[line_arr[i],'y'],0), (df.loc[line_arr[i+1],'x'],df.loc[line_arr[i+1],'y'],0), 1.4, 1.4 )
		cilinder = pymesh.form_mesh(cilinder.vertices*np.array((1,1,0.5)), cilinder.faces);
		cil_list.append(cilinder)

	cil_mesh = pymesh.merge_meshes(cil_list)

	print('lines'+Cosntellation)

	mesh = pymesh.boolean(mesh,cil_mesh,operation='union')

	DOWN = pymesh.meshutils.generate_box_mesh((-400,-400,-400),(400,400,0))
	mesh = pymesh.boolean(mesh,DOWN,'difference')

	pymesh.meshio.save_mesh('./'+Cosntellation+'.stl',mesh)
	
	print('completed'+Cosntellation)


const_lines_file = open("./constellationship.fab", "r")
Clines = const_lines_file.readlines()


#for line in Clines:
#	make(line.split()[0])

make('UMi')