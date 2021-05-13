#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import sys
from xml.etree import ElementTree as ET
from datetime import datetime,timedelta
import math
import subprocess
import shutil
import glob
import gdal
import os
from matplotlib import pyplot as plt


"""DEFINO ALGUNAS FUNCIONES"""

def rebin(a, shape):
	sh = (shape[0], a.shape[0]//shape[0], shape[1], a.shape[1]//shape[1])
	return a.reshape(sh).mean(-1).mean(1)

#Defino una s genérica: i es la coordenada a la que me refiero (0 para X, 1 para Y, 2 para Z) 

def orbitV(t, master, i):
	m = master[:,i]
	ss = m[0] + m[1]*t + m[2]*t**2 + m[3]*t**3 + m[4]*t**4 + m[5]*t**5 + m[6]*t**6 + m[7]*t**7 + m[8]*t**8 
	return ss

def timeSincro(t0, x, y, z, master, maxit=100, tolerance=1e-8):

	delta = tolerance
	
	tt = np.zeros(len(master)-2)
	
	while np.min(abs(delta)) >= tolerance and maxit > 0:
		t0 = np.min(t0)
		tt[0] = t0
		for i in range(1, len(master)-2):
			globals()['t%d'%i] = eval('t%d'%(i-1)) * t0
			tt[i] = eval('t%d'%i)
		
		sxm = master[0,0] + tt @ master[1:-1,0]
		sym = master[0,1] + tt @ master[1:-1,1]
		szm = master[0,2] + tt @ master[1:-1,2]
		
		vxm = master[0,3] + tt @ master[1:-1,3]
		vym = master[0,4] + tt @ master[1:-1,4]
		vzm = master[0,5] + tt @ master[1:-1,5]
 
		axm = master[0,6] + tt @ master[1:-1,6]
		aym = master[0,7] + tt @ master[1:-1,7]
		azm = master[0,8] + tt @ master[1:-1,8]

		delta = -(vxm*(x - sxm) + vym*(y - sym) + vzm*(z - szm)) / (axm*(x - sxm) + aym*(y - sym) + azm*(z - szm) - vxm**2 - vym**2 - vzm**2)
		
		#print(f'delta min {np.min(delta)}')
		#print(f'delta max {np.max(delta)}')
				
		t0 += delta
		maxit -= 1
	
	return t0


def rangeMat0(x, y, z, t0, master):
	t = timeSincro(t0, x, y, z, master)
	range_0 = np.sqrt((x - orbitV(t, master, 0))**2 + (y - orbitV(t, master, 1))**2 + (z - orbitV(t, master, 2))**2)
	return range_0 


def rangeMat(x, y, z, t, master):
	#Leo el shape que me da fila y columna
	ranget = np.zeros((x.shape), float)
	t0 = timeSincro(t, x[0,:], y[0,:], z[0,:], master)
	msx, msy, msz = master[:,0], master[:,1], master[:,2]
	"""
	El IDL hace el for para cada FILA (segunda dimensión). 
	En python definimos a X, Y, Z como arreglos de Fila por Columna 
	(dimensiones 2, 1 de IDL). Si en IDL se itera por siz(2), 
	tercer elemento de la variable SIZ, aca tenemos que iterar por siz(0)
	"""
	# ~ for i in range(3519):
	for i in range(x.shape[0]):	
		
		print(f'ciclo {i}')
		
		xt, yt, zt = x[i,:], y[i,:], z[i,:]
		t0 = timeSincro(t0, xt, yt, zt, master)
		#print(f't0 min {np.min(t0)}')
		#print(f't0 max {np.max(t0)}')
		
		for j in range(1, master.shape[0]-1):
			globals()['t'+str(j)] = eval('t' + str(j - 1)) * t0

		sxm = msx[0] + msx[1]*t0 + msx[2]*t1 + msx[3]*t2 + msx[4]*t3 + msx[5]*t4 + msx[6]*t5 + msx[7]*t6 + msx[8]*t7 #+ msx[9]*t8 + msx[10]*t9 + msx[11]*t10 + msx[12]*t11 + msx[13]*t12 + msx[14]*t13 #+ msx[15]*t14 + msx[16]*t15  + msx[17]*t16 + msx[18]*t17
		sym = msy[0] + msy[1]*t0 + msy[2]*t1 + msy[3]*t2 + msy[4]*t3 + msy[5]*t4 + msy[6]*t5 + msy[7]*t6 + msy[8]*t7 #+ msy[9]*t8 + msy[10]*t9 + msy[11]*t10 + msy[12]*t11 + msy[13]*t12 + msy[14]*t13 #+ msy[15]*t14 + msy[16]*t15  + msy[17]*t16 + msy[18]*t17
		szm = msz[0] + msz[1]*t0 + msz[2]*t1 + msz[3]*t2 + msz[4]*t3 + msz[5]*t4 + msz[6]*t5 + msz[7]*t6 + msz[8]*t7 #+ msz[9]*t8 + msz[10]*t9 + msz[11]*t10 + msz[12]*t11 + msz[13]*t12 + msz[14]*t13 #+ msz[15]*t14 + msz[16]*t15  + msz[17]*t16 + msz[18]*t17

		ranget[i,:] = np.sqrt((xt - sxm)**2 + (yt - sym)**2 + (zt - szm)**2)
		#print(ranget[i,:])
	return ranget


def coordEllipsoid(t, ro1, master, x0, y0, z0, maxit=100, tolerance=1e-8):

	a = 6378137
	b = 6356752.314

	x = x0.copy()
	y = y0.copy()
	z = z0.copy()

	aa = np.ones((3,3))
	f = np.ones((3,1))

	sxm = orbitV(t, master, 0)
	sym = orbitV(t, master, 1)
	szm = orbitV(t, master, 2)

	vxm = orbitV(t, master, 3)
	vym = orbitV(t, master, 4)
	vzm = orbitV(t, master, 5)

	vmod = np.sqrt(vxm**2 + vym**2 + vzm**2)
	#print(vmod.shape)

	vxm = vxm/vmod
	vym = vym/vmod
	vzm = vzm/vmod

	delta = tolerance
	while np.min(abs(delta)) >= tolerance and maxit > 0:

		ro1_temp = np.sqrt((x - sxm)**2 + (y - sym)**2 + (z - szm)**2)
		rr_app = np.sqrt(x**2 + y**2 + (a * z/b)**2)

		aa[0,0] = (x - sxm)/ro1_temp
		aa[0,1] = (y - sym)/ro1_temp
		aa[0,2] = (z - szm)/ro1_temp

		aa[1,0] = vxm
		aa[1,1] = vym
		aa[1,2] = vzm

		aa[2,0] = x/rr_app
		aa[2,1] = y/rr_app
		aa[2,2] = (a/b)**2 * z/rr_app
		#print(np.linalg.matrix_rank(aa))

		f[0] = ro1_temp - ro1
		f[1] = (x - sxm)*vxm + (y - sym)*vym + (z - szm)*vzm
		f[2] = rr_app - a

		aai = np.linalg.inv(aa)
		delta = -aai @ f
		x = x + delta[0]
		y = y + delta[1]
		z = z + delta[2]
		  
		maxit -= 1

	return np.hstack((x, y, z))

def crossProduct(a, b):
	
	c = np.zeros(a.shape)
	c[0] = a[1]*b[2] - a[2]*b[1]
	c[1] = -(a[0]*b[2] - a[2]*b[0])
	c[2] = a[0]*b[1] - a[1]*b[0]

	return c

"""orbitP"""

def orbitP(t, x, y, z):
	#Declaración de las variables
	N_stat = len(t)
	
	rvx1, rvy1, rvz1 = np.zeros(N_stat), np.zeros(N_stat), np.zeros(N_stat)
	rax1, ray1, raz1 = np.zeros(N_stat), np.zeros(N_stat), np.zeros(N_stat)
	
	master = np.zeros((N_stat, 9))
	
	par = [0, 0, 0, 0, 0, 0, N_stat, t.copy(), x, y, z]
	
	#Parámetros de la escena procesada:
	
	par[0] = t[0]
	par[1] = t[0]
	par[2] = t[1] - t[0]
	par[5] = np.sum(t)/len(t) #t_mid_m
	par[3] = par[0] - par[5] #tstart
	
	t -= par[5]
	
	a = np.array([t[i]**j for i in range(N_stat) for j in range(N_stat)])
	a = a.reshape(N_stat,N_stat)
	
	av = np.linalg.inv(a)
	
	rx1= (x.T @ av.T).T
	ry1= (y.T @ av.T).T
	rz1= (z.T @ av.T).T
	
	
	for i in range(1, N_stat):
		rvx1[i-1] = i*rx1[i]
		rvy1[i-1] = i*ry1[i]
		rvz1[i-1] = i*rz1[i]
		
	for i in range(1, N_stat):
		rax1[i-1] = i*rvx1[i]
		ray1[i-1] = i*rvy1[i]
		raz1[i-1] = i*rvz1[i]
		
	master[:,0] = rx1
	master[:,1] = ry1
	master[:,2] = rz1
	master[:,3] = rvx1
	master[:,4] = rvy1
	master[:,5] = rvz1
	master[:,6] = rax1
	master[:,7] = ray1
	master[:,8] = raz1
	
	#print(rx1)
	
	return master, par

def deltaOrb(X, Y, Z, Omaster, Oslave, DtB0, DtAngPar):

	DimXSar, DimYSar = X.shape[1], X.shape[0]
	ShAz = 10
	#Imprimir la cantidad de elementos o vectores orbitales que tiene la slave
	pOrbit = len(Oslave[:,0])
	step = int((DimYSar - ShAz*2) / (pOrbit - 1))
	indi = np.arange(pOrbit)*step
	
	BperpV = np.zeros((pOrbit, 3))
	VrcenV = np.zeros((pOrbit, 3))
	MidRange = int(DimXSar/2)
	xm = X[:, MidRange]
	ym = Y[:, MidRange]
	zm = Z[:, MidRange]
	RgMid = np.mean(rangeMat0(xm, ym, zm, 0, Omaster))
		
	################ Calcular los versores de los desplazamientos ##################
	# En este for lo que hace es ir moviendose en Azimut  y obtener el valor XYZ para el rango medio
	for i in range(pOrbit):
		Xmr = X[indi[i] + ShAz, MidRange]
		Ymr = Y[indi[i] + ShAz, MidRange]
		Zmr = Z[indi[i] + ShAz, MidRange]
		tcenref = timeSincro(0, Xmr, Ymr, Zmr, Omaster)
		ElliPoi = coordEllipsoid(tcenref, RgMid, Omaster, Xmr, Ymr, Zmr) #  Calcola corrispondenti punti su ellissoide
		Versref = np.array([orbitV(tcenref, Omaster, 3), orbitV(tcenref, Omaster, 4), orbitV(tcenref, Omaster, 5)])
		Versref = Versref/np.sqrt(np.sum(Versref**2))
		satcen = np.array([orbitV(tcenref, Omaster, 0), orbitV(tcenref, Omaster, 1), orbitV(tcenref, Omaster, 2)])
		rcen = ElliPoi - satcen
		VrcenV[i,:] = rcen/np.sqrt(np.sum(rcen**2))
		BperpV[i,:] = crossProduct(Versref, VrcenV[i,:])
	
	Tsn0 = timeSincro(0, X[indi+ShAz,MidRange], Y[indi+ShAz,MidRange], Z[indi+ShAz,MidRange], Oslave)
	Tsn0m = Tsn0 - np.mean(Tsn0)
	
	Sxs0 = orbitV(Tsn0, Oslave, 0)
	Sys0 = orbitV(Tsn0, Oslave, 1)
	Szs0 = orbitV(Tsn0, Oslave, 2)
	
	# Mueve los vectores de órbita de la Slave según las perturbaciones en Bperp y Bpar
	Sxsn0 = Sxs0 + DtB0*BperpV[:,0] + Tsn0m*DtAngPar*VrcenV[:,0]
	Sysn0 = Sys0 + DtB0*BperpV[:,1] + Tsn0m*DtAngPar*VrcenV[:,1]
	Szsn0 = Szs0 + DtB0*BperpV[:,2] + Tsn0m*DtAngPar*VrcenV[:,2]
	
	# ~ print(f'Tsn0:\n{Tsn0}')
	# ~ print(f'Sxsn0:\n{Sxsn0}')
	# ~ print(f'Sysn0:\n{Sysn0}')
	# ~ print(f'Szsn0:\n{Szsn0}')
	
	OslaveN, parsn = orbitP(Tsn0, Sxsn0, Sysn0, Szsn0)
	
	# ~ print(f'OslaveN:\n{OslaveN}')
	# ~ print(f'parsn:\n{parsn}')
	
	# ~ sys.exit(0)

	return OslaveN,parsn


def ll2ECEF(lat,lon,alt):
    # WGS84 ellipsoid constants:
    a = 6378137
    e = 8.1819190842622e-02

    # intermediate calculation
    # (prime vertical radius of curvature)
    N = a / np.sqrt(1 - e**2 * np.sin(lat)**2)
    X = (N+alt) * np.cos(lat) * np.cos(lon)
    Y = (N+alt) * np.cos(lat) * np.sin(lon)
    Z = ((1 - e**2) * N + alt) * np.sin(lat)

    return X,Y,Z
    
def ecef2lla(x, y, z):
	# x, y and z are scalars or vectors in meters

	a=6378137
	a_sq=a**2
	e = 8.181919084261345e-2
	e_sq = 6.69437999014e-3

	f = 1/298.257223563
	b = a*(1-f)

	# calculations:
	r = np.sqrt(x**2 + y**2)
	ep_sq  = (a**2-b**2)/b**2
	ee = (a**2-b**2)
	f = (54*b**2)*(z**2)
	g = r**2 + (1 - e_sq)*(z**2) - e_sq*ee*2
	c = (e_sq**2)*f*r**2/(g**3)
	s = (1 + c + np.sqrt(c**2 + 2*c))**(1/3.)
	p = f/(3.*(g**2)*(s + (1./s) + 1)**2)
	q = np.sqrt(1 + 2*p*e_sq**2)
	r_0 = -(p*e_sq*r)/(1+q) + np.sqrt(0.5*(a**2)*(1+(1./q)) - p*(z**2)*(1-e_sq)/(q*(1+q)) - 0.5*p*(r**2))
	u = np.sqrt((r - e_sq*r_0)**2 + z**2)
	v = np.sqrt((r - e_sq*r_0)**2 + (1 - e_sq)*z**2)
	z_0 = (b**2)*z/(a*v)
	h = u*(1 - b**2/(a*v))
	phi = np.arctan((z + ep_sq*z_0)/r)
	lambd = np.arctan2(y, x)


	return phi*180/np.pi, lambd*180/np.pi, h
   
  
	
def OrbitalCorrector(x,y,z,coh,fas,ref,sec,rlook=1,alook=1,gamma=0.25):
	
	rgref = rangeMat(x, y, z, 0, ref)
	rgsec = rangeMat(x, y, z, 0, sec)
	
	sint1 = (rgref - rgsec)*4*np.pi/0.23513133960784313
	fas2 = np.arctan2(np.sin(-fas + sint1), np.cos(-fas + sint1))

	"""Estimación del delta Bperp"""
	print("Estimando delta Línea de base perpendicular")

	gradv = np.zeros(2)
	DtAngPar = 0
	SemIntDtB = 10
	spb = [-SemIntDtB, SemIntDtB]

	## Calcula el gradiente a lo largo del rango moviendo la orbita
	## Lo hace para el intervalo B_|_ +/- SemIntDtB

	for i in range(2):
		
		DtB0 = spb[i]
		secn0, p = deltaOrb(x, y, z, ref, sec, DtB0, DtAngPar)
		
		# Calcula el nuevo diferencial con la nueva orbita esclava
		
		rgsecn = rangeMat(x, y, z, 0, secn0)
		sint1n0 = ((rgref - rgsecn)*4*np.pi/0.23513133960784313)
		Diffp = np.arctan2(np.sin(fas2 - sint1n0), np.cos(fas2 - sint1n0))
		# ~ Diffp = np.arctan2(rebin(np.sin(fas2 - sint1n0), (300, 500)), rebin(np.cos(fas2 - sint1n0), (300, 500)))
		# ~ Diffp = np.arctan2(rebin(np.sin(Diffp), (300, 500)), rebin(np.cos(Diffp), (300, 500)))
		
		#Calculo del gradiente en dirección del rango
		
		Diffp_shift = np.roll(Diffp, 1, axis=1) # Lo muevo por columnas.
		gradd = np.arctan2(np.sin((Diffp - Diffp_shift)), np.cos((Diffp - Diffp_shift)))
		# ~ coh > gamma = np.where(~np.isnan(gradd))
		grad = gradd[(coh > gamma) & (~np.isnan(gradd))]
		
		#Usa el histograma para el cálculo del gradiente medio:
		binsiz = 0.2
		bins_n = int(np.ceil((np.max(grad) - np.min(grad))/binsiz))
		
		y_hist, x_hist, _ = plt.hist(grad.ravel(), bins=bins_n)
		# ~ plt.show()
		#max_i = np.where(y_hist == np.max(y_hist))[0][0]
		res = np.poly1d(np.polyfit(np.arange(len(y_hist)-2), y_hist[1:-1], 4))
		tim = np.arange(1000) / 1000 * (len(y_hist)-2)
		histif1 = res(tim)
		
		
		gradv[i]= np.mean((np.where(histif1 == np.max(histif1))[0][0]*(len(y_hist)-2) / 1000. + 1.)*binsiz + np.min(grad))
	
	a0 = np.mean(gradv)
	a1 = (gradv[1]-gradv[0])/(spb[1]-spb[0])
	DtBv1=-a0/a1
	
	gradv1=np.zeros(6)
	DtBvv=np.zeros(6)
	
	## Calcola il gradiente di fase lungo il range variando la baseline | 
	## di +/- 5m attorno più la prima stima DtBv1
	# window,1,retain=2,title='Correzione Baseline Perpendicolare',xs=DimYsar/scl,ys=DimXsar/scl

	for i in range(6):
		DtB0 = DtBv1 + i*2 - 5
		DtBvv[i] = DtB0
		secn0, p = deltaOrb(x, y, z, ref, sec, DtB0, DtAngPar)
		
		# Calcola nuovo differenziale con la nuova orbita slave
		rgsecn = rangeMat(x,y,z,0,secn0)
		sint1n0 = (rgref-rgsecn)*4*np.pi/0.23513133960784313
		Diffp = np.arctan2(np.sin(fas2-sint1n0), np.cos(fas2-sint1n0))
		# ~ Diffp = np.arctan2(rebin(np.sin(fas2-sint1n0),(300,500)), rebin(np.cos(fas2-sint1n0),(300,500)))
		
		# Grdiente della fase lungo il range
		Diffp_shift = np.roll(Diffp, 1, axis=1) # Lo muevo por columnas.
		gradd = np.arctan2(np.sin((Diffp-Diffp_shift)),np.cos((Diffp-Diffp_shift)))
		grad = gradd[(coh > gamma) & (~np.isnan(gradd))]
		gradv1[i] = np.mean(grad)


	#Calcola la retta di regressione dei gradienti gradv1 al variare
	#della correzione di B_|_. Lo 0 ci dà la correzione DtB cercata

	# Recordar que np.polyfit da los coeficientes en orden decreciente
	resec = np.polyfit(DtBvv, gradv1, deg=1)
	DtB = -resec[1]/resec[0]
	# ~ resec = np.polynomial.polynomial.Polynomial.fit(DtBvv,gradv1,deg=1)
	# ~ DtB = -resec.coef[1]/resec.coef[0]
	
	print(f'DtBvv: {DtBvv}, gradv1: {gradv1}')
	print(f'resec: {resec}')
	print(f'DtB: {DtB}')
	# ~ DtB = 0
	# Acá podemos aplicar DtB al de baja y al de alta por separado

	secn0, p = deltaOrb(x, y, z, ref, sec, DtB, DtAngPar)
	rgsecn = rangeMat(x,y,z,0,secn0)
	sint1n0 = (rgref-rgsecn)*4*np.pi/0.23513133960784313
	Diffp1 = np.arctan2(np.sin(fas2-sint1n0), np.cos(fas2-sint1n0))


	print("Estimando delta Línea de base paralela")

	############################################# PARALELA #################################################################

	################ Ricerca su parallela +/- .0005*65 metri/sec #############

	gradv= np.zeros(2)
	SemIntDtAngPar=0.0005*65
	spb=[-SemIntDtAngPar,SemIntDtAngPar]

	## Calcola il gradiente di fase lungo l'Azimuth variando di +/- .0005 m/s 
	## la velocità di variazione della B|| lungo la direzione di volo

	for i in range(2):
		
		# Mover la orbita slave original con el DtB ya estimado y el DtAngPar
		DtAngPar0=spb[i]
		secn0,p = deltaOrb(x,y,z,ref,sec,DtB,DtAngPar0)
		
		# Calcola nuovo differenziale con la nuova orbita slave
		rgsecn=rangeMat(x,y,z,0,secn0)
		sint1n0 = (rgref-rgsecn)*4*np.pi/0.23513133960784313
		Diffp =np.arctan2(np.sin(fas2-sint1n0), np.cos(fas2-sint1n0))
		
		#Calculo del gradiente en direcciión de azimut:
		Diffp_shift = np.roll(Diffp,1, axis = 0) # Lo muevo por filas
		gradd=np.arctan2(np.sin((Diffp-Diffp_shift)),np.cos((Diffp-Diffp_shift)))
		grad = gradd[(coh > gamma) & (~np.isnan(gradd))]
		gradv[i] = np.mean(grad)


	# ;| Calcola la retta di regressione dei gradienti gradv1 al variare
	# ;| della correzione della velocit� di variazione della B||.
	# ;| Lo 0 della retta ci d� una prima stima della correzione DtAngPar	

	a0 = np.mean(gradv)
	a1 = (gradv[1]-gradv[0])/(spb[1]-spb[0])
	DtAngPar=-a0/a1

	gradv1P=np.ones(6)
	DtAngParV=np.ones(6)

	## Calcola il gradiente di fase lungo l'Azimuth variando di
	## +/-00025m/s + DtAngPar con passo .0001m/s
	## la velocit� di variazione della B|| lungo la direzione di volo

	for i in range(6):
		# Sposta l'orbita slave lungo la B_|_ di DtB metri
		# ed incrementa di DtAngPar0 m/s la velocit� di variazione di B||
		
		DtAngPar0 = DtAngPar + (i*0.0001 - 0.00025)*65
		DtAngParV[i] = DtAngPar0
		secn0,p = deltaOrb(x,y,z,ref,sec,DtB,DtAngPar0)
		
		# Calcola nuovo differenziale con la nuova orbita slave
		rgsecn=rangeMat(x,y,z,0,secn0)
		sint1n0 = (rgref-rgsecn)*4*np.pi/0.23513133960784313
		Diffp =np.arctan2(np.sin(fas2-sint1n0), np.cos(fas2-sint1n0))
		
		#Calculo del gradiente en direcciión de azimut:
		Diffp_shift = np.roll(Diffp,1, axis = 0) # Lo muevo por filas
		gradd=np.arctan2(np.sin((Diffp-Diffp_shift)),np.cos((Diffp-Diffp_shift)))
		grad = gradd[(coh > gamma) & (~np.isnan(gradd))]
		gradv1P[i] = np.mean(grad)

	# Recordar que np.polyfit da los coeficientes en orden decreciente	
	resecP = np.polyfit(DtAngParV,gradv1P, deg=1)
	DtAngPar = -resecP[1]/resecP[0]

	# Calcola e visualizza il differenziale con la B_|_ e  DtAngPar corrette

	secn0,p = deltaOrb(x,y,z,ref,sec,DtB,DtAngPar)
	rgsecn=rangeMat(x,y,z,0,secn0)
	sint1n0 = (rgref-rgsecn)*4*np.pi/0.23513133960784313
	Diffp1 =np.arctan2(np.sin(fas2-sint1n0), np.cos(fas2-sint1n0))
	
	return Diffp1,DtB,DtAngPar
