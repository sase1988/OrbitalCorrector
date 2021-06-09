#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt



class OrbitalCorrector:
	def __init__(self, 
	             x, 
	             y, 
	             hgt, 
	             coherence, 
	             phase, 
	             reference, 
	             secondary, 
	             rangeLooks=1, 
	             azimuthLooks=1, 
	             gamma=0.25,
	             wavelength=0.23513133960784313
	             ):
		
		self.x = x
		self.y = y
		self.z = hgt
		self.coh = coherence
		self.phase = phase
		self.ref = reference
		self.sec = secondary
		self.rlook = rangeLooks
		self.alook = azimuthLooks
		self.gamma = gamma
		self.ml = 0 #If rlooks = 0 and alooks = 0
		self.wavelength = wavelength
		self.a = 6378137 #semi-major axis of the WGS84 ellipsoid
		self.b = 6356752.314 #semi-minor axis of the WGS84 ellipsoid
		self.phase_hr = None
		
		if self.rlook != 0 or self.alook != 0:
			self.ml = 1
		
	def rebin(self, array):
		rl=self.rlook
		al=self.alook
		high_res_rows = (array.shape[0]//al) * al
		high_res_cols = (array.shape[1]//rl) * rl
		low_res_rows = high_res_rows//al
		low_res_cols = high_res_cols//rl
		new_arr = array[:high_res_rows, :high_res_cols]
		sh = (low_res_rows, al, low_res_cols, rl)
		return new_arr.reshape(sh).mean(-1).mean(1)

	def crossProduct(self, a, b):
		
		c = np.zeros(a.shape)
		c[0] = a[1]*b[2] - a[2]*b[1]
		c[1] = -(a[0]*b[2] - a[2]*b[0])
		c[2] = a[0]*b[1] - a[1]*b[0]
		return c

	def orbitV(self, t, array, i):
		m = array[:,i]
		ss = m[0] + m[1]*t + m[2]*t**2 + m[3]*t**3 + m[4]*t**4 + m[5]*t**5 + m[6]*t**6 + m[7]*t**7 + m[8]*t**8 
		return ss

	def orbitP(self, t, x, y, z):
		N_stat = len(t)
		array = np.zeros((N_stat, 9))
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
		
		array[:,0] = (x.T @ av.T).T #rx1
		array[:,1] = (y.T @ av.T).T #ry1
		array[:,2] = (z.T @ av.T).T #rz1
		array[:-1,3] = np.arange(1,9) * array[1:,0] #rvx1
		array[:-1,4] = np.arange(1,9) * array[1:,1] #rvy1
		array[:-1,5] = np.arange(1,9) * array[1:,2] #rvz1
		array[:-1,6] = np.arange(1,9) * array[1:,3] #rax1
		array[:-1,7] = np.arange(1,9) * array[1:,4] #ray1
		array[:-1,8] = np.arange(1,9) * array[1:,5] #raz1

		return array, par


	def timeSincro(self, t0, x, y, z, array, maxit=100, tolerance=1e-8):

		delta = tolerance
		tt = np.zeros(len(array)-2)

		while np.min(abs(delta)) >= tolerance and maxit > 0:
			t0 = np.min(t0)
			tt[0] = t0
			for i in range(1, len(array)-2):
				globals()['t%d'%i] = eval('t%d'%(i-1)) * t0
				tt[i] = eval('t%d'%i)

			sxm = array[0,0] + tt @ array[1:-1,0]
			sym = array[0,1] + tt @ array[1:-1,1]
			szm = array[0,2] + tt @ array[1:-1,2]
			
			vxm = array[0,3] + tt @ array[1:-1,3]
			vym = array[0,4] + tt @ array[1:-1,4]
			vzm = array[0,5] + tt @ array[1:-1,5]
	 
			axm = array[0,6] + tt @ array[1:-1,6]
			aym = array[0,7] + tt @ array[1:-1,7]
			azm = array[0,8] + tt @ array[1:-1,8]

			delta = -(vxm*(x - sxm) + vym*(y - sym) + vzm*(z - szm)) / (axm*(x - sxm) + aym*(y - sym) + azm*(z - szm) - vxm**2 - vym**2 - vzm**2)
			
			t0 += delta
			maxit -= 1
		
		return t0

	def rangeMat0(self, x, y, z, t0, array):
		t = self.timeSincro(t0, x, y, z, array)
		range_0 = np.sqrt(
		                  (x - self.orbitV(t, array, 0))**2 + 
		                  (y - self.orbitV(t, array, 1))**2 + 
		                  (z - self.orbitV(t, array, 2))**2
		                 )
		return range_0 

	def rangeMat(self, x, y, z, t, array):
		
		#Leo el shape que me da fila y columna
		ranget = np.zeros((x.shape), float)
		t0 = self.timeSincro(t, x[0,:], y[0,:], z[0,:], array)
		msx, msy, msz = array[:,0], array[:,1], array[:,2]

		for i in range(x.shape[0]):	
			
			xt, yt, zt = x[i,:], y[i,:], z[i,:]
			t0 = self.timeSincro(t0, xt, yt, zt, array)
			
			for j in range(1, array.shape[0]-1):
				globals()['t'+str(j)] = eval('t' + str(j - 1)) * t0

			sxm = msx[0] + msx[1]*t0 + msx[2]*t1 + msx[3]*t2 + msx[4]*t3 + msx[5]*t4 + msx[6]*t5 + msx[7]*t6 + msx[8]*t7
			sym = msy[0] + msy[1]*t0 + msy[2]*t1 + msy[3]*t2 + msy[4]*t3 + msy[5]*t4 + msy[6]*t5 + msy[7]*t6 + msy[8]*t7
			szm = msz[0] + msz[1]*t0 + msz[2]*t1 + msz[3]*t2 + msz[4]*t3 + msz[5]*t4 + msz[6]*t5 + msz[7]*t6 + msz[8]*t7
			
			ranget[i,:] = np.sqrt((xt - sxm)**2 + (yt - sym)**2 + (zt - szm)**2)
			
		return ranget

	def coordEllipsoid(self, t, ro1, array, x, y, z, maxit=100, tolerance=1e-8):

		a = self.a
		b = self.b

		aa = np.ones((3,3))
		f = np.ones((3,1))

		sxm = self.orbitV(t, array, 0)
		sym = self.orbitV(t, array, 1)
		szm = self.orbitV(t, array, 2)

		vxm = self.orbitV(t, array, 3)
		vym = self.orbitV(t, array, 4)
		vzm = self.orbitV(t, array, 5)

		vmod = np.sqrt(vxm**2 + vym**2 + vzm**2)

		vxm = vxm/vmod
		vym = vym/vmod
		vzm = vzm/vmod

		delta = tolerance
		
		aa[1,0] = vxm 
		aa[1,1] = vym
		aa[1,2] = vzm
		
		while np.min(abs(delta)) >= tolerance and maxit > 0:

			ro1_temp = np.sqrt((x - sxm)**2 + (y - sym)**2 + (z - szm)**2)
			rr_app = np.sqrt(x**2 + y**2 + (a * z/b)**2)

			aa[0,0] = (x - sxm)/ro1_temp
			aa[0,1] = (y - sym)/ro1_temp
			aa[0,2] = (z - szm)/ro1_temp
			
			aa[2,0] = x/rr_app
			aa[2,1] = y/rr_app
			aa[2,2] = (a/b)**2 * z/rr_app

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

	def deltaOrb(self, Oref, Osec, DtB0, DtAngPar, x, y, z):

		DimYSar, DimXSar = x.shape
		ShAz = 10
		pOrbit = len(Osec[:,0])
		step = int((DimYSar - ShAz*2) / (pOrbit - 1))
		indi = np.arange(pOrbit)*step
		
		BperpV = np.zeros((pOrbit, 3))
		VrcenV = np.zeros((pOrbit, 3))
		MidRange = int(DimXSar/2)
		xm = x[:, MidRange]
		ym = y[:, MidRange]
		zm = z[:, MidRange]
		RgMid = np.mean(self.rangeMat0(xm, ym, zm, 0, Oref))
			
		################ Calcular los versores de los desplazamientos ##################
		# En este for lo que hace es ir moviendose en Azimut  y obtener el valor XYZ para el rango medio
		for i in range(pOrbit):
			Xmr = x[indi[i] + ShAz, MidRange]
			Ymr = y[indi[i] + ShAz, MidRange]
			Zmr = z[indi[i] + ShAz, MidRange]
			tcenref = self.timeSincro(0, Xmr, Ymr, Zmr, Oref)
			ElliPoi = self.coordEllipsoid(tcenref, RgMid, Oref, Xmr, Ymr, Zmr) #  Calcola corrispondenti punti su ellissoide
			Versref = np.array([self.orbitV(tcenref, Oref, 3), self.orbitV(tcenref, Oref, 4), self.orbitV(tcenref, Oref, 5)])
			Versref = Versref/np.sqrt(np.sum(Versref**2))
			satcen = np.array([self.orbitV(tcenref, Oref, 0), self.orbitV(tcenref, Oref, 1), self.orbitV(tcenref, Oref, 2)])
			rcen = ElliPoi - satcen
			VrcenV[i,:] = rcen/np.sqrt(np.sum(rcen**2))
			BperpV[i,:] = self.crossProduct(Versref, VrcenV[i,:])
		
		Tsn0 = self.timeSincro(0, x[indi+ShAz,MidRange], y[indi+ShAz,MidRange], z[indi+ShAz,MidRange], Osec)
		Tsn0m = Tsn0 - np.mean(Tsn0)
		
		Sxs0 = self.orbitV(Tsn0, Osec, 0)
		Sys0 = self.orbitV(Tsn0, Osec, 1)
		Szs0 = self.orbitV(Tsn0, Osec, 2)
		
		# Mueve los vectores de órbita de la Slave según las perturbaciones en Bperp y Bpar
		Sxsn0 = Sxs0 + DtB0*BperpV[:,0] + Tsn0m*DtAngPar*VrcenV[:,0]
		Sysn0 = Sys0 + DtB0*BperpV[:,1] + Tsn0m*DtAngPar*VrcenV[:,1]
		Szsn0 = Szs0 + DtB0*BperpV[:,2] + Tsn0m*DtAngPar*VrcenV[:,2]
		
		OsecN, parsn = self.orbitP(Tsn0, Sxsn0, Sysn0, Szsn0)
		
		return OsecN, parsn
		


class mainCorrector(OrbitalCorrector):
	def configureParams(self):
		self.constant = 4 * np.pi / self.wavelength
	
	def newDifferential(self, phase, rgref, sec, DtB, DtAngPar, x, y, z):
		
		secn0, p = self.deltaOrb(self.ref, sec, DtB, DtAngPar, x, y, z)
		rgsecn = self.rangeMat(x, y, z, 0, secn0)
		sint1n0 = (rgref - rgsecn) * self.constant
		Diffp = np.arctan2(np.sin(phase - sint1n0), np.cos(phase - sint1n0))
		
		return Diffp

	def grossEstimation(self, coh, phase, rgref, spb, baseline, x, y, z, DtB0=None, DtAngPar0=None):
		
		gradv= np.zeros(2)
		
		for i in range(2):
			if baseline=='perpendicular':
				DtB0 = spb[i]
				axis = 1 #Columnas
			if baseline=='parallel':
				DtAngPar0 = spb[i]
				axis = 0 #Filas
			
			# Calcula el nuevo diferencial con la nueva orbita esclava
			Diffp = self.newDifferential(phase, rgref, self.sec, DtB0, DtAngPar0, x, y, z)
			#Calculo del gradiente en dirección:
			#rango si es perpendicular, azimut si es paralela
			Diffp_shift = np.roll(Diffp, 1, axis=axis)
			gradd = np.arctan2(np.sin((Diffp - Diffp_shift)), np.cos((Diffp - Diffp_shift)))
			grad = gradd[(coh > self.gamma) & (~np.isnan(gradd))]
			
			if baseline=='perpendicular':
				#Usa el histograma para el cálculo del gradiente medio:
				binsiz = 0.2
				bins_n = int(np.ceil((np.max(grad) - np.min(grad))/binsiz))
				y_hist, x_hist, _ = plt.hist(grad.ravel(), bins=bins_n)
				res = np.poly1d(np.polyfit(np.arange(len(y_hist)-2), y_hist[1:-1], 4))
				tim = np.arange(1000) / 1000 * (len(y_hist)-2)
				histif1 = res(tim)
				gradv[i] = np.mean((np.where(histif1 == np.max(histif1))[0][0]*(len(y_hist)-2) / 1000. + 1.)*binsiz + np.min(grad))
	
			else:
				gradv[i] = np.mean(grad)
				
		a0 = np.mean(gradv)
		a1 = (gradv[1] - gradv[0])/(spb[1] - spb[0])

		return -a0/a1	
	
	####################################################################
	
	def fineEstimation(self, coh, phase, rgref, DtBv1, x, y, z, DtAngPar, baseline):
		
		gradv1 = np.zeros(6)
		DtBvv = np.zeros(6)

		for i in range(6):
			if baseline=='perpendicular':
				DtB0 = DtBv1 + i*2 - 5
				DtBvv[i] = DtB0
				axis = 1 #Columnas
				DtAngPar0 = DtAngPar
			if baseline=='parallel':
				DtB0 = DtBv1
				DtAngPar0 = DtAngPar + (i*0.0001 - 0.00025)*65
				DtBvv[i] = DtAngPar0
				axis = 0 #Filas

			# Calcula el nuevo diferencial con la nueva orbita esclava
			Diffp = self.newDifferential(phase, rgref, self.sec, DtB0, DtAngPar0, x, y, z)

			#Calculo del gradiente en dirección:
			#rango si es perpendicular, azimut si es paralela
			Diffp_shift = np.roll(Diffp, 1, axis=axis) 
			gradd = np.arctan2(np.sin((Diffp - Diffp_shift)),np.cos((Diffp - Diffp_shift)))
			grad = gradd[(coh > self.gamma) & (~np.isnan(gradd))]
			gradv1[i] = np.mean(grad)


		#Calcola la retta di regressione dei gradienti gradv1 al variare
		#della correzione di B_|_. Lo 0 ci dà la correzione DtB cercata

		# Recordar que np.polyfit da los coeficientes en orden decreciente
		resec = np.polyfit(DtBvv, gradv1, deg=1)
		return -resec[1]/resec[0]
		
	####################################################################	

	def runCorrector(self):
		self.configureParams()
		
		if self.ml:
			x = self.rebin(self.x.copy())
			y = self.rebin(self.y.copy())
			z = self.rebin(self.z.copy())
			coh = self.rebin(self.coh.copy())
			phase = self.rebin(self.phase.copy())
			rgref_hr = self.rangeMat(self.x, self.y, self.z, 0, self.ref)
			rgsec_hr = self.rangeMat(self.x, self.y, self.z, 0, self.sec)
			sint1_hr = (rgref_hr - rgsec_hr) * self.constant
			phase2_hr = np.arctan2(np.sin(-self.phase + sint1_hr), 
		                           np.cos(-self.phase + sint1_hr)
		                          )
			
		else:
			x = self.x.copy()
			y = self.y.copy()
			z = self.z.copy()
			coh = self.coh.copy()
			phase = self.phase.copy()
			
		rgref = self.rangeMat(x, y, z, 0, self.ref)
		rgsec = self.rangeMat(x, y, z, 0, self.sec)
		sint1 = (rgref - rgsec) * self.constant
		
		phase2 = np.arctan2(np.sin(-phase + sint1), 
		                    np.cos(-phase + sint1)
		                   )

		"""Estimación del delta Bperp"""
		print("Estimando delta Línea de base perpendicular")
		baseline = 'perpendicular'

		DtAngPar = 0
		SemIntDtB = 10
		spb = [-SemIntDtB, SemIntDtB]

		## Calcula el gradiente a lo largo del rango moviendo la orbita
		## Lo hace para el intervalo B_|_ +/- SemIntDtB
		DtBv1 = self.grossEstimation(coh, phase2, rgref, spb, baseline, x, y, z, DtAngPar0=DtAngPar)
		
		## Calcula el gradiente de fase a lo largo del rango variando la línea
		## base perpendicular en +/- 5m alrededor más la primera estimación DtBv1 
		
		DtB = self.fineEstimation(coh, phase2, rgref, DtBv1, x, y, z, DtAngPar, baseline)
		
		# Calcular el nuevo diferencial con la nueva órbita esclava
		# Acá podemos aplicar DtB al de baja y al de alta por separado

		Diffp1 = self.newDifferential(phase2, rgref, self.sec, DtB, DtAngPar, x, y, z)
		
		########################### PARALELA ###########################
		########## Ricerca su parallela +/- .0005*65 metri/sec #########

		print("Estimando delta Línea de base paralela")
		baseline = 'parallel'

		SemIntDtAngPar = 0.0005*65
		spb = [-SemIntDtAngPar, SemIntDtAngPar]

		## Calcola il gradiente di fase lungo l'Azimuth variando di +/- .0005 m/s 
		## la velocità di variazione della B|| lungo la direzione di volo

		DtAngPar = self.grossEstimation(coh, phase2, rgref, spb, baseline, x, y, z, DtB0=DtB)

		## Calcola il gradiente di fase lungo l'Azimuth variando di
		## +/-00025m/s + DtAngPar con passo .0001m/s
		## la velocit� di variazione della B|| lungo la direzione di volo

		DtAngPar = self.fineEstimation(coh, phase2, rgref, DtB, x, y, z, DtAngPar, baseline)
		
		# Calcola e visualizza il differenziale con la B_|_ e  DtAngPar corrette
		# Acá podemos aplicar DtB al de baja y al de alta por separado

		Diffp1 = self.newDifferential(phase2, rgref, self.sec, DtB, DtAngPar, x, y, z)
		
		if self.ml:
			Diffp1_hr = self.newDifferential(phase2_hr, 
			                                 rgref_hr, 
			                                 self.sec, 
			                                 DtB, 
			                                 DtAngPar, 
			                                 self.x, 
			                                 self.y, 
			                                 self.z
			                                )
			return Diffp1_hr, DtB, DtAngPar
			
		return Diffp1, DtB, DtAngPar
