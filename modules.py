#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt
import warnings
warnings.filterwarnings('ignore')

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
	             wavelength=0.23513133960784313,
	             ramp = False
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
		self.ramp = ramp
		self.ml = 0 #If rlooks = 0 and alooks = 0
		self.wavelength = wavelength
		self.a = 6378137 #semi-major axis of the WGS84 ellipsoid
		self.b = 6356752.314 #semi-minor axis of the WGS84 ellipsoid

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

	def orbitV(self, t, polynomials, i):
		m = polynomials[:,i]
		orbitValue = m[0] + m[1]*t + m[2]*t**2 + m[3]*t**3 + m[4]*t**4 + m[5]*t**5 + m[6]*t**6 + m[7]*t**7 + m[8]*t**8 
		return orbitValue

	def orbitP(self, t, x, y, z):
		N_stat = len(t)
		polynomials = np.zeros((N_stat, 9))
		par = [0, 0, 0, 0, 0, 0, N_stat, t.copy(), x, y, z]
		
		#Input scene parameters:
		par[0] = t[0]
		par[1] = t[0]
		par[2] = t[1] - t[0]
		par[5] = np.sum(t)/len(t) #t_mid_m
		par[3] = par[0] - par[5] #tstart
		
		t -= par[5]
		
		a = np.array([t[i]**j for i in range(N_stat) for j in range(N_stat)])
		a = a.reshape(N_stat,N_stat)
		av = np.linalg.inv(a)
		
		polynomials[:,0] = (x.T @ av.T).T #rx1
		polynomials[:,1] = (y.T @ av.T).T #ry1
		polynomials[:,2] = (z.T @ av.T).T #rz1
		polynomials[:-1,3] = np.arange(1,9) * polynomials[1:,0] #rvx1
		polynomials[:-1,4] = np.arange(1,9) * polynomials[1:,1] #rvy1
		polynomials[:-1,5] = np.arange(1,9) * polynomials[1:,2] #rvz1
		polynomials[:-1,6] = np.arange(1,9) * polynomials[1:,3] #rax1
		polynomials[:-1,7] = np.arange(1,9) * polynomials[1:,4] #ray1
		polynomials[:-1,8] = np.arange(1,9) * polynomials[1:,5] #raz1

		return polynomials, par


	def timeSincro(self, t0, x, y, z, polynomials, maxit=100, tolerance=1e-8):

		delta = tolerance
		tt = np.zeros(len(polynomials)-2)

		while np.min(abs(delta)) >= tolerance and maxit > 0:
			t0 = np.min(t0)
			tt[0] = t0
			for i in range(1, len(polynomials)-2):
				globals()['t%d'%i] = eval('t%d'%(i-1)) * t0
				tt[i] = eval('t%d'%i)

			sx = polynomials[0,0] + tt @ polynomials[1:-1,0]
			sy = polynomials[0,1] + tt @ polynomials[1:-1,1]
			sz = polynomials[0,2] + tt @ polynomials[1:-1,2]
																										 
			vx = polynomials[0,3] + tt @ polynomials[1:-1,3]
			vy =polynomials[0,4] + tt @ polynomials[1:-1,4]
			vz = polynomials[0,5] + tt @ polynomials[1:-1,5]
																										 
			ax = polynomials[0,6] + tt @ polynomials[1:-1,6]
			ay =polynomials[0,7] + tt @ polynomials[1:-1,7]
			az = polynomials[0,8] + tt @ polynomials[1:-1,8]

			delta = -(vx*(x - sx) + vy*(y - sy) + vz*(z - sz)) / (ax*(x - sx) + ay*(y - sy) + az*(z - sz) - vx**2 - vy**2 - vz**2)
			
			t0 += delta
			maxit -= 1
		
		return t0

	def rangeMat0(self, x, y, z, t0, polynomials):
		t = self.timeSincro(t0, x, y, z, polynomials)
		range_0 = np.sqrt(
		                  (x - self.orbitV(t, polynomials, 0))**2 + 
		                  (y - self.orbitV(t, polynomials, 1))**2 + 
		                  (z - self.orbitV(t, polynomials, 2))**2
		                 )
		return range_0 

	def rangeMat(self, x, y, z, t, polynomials):
		
		ranget = np.zeros((x.shape), float)
		t0 = self.timeSincro(t, x[0,:], y[0,:], z[0,:], polynomials)
		msx, msy, msz = polynomials[:,0], polynomials[:,1], polynomials[:,2]

		for i in range(x.shape[0]):	
			
			xt, yt, zt = x[i,:], y[i,:], z[i,:]
			t0 = self.timeSincro(t0, xt, yt, zt, polynomials)
			
			for j in range(1, polynomials.shape[0]-1):
				globals()['t'+str(j)] = eval('t' + str(j - 1)) * t0

			sx = msx[0] + msx[1]*t0 + msx[2]*t1 + msx[3]*t2 + msx[4]*t3 + msx[5]*t4 + msx[6]*t5 + msx[7]*t6 + msx[8]*t7
			sy = msy[0] + msy[1]*t0 + msy[2]*t1 + msy[3]*t2 + msy[4]*t3 + msy[5]*t4 + msy[6]*t5 + msy[7]*t6 + msy[8]*t7
			sz = msz[0] + msz[1]*t0 + msz[2]*t1 + msz[3]*t2 + msz[4]*t3 + msz[5]*t4 + msz[6]*t5 + msz[7]*t6 + msz[8]*t7
			
			ranget[i,:] = np.sqrt((xt - sx)**2 + (yt - sy)**2 + (zt - sz)**2)
			
		return ranget

	def coordEllipsoid(self, t, ro1, polynomials, x, y, z, maxit=100, tolerance=1e-8):

		a = self.a
		b = self.b

		aa = np.ones((3,3))
		f = np.ones((3,1))

		sx = self.orbitV(t, polynomials, 0)
		sy = self.orbitV(t, polynomials, 1)
		sz = self.orbitV(t, polynomials, 2)

		vx = self.orbitV(t, polynomials, 3)
		vy = self.orbitV(t, polynomials, 4)
		vz = self.orbitV(t, polynomials, 5)

		vmod = np.sqrt(vx**2 + vy**2 + vz**2)

		vx = vx/vmod
		vy = vy/vmod
		vz = vz/vmod

		delta = tolerance
		
		aa[1,0] = vx 
		aa[1,1] = vy
		aa[1,2] = vz
		
		while np.min(abs(delta)) >= tolerance and maxit > 0:

			ro1_temp = np.sqrt((x - sx)**2 + (y - sy)**2 + (z - sz)**2)
			rr_app = np.sqrt(x**2 + y**2 + (a * z/b)**2)

			aa[0,0] = (x - sx)/ro1_temp
			aa[0,1] = (y - sy)/ro1_temp
			aa[0,2] = (z - sz)/ro1_temp
			
			aa[2,0] = x/rr_app
			aa[2,1] = y/rr_app
			aa[2,2] = (a/b)**2 * z/rr_app

			f[0] = ro1_temp - ro1
			f[1] = (x - sx)*vx + (y - sy)*vy + (z - sz)*vz
			f[2] = rr_app - a

			aai = np.linalg.inv(aa)
			delta = -aai @ f
			x = x + delta[0]
			y = y + delta[1]
			z = z + delta[2]
			  
			maxit -= 1

		return np.hstack((x, y, z))

	def deltaOrb(self, refPoly, secPoly, deltaPerp, deltaPar, x, y, z):

		dimYSar, dimXSar = x.shape
		ShAz = 10
		pOrbit = len(secPoly[:,0])
		step = int((dimYSar - ShAz*2) / (pOrbit - 1))
		indi = np.arange(pOrbit)*step
		
		deltaPerpVect = np.zeros((pOrbit, 3))
		vrCenVect = np.zeros((pOrbit, 3))
		midRange = int(dimXSar/2)
		xm = x[:, midRange]
		ym = y[:, midRange]
		zm = z[:, midRange]
		rgMid = np.mean(self.rangeMat0(xm, ym, zm, 0, refPoly))
			
		# Move along azimuth and get the XYZ values por the mid-range
		for i in range(pOrbit):
			Xmr = x[indi[i] + ShAz, midRange]
			Ymr = y[indi[i] + ShAz, midRange]
			Zmr = z[indi[i] + ShAz, midRange]
			tcenref = self.timeSincro(0, Xmr, Ymr, Zmr, refPoly)
			elliPos = self.coordEllipsoid(tcenref, rgMid, refPoly, Xmr, Ymr, Zmr) #  Compute the coordinates of these XYZ points on the ellipsoid
			Versref = np.array([self.orbitV(tcenref, refPoly, 3), self.orbitV(tcenref, refPoly, 4), self.orbitV(tcenref, refPoly, 5)])
			Versref = Versref/np.sqrt(np.sum(Versref**2))
			satcen = np.array([self.orbitV(tcenref, refPoly, 0), self.orbitV(tcenref, refPoly, 1), self.orbitV(tcenref, refPoly, 2)])
			rcen = elliPos - satcen
			vrCenVect[i,:] = rcen/np.sqrt(np.sum(rcen**2))
			deltaPerpVect[i,:] = self.crossProduct(Versref, vrCenVect[i,:])
		
		Tsn0 = self.timeSincro(0, x[indi+ShAz,midRange], y[indi+ShAz,midRange], z[indi+ShAz,midRange], secPoly)
		Tsn0m = Tsn0 - np.mean(Tsn0)
		
		Sxs0 = self.orbitV(Tsn0, secPoly, 0)
		Sys0 = self.orbitV(Tsn0, secPoly, 1)
		Szs0 = self.orbitV(Tsn0, secPoly, 2)
		
		# Perturbate the slave polynomials according to deltaPerp and deltaPar
		Sxsn0 = Sxs0 + deltaPerp*deltaPerpVect[:,0] + Tsn0m*deltaPar*vrCenVect[:,0]
		Sysn0 = Sys0 + deltaPerp*deltaPerpVect[:,1] + Tsn0m*deltaPar*vrCenVect[:,1]
		Szsn0 = Szs0 + deltaPerp*deltaPerpVect[:,2] + Tsn0m*deltaPar*vrCenVect[:,2]
		
		secPolyN, parsn = self.orbitP(Tsn0, Sxsn0, Sysn0, Szsn0)
		
		return secPolyN, parsn
		
class mainCorrector(OrbitalCorrector):
	def configureParams(self):
		self.constant = 4 * np.pi / self.wavelength
	
	def newDifferential(self, phase, rgref, sec, deltaPerp, deltaPar, x, y, z):
		
		#Compute the new orbital polynomials
		secn0, p = self.deltaOrb(self.ref, sec, deltaPerp, deltaPar, x, y, z)
		#Compute the new range matrix
		rgsecn = self.rangeMat(x, y, z, 0, secn0)
		#Compute the new synthetic file
		sint = (rgref - rgsecn) * self.constant
		#Compute the new differential interferogram
		Diffp = np.arctan2(np.sin(phase - sint), np.cos(phase - sint))
		
		return Diffp

	def calcRamp(self,phase,phase_n):
		ramp = np.arctan2(np.sin(phase_n - phase), np.cos(phase_n - phase))
		return ramp

	def coarseEstimation(self, coh, phase, rgref, deltaList, baseline, x, y, z, deltaPerp0=None, deltaPar0=None):
		
		gradVector= np.zeros(2)
		
		for i in range(2):
			if baseline=='perpendicular':
				deltaPerp0 = deltaList[i]
				axis = 1 #Columns
			if baseline=='parallel':
				deltaPar0 = deltaList[i]
				axis = 0 #Rows
			
			# Compute the new differential
			newPhase = self.newDifferential(phase, rgref, self.sec, deltaPerp0, deltaPar0, x, y, z)
			
			#Compute the phase gradient:
			newPhase_shift = np.roll(newPhase, 1, axis=axis)
			gradd = np.arctan2(np.sin((newPhase - newPhase_shift)), np.cos((newPhase - newPhase_shift)))
			grad = gradd[(coh > self.gamma) & (~np.isnan(gradd))]
			
			if baseline=='perpendicular':
				#Use the histogram to get the maximum gradient:
				binsiz = 0.2
				bins_n = int(np.ceil((np.max(grad) - np.min(grad))/binsiz))
				y_hist, x_hist, _ = plt.hist(grad.ravel(), bins=bins_n)
				polyFit = np.poly1d(np.polyfit(np.arange(len(y_hist)-2), y_hist[1:-1], 4))
				xValues = np.arange(1000) / 1000 * (len(y_hist)-2)
				fittedHist = polyFit(xValues)
				gradVector[i] = np.mean((np.where(fittedHist == np.max(fittedHist))[0][0]*(len(y_hist)-2) / 1000. + 1.)*binsiz + np.min(grad))
				plt.close()
				
	
			else:
				#If processing the parallel baseline, directly compute the mean of the gradient
				gradVector[i] = np.mean(grad)
		
		#Compute the coefficients of the secant line
		a0 = np.mean(gradVector)
		a1 = (gradVector[1] - gradVector[0])/(deltaList[1] - deltaList[0])

		return -a0/a1	
	
	####################################################################
	
	def fineEstimation(self, coh, phase, rgref, deltaPerp1, x, y, z, deltaPar1, baseline):
		
		gradVector = np.zeros(6)
		deltaList = np.zeros(6)

		for i in range(6):
			if baseline=='perpendicular':
				deltaPerp0 = deltaPerp1 + i*2 - 5
				deltaList[i] = deltaPerp0
				axis = 1 #Columns
				deltaPar0 = deltaPar1
			if baseline=='parallel':
				deltaPerp0 = deltaPerp1
				deltaPar0 = deltaPar1 + (i*0.0001 - 0.00025)*65
				deltaList[i] = deltaPar0
				axis = 0 #Rows

			#Compute the new differential
			newPhase = self.newDifferential(phase, rgref, self.sec, deltaPerp0, deltaPar0, x, y, z)

			#Compute the phase gradient:
			newPhase_shift = np.roll(newPhase, 1, axis=axis) 
			gradd = np.arctan2(np.sin((newPhase - newPhase_shift)),np.cos((newPhase - newPhase_shift)))
			grad = gradd[(coh > self.gamma) & (~np.isnan(gradd))]
			gradVector[i] = np.mean(grad)
			
		#Compute the 1D function that fits the gradient-DeltaB points, and its coefficients
		resec = np.polyfit(deltaList, gradVector, deg=1)
		return -resec[1]/resec[0]
		
	####################################################################	

	def runCorrector(self):
		self.configureParams()
		
		if self.ml:
			x = self.rebin(self.x.copy())
			y = self.rebin(self.y.copy())
			z = self.rebin(self.z.copy())
			coh = self.rebin(self.coh.copy())
			phase_hr =  self.phase.copy()
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
		sint = (rgref - rgsec) * self.constant
		
		phase2 = np.arctan2(np.sin(-phase + sint), 
		                    np.cos(-phase + sint)
		                   )

			
		########################### PERPENDICULAR ###########################
		
		print("Computing deltaPerp")
		baseline = 'perpendicular'

		deltaPar = 0
		delta = 10
		deltaList = [-delta, delta]

		coarseDeltaPerp = self.coarseEstimation(coh, phase2, rgref, deltaList, baseline, x, y, z, deltaPar0=deltaPar)
		fineDeltaPerp = self.fineEstimation(coh, phase2, rgref, coarseDeltaPerp, x, y, z, deltaPar, baseline)

		#Compute the new interferogram with the fine correction of Bperp
		#newPhase_perp = self.newDifferential(phase2, rgref, self.sec, fineDeltaPerp, deltaPar, x, y, z)
			
		########################### PARALLEL ##########################
		
		print("Computing deltaPar")
		baseline = 'parallel'

		delta = 0.0005*65
		deltaList = [-delta, delta]

		coarseDeltaPar = self.coarseEstimation(coh, phase2, rgref, deltaList, baseline, x, y, z, deltaPerp0=fineDeltaPerp)
		fineDeltaPar = self.fineEstimation(coh, phase2, rgref, fineDeltaPerp, x, y, z, coarseDeltaPar, baseline)
		
		correctedPhase = self.newDifferential(phase2, rgref, self.sec, fineDeltaPerp, fineDeltaPar, x, y, z)
		estRamp = self.calcRamp(phase, correctedPhase)
		
		if self.ml:
			correctedPhase = self.newDifferential(phase2_hr, 
			                                 rgref_hr, 
			                                 self.sec, 
			                                 fineDeltaPerp, 
			                                 fineDeltaPar, 
			                                 self.x, 
			                                 self.y, 
			                                 self.z
			                                )
			estRamp = self.calcRamp(phase_hr, correctedPhase)
			
		if self.ramp:
			return correctedPhase, fineDeltaPerp, fineDeltaPar, estRamp
		else:
			return correctedPhase, fineDeltaPerp, fineDeltaPar
