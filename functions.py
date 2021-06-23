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



def read_params_gmt(ref):

	prm = open(ref + '.PRM','r')
	
	lines = []
	for line in prm:
		if line.startswith('nrows'):
				lines.append(line.rstrip('\n'))
		if line.startswith('PRF'):
			lines.append(line.rstrip('\n'))
		if line.startswith('SC_clock_start'):
			lines.append(line.rstrip('\n'))
			
	prm.close()
			
	nrows = int(lines[0].split('=')[1])
	prf= float(lines[1].split('=')[1])
	SC_clock_start = float(lines[2].split('=')[1])
	
	s0_image = math.trunc((SC_clock_start - np.floor(SC_clock_start))*86400)
	dt_image = round(nrows/prf)
	
	led =  open(ref + '.LED', 'r')
	headers = next(led)

	s = []
	x = []
	y = []
	z = []
	vx = []
	vy = []
	vz = []
	
	for line in led:
			s.append(int(float(line.split(' ')[2])))
			x.append(float(line.split(' ')[3]))
			y.append(float(line.split(' ')[4]))
			z.append(float(line.split(' ')[5]))
			vx.append(float(line.split(' ')[6]))
			vy.append(float(line.split(' ')[7]))
			vz.append(float(line.split(' ')[8]))
					
	led.close()
	
	s = np.array(s)
			
	# Form a subset of 9 orbital positions separated every 2 seconds, starting from the image t0:

	nn = 9
	pos_i = np.where(s == s0_image)[0][0]

	s_orb = np.zeros(nn)
	x_orb = np.zeros(nn)
	y_orb = np.zeros(nn)
	z_orb = np.zeros(nn)

	for i in range(nn):
		s_orb[i] = s[pos_i + 2*i]
		x_orb[i] = x[pos_i + 2*i]
		y_orb[i] = y[pos_i + 2*i]
		z_orb[i] = z[pos_i + 2*i]

	t_orb = s_orb - np.mean(s_orb)
	orbitPolynomials, par = orbitP(t_orb, x_orb, y_orb, z_orb)

	return orbitPolynomials, par

def orbitP(t, x, y, z):
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
