#!/usr/bin/env python

import optparse
import itertools
import math
import os,sys
import glob
#import heapq
import numpy as np

def setupParserOptions():
	parser = optparse.OptionParser()
	parser.add_option("-f",dest="fpar",type="string",metavar="FILE",
		help="13pf_1_r1.par")
	parser.add_option("--grid",dest="grid",type="int",metavar="INT",default=6,
		help="generate a series of par files from -grid to grid")
	parser.add_option("--pf", dest="pf", type="int", metavar="int",
                help="protofilament number")
	options,args = parser.parse_args()

	if len(args) > 1:
		parser.error("Unknown commandline options: " +str(args))

	if len(sys.argv) < 2:
		parser.print_help()
		sys.exit()

	params={}

	for i in parser.option_list:
		if isinstance(i.dest,str):
			params[i.dest] = getattr(options,i.dest)
	return params

def shift40A(shx,shy,psi):
	MONOMER = 41.0
	shx_v1 = shx + MONOMER*math.cos(-psi/180.0*math.pi)
        shy_v1 = shy + MONOMER*math.sin(-psi/180.0*math.pi)
        shx_v2 = shx - MONOMER*math.cos(-psi/180.0*math.pi)
        shy_v2 = shy - MONOMER*math.sin(-psi/180.0*math.pi)
	if max(abs(shx_v1),abs(shy_v1)) < max(abs(shx_v2),abs(shy_v2)):
                shx = shx_v1
                shy = shy_v1
        else:
                shx = shx_v2
                shy = shy_v2
	return shx,shy

def closer2origin(shx,shy,psi):
	DIMER = 82.0
	shx_new1 = shx + DIMER*math.cos(-psi/180.0*math.pi)
	shy_new1 = shy + DIMER*math.sin(-psi/180.0*math.pi)		
	shx_new2 = shx - DIMER*math.cos(-psi/180.0*math.pi)
	shy_new2 = shy - DIMER*math.sin(-psi/180.0*math.pi)
	if max(abs(shx_new1),abs(shy_new1)) < max(abs(shx),abs(shy)):
		shx = shx_new1
		shy = shy_new1
	if max(abs(shx_new2),abs(shy_new2)) < max(abs(shx),abs(shy)):
		shx = shx_new2
		shy = shy_new2
	return shx,shy
	
def rotPF(params):
	PF = params['pf']
        MONOMER = 41.0
        PF = params['pf']
        RISE = MONOMER*3/PF
        TWIST = -360./PF
	FORMAT = "%7d%8.2f%8.2f%8.2f%10.2f%10.2f%8d%6d%9.1f%9.1f%8.2f%8.2f%10d%11.4f%8.2f%8.2f\n"
	grid = params["grid"]
	f1 = file(params["fpar"])
	ll1 = f1.readlines()
	l1 = [x for x in ll1 if x[0]!='C' and x!='\n']
	for pf in range(-grid,grid+1):
		foutn = file("%s_rotPF%dn"%(params["fpar"],-pf),"w")
		fouts = file("%s_rotPF%ds"%(params["fpar"],-pf),"w")
		for i in l1:
			l = i.split()
		        count = float(l[0])
		        psi = float(l[1])
		        theta = float(l[2])
		        phi = float(l[3])
		        shx = float(l[4])
		        shy = float(l[5])
		        mag = float(l[6])
		        micro = float(l[7])
		        df1 = float(l[8])
		        df2 = float(l[9])
		        angast = float(l[10])
			occ = float(l[11])
			logP = float(l[12])
			sigma = float(l[13])
			score = float(l[14])
			change = float(l[15])
			phi += TWIST*pf
			shx += pf*RISE*math.cos(-psi/180.0*math.pi)
			shy += pf*RISE*math.sin(-psi/180.0*math.pi)
			shx,shy = closer2origin(shx,shy,psi)
			foutn.write(FORMAT%(count,psi,theta,phi,shx,shy,mag,micro,df1,df2,angast,occ,logP,sigma,score,change))
			
			shx,shy = shift40A(shx,shy,psi)
			shx,shy = closer2origin(shx,shy,psi)
			fouts.write(FORMAT%(count,psi,theta,phi,shx,shy,mag,micro,df1,df2,angast,occ,logP,sigma,score,change))
		
	foutn.close()
	fouts.close()
	f1.close()

if __name__ == "__main__":
	params = setupParserOptions()
	rotPF(params)
