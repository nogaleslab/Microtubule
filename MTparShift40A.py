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
	parser.set_usage("%prog -f <EB3clK_5.par>")
	parser.add_option("-f",dest="fpar",type="string",metavar="FILE",
		help="13pf_1_r1.par")

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

def shift40A(params):
	FORMAT = "%7d%8.2f%8.2f%8.2f%10.2f%10.2f%8d%6d%9.1f%9.1f%8.2f%8.2f%10d%11.4f%8.2f%8.2f\n"
	MONOMER = 41.0
	DIMER = 82.0
	fout = file('%s-shift40A'%(params['fpar']),"w")
	f1 = file(params['fpar'])
        ll1 = f1.readlines()
        l1 = [x for x in ll1 if x[0]!='C' and x!='\n']
	for i in l1:
                t1 = i.split()
                count = float(t1[0])
                psi = float(t1[1])
                theta = float(t1[2])
                phi = float(t1[3])
                shx = float(t1[4])
                shy = float(t1[5])
                mag = float(t1[6])
                micro = float(t1[7])
                df1 = float(t1[8])
                df2 = float(t1[9])
                angast = float(t1[10])
                occ = float(t1[11])
                logP = float(t1[12])
                sigma = float(t1[13])
                score = float(t1[14])
                change = float(t1[15])
		
		shx_v1 = shx + MONOMER*math.cos(-psi/180.0*math.pi)
		shy_v1 = shy + MONOMER*math.sin(-psi/180.0*math.pi)		
		shx_v2 = shx - MONOMER*math.cos(-psi/180.0*math.pi)
		shy_v2 = shy - MONOMER*math.sin(-psi/180.0*math.pi)
		
		if max(abs(shx_v1),abs(shy_v1)) < max(abs(shx_v1),abs(shy_v2)):
                        shx = shx_v1
                        shy = shy_v1
		else:
                        shx = shx_v2
                        shy = shy_v2

		shx,shy = closer2origin(shx,shy,psi)
		fout.write(FORMAT%(count,psi,theta,phi,shx,shy,mag,micro,df1,df2,angast,occ,logP,sigma,score,change))
	
	f1.close()
	fout.close()

if __name__ == "__main__":
	params = setupParserOptions()
	shift40A(params)
