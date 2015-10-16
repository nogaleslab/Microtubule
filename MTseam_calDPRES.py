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
	parser.add_option("--pf", dest="pf", type="int", metavar="int",
                help="protofilament number")
        parser.add_option("--fv",dest="fv",type="choice", metavar="['v8','v9']",
                choices=['v8','v9'],default='v9', help="input frealign format, default v9")

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


def getMTlist(l1):
	l2 = [x for x in l1 if x[0]!='C' and x!='\n']
	MTlist_dup = [int(x.split()[7]) for x in l2]
	MTlist = list(set(MTlist_dup))
	MTlist.sort()
	return MTlist_dup,MTlist

def grepMT(fname,MT):
	f1 = file(fname)
	l1 = f1.readlines()
        l2 = [x for x in l1 if x[0]!='C' and x!='\n']
        l3 = [x for x in l2 if int(x.split()[7])==MT]
	f1.close()
        return l3

def grepPtcl(fname,ptcl):
	f1 = file(fname)
	l1 = f1.readlines()
        l2 = [x for x in l1 if x[0]!='C' and x!='\n']
        l3 = l2[ptcl-1]
	f1.close()
        return l3

def findclosest(mylist,value):
	for i in mylist:
		if abs(value-i) < 5.0 or abs(value-i) >355.0:
			return mylist.index(i)
	return None

def mymax(mylist):
	mylist2 = [i for i in mylist if abs(i) < 50]
	if len(mylist2) > 0:
		return max(mylist2)
	else:
		return 0

def pickRotDPRES(mylist):
	if max(mylist) < 50.0 and max(mylist) > 0.1:
		return mylist.index(max(mylist))
	else:
		return 6

def pickRotLPRES(mylist):
	if min(mylist) < 91.0:
		return mylist.index(min(mylist))
	else:
		return 6

def pickRotHPRES(mylist):
        if min(mylist) < 91.0:
                return mylist.index(max(mylist))
        else:
                return 6
	

def CutbyPercentage(mylist,perc):
        import copy
        n = len(mylist)
        mylist2 = copy.copy(mylist)
        mylist2.sort()
        cutoff = mylist2[int(perc*n)-1]
        return cutoff

def MTreduce(l_MT):
	n = len(l_MT)
	rot_MT = 0
	adpres_MT = -1
	dpres_MT = 1
	for i in l_MT:
		t1 = i.split()
		MT = int(t1[1])
		rot = int(t1[3])
		dpres = float(t1[4])
		adpres = abs(dpres)
		if adpres > adpres_MT:
			rot_MT = rot
			adpres_MT = adpres
			dpres_MT = dpres
	return ("%d\t%d\t%8.2f\n"%(MT,rot_MT,dpres_MT))

# pick the rot with most weights
def MTreduce2(l_MT,pf):
	n = len(l_MT)
	rot_MT = 0
	dpreslist = [0]*(pf+1)
	for i in l_MT:
		t1 = i.split()
		MT = int(t1[1])
		rot = int(t1[3])
		dpres = float(t1[4])
		rot2 = rot+6
		dpreslist[rot2] += dpres

	adpreslist = [abs(i) for i in dpreslist]	

	rot2_MT = adpreslist.index(max(adpreslist))
	dpres_MT = dpreslist[rot2_MT]
	rot_MT = rot2_MT-6

	return ("%d\t%d\t%8.2f\n"%(MT,rot_MT,dpres_MT))

def MTreduce_master(filename,pf):
	f1 = file(filename)
	l1 = f1.readlines()
	n1 = len(l1)
	MTlist_dup = [int(x.split()[1]) for x in l1]
	MTlist = list(set(MTlist_dup))
	MTlist.sort()
	
	lastMT = int(l1[-1].split()[1])
	MTparlist = [[] for i in range(lastMT+1)]

	for i in l1:
		MT = int(i.split()[1])
		MTparlist[MT].append(i)
	
	fout = file("%s-reduced"%filename,"w")
	for MT in MTlist:
		#print "working on MT %d\t\r"%MT,
		#l_MT = grepMT(params["fpar"],MT)
		l_MT = MTparlist[MT]
		fout.write(MTreduce2(l_MT,pf))
	f1.close()
	fout.close()
	
def mainloop(params):	
	f1 = file("%s_rotPF0n"%params["fpar"])
	ll1 = f1.readlines()
	l1 = [x for x in ll1 if x[0]!='C' and x!='\n']
	n = len(l1)
	MTlist_dup,MTlist = getMTlist(l1)
	
	if params['fv'] == 'v9':
		pres0list = [float(x.split()[14]) for x in l1]
		pres0 = min(np.array(pres0list))
	elif params['fv'] == 'v8':
		pres0list = [float(x.split()[11]) for x in l1]
		pres0 = max(np.array(pres0list))

	if params['pf'] == 12:
		RISE = 10.44
		TWIST = -29.86
		grid = 5
		start = -5
		end = 6
	elif params['pf'] == 13:
		RISE = 9.5
		TWIST = -27.71
		grid = 6
		start = -6
		end = 6
	elif params['pf'] == 14:
		RISE = 8.85
		TWIST = -25.76
		grid = 6
		start = -6
		end = 6
	else:
		print "please specify PF"
		sys.exit()
	
	pf = params['pf']
	
	# initialize phi_ref for each ptcl
	phi_ref = {}
	# phi_ref[im] = [0,27,54,...]
	for im in range(n):
		phi = float(l1[im].split()[3])
		phi_ref[im] = [999 for x in xrange(end-start+1)]
		
		for rot in range(start,end+1):
			rot2 = rot + grid
			# rot2 starts from 0, rot starts from -6
			phi_ref[im][rot2] = phi + -1*rot*TWIST
			phi_ref[im][rot2] %= 360.0
	
	# initialize presdict for each particle
	presdict = {}
	# presdict[im] = [[75,74],[75,74],[75,74]...]
	
	for im in range(n):
		presdict[im] = [0 for i in xrange(end-start+1)]
		for rot in range(start,end+1):
			rot2 = rot + grid
			#presdict[im][rot2] = [444,999]
			presdict[im][rot2] = [pres0,pres0]
	print "finish initialization"	

	# assign values
	for rot in range(start,end+1):
		rot2 = rot + grid
		f1 = file("%s_rotPF%dn"%(params["fpar"],rot))
		ll1 = f1.readlines()
		l1 = [x for x in ll1 if x[0]!='C' and x!='\n']
		f2 = file("%s_rotPF%ds"%(params["fpar"],rot))
		ll2 = f2.readlines()
		l2 = [x for x in ll2 if x[0]!='C' and x!='\n']
		
		for im in range(n):
			MT = MTlist_dup[im]
			
			t1 = l1[im].split()
			phi1 = float(t1[3])
			theta1 = float(t1[2])
			if params['fv'] == 'v9':
				pres1 = float(t1[14])
			elif params['fv'] == 'v8':
				pres1 = float(t1[11])

			t2 = l2[im].split()
			phi2 = float(t2[3])
			theta2 = float(t2[2])
			if params['fv'] == 'v9':
				pres2 = float(t2[14])
			elif params['fv'] == 'v8':
				pres2 = float(t2[11])	
			
			if theta1 > 110 or theta1 < 70 or theta2 > 110 or theta2 < 70:
				continue
	
			index1 = findclosest(phi_ref[im],phi1)
			if index1 != None: 
				#presdict[im][index1][0] = min(pres1,presdict[im][index1][0])
				presdict[im][index1][0] = pres1

			index2 = findclosest(phi_ref[im],phi2)
			if index2 != None:
				#presdict[im][index2][1] = min(pres2,presdict[im][index2][1])
				presdict[im][index2][1] = pres2
		f1.close()
		f2.close()
	
	# initialize dpreslist[im][rot2]
	dpreslist = [[0 for i in xrange(end-start+1)] for j in xrange(n)]
	absdpreslist = [[0 for i in xrange(end-start+1)] for j in xrange(n)]

	fout1 = file('perptcl_stats.txt',"w")
	for im in range(n):
		MT = MTlist_dup[im]
		for rot in range(start,end+1):
			rot2 = rot + grid
			pres1 = presdict[im][rot2][0]
			pres2 = presdict[im][rot2][1]
			dpres = pres1 - pres2
			if abs(dpres) > 50:
				dpres = 0
			dpreslist[im][rot2] = dpres
			if params['fv'] == 'v8':
				dpres *= -1
			absdpreslist[im][rot2] = abs(dpres)
			fout1.write("%d\tMT %d\t%d\t%8.2f\t%8.2f\t%8.2f%8.2f\n"%(im,MT,rot,phi_ref[im][rot2],dpres,pres1,pres2))
		fout1.write("\n")

	maxdpres = [max(x) for x in absdpreslist]
	maxdpres_pick = [pickRotDPRES(x) for x in absdpreslist]
	dpres_mean = np.mean(maxdpres)
	print "mean dpres = %.2f"%dpres_mean
	dpres_std = np.std(maxdpres)
	print "std dpres = %.2f"%dpres_std
	#cutoff = dpres_mean - 1.0*dpres_std
	cutoff = CutbyPercentage(maxdpres,0.2)
        print "discard 20 percent of MTs, dpres cutoff = %.2f"%cutoff
	fout1.close()

	fout2 = file('goodMT_cutoff%.1f.txt'%cutoff,"w")
        fout11 = file("pick_dpres_5c.txt","w")
        fout12 = file("pick_dpres_3c.txt","w")
        for im in range(n):
                MT = MTlist_dup[im]
                rotpick = maxdpres_pick[im]
                phi = phi_ref[im][rotpick]
                fout11.write("%d\t%d\t%.2f\t%d\t%.2f\n"%(im+1,MT,phi,rotpick-grid,dpreslist[im][rotpick]))
                fout12.write("%d\t%d\t%.2f\n"%(MT,rotpick-grid,dpreslist[im][rotpick]))
                if maxdpres[im] >= cutoff:
                        fout2.write("%d\t%.2f\t%.2f\n"%(MT,phi,dpreslist[im][rotpick]))
        fout11.close()
        fout12.close()
        fout2.close()


	MTreduce_master('pick_dpres_5c.txt',pf)

if __name__ == "__main__":
	params = setupParserOptions()
	mainloop(params)
