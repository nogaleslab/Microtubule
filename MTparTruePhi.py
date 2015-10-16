#!/usr/bin/env python
import optparse
import os,sys
import numpy as np
import math
import itertools

def setupParserOptions():
	parser = optparse.OptionParser()
	parser.set_usage("%prog -f <EB3clK_5.par>")
	parser.add_option("-f",dest="fpar",type="string",metavar="FILE",
		help="EB3clK_5.par")
	parser.add_option("--apix", dest="apix", type="float", metavar="FLOAT", default=1,
                help="pixel size in angstroms, should always be 1")
	parser.add_option("--pf", dest="pf", type="int", metavar="int", default=13,
                help="PF")
	parser.add_option("--break", action="store_true",dest="break",default=False,
		help="break")
	parser.add_option("--forceline", action="store_true",dest="forceline",default=False,
		help="forceline")
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

def GetRotfromfile(MT,philist):
	f9 = file('picklist2.txt')
	l9 = f9.readlines()
	rotdict = {}
	signdict = {}
	for i in l9:
		t9 = i.split()
		micro = int(t9[0])
		rotdict[micro] = float(t9[1])
		try:
			signdict[micro] = math.copysign(1,float(t9[-1]))
		except:
			signdict[micro] = 1
	try:
		rot_MT = phidict[MT]
		sign_MT = signdict[MT]
	except:
		rot_MT = mymedian(philist)
		sign_MT = 1
	return phi_MT,sign_MT
	#return phi_MT,1

def mymedian(mylist):
	import copy
	n = len(mylist)
	if n%2 == 1:
		return np.median(mylist)
	else:
		mylist2 = copy.copy(mylist)
		mylist2.sort()
		return np.median(mylist2[1:])

def getMTlist(params):
	f1 = file(params["fpar"])
	l1 = f1.readlines()
	l2 = [x for x in l1 if x[0]!='C' and x!='\n']
	MTlist_dup = [int(x.split()[7]) for x in l2]
	MTlist = list(set(MTlist_dup))
	MTlist.sort()
	f1.close()
	return MTlist

def findMinShxy(shxlist,shylist):
	shxlistn = np.array(shxlist)
	shylistn = np.array(shylist)
	shxylistn = abs(shxlistn)+abs(shylistn)
	# find the location when shx and shy cross each other	
	dshxylistn = abs(shxlistn - shylistn)
	#shxylist = list(shxylistn)
	dshxylist = list(dshxylistn)
	#return shxylist.index(min(shxylist))
	return dshxylist.index(min(dshxylist))

def findClosestShxy(shx,shy,shx_ref,shy_ref,psi,apix):
	SHIFTS = np.array([164./apix,-164./apix,123./apix,-123./apix,82./apix,-82./apix,41./apix,-41./apix])
	dx = abs(shx - shx_ref)
	dy = abs(shy - shy_ref)
	shx_cand = shx + SHIFTS*math.cos(-psi/180.0*math.pi)
	shy_cand = shy + SHIFTS*math.sin(-psi/180.0*math.pi)
	dx_cand = abs(shx_cand - shx_ref)
	dy_cand = abs(shy_cand - shy_ref)
	dtotal = dx_cand + dy_cand
	dtotal = list(dtotal)
	if min(dtotal) < (dx+dy):
		pick = dtotal.index(min(dtotal))
		shx = shx_cand[pick]
		shy = shy_cand[pick]
	return shx,shy

def findclosest(target, collection):
    return min((abs(target - i), i) for i in collection)[1]

# make phi follow a line
def fixPhilist(MT,philist,params):
	import copy
	apix = params['apix']
	if params['pf'] == 12:
		RISE = 10.44/apix
		TWIST = -29.86
	elif params['pf'] == 13:
		RISE = 9.5/apix
		TWIST = -27.71
	elif params['pf'] == 14:
		RISE = 8.85/apix
		TWIST = -25.76
	else:
		print "please specify PF"
		sys.exit()
	if os.path.isfile('picklist.txt'):
		print "will use picklist.txt\r",
		phi_MT,sign_MT = GetPhifromfile(MT,philist)
	else:
		phi_MT = mymedian(philist)
		sign_MT = 1
	n = len(philist)
	philist2 = copy.copy(philist)

	#if params['pf14']:
	loc = philist.index(findclosest(phi_MT,philist))
	philist2[loc] = phi_MT
	#print "MT %d: loc = %.2f"%(MT,findclosest(phi_MT,philist))		

	for i in range(loc+1,n):
		dphi = philist2[i-1] - philist2[i]
		philist2[i] += TWIST*round(dphi/TWIST)

	for i in range(loc-1,-1,-1):
		dphi = philist2[i+1] - philist2[i]
		philist2[i] += TWIST*round(dphi/TWIST)
	#else:
	#	for i in range(n):			
	#		dphi = phi_MT - philist2[i]
	#		philist2[i] += TWIST*round(dphi/TWIST)
	if params['forceline']:
		philist2 = removeOutliers(philist2,1.0)

	return philist2,sign_MT

def removeOutliers(mylist,tol):
	#tol = 2.0
	n = len(mylist)
	xi = np.arange(0,n)
	slope, intercept, r_value, p_value, std_err = stats.linregress(xi,mylist)
	line = slope*xi+intercept
	d = abs(mylist-line)
	mdev = np.median(d)
	s = d/mdev if mdev else [0 for i in d]
	
	# generate a clean list
	xx = []
	mylist_tmp = []
	for i in range(n):
		if s[i] > tol:
			continue
		else:
			mylist_tmp.append(mylist[i])
			xx.append(i)
	xx = np.array(xx)
	
	if len(xx) > 2:
		slope, intercept, r_value, p_value, std_err = stats.linregress(xx,mylist_tmp)
        	line = slope*xi+intercept
        d = abs(mylist-line)
	mdev = np.median(d)
	s = d/mdev if mdev else [0 for i in d]
	#print s

	for i in range(n):
		if s[i] > tol:
			mylist[i] = line[i]
	return mylist

def sumAbs(mylist):
	mylist2 = [abs(i) for i in mylist]
	return sum(mylist2)
	
def shift40A_MT(shxlist,shylist,psi,apix):
	MONOMER = 41.0/apix
	shxlist_v1 = [i + MONOMER*math.cos(-psi/180.0*math.pi) for i in shxlist]
        shylist_v1 = [i + MONOMER*math.sin(-psi/180.0*math.pi) for i in shylist]
        shxlist_v2 = [i - MONOMER*math.cos(-psi/180.0*math.pi) for i in shxlist]
        shylist_v2 = [i - MONOMER*math.sin(-psi/180.0*math.pi) for i in shylist]

	if max(sumAbs(shxlist_v1),sumAbs(shylist_v1)) < max(sumAbs(shxlist_v2),sumAbs(shylist_v2)):
	#if max(abs(shx_v1),abs(shy_v1)) < max(abs(shx_v2),abs(shy_v2)):
                shxlist = shxlist_v1
                shylist = shylist_v1
        else:
                shxlist = shxlist_v2
                shylist = shylist_v2
	return shxlist,shylist

def closer2origin(shx,shy,psi,apix):
	DIMER = 82.0/apix
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

def shift40A(shx,shy,psi,apix):
	MONOMER = 41.0/apix
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

def shift80A(shx,shy,psi,DIMER):
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

def convertSameLine(l_MT,rotdict,signdict,params):
	apix = params['apix']
	if params['pf'] == 12:
		RISE = 10.44/apix
		TWIST = -29.86
	elif params['pf'] == 13:
		RISE = 9.5/apix
		TWIST = -27.71
	elif params['pf'] == 14:
		RISE = 8.85/apix
		TWIST = -25.76
	else:
		print "please specify PF"
	MONOMER = 41.0/apix
	
	n = len(l_MT)
	
	#psilist = [float(x.split()[1]) for x in l_MT]
	#thetalist = [float(x.split()[2]) for x in l_MT]
	#philist = [float(x.split()[3]) for x in l_MT]	
	#shxlist = [float(x.split()[4]) for x in l_MT]
	#shylist = [float(x.split()[5]) for x in l_MT]
	MT = int(l_MT[0].split()[7])
	
	try:	
		rot_MT = rotdict[MT]
		sign_MT = signdict[MT]
	except:
		rot_MT = 0
		sign_MT = 1
	
	
	l_MT_new = []
	#FORMAT = "%7d%8.2f%8.2f%8.2f%8.2f%8.2f%7.0f.%6d%9.1f%9.1f%8.2f%7.2f%8.2f\n"
	FORMAT = "%7d%8.2f%8.2f%8.2f%10.2f%10.2f%8d%6d%9.1f%9.1f%8.2f%8.2f%10d%11.4f%8.2f%8.2f\n"
	for i in range(n):
		t1 = l_MT[i].split()
		count = int(t1[0])
                psi = float(t1[1])
                theta = float(t1[2])
                phi = float(t1[3])
                shx = float(t1[4])
                shy = float(t1[5])
		#psi = psilist[i]
		#theta = thetalist[i]
		#phi = philist[i]
		#shx = shxlist[i]
		#shy = shylist[i]
		phi -= TWIST*rot_MT
		shx -= rot_MT*RISE*math.cos(-psi/180.0*math.pi)
		shy -= rot_MT*RISE*math.sin(-psi/180.0*math.pi)
		if phi < 0:
			phi += 360.0
		if phi > 360.0:
			phi -= 360.0
		if sign_MT < 0:
			shx,shy = shift40A(shx,shy,psi,apix)
		shx,shy = closer2origin(shx,shy,psi,apix)
                mag = float(t1[6])
                micro = int(t1[7])
                df1 = float(t1[8])
                df2 = float(t1[9])
                angast = float(t1[10])
                occ = float(t1[11])
                logP = float(t1[12])
                sigma = float(t1[13])
                score = float(t1[14])
                change = float(t1[15])
                l_MT_new.append(FORMAT%(count,psi,theta,phi,shx,shy,mag,micro,df1,df2,angast,occ,logP,sigma,score,change))
                #pres = float(t1[11])
                #dpres = float(t1[12])
		#l_MT_new.append(FORMAT%(count,psi,theta,phi,shx,shy,mag,micro,df1,df2,angast,pres,dpres))
	return l_MT_new


def mainloop(params):
	f1 = file(params["fpar"])
	ll1 = f1.readlines()
	l1 = [x for x in ll1 if x[0]!='C' and x!='\n']
	n = len(l1)
	MTlist = getMTlist(params)
	lastMT = int(l1[-1].split()[7])
	#print lastMT
	MTparlist = [[] for i in range(lastMT+1)]
	for i in l1:
		MT = int(i.split()[7])
		MTparlist[MT].append(i)
	
	fout = file("%s-samePhi2"%params["fpar"],"w")
	SEG = 7

	f9 = file('picklist2.txt')
	l9 = f9.readlines()
	rotdict = {}
	signdict = {}
	for i in l9:
		t9 = i.split()
		micro = int(t9[0])
		rotdict[micro] = float(t9[1])
		try:
			signdict[micro] = math.copysign(1,float(t9[-1]))
		except:
			signdict[micro] = 1
	f9.close()	

	for MT in MTlist:
		print "working on MT %d\t\r"%MT,
		#l_MT = grepMT(params["fpar"],MT)
		l_MT = MTparlist[MT]
		fout.writelines(convertSameLine(l_MT,rotdict,signdict,params))
			
	f1.close()
	fout.close()

if __name__ == "__main__":
	params = setupParserOptions()
	mainloop(params)
