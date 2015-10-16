#!/usr/bin/env python
import optparse
import os,sys
import numpy as np
import math
import itertools

def setupParserOptions():
	parser = optparse.OptionParser()
	parser.add_option("-f",dest="fpar",type="string",metavar="FILE",
		help="13pf_1_r1.par")
	parser.add_option("--apix", dest="apix", type="float", metavar="FLOAT",
                help="pixel size in angstroms")
	parser.add_option("--pf", dest="pf", type="int", metavar="int",
                help="protofilament number")
	parser.add_option("--fv",dest="fv",type="choice", metavar="['v8','v9']",
                choices=['v8','v9'],default='v9', help="input frealign format, default v9, output is always v9")

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

def GetPhifromfile(MT,philist):
	f9 = file('philist.txt')
	l9 = f9.readlines()
	phidict = {}
	signdict = {}
	for i in l9:
		t9 = i.split()
		micro = int(t9[0])
		phidict[micro] = float(t9[1])
		try:
			signdict[micro] = math.copysign(1,float(t9[-1]))
		except:
			signdict[micro] = 1
	try:
		phi_MT = phidict[MT]
		sign_MT = signdict[MT]
	except:
		phi_MT = mymedian(philist)
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
	MTlist_dup = [abs(int(x.split()[7])) for x in l2]
	MTlist = list(set(MTlist_dup))
	MTlist.sort()
	f1.close()
	return MTlist


def findClosestShxy(shx,shy,shx_ref,shy_ref,psi):
	SHIFTS = np.array([164.,-164.,123.,-123.,82.,-82.,41.,-41.])
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
def unifyPhi(params,MT,philist):
	import copy
	MONOMER = 41.0
        PF = params['pf']
        RISE = MONOMER*3/PF
        TWIST = -360./PF
	if os.path.isfile('philist.txt'):
		#print "will use philist.txt\r",
		phi_MT,sign_MT = GetPhifromfile(MT,philist)
	else:
		#print "no philist.txt file specified, I will use median PHI"		
		phi_MT = mymedian(philist)
		sign_MT = 1
	n = len(philist)
	philist2 = copy.copy(philist)

	loc = philist.index(findclosest(phi_MT,philist))
	philist2[loc] = phi_MT
	#print "MT %d: loc = %.2f"%(MT,findclosest(phi_MT,philist))		

	for i in range(loc+1,n):
		dphi = philist2[i-1] - philist2[i]
		philist2[i] += TWIST*round(dphi/TWIST)

	for i in range(loc-1,-1,-1):
		dphi = philist2[i+1] - philist2[i]
		philist2[i] += TWIST*round(dphi/TWIST)

	return philist2,sign_MT

def sumAbs(mylist):
	mylist2 = [abs(i) for i in mylist]
	return sum(mylist2)
	
def shift40A_MT(shxlist,shylist,psi):
	MONOMER = 41.0
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

def followline(shxlist,shylist,psilist,loc):
	n = len(shxlist)
	for i in range(loc+1,n):
		shxlist[i],shylist[i] = findClosestShxy(shxlist[i],shylist[i],shxlist[i-1],shylist[i-1],psilist[i])
	for i in range(loc-1,-1,-1):
		shxlist[i],shylist[i] = findClosestShxy(shxlist[i],shylist[i],shxlist[i+1],shylist[i+1],psilist[i])
	return shxlist,shylist

def unifyAll(params,l_MT,MT):
	MONOMER = 41.0
        PF = params['pf']
        RISE = MONOMER*3/PF
        TWIST = -360./PF
	
	n = len(l_MT)
	
	psilist = [float(x.split()[1]) for x in l_MT]
	thetalist = [float(x.split()[2]) for x in l_MT]
	philist = [float(x.split()[3]) for x in l_MT]	
	shxlist = [float(x.split()[4]) for x in l_MT]
	shylist = [float(x.split()[5]) for x in l_MT]
	scorelist = [float(x.split()[14]) for x in l_MT]
	#MT = int(l_MT[0].split()[7])

	philist2,sign_MT = unifyPhi(params,MT,philist)

	for i in range(n):
		dphi = philist2[i] - philist[i]
		if abs(dphi) > 15.0:
			#scale = (phi_new-philist[i])/TWIST
			scale = round(dphi/TWIST)
			shxlist[i] += scale*RISE*math.cos(-psilist[i]/180.0*math.pi)
			shylist[i] += scale*RISE*math.sin(-psilist[i]/180.0*math.pi)
		philist[i] = philist2[i]
		
	#start from the particle with the best score
	if params['fv'] == 'v8':
		loc = scorelist.index(min(scorelist))
	else:
		loc = scorelist.index(max(scorelist))

	shxlist,shylist = followline(shxlist,shylist,psilist,loc)

	psi_MT = np.median(psilist)
	theta_MT = np.median(thetalist)

	if sign_MT < 0:
		shxlist,shylist = shift40A_MT(shxlist,shylist,psi_MT)

	l_MT_new = []
	#FORMAT = "%7d%8.2f%8.2f%8.2f%10.2f%10.2f%8d%6d%9.1f%9.1f%8.2f%8.2f%10d%11.4f%8.2f\n"
	FORMAT = "%7d%8.2f%8.2f%8.2f%10.2f%10.2f%8d%6d%9.1f%9.1f%8.2f%8.2f%10d%11.4f%8.2f%8.2f\n"
	for i in range(n):
		t1 = l_MT[i].split()
		count = int(t1[0])
                #psi = float(t1[1])
                #theta = float(t1[2])
                #phi = float(t1[3])
                #shx = float(t1[4])
                #shy = float(t1[5])
		psi = psilist[i]
		if abs(psi-psi_MT) < 30.0 or abs(psi-psi_MT) > 330:
			psi = psilist[i]
		else:
                        psi = psi_MT
		theta = thetalist[i]
		if abs(theta-theta_MT) > 5.0:
			theta = theta_MT
		phi = philist[i]
		shx = shxlist[i]
		shy = shylist[i]
		shx,shy = closer2origin(shx,shy,psi)
                mag = float(t1[6])
                #micro = int(t1[7])
		micro = MT
                df1 = float(t1[8])
                df2 = float(t1[9])
                angast = float(t1[10])
		occ = float(t1[11])
               	logP = float(t1[12])
                sigma = float(t1[13])
                score = float(t1[14])
                change = float(t1[15])
		l_MT_new.append(FORMAT%(count,psi,theta,phi,shx,shy,mag,micro,df1,df2,angast,occ,logP,sigma,score,change))
	return l_MT_new

def fv8tov9(l_MT,params):
	FORMAT = "%7d%8.2f%8.2f%8.2f%10.2f%10.2f%8d%6d%9.1f%9.1f%8.2f%8.2f%10d%11.4f%8.2f%8.2f\n"
	l_MT_v9 = []
	apix = params['apix']
        for i in l_MT:
                t1 = i.split()
                count = int(t1[0])
                psi = float(t1[1])
                theta = float(t1[2])
                phi = float(t1[3])
                shx = float(t1[4])*apix
                shy = float(t1[5])*apix
		mag = float(t1[6])
                micro = int(t1[7])
                df1 = float(t1[8])
                df2 = float(t1[9])
                angast = float(t1[10])
                occ = 100.0
                logP = -22000.0
                sigma = 0.6
                score = 20.0
                change = 0.0
		l_MT_v9.append(FORMAT%(count,psi,theta,phi,shx,shy,mag,micro,df1,df2,angast,occ,logP,sigma,score,change))
	return l_MT_v9


def mainloop(params):

	f1 = file(params["fpar"])
	ll1 = f1.readlines()
	l1 = [x for x in ll1 if x[0]!='C' and x!='\n']
	n = len(l1)
	MTlist = getMTlist(params)
	lastMT = int(l1[-1].split()[7])
	MTparlist = [[] for i in range(lastMT+1)]
	for i in l1:
		MT = int(i.split()[7])
		MTparlist[MT].append(i)

	if os.path.isfile('philist.txt'):
		print "will use philist.txt"
	else:
		print "no philist.txt file specified, I will use median PHI"
	
	fout = file("%s-uniPhi"%params["fpar"],"w")
	MTnew = 0
	for MT in MTlist:
		print "working on MT %d\t\r"%MT,
		l_MT = MTparlist[MT]
		if params['fv'] == 'v8':
			l_MT = fv8tov9(l_MT,params)
		fout.writelines(unifyAll(params,l_MT,MT))
			
	f1.close()
	fout.close()

if __name__ == "__main__":
	params = setupParserOptions()
	mainloop(params)
