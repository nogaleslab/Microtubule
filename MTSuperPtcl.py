#!/usr/bin/env python
from EMAN2 import *
import math
import optparse
import os,sys
import glob
import numpy as np

def setupParserOptions():
        parser = optparse.OptionParser()
        parser.add_option("-f",dest="fpar",type="string",metavar="FILE",
                help="13pf_1_r1.par")
	parser.add_option("-s",dest="stack",type="string",metavar="FILE",
		help="start.hed")
	parser.add_option("--outmrc", action="store_true",dest="outmrc",default=False,
		help="output to MRC format")
	parser.add_option("--apix", dest="apix", type="float", metavar="FLOAT",
                help="pixel size in angstroms")

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
	return MTlist

def grepMT(l1,MT):
	#f1 = file(fname)
	#l1 = f1.readlines()
        l2 = [x for x in l1 if x[0]!='C' and x!='\n']
        l3 = [x for x in l2 if int(x.split()[7])==MT]
	#f1.close()
        return l3

def mymedian(mylist):
	import copy
	n = len(mylist)
	if n%2 == 1:
		return np.median(mylist)
	else:
		mylist2 = copy.copy(mylist)
		mylist2.sort()
		return np.median(mylist2[1:])

def MakeMTavg(l_MT,stack,nx,apix):
	start = int(l_MT[0].split()[0])
	end = int(l_MT[-1].split()[0])
	MTstack = EMData.read_images(stack,range(start-1,end))
	
	MTavg = EMData(nx,nx)
	MTavg.to_zero()
	n = len(l_MT)
	
	philist = []
	thetalist = []
	logPlist = []
	sigmalist = []
	scorelist = []
	changelist = []
	
	for i in range(n):
		img = MTstack[i]
		t1 = l_MT[i].split()
		psi = float(t1[1])
		theta = float(t1[2])
		phi = float(t1[3])
		shx = float(t1[4])
		shy = float(t1[5])
		mag = float(t1[6])
		micro = int(t1[7])
		df1 = float(t1[8])
		df2 = float(t1[9])
		angast = float(t1[10])
                #occ = float(t1[11])
                logP = float(t1[12])
                sigma = float(t1[13])
                score = float(t1[14])
                change = float(t1[15])
		#
		philist.append(phi)
		thetalist.append(theta)
		logPlist.append(logP)
		sigmalist.append(sigma)
		scorelist.append(score)
		changelist.append(change)
		# Frealign applies the shifts first and then the rotations.
		t1 = Transform()
		t1.set_trans(-shx/apix,-shy/apix)
		img2 = img.process("xform",{"transform":t1})
		t2 = Transform()
		t2.set_rotation({"type":"2d","alpha":-psi})
		img3 = img2.process("xform",{"transform":t2})
		MTavg.add(img3)
	
	phi_MT = mymedian(philist)
	theta_MT = mymedian(thetalist)
	logP_MT = np.median(logPlist)
	sigma_MT = np.median(sigmalist)
	score_MT = np.median(scorelist)
	change_MT = np.median(changelist)
	#FORMAT = "%7d%8.2f%8.2f%8.2f%10.2f%10.2f%8d%6d%9.1f%9.1f%8.2f%8.2f%10d%11.4f%8.2f%8.2f\n"
	FORMAT2 = "%8.2f%8.2f%8.2f%10.2f%10.2f%8d%6d%9.1f%9.1f%8.2f%8.2f%10d%11.4f%8.2f%8.2f\n"
	l_MTavg = FORMAT2%(0,theta_MT,phi_MT,0,0,mag,micro,df1,df2,angast,100,logP_MT,sigma_MT,score_MT,change_MT)
	return MTavg,l_MTavg

def mainloop(params):
	apix = params['apix']
	f1 = file(params["fpar"])
        ll1 = f1.readlines()
        l1 = [x for x in ll1 if x[0]!='C' and x!='\n']
	
	MTlist = getMTlist(l1)
	#print MTlist
	#nMT = len(MTlist)

	lastMT = int(l1[-1].split()[7])
	#print lastMT
	MTparlist = [[] for i in range(lastMT+1)]
	for i in l1:
		MT = int(i.split()[7])
		MTparlist[MT].append(i)
	
	# first calculate the num of particles
	count = 0
	SEG = 7
	for MT in MTlist:
		#print "pre-working on MT %d\t\r"%MT,
		#l_MT = grepMT(l1,MT)
		l_MT = MTparlist[MT]
		nptcl = len(l_MT)
		if nptcl < 3:
			continue
		
		nSEG = nptcl/SEG
		if nptcl%SEG > SEG-2 or nSEG == 0:
			nSEG += 1
		count += nSEG
	
	stack = params["stack"]
	im = EMData(stack,0)
	nx = im.get_xsize()
	del im
	print "# of MTs: %d"%count
	if params['outmrc']:
		print "Allocating space for MTSuperPtcl...\n"
		MTavgstack = EMData(nx,nx,count)
		MTavgstack.write_image("MTSuperPtcl.mrc")
		print "Done allocating"

	f2 = file("%s_MTSuperPtcl"%params["fpar"],"w")
	#FORMAT = "%7d%8.2f%8.2f%8.2f%8.2f%8.2f%8.0f%6d%9.1f%9.1f%8.2f%7.2f%8.2f\n"
	FORMAT = "%7d%8.2f%8.2f%8.2f%10.2f%10.2f%8d%6d%9.1f%9.1f%8.2f%8.2f%10d%11.4f%8.2f%8.2f\n"
	
	# reset count = 0
	count = 0
	for MT in MTlist:
		print "working on MT %d\t\r"%MT,
		#l_MT = grepMT(l1,MT)
		l_MT = MTparlist[MT]
		nptcl = len(l_MT)
		if nptcl < 3:
			continue
		tmp0 = l_MT[0].split()
		mag = float(tmp0[6])
		#df1 = float(tmp0[8])
		#df2 = float(tmp0[9])
		angast = float(tmp0[10])
		pres = float(tmp0[11])
		
		nSEG = nptcl/SEG
		ttt2 = nptcl%SEG
		if nptcl%SEG > SEG-2 or nSEG == 0:
			nSEG += 1
		for i in range(nSEG):
			start = i*SEG
			#end = (i+1)*SEG
			if i == nSEG-1:
				end = nptcl
			else:
				end = (i+1)*SEG
			l_MT_SEG = l_MT[start:end]
			MTavg,l_MTavg = MakeMTavg(l_MT_SEG,stack,nx,apix)
			if params['outmrc']:		    
				region = Region(0,0,count,nx,nx,1)
				MTavg.write_image("MTSuperPtcl.mrc",0,EMUtil.get_image_ext_type("mrc"), False, region, EMUtil.EMDataType.EM_FLOAT, True)
			else:
				MTavg.write_image("MTSuperPtcl.hed",-1)
			f2.write("%7d%s"%(count+1,l_MTavg))
			count += 1

	f1.close()
	f2.close()


if __name__ == "__main__":
        params = setupParserOptions()
	if params['outmrc']:
		try:
			os.remove("MTSuperPtcl.mrc")
		except:
			pass
	else:
		try:
                        os.remove("MTSuperPtcl.hed")
			os.remove("MTSuperPtcl.img")
                except:
                        pass
	mainloop(params)
