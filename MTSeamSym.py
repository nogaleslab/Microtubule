#!/usr/bin/env python

import os,sys
import subprocess
import optparse
import math
import time
from itertools import product

#==========================
def setupParserOptions():
	parser = optparse.OptionParser()
	#parser.set_usage("%prog -s <stack>")
	parser.add_option("-v",dest="volume",type="string",metavar="FILE",
		help="input volume")
	parser.add_option("--savepf",dest="savepf",action="store_true",default=False,
                help="output pf.mrc")
	parser.add_option("--apix",dest="apix",type="float",metavar="FLOAT",
		help="pixel size in angstroms")
	parser.add_option("--rise",dest="rise", type="float",metavar="FLOAT",
		help="helical rise, in angstrom")
	parser.add_option("--twist",dest="twist", type="float",metavar="FLOAT",
		help="helical twist, in degrees")
	parser.add_option("--orad", dest="orad", type="int", metavar="INT",default=215,
		help="outer radius of the 3D mask in x-y plane, in angstrom")
	parser.add_option("--irad", dest="irad", type="int", metavar="INT", default=50,
		help="inner radius of the 3D mask in x-y plane, in angstrom")
	parser.add_option("--zrad", dest="zrad", type="int", metavar="INT",
		help="radius of the 3D mask in z direction, in angstrom")
	parser.add_option("--decor",dest="decor",type="choice", metavar="['kinesin','none']",
		choices=['kinesin','none'], default = 'kinesin',
		help="MT bound with kinesin or nothing")

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

#==========================
def checkConflicts(params):
	if not params['volume']:
		print "Specify a volume"
		sys.exit()
	if not os.path.isfile(params['volume']):
		print "the specified volume '%s' does not exist"%params['volume']
		sys.exit()
	# read volume and get size
	print "reading file: %s"%params['volume']
	params['vol'] = EMData()
	params['vol'].read_image(params['volume'])
	params['nx'] = params['vol'].get_xsize()

	# get helical parameters
	params['pf']=int(round(360.0/abs(params['twist'])))

def circularMask2D(nx):
	falloff = 30.0					# cosine falloff on the edge
	smask2D = EMData(nx,nx)
	smask2D.to_one()
	rad = nx/2-falloff

	for i,j in product(range(nx),range(nx)):
			dx = abs(i-nx/2)
			dy = abs(j-nx/2)

			r2 = dx**2 + dy**2
			if r2 > rad**2:
				wt = 0.5*(1 + math.cos(math.pi*min(1,(math.sqrt(r2)-rad)/falloff)))
			else:
				wt = 1
			smask2D.set(i,j,wt)
	return smask2D

#==========================
def createMask(params):
	nx = params['nx']
	apix = params['apix']

	img = EMData(nx,nx)
	img.to_zero()
	#add 3 degrees to overlap with the neighboring density
	overlap=3*math.pi/180.0
	alpha = math.pi/2 - math.pi/params['pf'] - overlap
	for x,y in ((x,y) for x in range(0,nx) for y in range(nx/2,nx)):
		dx = abs(x-nx/2)
		dy = abs(y-nx/2)
		# if above the line y = tan(alpha)*x
		inner = dx*math.tan(alpha)
		outer = dx*math.tan(alpha-overlap)
		if dy >= inner:
			img.set(x,y,1)
		elif dy >= outer:
			pos = (inner-dy)/(inner-outer)
			img.set(x,y,1-pos)

	smask2D = circularMask2D(nx)
	img.mult(smask2D)
	wedge = EMData(nx,nx,nx)
	twist = params['twist']
	rise = params['rise']
	alpha = 360+(params['pf']*twist)
	for z in range(nx):
		l = params['pf']*rise
		rot = alpha/l*apix
		finalrot = ((z-nx/2)*rot)/3
		t=Transform()
		t.set_rotation({"type":"2d","alpha":-finalrot})
		newslice=img.process("xform",{"transform":t})
		wedge.insert_clip(newslice,(0,0,z))

	if params["decor"] == "kinesin":
		print "decor = kinesin"
		ymsk = int(148/apix)
		xmsk = int(49/apix)
		mskrad = int(20/apix)
		if params['pf'] == 12:
			ymsk = int(142/apix)
			xmsk = int(42/apix)
			mskrad = int(16/apix)
	
		# see if mask is near the edge:
		edge=ymsk*math.atan(math.pi/params['pf'])
		if (abs(xmsk)+mskrad)>=edge:
			# distance for corresponding positive mask
			edge = int(2*edge)
			xmsk2 = int(math.copysign(edge-abs(xmsk),xmsk)*-1)
			# take max of 1 mask
			avgr = Averagers.get("minmax",{"max":1})
			avgr.add_image_list([wedge,kinesinMask(nx,mskrad,xmsk2,ymsk,pos=True)])
			wedge=avgr.finish()
		# multiply 0 mask
		wedge *= kinesinMask(nx,mskrad,xmsk,ymsk)
	
	# odd-numbered protofilaments are off by 1/2 twist
	if params['pf']%2==1:
		t = Transform({"type":"spider","psi":twist/2})
		wedge.process_inplace("xform",{"transform":t})
	
	if params['pf'] == 12:
		print "pf = 12"
		# apply additional rotation so that the wedge covers the groove
		t2 = Transform({"type":"spider","psi":-17.0})
		wedge.process_inplace("xform",{"transform":t2})

	return wedge

#===========================
def regenerateFromPF(params,wedgemask):
	"""
	mask out one away from seam, regenerate microtubule with seam 
	"""
	import shutil,subprocess

	# convert rise to pixels
	nx = params['nx'] 
	rise = params['rise']/params['apix']
	twist = params['twist']

	if params['savepf']:
		pfvol = params['vol']*wedgemask
		pfvol.write_image("pf.mrc")
		sys.exit()
	
	sumvol = EMData(nx,nx,nx)
	sumvol.to_zero()
	pfoffset=int(params['pf']/2)

	start_time = time.time()

	for pnum in range(-pfoffset,params['pf']-pfoffset):
		print "preparing copy %i"%pnum
		ang = twist*pnum
		trans = -(rise*pnum)
		ang*=-1
		trans*=-1
		t = Transform({"type":"spider","psi":ang})
		t.set_trans(0,0,trans)
		volcopy = params['vol'].process("xform",{"transform":t})
		seammaskcopy = wedgemask.process("xform",{"transform":t})
		
		sumvol = sumvol*(1-seammaskcopy)+volcopy*seammaskcopy

	print "Seamed MT regenerated in %.2f minutes"%((time.time()-start_time)/60.0)

	params['vol']=sumvol.process("normalize")
	

#===========================
def kinesinMask(nx,rad,cx,cy,pos=False):
	# soft edge cylinder mask for kinesin position
	img = EMData(nx,nx)
	img.to_one()
	if pos is True:
		img.to_zero()

	# outer radius
	orad = (rad+rad*.5)

	if abs(cy) > (nx/2-orad) : cy = int((cy/abs(cy))*(nx/2-orad))
	if abs(cx) > (nx/2-orad) : cx = int((cx/abs(cx))*(nx/2-orad))
	for x,y in ((x,y) for x in range(-nx/2,nx/2) for y in range(-nx/2,nx/2)):
		r2 = x**2+y**2
		if r2 < (orad*orad):
			if r2 < rad*rad:
				val=1
			else:
				diff=orad**2-rad**2
				val=1-((r2-rad*rad)/(diff))
			if pos is True:
				img.set(nx/2-x+cx,nx/2+y+cy,val)
			else:
				img.set(nx/2+x+cx,nx/2+y+cy,1-val)
		
	cylmask = EMData(nx,nx,nx)
	twist = params['twist']
	rise = params['rise']
	alpha = 360+(params['pf']*twist)
	for z in range(nx):
		l = params['pf']*rise
		rot = alpha/l*params['apix']
		finalrot = (z-nx/2)*rot
		t=Transform()
		t.set_rotation({"type":"2d","alpha":-finalrot/3})
		newslice=img.process("xform",{"transform":t})
		cylmask.insert_clip(newslice,(0,0,z))
	return cylmask

#==========================
def edgeMask(params):
	"""
	create a 3D cylinder mask to remove edges and artifacts from symmetrization
	"""
	nx = params['nx']
	nxm = int(nx+(nx*0.3))
	nym = int(nx-(nx*0.3))

	mask2d = EMData(nxm,nym)
	mask2d.to_one()
	mask2d.process_inplace("mask.decayedge2d",{"width":nx*0.1})
	mask2d.clip_inplace(Region(int(nx*0.3)/2,-int(nx*0.3)/2,nx,nx))
	mask=EMData(nx,nx,nx)
	for i in xrange(nx):
		mask.insert_clip(mask2d,(0,0,i))
	
	t = Transform({"type":"spider","theta":90.0,"phi":90.0})
	mask.process_inplace("xform",{"transform":t})

	irad = int(nx/2*0.85)
	orad = (nx/2)-2
	falloff = orad - irad
	for i in xrange(nx):
		slice2d = mask.get_clip(Region(0,0,i,nx,nx,1))
		slice2d.process_inplace("mask.gaussian",{"inner_radius":irad,"outer_radius":falloff})
		mask.insert_clip(slice2d,[0,0,i])
	params['vol']*=mask

# Generate 2D slices to be inserted into mask3D volume
def createMask2D(params):
	from itertools import product
	apix = params['apix']
	orad = float(params['orad'])/apix
	irad = float(params['irad'])/apix
	nx = params['nx']
	falloff_r = 30			# use steeper falloff
	mask2D = EMData(nx,nx)
	mask2D.to_one()
	for x,y in product(range(nx),range(nx)):
		dx = abs(x-nx/2)
		dy = abs(y-nx/2)
		r2 = dx**2+dy**2
		if r2 > orad*orad:
			wt1 = 0.5*(1 + math.cos(math.pi*min(1,(math.sqrt(r2)-orad)/falloff_r)))
			mask2D.set(x,y,wt1)
		elif r2 < irad*irad:
			wt2 = 0.5*(1 + math.cos(math.pi*min(1,(irad-math.sqrt(r2))/falloff_r)))
			mask2D.set(x,y,wt2)
	return mask2D

def createMask3D(params,mask2D):
	apix = params['apix']
	orad = float(params['orad'])/apix
	irad = float(params['irad'])/apix
	nx = mask2D.get_xsize()
	if params['zrad']:
		zrad = float(params['zrad'])/apix
	else:
		zrad = int(nx/2*0.6)
	mask3D = EMData(nx,nx,nx)
	falloff_z = 30.0
	# now apply soft mask
	for z in range(nx):
		img = EMData(nx,nx)
		img = mask2D.copy()
		dz = abs(z-nx/2)
		if dz > zrad:
			wt3 = 0.5*(1 + math.cos(math.pi*min(1,(dz-zrad)/falloff_z)))
			img.mult(wt3)
		mask3D.insert_clip(img,(0,0,z))
	
	return mask3D

#==========================
def getEMANPath():
	### get the eman2 directory	
	emanpath = subprocess.Popen("env | grep EMAN2DIR", shell=True, stdout=subprocess.PIPE).stdout.read().strip()
	if emanpath:
		emanpath = emanpath.replace("EMAN2DIR=","")
	if os.path.exists(emanpath):
		return emanpath
	print "EMAN2 was not found, make sure it is in your path"
	sys.exit()

#==========================
if __name__ == "__main__":
	getEMANPath()
	from EMAN2 import *
	from sparx import *
	params=setupParserOptions()
	checkConflicts(params)
	
	wedgemask = createMask(params)
	wedgemask.write_image("wedgemask.mrc")
	
	regenerateFromPF(params,wedgemask)
	
	mask2D = createMask2D(params)
	mask3D = createMask3D(params,mask2D)
	#mask3D.write_image('mask3D.mrc')
	params['vol']*= mask3D
	params['vol'].write_image("output.mrc")
