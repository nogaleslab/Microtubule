#!/usr/bin/env python

import optparse
import itertools
import math
import os,sys
import glob

def setupParserOptions():
	parser = optparse.OptionParser()
	parser.set_usage("%prog -f <EB3clK_5.par>")
	parser.add_option("--f1",dest="fpar1",type="string",metavar="FILE",
		help="13pf_2.par")
	parser.add_option("--f2",dest="fpar2",type="string",metavar="FILE",
		help="13pf_4.par")

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

def lowpres(params):
	FORMAT = "%7d%8.2f%8.2f%8.2f%10.2f%10.2f%8d%6d%9.1f%9.1f%8.2f%8.2f%10d%11.4f%8.2f%8.2f\n"
	count = 0

	f1 = file(params['fpar1'])
	ll1 = f1.readlines()
	f2 = file(params['fpar2'])
	ll2 = f2.readlines()
	l1 = [x for x in ll1 if x[0]!='C' and x!='\n']
	l2 = [x for x in ll2 if x[0]!='C' and x!='\n']
	fout = file('%s-lowpres'%params['fpar1'],"w")
	n1 = len(l1)
	n2 = len(l2)
	if n1 != n2:
		print "fatal error, n1 != n2"
		sys.exit()

	for i in range(n1):
		t1 = l1[i].split()
		t2 = l2[i].split()
		score1 = float(t1[-2])
		score2 = float(t2[-2])
		if (score2-score1) > 1.5:
			fout.write(l2[i])
			count += 1
		else:
			fout.write(l1[i])
	
	print "%s lines changed"%count

	for i in range(n1):
                t1 = l1[i].split()
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
                score1 = float(t1[14])
		t2 = l2[i].split()
		score2 = float(t2[14])
	
	f1.close()
	f2.close()
	fout.close()

if __name__ == "__main__":
	params = setupParserOptions()
	lowpres(params)

