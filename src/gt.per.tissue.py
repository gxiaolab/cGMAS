#!/usr/bin/python

import sys
import argparse
import glob
from time import strftime
import os
import time
import gzip

###########
# Use GT info to validate method to call homo snp from RNA-Seq
#	- min coverage: 10 or 20
#	- AR: [0, 0.1), (0.9, 1]
#	=> Use GT 0/0 or 1/1  to validate!
###########

parser = argparse.ArgumentParser(description='Script descriptions here')
parser.add_argument('-i', metavar='infdir', required=True, help='Input SNV count files dir')
parser.add_argument('-a', metavar='asasf', required=True, help='Input asas file,asas events that are in at least X samples, sep by comma')
parser.add_argument('-r', metavar='refdir', required=True, help='ref GT files dir')
parser.add_argument('-t', metavar='total', required=True, help='total coverage to tell homo dbSNPs')
parser.add_argument('-m', metavar='mono', required=True, help='mono coverage to tell homo dbSNPs')
parser.add_argument('-u', metavar='upperAR', required=True, help='upper-bound allelic ratio for hetero dbSNPs (<=)')
parser.add_argument('-l', metavar='lowerAR', required=True, help='lower-bound allelic ratio for hetero dbSNPs (>=)')
parser.add_argument('-o', metavar='outf', required=True, help='Output file')

opts = parser.parse_args()
print 'Infs dir: %s' % opts.i
print 'asas f: %s' % opts.a
print 'refs dir: %s' % opts.r
print 'total: %s' % opts.t
print 'mono: %s' % opts.m
print 'highAR: %s' % opts.u
print 'lowAR: %s' % opts.l
print 'Outf: %s' % opts.o

TOTO = int(opts.t)
MONO = int(opts.m)
UPAR = float(opts.u)
LOAR = float(opts.l)

def main(argv):
	print 'job starts', strftime('%a, %d %b %Y %I:%M:%S')
	start_time = time.time()

	asasf,asasf2 = opts.a.split(',')

	asasall = set()
	with open(asasf2) as f:
		for l in f:
			if not l.startswith('case'):
				l = l.split('\t')
				asasall.add(tuple(l[0].split('.')))

	asas = set()
	with open(asasf) as f:
		for l in f:
			chrm,pos,strd,g,exon,typ = l.strip().split('\t')
			if (chrm,pos,strd,exon,typ) in asasall:
				snp = (chrm,pos)
				asas.add(snp)

	#collect all gt snps
	ref = {}
	for ff in glob.glob(opts.r + '/*'):
		print ff
		with gzip.open(ff) as f:
			for l in f:
				l = l.split('\t')
				snp = ('chr'+l[4], l[5])
				nts = '{}>{}'.format(l[-2],l[-1][1])
				ref[snp] = (l[2],nts)

	#check RNA-Seq read counts for homo snps
	out = open(opts.o,'w')
	out.write('chrm\tpos\tgt\tsource\talleles\tcounts\n')
	for ff in glob.glob(opts.i+'*'):
		print ff
		with open(ff) as f:
			for l in f:
				l = l.strip().split()
				chrm = l[0]
				snp = (chrm, l[1])
				counts = l[-1]
				nt = l[2]
				#gt snps
				if snp in ref:
					if snp in asas: out.write('{}\t{}\t{}\tASAS-GT\t{}\t{}\n'.format(chrm,l[1],ref[snp][0],ref[snp][1],counts))
					else: out.write('{}\t{}\t{}\tGT\t{}\t{}\n'.format(chrm,l[1],ref[snp][0],ref[snp][1],counts))
				#snv that are asas
				elif snp in asas:
					out.write('{}\t{}\t0/1\tASAS-snv\t{}\t{}\n'.format(chrm,l[1],nt,counts))
				#snv from rna-seq
				else:
					a,b,c = map(int,counts.split(':'))
					ar = 1.*a/(a+b+c)
					if min(a,b) >= MONO:
						if ar >= LOAR and ar <= UPAR:
							out.write('{}\t{}\t0/1\tsnv\t{}\t{}\n'.format(chrm,l[1],nt,counts))
						elif sum([a,b]) >= TOTO:
							if ar < LOAR:
								out.write('{}\t{}\t1/1\tsnv\t{}\t{}\n'.format(chrm,l[1],nt,counts))
							elif ar > UPAR:
								out.write('{}\t{}\t0/0\tsnv\t{}\t{}\n'.format(chrm,l[1],nt,counts))
							else: print 'what',l; sys.exit()
					elif sum([a,b]) >= TOTO:
						if ar < LOAR:
							out.write('{}\t{}\t1/1\tsnv\t{}\t{}\n'.format(chrm,l[1],nt,counts))
						elif ar > UPAR:
							out.write('{}\t{}\t0/0\tsnv\t{}\t{}\n'.format(chrm,l[1],nt,counts))
#						else: out.write('{}\t{}\tna\tno1\t{}\t{}\n'.format(chrm,l[1],nt,counts))
#					else: out.write('{}\t{}\tna\tno2\t{}\t{}\n'.format(chrm,l[1],nt,counts))

	out.close()

	print("--- %s seconds ---" % (time.time() - start_time))
	print 'DONE!', strftime('%a, %d %b %Y %I:%M:%S')

if __name__ == '__main__':
	main(sys.argv[1:])


