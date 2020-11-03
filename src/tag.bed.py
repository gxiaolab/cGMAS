#!/usr/bin/python

import sys
import argparse
import glob
from time import strftime
import os
import time
from collections import defaultdict

###########
# get tag snv bed
# list of tag snps (bed format)
#	*info column: refcount:varcount[:othercount]|chrm:exon_st1:exon_end:strd|genotypeSource
#----------
# NEW: add genotype source info for tag snv!
# NEW: 6/8/2018: Union of GMAS, get all het & hom snps for tag
###########

parser = argparse.ArgumentParser(description='Script descriptions here')
parser.add_argument('-i', metavar='inf', required=True, help='Input candid file')
parser.add_argument('-a', metavar='asasf', required=True, help='Input file w/ asas events that are in at least X samples')
parser.add_argument('-r', metavar='ref', required=True, help='alleles count ref file')
parser.add_argument('-c', metavar='cov', required=True, help='tag snv cov')
parser.add_argument('-o', metavar='outf', required=True, help='Output file')

opts = parser.parse_args()
print 'Inf: %s' % opts.i
print 'cov: %s' % opts.c
print 'asasf: %s' % opts.a
print 'ref: %s' % opts.r
print 'Outf: %s' % opts.o

COV = int(opts.c)

def main(argv):
	print 'job starts', strftime('%a, %d %b %Y %I:%M:%S')
	start_time = time.time()

	event = defaultdict(set)
	with open(opts.a) as f:
		for l in f:
			if not l.startswith('case'):
				l = l.split('\t')
				chrm,pos,strd,exon,typ = l[0].split('.')
				event[(chrm,pos)].add('{}:{}:{}'.format(chrm,exon,strd))

	GT = {}
	with open(opts.i) as f:
		for l in f:
			l = l.strip().split('\t')
			if (l[0],l[2]) in event:
				GT[(l[0],l[2])] = l[-1].split('|')

	out = open(opts.o,'w')
	for ff in glob.glob(opts.r+'/chr*snp'):
		with open(ff) as f:
			for l in f:
				chrm,pos,nt,rs,counts = l.strip().split()
				if (chrm,pos) in event and (chrm,pos) in GT:
					gt,source,aa = GT[(chrm,pos)]
					if gt == '0/0' or gt == '1/1' or sum(map(int,counts.split(':'))[:-1]) < COV: continue
					for exon in event[(chrm,pos)]:
						out.write('{}\t{}\t{}\t{}|{}|{}\t0\t{}\n'.format(chrm,int(pos)-1,pos,counts,exon,source,exon[-1]))

	out.close()

	print("--- %s seconds ---" % (time.time() - start_time))
	print 'DONE!', strftime('%a, %d %b %Y %I:%M:%S')

if __name__ == '__main__':
	main(sys.argv[1:])


