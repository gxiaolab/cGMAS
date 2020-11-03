#!/usr/bin/python

import sys
import argparse
import glob
from time import strftime
import os
import time
from collections import defaultdict

sys.path.append('./lib')
import GenomeFetch as gf
gf = gf.GenomeFetch('hg19')

###########
# combine inferred gt of snv from multiple tissues per individual
#	- GT: trust as is
#	- snv: if 0/1 in X tissues => 0/1
#		if 0/0 or 1/1 in X tissues and the rest are 1/1 or 0/0 => 0/1
#		other => 0/0 or 1/1 by majority
# output: bed
#	info column: gt|source|nt
#		gt => 0/0; 0/1; 1/1
#		source => GT; snv #(snv means that the gt is inferred from snv counts)
#	strd not needed!
###########

parser = argparse.ArgumentParser(description='Script descriptions here')
parser.add_argument('-i', metavar='ind', required=True, help='Individual of interest')
parser.add_argument('-d', metavar='indir', required=True, help='Input file dir')
parser.add_argument('-s', metavar='suff', required=True, help='file suffix')
parser.add_argument('-m', metavar='min', required=True, help='minT to decide heterozygous')
parser.add_argument('-o', metavar='outf', required=True, help='Output file')

opts = parser.parse_args()
print 'Indiv: %s' % opts.i
print 'Indir: %s' % opts.d
print 'suffix: %s' % opts.s
print 'min: %s' % opts.m
print 'Outf: %s' % opts.o

minT = int(opts.m)

def main(argv):
	print 'job starts', strftime('%a, %d %b %Y %I:%M:%S')
	start_time = time.time()

	zero = lambda: defaultdict(int)
	res = defaultdict(zero)
	info = {}
	for ff in glob.glob('{}/*/{}.*{}'.format(opts.d,opts.i,opts.s)):
		print ff
		with open(ff) as f:
			for l in f:
				chrm,pos,gt,sor,nt,counts = l.strip().split('\t')
				snv = (chrm,pos)
				res[snv][gt] += 1
				info[snv] = '{}|{}'.format(sor,nt)

	xout = open(opts.o+'.excluded','w')
	with open(opts.o,'w') as out:
		for (chrm,pos) in res.iterkeys():
			if 'GT' in info[(chrm,pos)]:
				if len(res[(chrm,pos)]) != 1: print (chrm,pos),res[(chrm,pos)]; sys.exit()
				out.write('{}\t{}\t{}\t{}|{}\n'.format(chrm,int(pos)-1,pos,res[(chrm,pos)].keys()[0],info[(chrm,pos)]))
			else:
				rr, ht, vv = 0, 0, 0
				try: rr = res[(chrm,pos)]['0/0']
				except KeyError: pass
				try: ht += res[(chrm,pos)]['0/1']
				except KeyError: pass
				try: ht += res[(chrm,pos)]['1/0']
				except KeyError: pass
				try: vv = res[(chrm,pos)]['1/1']
				except KeyError: pass
				if ht >= minT: #not sure whether this is needed =>  or (rr >= minT and vv >= minT):
					#print 'ht',(chrm,pos),info[(chrm,pos)],rr,ht,vv
					out.write('{}\t{}\t{}\t0/1|{}\n'.format(chrm,int(pos)-1,pos,info[(chrm,pos)]))
				elif rr > vv and rr >= minT:
					out.write('{}\t{}\t{}\t0/0|{}\n'.format(chrm,int(pos)-1,pos,info[(chrm,pos)]))
				elif rr < vv and vv >= minT:
					out.write('{}\t{}\t{}\t1/1|{}\n'.format(chrm,int(pos)-1,pos,info[(chrm,pos)]))
				else:
					xout.write('{}.{}\t{}\t{}\n'.format(chrm,pos,info[(chrm,pos)],res[(chrm,pos)].items()))

	xout.close()

	print("--- %s seconds ---" % (time.time() - start_time))
	print 'DONE!', strftime('%a, %d %b %Y %I:%M:%S')

if __name__ == '__main__':
	main(sys.argv[1:])


