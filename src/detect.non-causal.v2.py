#!/usr/bin/python

import sys
import argparse
import glob
from time import strftime
import os
import time
import numpy as np
from scipy import stats
from collections import defaultdict

sys.path.append('./lib')
import GenomeFetch as gf
gf = gf.GenomeFetch('hg19')

###########
# detect non-causal snvs leveraging the different genotypes
# ONLY WORKs AFTER get.causal.py because we binarize the Si in this script - either Si->1 or Si->0
# => The 1st round of filtering guarentees that the candidates have Si->1. This script will look at 
#	the Si of each individual to say whether the Si is ->1 or ->0.
#	Or i.e. Si>0.5 or <0.5
# V2: For nComp of 2
#	=> the cutoff for Si > or < determined by the midpoint of Si peak 1 and peak 2
###########

parser = argparse.ArgumentParser(description='Script descriptions here')
parser.add_argument('-i', metavar='inf', required=True, help='Input file prefix')
parser.add_argument('-r', metavar='ref', help='ref file with SNV tested')
parser.add_argument('-s', metavar='suff', required=True, help='suff')
parser.add_argument('-o', metavar='outf', required=True, help='Output file')

opts = parser.parse_args()
print 'Inf: %s' % opts.i
print 'suff: %s' % opts.s
print 'ref for testing SNV: %s' % opts.r
print 'Outf: %s' % opts.o

def main(argv):
	print 'job starts', strftime('%a, %d %b %Y %I:%M:%S')
	start_time = time.time()

	#store the snv to be tested
	testing = defaultdict(set) #(chrm.pos.strd,exon): set of lines
	thresh = {} #(candid,exon,tag): si thresh
#	means = defaultdict(list) #(chrm.pos.strd,exon): list of peakSiMean
#	stdev = defaultdict(list) #(chrm.pos.strd,exon): list of peakSiStdev
	if opts.r:
		with open(opts.r) as f:
			for ll in f:
				ll = ll.strip()
				l = ll.split()
				testing[(l[0],l[3])].add(ll)
				try:
					thresh[(l[0],l[3],l[4])] = sum(map(float,l[13].split('|')))/2
				except ValueError: pass
#				means[(l[0],l[3])].append(float(l[13]))
#				stdev[(l[0],l[3])].append(float(l[14]))

	#data[case]: [het-good,het-bad,hom-good,hom-bad]; good: si->1, bad: si->0 
	zero = lambda:{'hetgood':0, 'hetbad':0, 'homgood':0, 'hombad':0}
	data = defaultdict(zero)
#	hetgood,hetbad,homgood,hombad = 0, 0, 0, 0
	for ff in glob.glob('{}/*2/{}*'.format(opts.i,opts.s)):
		with open(ff) as f:
			for l in f:
				if l.startswith('causalCandidate'): continue
				l = l.strip().split('\t')
				if (l[0],l[2],l[4]) not in thresh: continue
				gt,source,nt = l[1].split('|')
#				data[(l[0],l[3])][gt].append(float(l[-1]))
				si = float(l[-1])
				if si > thresh[(l[0],l[2],l[4])]:
					if gt in ('0/1','1/0'): #het
						data[(l[0],l[2])]['hetgood'] += 1
					else:
						data[(l[0],l[2])]['homgood'] += 1
				else:
					if gt in ('0/1','1/0'): #het
						data[(l[0],l[2])]['hetbad'] += 1
					else:
						data[(l[0],l[2])]['hombad'] += 1

	#test bias and write output
	with open(opts.o, 'w') as out:
		out.write('causalCandidate\tsource\tnt\texon\ttagSNV\ttotalIndiv\tRR\tRV\tVV\tpeakSi\tzScore\tpValue\tnComp\tpeakSiMean\tpeakSiStdev\tpeakSiN\tpeakSiR\ttrx\tdistType\tdist\tp0\tp1\tfisherP\n')
		for case in data.iterkeys():
			hetgood = data[case]['hetgood']
			homgood = data[case]['homgood']
			hetbad = data[case]['hetbad']
			hombad = data[case]['hombad']
			oddsratio, pvalue = stats.fisher_exact([[hetgood, hetbad], [homgood, hombad]])
			#if pvalue <= 0.05: continue
			#exicute below if no bias based on gt => GOOD!
			for x in testing[case]:
				#out.write(x)
				out.write('{}\t{}\n'.format(x,pvalue))

	#generate 2x2 
	#res[case]: [het-good,het-bad,hom-good,hom-bad]; good: si->1, bad: si->0 
#	for case in data.iterkeys():
#		mu = np.mean(means[case])
#		sigma = np.mean(stdev[case])
#		norm = np.random.normal(mu, sigma, 1000)
#		het = data[case]['0/1'] + data[case]['1/0']
#		hom = data[case]['0/0'] + data[case]['1/1']
#		ht_t, ht_p = stats.ttest_ind(het, norm, equal_var = False)
#		hm_t, hm_p = stats.ttest_ind(hom, norm, equal_var = False)
#		if ht_p < 0.05 and hm_p < 0.05:
#			print 'uh-oh1',case
#			print 'mu:', mu
#			print 'het:',het,ht_p
#			print 'hom:',hom,hm_p
#			sys.exit()
#		elif ht_p < 0.05: #ht is deviate from mu Si => bad group
#			
#		elif hm_p < 0.05: #hm is deviate from mu Si => bad group
#		else: #both are insignificant...
#			print 'uh-oh2',case
#			print 'mu:', mu
#			print 'het:',het,ht_p
#			print 'hom:',hom,hm_p

	print("--- %s seconds ---" % (time.time() - start_time))
	print 'DONE!', strftime('%a, %d %b %Y %I:%M:%S')

if __name__ == '__main__':
	main(sys.argv[1:])


