#!/usr/bin/python

import sys
import argparse
import glob
from time import strftime
import os
import time
import numpy as np
from sklearn.mixture import GMM
from collections import defaultdict
import scipy as sp
import math
import smartFunctions

###########
# Need to find the snv with Si number close to 1, which is the predicted causal snv.
# => use GMM to decompose the Si distribution, and check where the max mean value is among different components.
#    (http://www.astroml.org/book_figures/chapter4/fig_GMM_1D.html#)
#	1) Need to check how many component there are based on BIC. (http://scikit-learn.org/stable/auto_examples/mixture/plot_gmm_selection.html)
#	2) Which component is the max peak?
#
# explore resulting si => want to find the following cases:
#	1) SNVs w/ all 0/0, 0/1, and 1/1 GT, so we can plot GT vs. PSI plots.
#	2) SNVs that are found in many individuals, so we have enough data points to plot Si distribution.
#	3) tissue differences (good to have but not necessary)
# Output 1 file per tissue:
#	snv, exon, source|nt, n(indiv), n(0/0), n(0/1), n(1/1), mean(Si), median(Si), min(Si), max(Si)
#
# *** Consider ALL candidate-tag pairs. Real causal SNP should have sig high Si (->1) w/ all tags!
###########
#	0) characteristics of non-causal snps
#		1> Have >=2 peaks - 1 @ X=0 and 1 @ X=1 => test GMM_min_mean vs. X=0 & GMM_max_mean vs. X=1
#		2> The heights of peak_X=0 & peak_X=1 is similar => the height is "p(X)" in GMM term
# 	1) add component means (X value of peaks) in output
#	2) add a column for p-value & z-score from testing agains x=0
###########
# Need additional filter: randomly permute data points to be evenly distributed to x-axis, repeat this 100x times to estimate the variance of bg,
# and then compare it with peak height to see if it's still sig
###########

parser = argparse.ArgumentParser(description='Script descriptions here')
parser.add_argument('-i', metavar='inf', required=True, help='Input file prefix')
parser.add_argument('-r', metavar='ref', help='ref file with SNVs to be filtered out')
parser.add_argument('-m', metavar='min', help='min data points (individuals) in the Si distri')
parser.add_argument('-n', metavar='N', required=True, help='Number of GMM components fitted')
parser.add_argument('-o', metavar='outf', required=True, help='Output file')
parser.add_argument('-b', metavar='bin', required=True, help='bin for randomization')
parser.add_argument('-p', metavar='het', required=True, help='% heterozygous individuals')

opts = parser.parse_args()
print 'Inf: %s' % opts.i
print 'min indiv: %s' % opts.m
print 'N of GMMs: %s' % opts.n
print 'ref to be filtered out: %s' % opts.r
print 'Outf: %s' % opts.o
print 'bin size: %s' % opts.b
print 'percent het: %s' % opts.p

MIN = int(opts.m)
NNN = int(opts.n)
BIN = int(opts.b)
percHET = float(opts.p)

def main(argv):
	print 'job starts', strftime('%a, %d %b %Y %I:%M:%S')
	start_time = time.time()

	res = defaultdict(dict)
	for ff in glob.glob(opts.i +"*"):
		with open(ff) as f:
			for l in f:
				if not l.startswith('causalCandidate'):
					l = l.strip().split('\t')
					gt,source,nt = l[1].split('|')
					if 'info' not in res[(l[0],l[2],l[4])]: res[(l[0],l[2],l[4])] = {'info': '{}\t{}'.format(source,nt),
												'0/0': 0, '0/1': 0, '1/0': 0, '1/1': 0, 'si': []}
					try:
						res[(l[0],l[2],l[4])][gt] += 1
						res[(l[0],l[2],l[4])]['si'].append(float(l[-1]))
					except KeyError:
						print l
						pass

	no = set() #chrm.pos.strd
	if opts.r:
		with open(opts.r) as f:
			for l in f:
				snv = l.strip().split()[0]
				no.add(snv)

	with open(opts.o,'w') as out:
		out.write('causalCandidate\tsource\tnt\texon\ttagSNV\ttotalIndiv\tRR\tRV\tVV\tpeakSi\tz0;1\tp0;1\tnComp\tpeakSiMeans\tpeakSiStdevs\tpeakSiNs\tpeakSiRs\n')
		for (snv, exon, tag) in res.iterkeys():
			if len(res[(snv, exon, tag)]['si']) >= MIN and snv not in no:

				het = res[(snv, exon, tag)]['0/1'] + res[(snv, exon, tag)]['1/0']
				tot = res[(snv, exon, tag)]['0/0'] + het + res[(snv, exon, tag)]['1/1']

				# need to pass the min % heterogzyous threshold
				if 100.*het/tot <= percHET: continue

				# Learn the best-fit GMM models
				#  Here we'll use GMM in the standard way: the fit() method
				#  uses an Expectation-Maximization approach to find the best
				#  mixture of Gaussians for the data

				# fit models with 1 to NNN components
				si = np.array(res[(snv, exon, tag)]['si'])
				X = np.reshape(si,(len(si),1))
				#N = np.arange(1, 11)
				N = np.arange(1, NNN+1)
				models = [None for i in range(len(N))]

				for i in range(len(N)):
					models[i] = GMM(N[i]).fit(X)

				# compute the BIC
				BIC = [m.bic(X) for m in models]

				# best-fit mixture
				M_best = models[np.argmin(BIC)]

				# the largest peak in the distri
				#  score: returns an array of Log probabilities of each data point in X
				#  np.exp(logProb) => gives p(x) -> y-axis of the density or histogram
				#  X: si numbers
				logp = M_best.score(X)
				p = np.exp(logp)
				maxP = np.argmax(p)
				peakX = X[maxP][0]
				nComp = M_best.n_components
				labs =  M_best.predict(X)
				peakComp = labs[maxP]

				#randomly distribute the data points
				#discretize the x-axis to bins: https://stackoverflow.com/questions/31730028/how-can-i-generate-a-random-sample-of-bin-counts-given-a-sequence-of-bin-probabi
				nbin = BIN
				targBINindex = int(peakX*nbin)-1
				bgprob = [0] + [1./nbin for i in range(nbin)]
				cdf=np.cumsum(bgprob)

				thresh=[]
				trials=500
				for i in range(trials):
					thresh.append(int(math.ceil(smartFunctions.distributeBins(nbin,len(si),bgprob,targBINindex)[1])))

				#do for all components
				# To decide whether the peakX is close to Si = 0 or 1
				#  calculate the z-score of 0 or 1 to the peak component's mean
				zs, ps, means, stdevs, pns, ratios = [], [], [], [], [], []
				for i in range(nComp):
					peakM = X[labs == i]
					peakMn = len(peakM)
					mean = np.mean(peakM)
					stdev = np.std(peakM, dtype=np.float64)
					z0 = (0 - mean)/stdev
					p0 = sp.stats.norm.sf(abs(z0))*2 #twosided
					z1 = (1 - mean)/stdev
					p1 = sp.stats.norm.sf(abs(z1))*2 #twosided
					means.append(str(mean))
					stdevs.append(str(stdev))
					pns.append(str(peakMn))
					ratios.append(str(1.*peakMn/tot))

					### rm.bg ###
					thispeak = defaultdict(list) #len(p[labs == i])
					try:
						ind = int(X[labs == i][np.argmin(abs(peakM-mean))][0]*nbin)
					except ValueError:
						zs.append('NA;NA')
						ps.append('NA;NA')
						continue
					for xxx in X[labs == i]:
						thispeak[int(xxx[0]*nbin)].append(xxx[0])

					#real peak needs to be > than mean + 2.58 times the stdev (99% confidence intervals assuming 500 trials give us normal distri)
					testingp = len(thispeak[ind])
					if testingp <= np.mean(thresh) + 2.58*np.std(thresh, dtype=np.float64, axis=0):
						zs.append('NA;NA')
						ps.append('NA;NA')
						continue
					#############

					zs.append('{};{}'.format(z0,z1))
					ps.append('{};{}'.format(p0,p1))

				out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(snv,res[(snv,exon,tag)]['info'],exon,tag,
						tot,res[(snv,exon,tag)]['0/0'],het,res[(snv,exon,tag)]['1/1'],peakX,
						'|'.join(zs),'|'.join(ps),nComp,'|'.join(means),'|'.join(stdevs),'|'.join(pns),'|'.join(ratios)))

	print("--- %s seconds ---" % (time.time() - start_time))
	print 'DONE!', strftime('%a, %d %b %Y %I:%M:%S')

if __name__ == '__main__':
	main(sys.argv[1:])


