#!/usr/bin/python

import sys
import argparse
import glob
from time import strftime
import os
import time
from collections import defaultdict

sys.path.append('./lib')

###########
# If 1 candid causal snv has multiple tags, need to make sure all the candid-tag pairs are predicted to be causal!!
#	=> real causal snv doesn't matter which tag it is. all tags show show the same results!
# Filter out also the snv & exon not on the same strd! in this case, the snv usually have both +/-! so it's okay to filter out one.
# output: combine all chrm and generate 1 file for each tissue
#####
#12/21/2017: add anno => putative causal snp dist in neigh exons & neigh introns
#		+ strand:
#			===[* upE *]* === flankIup === *[* targetE *]* === flankIdn === *[* dnE *]===
#			    +     - +                  - +         - +                  - +     -
#
#				1) causal snp in target exon: distType = targE
#						a) exonic causal snp is upstream of targetE end: dist < 0
#						b) exonic causal snp is downstream of targetE start: dist > 0
#				2) causal snp in upstream intron: distType = flankI
#						a) intronic causal snp is upstream of targetE start: dist < 0
#						b) intronic causal snp is downstream of upE end: dist > 0
#				3) causal snp in upstream exon: distType = upE
#						a) exonic causal snp is upstream of upE end: dist < 0
#						b) exonic causal snp is downstream of upE start: dist > 0
#				4) causal snp in downstream intron: distType = flankI
#						a) intronic causal snp is upstream of dnE start: dist < 0
#						b) intronic causal snp is downstream of targetE end: dist > 0
#				5) causal snp in downstream exon: distType = dnE
#						a) exonic causal snp is upstream of dnE end: dist < 0
#						b) exonic causal snp is downstream of dnE start: dist > 0
###########

parser = argparse.ArgumentParser(description='Script descriptions here')
parser.add_argument('-i', metavar='annoI', required=True, help='intron anno bed')
parser.add_argument('-e', metavar='annoE', required=True, help='exon anno bed')
parser.add_argument('-r', metavar='causalf', required=True, help='ref causal si file dir') #peak.si-minI10.minT2/Artery-Aorta/chr16.peak.si.txt
parser.add_argument('-o', metavar='outf', required=True, help='Output file')
parser.add_argument('-t', metavar='tissue', required=True, help='tissue of interest')
parser.add_argument('-s', metavar='si', required=True, help='min Si')
parser.add_argument('-p', metavar='pval', required=True, help='min pval; pval is testing whether si is diff from 1')
parser.add_argument('-n', metavar='minPt', required=True, help='min data points (indiv) per causal-exon-tag pair')
parser.add_argument('-m', metavar='major', required=True, help='min membership ratio of the major component')

opts = parser.parse_args()
print 'intron bed anno: %s' % opts.i
print 'exon bed anno: %s' % opts.e
print 'Outf: %s' % opts.o
print 'causal ref: %s' % opts.r
print 'tissue: %s' % opts.t
print 'min pval: %s' % opts.p
print 'min si: %s' % opts.s
print 'min points: %s' % opts.n
print 'min membership ratio of the major component: %s' % opts.m

SI = float(opts.s)
PV = float(opts.p)
MEM = float(opts.m)
N = int(opts.n)

#use dirs to store all candid-tag pairs
mem = defaultdict(dict) #membership ratio
si = defaultdict(dict) #peak si
nn = defaultdict(dict) #total number of si
pv = defaultdict(dict) #pvals
li = defaultdict(dict) #entire lines
res = defaultdict(dict) #entire lines for final results
intronup = defaultdict(set) #neighboring upstream intron
introndn = defaultdict(set) #neighboring downstream intron
exonup = defaultdict(set) #neighboring upstream exon
exondn = defaultdict(set) #neighboring downstream exon

# When AS region is exon,
#	annoTarget: intron bed
#	annoNeigh: exon bed
# When AS region is intron,
#	annoTarget: exon bed
#	annoNeigh: intron bed
def calcDist(annoTarget,annoNeigh,r1,r2,r3,out):
	#The comments here is assuming the AS region is exon
	#get introns
	setup = lambda: {'iup':'NA', 'idn':'NA', 'eup':'NA', 'edn':'NA'}
	setup2=lambda: defaultdict(setup)
	anno = defaultdict(setup2) #anno[(cand,exon,tag)][trx]:{'iup':coord, 'idn':coord, 'eup':coord, 'edn':coord}
	with open(annoTarget) as f:
		for l in f:
			chrm,st0,end,info,x,strd = l.strip().split('\t')
			#g,trx,x = info.split('|')
			trx = info.split('_')[0]
			if strd == '+':
				if (chrm,end,strd) in intronup:
					for (cand,exon,tag) in intronup[(chrm,end,strd)]:
						anno[(cand,exon,tag)][trx]['iup'] = (st0,end)
						exonup[(chrm,st0,strd)].add((cand,exon,tag))
				if (chrm,st0,strd) in introndn:
					for (cand,exon,tag) in introndn[(chrm,st0,strd)]:
						anno[(cand,exon,tag)][trx]['idn'] = (st0,end)
						exondn[(chrm,end,strd)].add((cand,exon,tag))
			else:
				if (chrm,st0,strd) in intronup:
					for (cand,exon,tag) in intronup[(chrm,st0,strd)]:
						anno[(cand,exon,tag)][trx]['iup'] = (st0,end)
						exonup[(chrm,end,strd)].add((cand,exon,tag))
				if (chrm,end,strd) in introndn:
					for (cand,exon,tag) in introndn[(chrm,end,strd)]:
						anno[(cand,exon,tag)][trx]['idn'] = (st0,end)
						exondn[(chrm,st0,strd)].add((cand,exon,tag))

	#get exons
	with open(annoNeigh) as f:
		for l in f:
			chrm,st0,end,info,x,strd = l.strip().split('\t')
			#g,trx,x = info.split('|')
			trx = info.split('_')[0]
			if strd == '+':
				if (chrm,end,strd) in exonup:
					for (cand,exon,tag) in exonup[(chrm,end,strd)]:
						if trx in anno[(cand,exon,tag)]: anno[(cand,exon,tag)][trx]['eup'] = (st0,end)
				if (chrm,st0,strd) in exondn:
					for (cand,exon,tag) in exondn[(chrm,st0,strd)]:
						if trx in anno[(cand,exon,tag)]: anno[(cand,exon,tag)][trx]['edn'] = (st0,end)
			else:
				if (chrm,st0,strd) in exonup:
					for (cand,exon,tag) in exonup[(chrm,st0,strd)]:
						if trx in anno[(cand,exon,tag)]: anno[(cand,exon,tag)][trx]['eup'] = (st0,end)
				if (chrm,end,strd) in exondn:
					for (cand,exon,tag) in exondn[(chrm,end,strd)]:
						if trx in anno[(cand,exon,tag)]: anno[(cand,exon,tag)][trx]['edn'] = (st0,end)

	#calc dist and write output
	#+ strand:
	#		===[* upE *]* === flankIup === *[* targetE *]* === flankIdn === *[* dnE *]===
	#			+     - +                  - +         - +                  - +     -
	#anno[(cand,exon,tag)][trx]:{'iup':coord, 'idn':coord, 'eup':coord, 'edn':coord}
	for (cand,exon) in res.iterkeys():
		for tag,l in res[(cand,exon)].iteritems():
			chrm,pos,strd = cand.split('.')
			pos = int(pos)
			targst1, targend = map(int,exon.split(':')[1:-1])
			#causal snp in targE
			if pos >= targst1 and pos <= targend:
				dist = min(pos-targst1+1, pos-targend-1, key=abs)
				if strd == '-': dist = -1*dist
				for trx in anno[(cand,exon,tag)].iterkeys():
					if 'NA' not in anno[(cand,exon,tag)][trx].values():
						out.write('{}\t{}\t{}\t{}\n'.format(l,trx,r1,dist))
			#causal snp is upstream of target exon
			elif pos < targst1:
				for trx in anno[(cand,exon,tag)].iterkeys():
					if 'NA' not in anno[(cand,exon,tag)][trx].values():
						if strd == '+':
							ist0,iend = map(int,anno[(cand,exon,tag)][trx]['iup'])
							est0,eend = map(int,anno[(cand,exon,tag)][trx]['eup'])
							side = 'up'
						else:
							ist0,iend = map(int,anno[(cand,exon,tag)][trx]['idn'])
							est0,eend = map(int,anno[(cand,exon,tag)][trx]['edn'])
							side = 'down'
						#causal snp is in upstream intron
						if pos > ist0 and pos <= iend:
							dist = min(pos-ist0, pos-iend-1, key=abs)
							if strd == '-': dist = -1*dist
							out.write('{}\t{}\tflank{}{}\t{}\n'.format(l,trx,r2,side,dist))
						#causal snp is in upstream exon
						elif pos > est0 and pos <= eend:
							dist = min(pos-est0, pos-eend-1, key=abs)
							if strd == '-': dist = -1*dist
							out.write('{}\t{}\t{}{}\t{}\n'.format(l,trx,side,r3,dist))
						#not in range of interest in this trx!
						#else: print 'filter out {}stream:'.format(side),cand,exon,tag,trx
			#causal snp is downstream of target exon
			else:
				for trx in anno[(cand,exon,tag)].iterkeys():
					if 'NA' not in anno[(cand,exon,tag)][trx].values():
						if strd == '+':
							ist0,iend = map(int,anno[(cand,exon,tag)][trx]['idn'])
							est0,eend = map(int,anno[(cand,exon,tag)][trx]['edn'])
							side = 'down'
						else:
							ist0,iend = map(int,anno[(cand,exon,tag)][trx]['iup'])
							est0,eend = map(int,anno[(cand,exon,tag)][trx]['eup'])
							side = 'up'
						#causal snp is in downstream intron
						if pos > ist0 and pos <= iend:
							dist = min(pos-ist0, pos-iend-1, key=abs)
							if strd == '-': dist = -1*dist
							out.write('{}\t{}\tflank{}{}\t{}\n'.format(l,trx,r2,side,dist))
						#causal snp is in downstream exon
						elif pos > est0 and pos <= eend:
							dist = min(pos-est0, pos-eend-1, key=abs)
							if strd == '-': dist = -1*dist
							out.write('{}\t{}\t{}{}\t{}\n'.format(l,trx,side,r3,dist))
						#not in range of interest in this trx!
						#else: print 'filter out {}stream:'.format(side),cand,exon,tag,trx

def main(argv):
	print 'job starts', strftime('%a, %d %b %Y %I:%M:%S')
	start_time = time.time()

	#get causality info
	for ff in glob.glob('{}/{}/*.txt'.format(opts.r,opts.t)): #ff: 1 chrm in 1 tissue at a time
		with open(ff) as f:
			for ll in f:
				if not ll.startswith('causalCandidate'):
					l = ll.split('\t')
					if int(l[5]) == int(l[6]) or int(l[5]) == int(l[8]): continue
					cstrd = l[0][-1]
					estrd = l[3][-1]
					tstrd = l[4][-1]
					if len(set([cstrd,estrd,tstrd])) == 1: #and l[5] not in l[6:9]: #strands agree AND not everyone has the same GT!! ==> 1/9/2017: it is okay to have same gt!!!
						if int(l[-12]) >= N: mem[(l[0],l[3])][l[4]] = float(l[-1].split('|')[0]) #mem ratio of the major component
						try:
							si[(l[0],l[3])][l[4]] = float(l[-4].split('|')[0]) #si mean of the major component
							nn[(l[0],l[3])][l[4]] = int(l[-12])
							pv[(l[0],l[3])][l[4]] = float(l[-6].split('|')[0].split(';')[-1]) #pval of the major component, testing diff from Si = 1
							li[(l[0],l[3])][l[4]] = ll
						except ValueError: continue #this means the pval of the major component is 'NA' => meaning it didn't pass rm.bg in the previous step

		for (cand,exon) in pv.iterkeys():
			LEN = len(pv[(cand,exon)])
			LEN2 = len([x for x in nn[(cand,exon)].values() if x >= N])
			if len([x for x in pv[(cand,exon)].values() if x <= PV]) == 0: #no tag snvs have sig pv (si away from 1)
				if len([x for x in si[(cand,exon)].values() if x >= SI]) == LEN: #all tag snvs pass si thresh
					if len([x for x in mem[(cand,exon)].values() if x >= MEM]) == LEN2: #all tag snvs that have enough indiv (N) pass membership ratio thresh
						for tag in li[(cand,exon)].iterkeys():
							if nn[(cand,exon)][tag] >= N: #only print out the entry that the tag snv has enough individuals.
								#don't need to require it for all tags of a cand-exon pair because not all indiv have the same tag.
								chrm,st1,end,strd = exon.split(':')
								if strd == '+':
									intronup[(chrm,str(int(st1)-1),strd)].add((cand,exon,tag))
									introndn[(chrm,end,strd)].add((cand,exon,tag))
								else:
									introndn[(chrm,str(int(st1)-1),strd)].add((cand,exon,tag))
									intronup[(chrm,end,strd)].add((cand,exon,tag))
								res[(cand,exon)][tag] = li[(cand,exon)][tag].strip()

	out = open(opts.o,'w')
	out.write('causalCandidate\tsource\tnt\texon\ttagSNV\ttotalIndiv\tRR\tRV\tVV\tpeakSi\tzScore\tpValue\tnComp\tpeakSiMean\tpeakSiStdev\tpeakSiN\tpeakSiR\ttrx\tdistType\tdist\n')
	#AS region is exon
	calcDist(opts.i,opts.e,'targE','I','E',out)

	#AS region is intron
	calcDist(opts.e,opts.i,'targI','E','I',out)

	out.close()

	print("--- %s seconds ---" % (time.time() - start_time))
	print 'DONE!', strftime('%a, %d %b %Y %I:%M:%S')

if __name__ == '__main__':
	main(sys.argv[1:])

 
