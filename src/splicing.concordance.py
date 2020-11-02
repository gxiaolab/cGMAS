#!/usr/bin/python

import sys
import argparse
import glob
from time import strftime
import os
import time


#sys.path.append('/home/estherhsiao/lib')
sys.path.append('/u/home/s/s8600192/lib')
import GenomeFetch as gf
gf = gf.GenomeFetch('hg19')
import subprocessFxn

###########
#- inputs:
#	- list of tag snps (bed format)
#		*info column: refcount:varcount[:othercount]|chrm:exon_st1:exon_end:strd|genotypeSource
#	- list of gt snps (bed format)
#		*info column: gt => 0/0; 0/1; 1/1
#- output:
#	- Si for each input snps
#
# To reduce some pairwise comparisons b/t SNPs:
#	- Tag SNP + candidate SNPs w/in 500bp near exon-intron boundaries
#	- if "INF" for -s option, use anno => only search for snp pairs w/in the same gene
#
# * In /u/home/s/s8600192/gxxiao3/GTEx/genotype/54045/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v6.p1.c1.GRU/GenotypeFiles/wes11.pass.snp.vcf.gz, for example,
#	there are complicated GT such like chr1.16388709. In this case, filter out in splicing.concordance.py!
#----------
# NEW: use genotype source info from tag snv.
#	=> if source is from GT, max d is as expected @ 0.5
#	=> if source is from db, max d depends on AR thresholds used to determine trustable heterozygous ASAS tag snvs.
###########

parser = argparse.ArgumentParser(description='Script descriptions here')
parser.add_argument('-i', metavar='inf', required=True, help='Input causal candidate bed file')
parser.add_argument('-d', metavar='indiv', required=True, help='Input individual ID')
parser.add_argument('-t', metavar='tag', required=True, help='tag snp bed file')
parser.add_argument('-m', metavar='maxD', required=True, help='max d for RNA-seq defined tag snvs')
parser.add_argument('-o', metavar='outf', required=True, help='Output file')
parser.add_argument('-a', metavar='anno', required=True, help='gene anno bed file')
parser.add_argument('-s', metavar='search', required=True, help='max dist in nt from candidate causal snp to the AS exon to be tested; input "INF" to test all possible snp pairs w/in the same gene')

opts = parser.parse_args()
print 'Inf candidate causal snp: %s' % opts.i
print 'indiv id: %s' % opts.d
print 'tag snp: %s' % opts.t
print 'max d for RNA-seq defined tag snvs: %s' % opts.m
print 'anno: %s' % opts.a
print 'search: %s' % opts.s
print 'Outf: %s' % opts.o

maxD = float(opts.m)

try:
	maxdist = int(opts.s)
except ValueError:
	maxdist = 'INF'

def snpInGene(inf,anno):
	res = {}
	stdout, stderr = subprocessFxn.run_command(['/u/nobackup/gxxiao2/apps/bin/bedtools2/bin/intersectBed', '-wo', '-a', inf, '-b', anno])
	if stderr: print 1,inf; sys.exit()
	for l in stdout.strip().split('\n'):
		l = l.split('\t')
		g = l[-4].split('|')[0]
		if g not in res: res[g] = {}
		res[g][tuple(l[:-7]+[l[-2]])] = 1 #with strd ==> check for same strd later in get.causal.py!
	return res

def main(argv):
	print 'job starts', strftime('%a, %d %b %Y %I:%M:%S')
	start_time = time.time()

	#always search snp pairs w/in genes
	#candidate causal snp overlap w/ gene anno
	candidate = snpInGene(opts.i,opts.a) #candidate[gene]:snp-bed-file
	#tag snp overlap w/ gene anno
	tag = snpInGene(opts.t,opts.a)

	out = open(opts.o,'w')
	out.write('causalCandidate\tgenotype\texon\tindiv|tagGTsource\ttagSNP\tcoverage\tref\tvar\tallelicR\tdi\tdi2\tSi\n')
	if maxdist == 'INF': 
		for g in set(candidate.keys()).intersection(set(tag.keys())): #intersection means candidate and tag are in the same gene
			for cl in candidate[g].iterkeys():
				csnp = '{}.{}.{}'.format(cl[0],cl[2],cl[-1]) #with strd ==> check for same strd later in get.causal.py!
				gt = cl[3]
				try:
					gtreal = gt.split('|')[0]
				except IndexError: print 'index err',gt; continue
				for tl in tag[g].iterkeys():
					tsnp = '{}.{}.{}'.format(tl[0],tl[2],tl[-1]) #with strd ==> check for same strd later in get.causal.py!
					counts,exon,source = tl[3].split('|')
					chrm,st1,end,strd = exon.split(':')
					if cl[-1] == tl[-1] and cl[-1] == strd:
						counts = map(int,counts.split(':'))
						ref,var = counts[:2]
						tot = sum(counts)
						ri = 1.*ref/tot
						di = abs(0.5-ri)
						if tsnp == csnp or gtreal in ('0/1', '1/0'):
							if source == 'GT': si = di**2/0.5**2
							else: si = di**2/(maxD)**2
						elif gtreal in ('0/0', '1/1'):
							if source == 'GT': si = 1-di**2/0.5**2
							else: si = 1-di**2/(maxD)**2
						#1/2/2018: haven't tried the following 2 lines yet...
						#if tsnp == csnp or gtreal[0] != gtreal[-1]: si = di**2/0.5**2
						#elif gtreal[0] == gtreal[-1]: si = 1-di**2/0.5**2
						else: print 'other gt',csnp,gt; continue
						out.write('{}\t{}\t{}\t{}|{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(csnp,gt,exon,opts.d,source,tsnp,tot,ref,var,ri,di,di**2,si))
	else:
		for g in set(candidate.keys()).intersection(set(tag.keys())): #intersection means candidate and tag are in the same gene
			for cl in candidate[g].iterkeys():
				csnp = '{}.{}.{}'.format(cl[0],cl[2],cl[-1]) #with strd ==> check for same strd later in get.causal.py!
				gt = cl[3]
				try:
					gtreal = gt.split('|')[0]
				except IndexError: print 'index err',gt; continue
				for tl in tag[g].iterkeys():
					tsnp = '{}.{}.{}'.format(tl[0],tl[2],tl[-1]) #with strd ==> check for same strd later in get.causal.py!  
					counts,exon,source = tl[3].split('|')
					chrm,st1,end,strd = exon.split(':')
					if cl[-1] == tl[-1] and cl[-1] == strd:
						st1,end,pos = int(st1), int(end), int(cl[2])
						if abs(st1-pos) <= maxdist or abs(end-pos) <= maxdist or pos >= st1 and pos <= end:
							counts = map(int,counts.split(':'))
							ref,var = counts[:2]
							tot = sum(counts)
							ri = 1.*ref/tot
							di = abs(0.5-ri)
							if tsnp == csnp or gtreal in ('0/1', '1/0'): 
								if source == 'GT': si = di**2/0.5**2
								else: si = di**2/(maxD)**2
							elif gtreal in ('0/0', '1/1'):
								if source == 'GT': si = 1-di**2/0.5**2
								else: si = 1-di**2/(maxD)**2
							#1/2/2018: haven't tried the following 2 lines yet...
							#if tsnp == csnp or gtreal[0] != gtreal[-1]: si = di**2/0.5**2
							#elif gtreal[0] == gtreal[-1]: si = 1-di**2/0.5**2
							else: print 'other gt',csnp,gt; continue
							out.write('{}\t{}\t{}\t{}|{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(csnp,gt,exon,opts.d,source,tsnp,tot,ref,var,ri,di,di**2,si))

	print("--- %s seconds ---" % (time.time() - start_time))
	print 'DONE!', strftime('%a, %d %b %Y %I:%M:%S')

if __name__ == '__main__':
	main(sys.argv[1:])


#			cbed = '\n'.join(['\t'.join(x) for x in candidate[g].keys()])
#			tbed = '\n'.join(['\t'.join(x) for x in tag[g].keys()])
#
#			stdout, stderr = subprocessFxn.run_command(['less', cbed])
#			if stderr: print stderr, 'less cbed'; sys.exit()
#
#			stdout, stderr = subprocessFxn.run_command(['sort', '-k1,1','-k2,2n'], stdout)
#			if stderr: print stderr, 'sort cbed'; sys.exit()
#
#			stdout2, stderr = subprocessFxn.run_command(['less', tbed])
#			if stderr: print stderr, 'less tbed'; sys.exit()
#
#			stdout2, stderr = subprocessFxn.run_command(['sort', '-k1,1','-k2,2n'], stdout2)
#			if stderr: print stderr, 'sort tbed'; sys.exit()
#			with open('{}.ttemp.bed'.format(opts.d),'w') as tout: tout.write(stdout2)
#
#			stdout, stderr = subprocessFxn.run_command(['closestBed','-d','-a','{}.ttemp.bed'.format(opts.d),'-b','-'],stdout2)
#			if stderr: print stderr, 'closetBed'; sys.exit()

#			for l in stdout.split('\n'):
#				l = l.split('\t')
#				if abs(int(l[-1])) <= maxdist:

