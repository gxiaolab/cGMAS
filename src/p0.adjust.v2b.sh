#!/bin/bash

method=fdr
minT=minT2

for name in detail.peak.si-minI10 #repr.si.detail.peak.si-minI10 #detail.peak.si-minI10
do
	indir=$name.v2.minT2
	indir0=$name.v2b.minT2.filtered2

	for pthresh in 0.05 0.1
	do
		outdir0=$name.v2b.minT2.filtered2.p0.$method.$pthresh
		mkdir $outdir0
		for pval in 0.05,0.05
		do
			for n in 40 
			do
				for major in 0.9
				do
					for gtr in 0.95
					do
						indir1=$indir0/$indir.pval$pval.n$n.major$major.gtr$gtr.causal
						outdir=$outdir0/$indir.pval$pval.n$n.major$major.gtr$gtr.causal
						mkdir $outdir

						outf=$outdir
						R CMD BATCH --no-save '--args indir="'$indir1'" padjmethod="'$method'" pthresh='$pthresh' outdir="'$outf'"' p.adjust.R log/$outdir0.$indir.pval$pval.n$n.major$major.gtr$gtr.causal.$method.$pthresh
					done
				done
			done
		done
	done
done

