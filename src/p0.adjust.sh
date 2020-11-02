#!/bin/bash

method=fdr
minT=minT2

for indir in detail.peak.si-minI10.v2.minT2 repr.si.detail.peak.si-minI10.v2.minT2 
do
	indir0=$indir.filtered2

	for pthresh in 0.05 0.1
	do
		outdir0=$indir.filtered2.p0.$method.$pthresh
		mkdir $outdir0
		for pval in 0.1,0.01 #0.1
		do
			for n in 40
			do
				indir1=$indir0/$indir.pval$pval.n$n.causal
				outdir=$outdir0/$indir.pval$pval.n$n.causal
				mkdir $outdir

				outf=$outdir
				R CMD BATCH --no-save '--args indir="'$indir1'" padjmethod="'$method'" pthresh='$pthresh' outdir="'$outf'"' p.adjust.R log/$outdir0.$indir.pval$pval.n$n.causal.$method.$pthresh
			done
		done
	done
done

