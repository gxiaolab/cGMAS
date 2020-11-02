##R CMD BATCH --no-save --no-restore '--args indir=xxx padjmethod=fdr pthresh=ddd outdir=xxx' fisherP.adjust.R out
library(stringr)

##First read in the arguments listed at the command line
args=(commandArgs(TRUE))
print(args)

for(i in 1:length(args)){
	eval(parse(text=args[[i]]))
}

files = Sys.glob(paste(indir,"*txt",sep='/'))
print(files)

dataFiles <- lapply(files, read.table, header=TRUE, sep='\t')
dat = do.call(rbind, dataFiles)
fisherPadj=p.adjust(dat$fisherP,padjmethod)
dat = cbind(dat,fisherPadj)

head(dat)

i=1
for(ff in files){
	print(ff)
	f=read.table(ff,sep="\t",header=T)
	outf=paste(outdir,str_split_fixed(ff, "/", 3)[,3],sep="/")
	n=length(f$fisherP)
	if(n>0){
		d=dat[i:c(i+n-1),]
		write.table(d[which(d$fisherPadj>=pthresh),],outf,quote=F,sep="\t",row.names = F)
		i=i+n
	}
}
print('DONE')
