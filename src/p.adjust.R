##R CMD BATCH --no-save --no-restore '--args indir=xxx padjmethod=fdr pthresh=ddd outdir=xxx' p.adjust.R out
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
p0adj=p.adjust(dat$p0,padjmethod)
dat = cbind(dat,p0adj)

head(dat)

i=1
for(ff in files){
	print(ff)
	f=read.table(ff,sep="\t",header=T)
	outf=paste(outdir,str_split_fixed(ff, "/", 3)[,3],sep="/")
	n=length(f$p0)
	if(n>0){
		print(paste(i,n,i+n-1))
		d=dat[i:(i+n-1),]
		write.table(d[which(d$p0adj<=pthresh),],outf,quote=F,sep="\t",row.names = F)
		i=i+n
	}
}
print('DONE')
