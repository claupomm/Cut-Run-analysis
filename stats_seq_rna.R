###################################################
## sequencing statistics
###################################################



# system('mkdir tables')
# system('mkdir plots')

dirElem  = strsplit(getwd(),'/')[[1]]
expname  = paste('_',dirElem[length(dirElem)],sep='');


## star stats files
files = list.files(path="star",pattern=".Log.final.out$", full.names=T)
print(labs <- gsub("star/", "", gsub(".Log.final.out", "", files, perl=TRUE)))

cols = c("TotalReads", "Multihits", "TooManyLoci", "Splices", "TooManyMM", "TooShort", "UniqueHits")

i = 1
stats = read.csv(files[i], sep="\t",as.is=T)
if (length(labs)>1) {
  for (i in 2:length(labs)) {
	  stats = merge(stats, read.csv(files[i], sep="\t",as.is=T), by="Started.job.on..")
  }
}
colnames(stats) = c("category", labs)
stats = stats[c(grep("Number of input reads", stats$category), grep("Number of reads mapped to multiple loci", stats$category), grep("Number of reads mapped to too many loci", stats$category), grep("Number of splices: Total", stats$category), grep("% of reads unmapped: too many mismatches", stats$category), grep("% of reads unmapped: too short", stats$category), grep("Uniquely mapped reads %", stats$category)),]
rownames(stats) = cols
stats = t(stats[,-1])
stats[,-1] = gsub("%","",stats[,-1])
stats = data.frame(t(matrix(as.numeric(stats), length(labs), length(cols), dimnames=list(labs, cols))))
stats = rbind(stats, short=stats["TooShort",]*stats["TotalReads",]/100)
stats = as.matrix(rbind(stats, unique=stats["UniqueHits",]*stats["TotalReads",]/100))
stats = rbind(stats, sumOther=(stats["TotalReads",]-stats["unique",]))

# unique hits and total reads
tab = stats[match(c("TotalReads", "unique", "sumOther"),dimnames(stats)[[1]]),]/10^6
pdf(paste('plots/ReadStats',expname,'.pdf',sep=''), height=5,width=7)
par(mar=c(4,9,2,1))
bp=barplot(tab[-1,],main='Mapping results', las=1, border=NA, horiz=T, xlab='Reads in Mio',ylab='',xlim=c(0,max(tab[1,])*1.3), col=c("yellowgreen", "gray88"))
legend("topright", dimnames(tab)[[1]][-1], bty='n', col=c("yellowgreen", "gray88"),pch=15)
mtext(side=2, las=1, line=-2, at=bp, text=paste(round(tab[2,]/tab[1,]*100),"%",sep="") )
dev.off()

# splice sites/reads
tab = stats[match(c("TotalReads", "unique", "Splices"),dimnames(stats)[[1]]),]/10^6
pdf(paste('plots/ReadStats_splices',expname,'.pdf',sep=''), height=5,width=7)
par(mar=c(4,9,2,1))
bp=barplot(tab[-1,],main='Mapping results', las=1, border=NA, horiz=T, beside=T, xlab='Reads in Mio',ylab='',xlim=c(0,max(tab[1,])*1.3), col=c("yellowgreen", "SkyBlue2"))
legend("topright", dimnames(tab)[[1]][-1], bty='n', col=c("yellowgreen", "SkyBlue2"),pch=15)
# mtext(side=2, las=1, line=-2, at=bp, text=paste(round(tab[2,]/tab[1,]*100),"%",sep="") )
dev.off()





sessionInfo()

