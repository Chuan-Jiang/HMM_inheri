chr <- read.table("chrsize.txt")
seg <- read.table("state_block.txt")
his <- read.table("1m_win.txt")

his_spl <- split(his,his$V1)
seg_spl <- split(seg,seg$V1)

chr_n <- chr$V1
chr<- chr$V2
names(chr) <- chr_n


pdf("visualization.pdf",width=35,height=5)
for ( n in names(chr)){
	len <- chr[n]
	plot(c(0,len), c(-15, 200), type = "n" ,xlab=n,ylab="Counts",axes=F)
	axis(2,seq(0,200,50))
	his_chr <- his_spl[[n]]
	seg_chr <- seg_spl[[n]]
	apply(his_chr,1,function(x){rect(x[2],x[3],x[4],x[5],col=x[6],border=NA)})
	abline(h=-3,lty=2,col="grey60",lwd=3)
	apply(seg_chr,1,function(x){segments(as.numeric(x[3]),y0=-10,x1=as.numeric(x[2]),col=x[4],lend="butt",lwd=15)})
}
dev.off()
