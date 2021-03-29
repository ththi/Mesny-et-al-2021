


###

tab=read.table("count_table_tran_its1.txt",header=T)
tab2=read.table("count_table_tran_its2.txt",header=T)


# remove low abun

high_read=colSums(tab)>=1000
tab=tab[,high_read]

high_read2=colSums(tab2)>=1000
tab2=tab2[,high_read2]


tree_ord=read.table("nam_tree_order.txt",header=F)


tab_ra=sweep(tab,2,colSums(tab),"/")
tab_ra2=sweep(tab2,2,colSums(tab2),"/")

row_ord=intersect(rownames(tab),rownames(tab2))


xx=match(tree_ord[,1],row_ord)


col_ro=grep("Root.*(P|S).ITS",colnames(tab))
col_ro2=grep("Root.*(P|S).ITS",colnames(tab2))

col_soi2=grep("Soil.*(P|S).ITS",colnames(tab2))
col_soi=grep("Soil.*(P|S).ITS",colnames(tab))


tab_c=tab_ra>0.0001

tab_c2=tab_ra2>0.0001


### make log boxplots with ra, from present samples !

tab_ra_rem=tab_ra
tab_ra_rem2=tab_ra2

tab_ra_rem[tab_ra_rem==0]<-NA
tab_ra_rem2[tab_ra_rem2==0]<-NA



##fc root soil

fc_vals=log(rowMeans(tab_ra_rem[row_ord[xx],col_ro],na.rm=T) / rowMeans(tab_ra_rem[row_ord[xx],col_soi],na.rm=T ))

p_list=matrix(1,nrow=length(xx),ncol=3,dimnames=list(rownames(tab_ra_rem[row_ord[xx],]),c("pval","pvalad","cex_val")))

tab_ra_rem_x=tab_ra_rem[row_ord[xx],]
for(i in 1:nrow(tab_ra_rem_x)){

	if ( all(is.na(tab_ra_rem_x[i,col_ro]))==T  ){next}

	test_res=wilcox.test(as.numeric(tab_ra_rem_x[i,col_ro]),as.numeric(tab_ra_rem_x[i,col_soi]))
	p_list[i,1]=test_res$p.value
		
}
p_list[,2]=p.adjust(p_list[,1],method="fdr")
zz=p_list[,2]>0.05
p_list[zz,3]=0.5


##### first panel

dev.new(height=7.058823,width=6.1)
par(mfrow=c(5,1),mar=c(0.01,4,4,4))
barplot(rowSums(tab_c[row_ord[xx],col_ro])/length(col_ro),names=F)
boxplot(t(log2(tab_ra_rem[row_ord[xx],col_ro])),las=2,col="light blue",outline=F,names=F)
plot(fc_vals,pch=19,cex=p_list[,3],xaxt="n")
axis(1, at=1:41, labels=tree_ord[,1],las=2)
segments(seq(1,41,1),fc_vals,seq(1,41,1),0)
segments(0,0,41,0)


## global fungal part (from first panel), need gplots package (for alternative, use standard heatmap)

dev.new()
library("gplots")
tab_glob=read.table("glob_root_count.txt",header=T)
tab_glob_ra=sweep(tab_glob[-42,],2,colSums(tab_glob[42,]),"/")
my_col <- colorRampPalette(c("white", "black"))(n = 100)
heatmap.2(as.matrix(tab_glob_ra[,]),tracecol=F,Rowv=F,col=my_col,margin=c(10,20))



##### second panel

dev.new(height=4.2,width=7.2)
par(mfrow=c(1,2))

plot(rowSums(tab_c[row_ord[xx],col_ro])/length(col_ro),rowSums(tab_c2[row_ord[xx],col_ro2])/length(col_ro2),xlab="ITS1",ylab="ITS2")
plot(rowMeans(tab_ra_rem[row_ord[xx],col_ro],na.rm=T), rowMeans(tab_ra_rem2[row_ord[xx],col_ro2],na.rm=T),log="xy",xlab="ITS1",ylab="ITS2")

### for third panel, use read counts from "tax_count_both.txt"

