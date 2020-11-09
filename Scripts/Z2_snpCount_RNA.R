load("~/cloud/yeastAse_redo/Z1_out/posmap52.RData")
genN <- read.table("/home/xinwen/mnt/nnp/0_3_1_assembly_crossMap/yps128_5_snp_gene_nofilter.txt",header = F)

########################## get yps genome mapping expression #############
setwd("~/mnt/nnp/2_0_2_express_hisat2/hst_yps128_5/")

fileNum <- sapply(strsplit(x=list.files(path=".", pattern = "*.snpCount"),split = "_"),"[",1)
fileLab <- rep(c("A","H"),times=10)
smpID <- paste0(fileNum,fileLab)
file <- paste0(fileNum,"_", fileLab,".snpCount")

exp_yps <- read.table(file[1],header=T,stringsAsFactors = F)
colnames(exp_yps)[-(1:2)] <- paste0("y",smpID[1],"_",c("rf","rc","af","ac")) # coverage,ref,refcount, alter fer, alter count


for (i in 2:20) {
    print(file[i])
    print(smpID[i])
    tmp <- read.table(file[i],header=T,stringsAsFactors = F)
    colnames(tmp)[-c(1,2)] <- paste0("y",smpID[i],"_",c("rf","rc","af","ac"))
    exp_yps <- merge.data.frame(exp_yps,tmp,by.x=c(1,2),by.y=c(1,2),all=F)  # if all=T, then has 87381 snps here; if all=F:42647
}


colnames(exp_yps)[c(1,2)] <- c("ypsChrom","ypsPosit")
###################### get rm genome mapping expression #############
setwd("~/mnt/nnp/2_0_2_express_hisat2/hst_rm11_B/")
fileNum <- sapply(strsplit(x=list.files(path=".", pattern = "*.snpCount"),split = "_"),"[",1)
fileLab <- rep(c("A","H"),times=10)
smpID <- paste0(fileNum,fileLab)
file <- paste0(fileNum,"_", fileLab,".snpCount")

exp_rm <- read.table(file[1],header=T,stringsAsFactors = F)
colnames(exp_rm)[-(1:2)] <- paste0("r",smpID[1],"_",c("rf","rc","af","ac")) # coverage,ref,refcount, alter fer, alter count

for (i in 2:20) {
    print(file[i])
    print(smpID[i])
    tmp <- read.table(file[i],header=T,stringsAsFactors = F)
    colnames(tmp)[-c(1,2)] <- paste0("r",smpID[i],"_",c("rf","rc","af","ac"))
    exp_rm <- merge.data.frame(exp_rm,tmp,by.x=c(1,2),by.y=c(1,2),all=F)  # if all=T, then has 87381 snps here; if all=F
}

colnames(exp_rm)[c(1,2)] <- c("rmChrom","rmPosit")

#### merge yps_rm data based on posion map ##########
setwd('~/cloud/yeastAse_redo/')


mer1 <- merge.data.frame(posmap52,exp_yps,by = c(1,2),sort = F,all = F)
mer2 <- merge.data.frame(mer1,exp_rm,by.x = c(3,4),by.y = c(1,2), sort = F, all = F )
expf <- mer2[,c(3,4,1,2,5:169)] # change order of columns 
#rm(list=ls()[-c(1,2,3)])
j <- sapply(expf,is.factor)
expf[j] <- lapply(expf[j],as.character)

x <- colSums(expf[,grep("^y.*[HA]_rc",names(expf),value=T)],na.rm=T)
y <- colSums(expf[,grep("^y.*[HA]_ac",names(expf),value=T)],na.rm=T)
xy <- data.frame(x,y)
ggplot(xy,aes(x=x,y=y)) + geom_point() + labs(title=" Total Expresssion Read counts in 20 hybrid samples using YPS as reference and correct SNPS" , x = "total read counts mapped to Yps128 allele" , y = "total read counts mapped to Rm11-1a allele" ) + scale_y_continuous(labels = scales::scientific) + scale_x_continuous(labels = scales::scientific) + geom_abline(slope = 1,intercept = 0) + theme(text=element_text(size=15))

x <- colSums(expf[,grep("^r.*[HA]_ac",names(expf),value=T)],na.rm=T)
y <- colSums(expf[,grep("^r.*[HA]_rc",names(expf),value=T)],na.rm=T)
xy <- data.frame(x,y)
ggplot(xy,aes(x=x,y=y)) + geom_point() + labs(title=" Total Expresssion Read counts in 20 hybrid samples using RM as reference and correct SNPS" , x = "total read counts mapped to Yps128 allele" , y = "total read counts mapped to Rm11-1a allele" ) + scale_y_continuous(labels = scales::scientific) + scale_x_continuous(labels = scales::scientific) + geom_abline(slope = 1,intercept = 0) + theme(text=element_text(size=15))

x <- colSums(expf[,grep("^y.*[HA]_rc",names(expf),value=T)],na.rm=T)
y <- colSums(expf[,grep("^r.*[HA]_rc",names(expf),value=T)],na.rm=T)
xy <- data.frame(x,y)
p3 <- ggplot(xy,aes(x=x,y=y)) + geom_point() + labs(title=" Total Expresssion Read counts in 20 hybrid samples using both reference and correct SNPS" , x = "total read counts mapped to Yps128 allele" , y = "total read counts mapped to Rm11-1a allele" ) + scale_y_continuous(labels = scales::scientific) + scale_x_continuous(labels = scales::scientific) + geom_abline(slope = 1,intercept = 0) + theme(text=element_text(size=15))
ggsave(p3, file="./Z2_out/sum_snp_refcounts.pdf")

yps_block <- expf[,1:4] %>% arrange(ypsChrom,ypsPosit)

overlap_gene <- c()
for(i in 2:nrow(genN)) {
    if(genN[i,2] == genN[i-1,2] ){
        print(unname(genN[i-1,]))
        print(unname(genN[i,]))
        overlap_gene <- c(overlap_gene,i-1,i)
    }
}
geneNu <- genN[-overlap_gene,]
names(geneNu) <- c("ypsChrom","ypsPosit","geneN")
expf_snp_gene <- merge.data.frame(expf,geneNu,by.x = c(1,2),by.y= c(1,2),sort=F,all=F)


exp_snp_gene <- expf_snp_gene[, c("geneN", "ypsChrom", "ypsPosit", "rmChrom", "rmPosit", grep("^y.*[AH]_rc",names(expf_snp_gene),value=T), grep("^r.*[AH]_rc",names(expf_snp_gene),value=T))]


save(exp_snp_gene, file="./Z2_out/exp_snp_gene.RData")

