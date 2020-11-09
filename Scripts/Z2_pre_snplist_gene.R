setwd("~/cloud/yeastAse_redo/")
load("./Z1_out/cal_evd_p_allinfo.RData")
yps_block <- cal_evd_p_allinfo[,1:4] %>% arrange(ypsChrom,ypsPosit) 

genN <- read.table("~/mnt/nnp/0_3_1_assembly_crossMap/yps128_5_snp_gene_nofilter.txt",header = F)
overlap_gene <- c()
for(i in 2:nrow(genN)) {
    if(genN[i,2] == genN[i-1,2] ){
        print(unname(genN[i-1,]))
        print(unname(genN[i,]))
        overlap_gene <- c(overlap_gene,i-1,i)
    }
}
geneNu <- genN[-overlap_gene,]


yps_rm_gene <- merge.data.frame(yps_block,geneNu,by.x = c(1,2),by.y= c(1,2),sort=F,all=F)
names(yps_rm_gene)[5] <- "geneN"

### 
write.table(unique(yps_rm_gene[,c(1,2,5)]),file="./Z2_out/yps128_5_snpls_evdp_gene",row.names = F,quote=F,col.names = F)

write.table(unique(yps_rm_gene[,c(3,4,5)]),file="./Z2_out/rm11_B_snpls_evdp_gene",row.names = F,quote=F,col.names = F)

write.table(unique(yps_rm_gene[,c(1,2)]),file="./Z2_out/yps128_5_snpls_evdp",row.names = F,quote=F,col.names = F)

write.table(unique(yps_rm_gene[,c(3,4)]),file="./Z2_out/rm11_B_snpls_evdp",row.names = F,quote=F,col.names = F)
