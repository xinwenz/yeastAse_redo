setwd("~/cloud/yeastAse_redo/")
load("./Z2_out/cal.RData")
######### filter snps with no evidence ######## 
evd_rcd_yA <- c()  #check yps_altnative hits match with rm 
evd_rcd_rA <- c()
for(i in 1:nrow(cal)) {
    if(cal[i,"ypsN.x"] == '.') {
        alt_ans <- "IST"
    } else if(cal[i,"rmN.x"] == '.'){
        alt_ans <- "DEL" 
    } else {
        alt_ans <- cal[i,"rmN.x"]
    }  
    
    ck1 <- sum(startsWith(as.character(cal[i,grep("^y.*af",names(cal),value=T)]), alt_ans),na.rm=T) 
    
    if(cal[i,"rmN.y"] == '.') {
        alt_ans <- "IST"
    } else if(cal[i,"ypsN.y"] == '.') {
        alt_ans <- "DEL" 
    } else {
        alt_ans <- cal[i,"ypsN.y"]
    }
    
    ck2 <- sum(startsWith(as.character(cal[i,grep("^r.*af",names(cal),value=T)]),alt_ans),na.rm = T)
    
    evd_rcd_yA <- c(evd_rcd_yA,ck1)
    evd_rcd_rA <- c(evd_rcd_rA,ck2)
} 

no_evd_rows <- which(evd_rcd_yA < 2 | evd_rcd_rA < 2 ) 

cal_evd <- cal[-no_evd_rows,grep("*_rc", names(cal), value=T)]

#save(cal_evd, file = "~/cloud/project/F1_snpfilter/3_5Bsnp/cal_evd.RData")
x <- colSums(cal_evd[,grep("^y.*[HA]_rc",names(cal),value=T)],na.rm=T)
y <- colSums(cal_evd[,grep("^r.*[HA]_rc",names(cal),value=T)],na.rm=T)
xy <- data.frame(x,y)
p1 <- ggplot(xy,aes(x=x,y=y)) + geom_point() + labs(title="DNA mapping on snps with evidence" , x = "total read counts mapped to Yps128 allele" , y = "total read counts mapped to Rm11-1a allele" ) + scale_y_continuous(labels = scales::scientific) + scale_x_continuous(labels = scales::scientific) + geom_abline(slope = 1,intercept = 0) + theme(text=element_text(size=15))

#ggsave(p1, filename = "Z3_out/evidence_snps_counts_sum.pdf")

##  second step is to remove p value too small positions### 
x1 <- rowSums(cal_evd[,c(1:20)[-13]],na.rm = TRUE)
y1 <- rowSums(cal_evd[,c(1:20)[-13]+20],na.rm = TRUE)
pval_df <- data.frame(yps=x1,ypsrm = x1+y1)
pres <- apply(pval_df, 1 , function(x) binom.test(x[1],x[2])$p.value)

#pres_rank <- rank(pres)

cal_evd_p <- cal_evd[pres > 0.05, ] 

##### get postions and group ### 
cal_evd_p_allinfo <- cal[-no_evd_rows,] %>% filter(pres > 0.05)
save(cal_evd_p_allinfo,file= "./Z1_out/cal_evd_p_allinfo.RData")

#### 
############ 
x <- colSums(cal_evd_p_allinfo[,grep("^y.*[HA]_rc",names(cal_evd_p_allinfo),value=T)],na.rm=T)
y <- colSums(cal_evd_p_allinfo[,grep("^r.*[HA]_rc",names(cal_evd_p_allinfo),value=T)],na.rm=T)
xy <- data.frame(x,y)
p2 <- ggplot(xy,aes(x=x,y=y)) + geom_point() + labs(title="DNA mapping on SNPs with evidence and P value > 0.05" , x = "total read counts mapped to Yps128 allele" , y = "total read counts mapped to Rm11-1a allele" ) + scale_y_continuous(labels = scales::scientific) + scale_x_continuous(labels = scales::scientific) + geom_abline(slope = 1,intercept = 0) + theme(text=element_text(size=15))
#ggsave(p2, filename = "Z3_out/evidence_p_snps_counts_sum.pdf")

### ########### 
mydisc <- function(chr1,pos1,chr2,pos2) {
    if(chr1 != chr2) {
        return(Inf)
    }else{
        return(abs(as.numeric(pos1)-as.numeric(pos2)))
    }
}


tf2groupName <- function(x) {
    ans <- vector(length=length(x))
    k <- 1 
    ans[1] <- 1 
    for(i in 2:length(x)) {
        if(x[i] == TRUE) {
            k <- k+1
            ans[i] <- k
        }else{
            ans[i] <- ans[i-1]
        }
    }
    return(ans)
}

######## add group name for DNA #######
yps_block <- cal_evd_p_allinfo[,1:4] %>% arrange(ypsChrom,ypsPosit)
yps_block_tf <- vector(length=nrow(yps_block))
yps_block_tf[1] <- TRUE  
for(i in 2:nrow(yps_block)) {
    yps_block_tf[i] <- 
        mydisc(yps_block[i-1,"ypsChrom"], yps_block[i-1,"ypsPosit"],
               yps_block[i,"ypsChrom"], yps_block[i,"ypsPosit"]) > 700000
} 

yps_rm_group <- cbind(yps_block,yps_block_tf,gN= tf2groupName(yps_block_tf))


write.table(unique(yps_rm_group[,c(1,2,6)]),file="~/cloud/yeastAse_redo/001_snpPos/Z3_out/yps128_5_snpls_evdp_group",row.names = F,quote=F,col.names = F)

write.table(unique(yps_rm_group[,c(3,4,6)]),file="~/cloud/yeastAse_redo/001_snpPos/Z3_out/rm11_B_snpls_evdp_group",row.names = F,quote=F,col.names = F)


   