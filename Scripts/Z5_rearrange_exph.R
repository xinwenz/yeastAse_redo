library(combinat)
setwd("~/cloud/yeastAse_redo/")
load("./Z2_out/exph.RData")


tmplt <- (1:20)[-c(7,19)]

# this function, matrix + orignal data. 
trsf_sum <- function(mat, orginEf) {
    yps_mat_col <- mat + 1 
    rm_mat_col <- yps_mat_col + 20 
    df <- data.frame(geneName = as.character(orginEf[,1,]))
   # print(yps_mat_col)
   # print(df)
    
    for(i in 1:ncol(yps_mat_col)) {
        yps_tep <- yps_mat_col[,i]
        yps_vname <- paste0("yps_",i)

        assign(x= yps_vname,value=rowSums(orginEf[,yps_tep,drop=F]))
        df <- cbind(df,get(yps_vname))
        names(df)[ncol(df)] <- yps_vname
    }

    for(i in 1:ncol(rm_mat_col)) {
        rm_tep <- rm_mat_col[,i]
        rm_vname <- paste0("rm_",i)
        assign(x=rm_vname, value= rowSums(orginEf[,rm_tep,drop=F]))
        df <- cbind(df,get(rm_vname))
        names(df)[ncol(df)] <- rm_vname
    }
    return(df)
}


chs_tmplt_181 <- combn(tmplt,m=1)
lv0_depth_100sample <- trsf_sum(chs_tmplt_181[,sample(ncol(chs_tmplt_181),size=100,replace = T),drop=F],exph)

chs_tmplt_182 <- combn(tmplt,m=2)
lv1_depth_100sample <- trsf_sum(chs_tmplt_182[,sample(ncol(chs_tmplt_182),size=100,replace = F)],exph)


chs_tmplt_184 <- combn(tmplt,m=4)
lv2_depth_100sample <- trsf_sum(chs_tmplt_184[,sample(ncol(chs_tmplt_184),size=100,replace = F)],exph)

chs_tmplt_186 <- combn(tmplt,m=6)
lv3_depth_100sample <- trsf_sum(chs_tmplt_186[,sample(ncol(chs_tmplt_186),size=100,replace = F)],exph)


chs_tmplt_188 <- combn(tmplt,m=8)
lv4_depth_100sample <- trsf_sum(chs_tmplt_188[,sample(ncol(chs_tmplt_188),size=100,replace = F)],exph)

save(list = ls(pattern="lv*_depth_100sample"), file="./Z5_out/depth_sample.RData")
######
t0<-colSums(lv0_depth_100sample[,-1])
plot(t0[1:100],t0[101:200])

t1<-colSums(lv1_depth_100sample[,-1])
plot(t1[1:100],t1[101:200])       

t2<-colSums(lv2_depth_100sample[,-1])
plot(t2[1:100],t2[101:200])

t3<-colSums(lv3_depth_100sample[,-1])
plot(t3[1:100],t3[101:200])

t4<-colSums(lv4_depth_100sample[,-1])
plot(t4[1:100],t4[101:200])
