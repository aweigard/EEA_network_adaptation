rm(list=ls())

### read in
abcd.dat<-read.csv("data/abcd_base_nback_postmedians_v3_new.csv")
hcp.dat<-read.csv("data/hcp_nback_postmedians_v3_new.csv")

### bootstrapping function to account for family clustering

boot.CI.cor <- function(x, y, data, clust, reps=1000){
  cs <- unique(data[,clust])
  samps <- rep(NA,reps)
  for(i in 1:reps){
    index <- sample(1:length(cs), length(cs), replace=TRUE)
    c1 <- cs[index]
    c1.t <- table(c1)
    sampdat <- data[data[,clust]%in%names(c1.t),]
    sampdat <- sampdat[rep(seq_len(nrow(sampdat)), c1.t[match(sampdat[,clust], names(c1.t ))]), ]
    samps[i] <- cor(sampdat[,x],sampdat[,y],use="complete")
  }
 samps
}

# robust to whether family or site (with families nested) is used in ABCD

s.abcd<-boot.CI.cor("zero.v.mean","two.v.mean",abcd.dat,"Family_ID")
quantile(s.abcd,probs=c(0.025,0.50,0.975))

s.abcd<-boot.CI.cor("zero.v.mean","two.v.mean",abcd.dat,"Site")
quantile(s.abcd,probs=c(0.025,0.50,0.975))

s.hcp<-boot.CI.cor("zero.v.mean","two.v.mean",hcp.dat,"Family_ID")
quantile(s.hcp,probs=c(0.025,0.50,0.975))

### partial network mean and difference scores:

covs<-c("Age","Age_squared" ,"Sex","Race","meanFD_nback","meanFD_nback_squared")
covs.hcp<-c("Age","Age_squared" ,"Sex","Race","fd_lr","fd_lr_squared","fd_rl","fd_rl_squared")

abcd.dat$FPN_0bk.resid<-resid(lm(as.formula(paste("FPN_0bk ~ ",
                                                                        paste(covs,collapse=' + '))),
                                                       data=abcd.dat))

hcp.dat$FPN_0bk.resid<-resid(lm(as.formula(paste("FPN_0bk ~ ",
                                                                       paste(covs.hcp,collapse=' + '))),
                                                      data=hcp.dat))

abcd.dat$FPN_2bk.resid<-resid(lm(as.formula(paste("FPN_2bk ~ ",
                                                  paste(covs,collapse=' + '))),
                                 data=abcd.dat))

hcp.dat$FPN_2bk.resid<-resid(lm(as.formula(paste("FPN_2bk ~ ",
                                                 paste(covs.hcp,collapse=' + '))),
                                data=hcp.dat))

abcd.dat$DAN_0bk.resid<-resid(lm(as.formula(paste("DAN_0bk ~ ",
                                                  paste(covs,collapse=' + '))),
                                 data=abcd.dat))

hcp.dat$DAN_0bk.resid<-resid(lm(as.formula(paste("DAN_0bk ~ ",
                                                 paste(covs.hcp,collapse=' + '))),
                                data=hcp.dat))

abcd.dat$DAN_2bk.resid<-resid(lm(as.formula(paste("DAN_2bk ~ ",
                                                  paste(covs,collapse=' + '))),
                                 data=abcd.dat))

hcp.dat$DAN_2bk.resid<-resid(lm(as.formula(paste("DAN_2bk ~ ",
                                                 paste(covs.hcp,collapse=' + '))),
                                data=hcp.dat))

abcd.dat$FPN_2_0_contrast.resid<-resid(lm(as.formula(paste("FPN_2_0_contrast ~ ",
                                                  paste(covs,collapse=' + '))),
                                 data=abcd.dat))

abcd.dat$DAN_2_0_contrast.resid<-resid(lm(as.formula(paste("DAN_2_0_contrast ~ ",
                                                           paste(covs,collapse=' + '))),
                                          data=abcd.dat))

hcp.dat$FPN_2_0_contrast.resid<-resid(lm(as.formula(paste("FPN_2_0_contrast ~ ",
                                                          paste(covs.hcp,collapse=' + '))),
                                         data=hcp.dat))

hcp.dat$DAN_2_0_contrast.resid<-resid(lm(as.formula(paste("DAN_2_0_contrast ~ ",
                                                 paste(covs.hcp,collapse=' + '))),
                                data=hcp.dat))


#############################################
# bootstrap correlations for ABCD table ######
#############################################

cor(abcd.dat$overall.v.mean.resid,abcd.dat$FPN_0bk.resid)
quantile(boot.CI.cor("overall.v.mean.resid","FPN_0bk.resid",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

cor(abcd.dat$overall.v.mean.resid,abcd.dat$FPN_2bk.resid)
quantile(boot.CI.cor("overall.v.mean.resid","FPN_2bk.resid",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

cor(abcd.dat$overall.v.mean.resid,abcd.dat$FPN_2_0_contrast.resid)
quantile(boot.CI.cor("overall.v.mean.resid","FPN_2_0_contrast.resid",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

cor(abcd.dat$overall.v.mean.resid,abcd.dat$DAN_0bk.resid)
quantile(boot.CI.cor("overall.v.mean.resid","DAN_0bk.resid",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

cor(abcd.dat$overall.v.mean.resid,abcd.dat$DAN_2bk.resid)
quantile(boot.CI.cor("overall.v.mean.resid","DAN_2bk.resid",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

cor(abcd.dat$overall.v.mean.resid,abcd.dat$DAN_2_0_contrast.resid)
quantile(boot.CI.cor("overall.v.mean.resid","DAN_2_0_contrast.resid",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

#############################################
# bootstrap correlations for HCP table ######
#############################################

cor(hcp.dat$overall.v.mean.resid,hcp.dat$FPN_0bk.resid)
quantile(boot.CI.cor("overall.v.mean.resid","FPN_0bk.resid",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$overall.v.mean.resid,hcp.dat$FPN_2bk.resid)
quantile(boot.CI.cor("overall.v.mean.resid","FPN_2bk.resid",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$overall.v.mean.resid,hcp.dat$FPN_2_0_contrast.resid)
quantile(boot.CI.cor("overall.v.mean.resid","FPN_2_0_contrast.resid",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$overall.v.mean.resid,hcp.dat$DAN_0bk.resid)
quantile(boot.CI.cor("overall.v.mean.resid","DAN_0bk.resid",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$overall.v.mean.resid,hcp.dat$DAN_2bk.resid)
quantile(boot.CI.cor("overall.v.mean.resid","DAN_2bk.resid",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$overall.v.mean.resid,hcp.dat$DAN_2_0_contrast.resid)
quantile(boot.CI.cor("overall.v.mean.resid","DAN_2_0_contrast.resid",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))


#############################################
# Supplemental correlation table ############
#############################################

cor(abcd.dat$FPN_0bk,abcd.dat$overall.v.mean)
quantile(boot.CI.cor("FPN_0bk","overall.v.mean",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

cor(abcd.dat$FPN_2bk,abcd.dat$overall.v.mean)
quantile(boot.CI.cor("FPN_2bk","overall.v.mean",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

cor(abcd.dat$FPN_2_0_contrast,abcd.dat$overall.v.mean)
quantile(boot.CI.cor("FPN_2_0_contrast","overall.v.mean",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

cor(abcd.dat$DAN_0bk,abcd.dat$overall.v.mean)
quantile(boot.CI.cor("DAN_0bk","overall.v.mean",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

cor(abcd.dat$DAN_2bk,abcd.dat$overall.v.mean)
quantile(boot.CI.cor("DAN_2bk","overall.v.mean",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

cor(abcd.dat$DAN_2_0_contrast,abcd.dat$overall.v.mean)
quantile(boot.CI.cor("DAN_2_0_contrast","overall.v.mean",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$FPN_0bk,hcp.dat$overall.v.mean)
quantile(boot.CI.cor("FPN_0bk","overall.v.mean",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$FPN_2bk,hcp.dat$overall.v.mean)
quantile(boot.CI.cor("FPN_2bk","overall.v.mean",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$FPN_2_0_contrast,hcp.dat$overall.v.mean)
quantile(boot.CI.cor("FPN_2_0_contrast","overall.v.mean",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$DAN_0bk,hcp.dat$overall.v.mean)
quantile(boot.CI.cor("DAN_0bk","overall.v.mean",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$DAN_2bk,hcp.dat$overall.v.mean)
quantile(boot.CI.cor("DAN_2bk","overall.v.mean",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$DAN_2_0_contrast,hcp.dat$overall.v.mean)
#[1] 0.4991206
quantile(boot.CI.cor("DAN_2_0_contrast","overall.v.mean",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))


