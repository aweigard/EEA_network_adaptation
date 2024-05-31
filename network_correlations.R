rm(list=ls())

### read in
abcd.dat<-read.csv("data/abcd_base_nback_postmedians_v5_asw.csv")
hcp.dat<-read.csv("data/hcp_nback_postmedians_v5_asw.csv")

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
    if(!inherits(x,"formula")){
      samps[i] <- cor(sampdat[,x],sampdat[,y],use="complete")
      }
    else { 
      pred<-predict(lm(x,sampdat))
      samps[i] <- cor(pred,sampdat[,y],use="complete")
      }
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

# contrasts

abcd.dat$FPN_2_0_contrast.resid<-resid(lm(as.formula(paste("FPN_2_0_contrast ~ ",
                                                           paste(covs,collapse=' + '))),
                                          data=abcd.dat))

abcd.dat$DAN_2_0_contrast.resid<-resid(lm(as.formula(paste("DAN_2_0_contrast ~ ",
                                                           paste(covs,collapse=' + '))),
                                          data=abcd.dat))

abcd.dat$VIS_2_0_contrast.resid<-resid(lm(as.formula(paste("VIS_2_0_contrast ~ ",
                                                           paste(covs,collapse=' + '))),
                                          data=abcd.dat))

abcd.dat$SMN_2_0_contrast.resid<-resid(lm(as.formula(paste("SMN_2_0_contrast ~ ",
                                                           paste(covs,collapse=' + '))),
                                          data=abcd.dat))

abcd.dat$VAN_2_0_contrast.resid<-resid(lm(as.formula(paste("VAN_2_0_contrast ~ ",
                                                           paste(covs,collapse=' + '))),
                                          data=abcd.dat))

abcd.dat$LIM_2_0_contrast.resid<-resid(lm(as.formula(paste("LIM_2_0_contrast ~ ",
                                                           paste(covs,collapse=' + '))),
                                          data=abcd.dat))

abcd.dat$DMN_2_0_contrast.resid<-resid(lm(as.formula(paste("DMN_2_0_contrast ~ ",
                                                           paste(covs,collapse=' + '))),
                                          data=abcd.dat))




hcp.dat$FPN_2_0_contrast.resid<-resid(lm(as.formula(paste("FPN_2_0_contrast ~ ",
                                                          paste(covs.hcp,collapse=' + '))),
                                         data=hcp.dat))

hcp.dat$DAN_2_0_contrast.resid<-resid(lm(as.formula(paste("DAN_2_0_contrast ~ ",
                                                          paste(covs.hcp,collapse=' + '))),
                                         data=hcp.dat))


hcp.dat$VIS_2_0_contrast.resid<-resid(lm(as.formula(paste("VIS_2_0_contrast ~ ",
                                                          paste(covs.hcp,collapse=' + '))),
                                         data=hcp.dat))

hcp.dat$SMN_2_0_contrast.resid<-resid(lm(as.formula(paste("SMN_2_0_contrast ~ ",
                                                          paste(covs.hcp,collapse=' + '))),
                                         data=hcp.dat))

hcp.dat$VAN_2_0_contrast.resid<-resid(lm(as.formula(paste("VAN_2_0_contrast ~ ",
                                                          paste(covs.hcp,collapse=' + '))),
                                         data=hcp.dat))

hcp.dat$LIM_2_0_contrast.resid<-resid(lm(as.formula(paste("LIM_2_0_contrast ~ ",
                                                          paste(covs.hcp,collapse=' + '))),
                                         data=hcp.dat))

hcp.dat$DMN_2_0_contrast.resid<-resid(lm(as.formula(paste("DMN_2_0_contrast ~ ",
                                                          paste(covs.hcp,collapse=' + '))),
                                         data=hcp.dat))




# single conditions for relevant networks 

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


abcd.dat$SMN_0bk.resid<-resid(lm(as.formula(paste("SMN_0bk ~ ",
                                                  paste(covs,collapse=' + '))),
                                 data=abcd.dat))

hcp.dat$SMN_0bk.resid<-resid(lm(as.formula(paste("SMN_0bk ~ ",
                                                 paste(covs.hcp,collapse=' + '))),
                                data=hcp.dat))

abcd.dat$SMN_2bk.resid<-resid(lm(as.formula(paste("SMN_2bk ~ ",
                                                  paste(covs,collapse=' + '))),
                                 data=abcd.dat))

hcp.dat$SMN_2bk.resid<-resid(lm(as.formula(paste("SMN_2bk ~ ",
                                                 paste(covs.hcp,collapse=' + '))),
                                data=hcp.dat))

#############################################
# bootstrap correlations for ABCD table ######
#############################################


cor(abcd.dat$overall.v.mean.resid,abcd.dat$FPN_2_0_contrast.resid)
quantile(boot.CI.cor("overall.v.mean.resid","FPN_2_0_contrast.resid",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

cor(abcd.dat$overall.v.mean.resid,abcd.dat$DAN_2_0_contrast.resid)
quantile(boot.CI.cor("overall.v.mean.resid","DAN_2_0_contrast.resid",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

cor(abcd.dat$overall.v.mean.resid,abcd.dat$VIS_2_0_contrast.resid)
quantile(boot.CI.cor("overall.v.mean.resid","VIS_2_0_contrast.resid",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

cor(abcd.dat$overall.v.mean.resid,abcd.dat$SMN_2_0_contrast.resid)
quantile(boot.CI.cor("overall.v.mean.resid","SMN_2_0_contrast.resid",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

cor(abcd.dat$overall.v.mean.resid,abcd.dat$VAN_2_0_contrast.resid)
quantile(boot.CI.cor("overall.v.mean.resid","VAN_2_0_contrast.resid",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

cor(abcd.dat$overall.v.mean.resid,abcd.dat$LIM_2_0_contrast.resid)
quantile(boot.CI.cor("overall.v.mean.resid","LIM_2_0_contrast.resid",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

cor(abcd.dat$overall.v.mean.resid,abcd.dat$DMN_2_0_contrast.resid)
quantile(boot.CI.cor("overall.v.mean.resid","DMN_2_0_contrast.resid",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))


#############################################
# bootstrap correlations for HCP table ######
#############################################


cor(hcp.dat$overall.v.mean.resid,hcp.dat$FPN_2_0_contrast.resid)
quantile(boot.CI.cor("overall.v.mean.resid","FPN_2_0_contrast.resid",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$overall.v.mean.resid,hcp.dat$DAN_2_0_contrast.resid)
quantile(boot.CI.cor("overall.v.mean.resid","DAN_2_0_contrast.resid",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$overall.v.mean.resid,hcp.dat$VIS_2_0_contrast.resid)
quantile(boot.CI.cor("overall.v.mean.resid","VIS_2_0_contrast.resid",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$overall.v.mean.resid,hcp.dat$SMN_2_0_contrast.resid)
quantile(boot.CI.cor("overall.v.mean.resid","SMN_2_0_contrast.resid",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$overall.v.mean.resid,hcp.dat$VAN_2_0_contrast.resid)
quantile(boot.CI.cor("overall.v.mean.resid","VAN_2_0_contrast.resid",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$overall.v.mean.resid,hcp.dat$LIM_2_0_contrast.resid)
quantile(boot.CI.cor("overall.v.mean.resid","LIM_2_0_contrast.resid",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$overall.v.mean.resid,hcp.dat$DMN_2_0_contrast.resid)
quantile(boot.CI.cor("overall.v.mean.resid","DMN_2_0_contrast.resid",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))


########################################
### contributions of FPN, DAN and SMN ##
########################################

three_network<-paste0(c("FPN","DAN","SMN"),"_2_0_contrast.resid")
all_network<-paste0(c("FPN","DAN","VIS","SMN","VAN","LIM","DMN"),"_2_0_contrast.resid")

three_formula<-as.formula(paste("overall.v.mean.resid ~ ",
                                paste(three_network,collapse=' + ')))
all_formula<-as.formula(paste("overall.v.mean.resid ~ ",
                                paste(all_network,collapse=' + ')))

three_lm_abcd<-lm(three_formula,data=abcd.dat)

all_lm_abcd<-lm(all_formula,data=abcd.dat)


three_lm_hcp<-lm(three_formula,data=hcp.dat)

all_lm_hcp<-lm(all_lm_abcd,data=hcp.dat)


abcd.dat$v.pred.three<-predict(three_lm_abcd)
cor(abcd.dat$v.pred.three,abcd.dat$overall.v.mean.resid)
quantile(boot.CI.cor(three_formula,"overall.v.mean.resid",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

abcd.dat$v.pred.all<-predict(all_lm_abcd)
cor(abcd.dat$v.pred.all,abcd.dat$overall.v.mean.resid)
quantile(boot.CI.cor(all_formula,"overall.v.mean.resid",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

hcp.dat$v.pred.three<-predict(three_lm_hcp)
cor(hcp.dat$v.pred.three,hcp.dat$overall.v.mean.resid)
quantile(boot.CI.cor(three_formula,"overall.v.mean.resid",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

hcp.dat$v.pred.all<-predict(all_lm_hcp)
cor(hcp.dat$v.pred.all,hcp.dat$overall.v.mean.resid)
quantile(boot.CI.cor(all_formula,"overall.v.mean.resid",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))


#################################
## single-condition estimates ###
#################################

# ABCD

cor(abcd.dat$overall.v.mean.resid,abcd.dat$FPN_0bk.resid)
quantile(boot.CI.cor("overall.v.mean.resid","FPN_0bk.resid",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

cor(abcd.dat$overall.v.mean.resid,abcd.dat$FPN_2bk.resid)
quantile(boot.CI.cor("overall.v.mean.resid","FPN_2bk.resid",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

cor(abcd.dat$overall.v.mean.resid,abcd.dat$DAN_0bk.resid)
quantile(boot.CI.cor("overall.v.mean.resid","DAN_0bk.resid",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

cor(abcd.dat$overall.v.mean.resid,abcd.dat$DAN_2bk.resid)
quantile(boot.CI.cor("overall.v.mean.resid","DAN_2bk.resid",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

cor(abcd.dat$overall.v.mean.resid,abcd.dat$SMN_0bk.resid)
quantile(boot.CI.cor("overall.v.mean.resid","SMN_0bk.resid",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

cor(abcd.dat$overall.v.mean.resid,abcd.dat$SMN_2bk.resid)
quantile(boot.CI.cor("overall.v.mean.resid","SMN_2bk.resid",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

# HCP 

cor(hcp.dat$overall.v.mean.resid,hcp.dat$FPN_0bk.resid)
quantile(boot.CI.cor("overall.v.mean.resid","FPN_0bk.resid",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$overall.v.mean.resid,hcp.dat$FPN_2bk.resid)
quantile(boot.CI.cor("overall.v.mean.resid","FPN_2bk.resid",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$overall.v.mean.resid,hcp.dat$DAN_0bk.resid)
quantile(boot.CI.cor("overall.v.mean.resid","DAN_0bk.resid",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$overall.v.mean.resid,hcp.dat$DAN_2bk.resid)
quantile(boot.CI.cor("overall.v.mean.resid","DAN_2bk.resid",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$overall.v.mean.resid,hcp.dat$SMN_0bk.resid)
quantile(boot.CI.cor("overall.v.mean.resid","SMN_0bk.resid",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$overall.v.mean.resid,hcp.dat$SMN_2bk.resid)
quantile(boot.CI.cor("overall.v.mean.resid","SMN_2bk.resid",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))


######################################
## 0-back correlations with 2-back ###
######################################

cor(abcd.dat$FPN_0bk.resid,abcd.dat$FPN_2bk.resid)
quantile(boot.CI.cor("FPN_0bk.resid","FPN_2bk.resid",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

cor(abcd.dat$DAN_0bk.resid,abcd.dat$DAN_2bk.resid)
quantile(boot.CI.cor("DAN_0bk.resid","DAN_2bk.resid",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

cor(abcd.dat$SMN_0bk.resid,abcd.dat$SMN_2bk.resid)
quantile(boot.CI.cor("SMN_0bk.resid","SMN_2bk.resid",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))



cor(hcp.dat$FPN_0bk.resid,hcp.dat$FPN_2bk.resid)
quantile(boot.CI.cor("FPN_0bk.resid","FPN_2bk.resid",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$DAN_0bk.resid,hcp.dat$DAN_2bk.resid)
quantile(boot.CI.cor("DAN_0bk.resid","DAN_2bk.resid",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$SMN_0bk.resid,hcp.dat$SMN_2bk.resid)
quantile(boot.CI.cor("SMN_0bk.resid","SMN_2bk.resid",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))


#############################################
# Supplemental correlation table ############
#############################################


cor(abcd.dat$FPN_2_0_contrast,abcd.dat$overall.v.mean)
quantile(boot.CI.cor("FPN_2_0_contrast","overall.v.mean",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))


cor(abcd.dat$DAN_2_0_contrast,abcd.dat$overall.v.mean)
quantile(boot.CI.cor("DAN_2_0_contrast","overall.v.mean",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

cor(abcd.dat$VIS_2_0_contrast,abcd.dat$overall.v.mean)
quantile(boot.CI.cor("VIS_2_0_contrast","overall.v.mean",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

cor(abcd.dat$SMN_2_0_contrast,abcd.dat$overall.v.mean)
quantile(boot.CI.cor("SMN_2_0_contrast","overall.v.mean",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

cor(abcd.dat$VAN_2_0_contrast,abcd.dat$overall.v.mean)
quantile(boot.CI.cor("VAN_2_0_contrast","overall.v.mean",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

cor(abcd.dat$LIM_2_0_contrast,abcd.dat$overall.v.mean)
quantile(boot.CI.cor("LIM_2_0_contrast","overall.v.mean",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

cor(abcd.dat$DMN_2_0_contrast,abcd.dat$overall.v.mean)
quantile(boot.CI.cor("DMN_2_0_contrast","overall.v.mean",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))



cor(hcp.dat$FPN_2_0_contrast,hcp.dat$overall.v.mean)
quantile(boot.CI.cor("FPN_2_0_contrast","overall.v.mean",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$DAN_2_0_contrast,hcp.dat$overall.v.mean)
quantile(boot.CI.cor("DAN_2_0_contrast","overall.v.mean",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$VIS_2_0_contrast,hcp.dat$overall.v.mean)
quantile(boot.CI.cor("VIS_2_0_contrast","overall.v.mean",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$SMN_2_0_contrast,hcp.dat$overall.v.mean)
quantile(boot.CI.cor("SMN_2_0_contrast","overall.v.mean",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$VAN_2_0_contrast,hcp.dat$overall.v.mean)
quantile(boot.CI.cor("VAN_2_0_contrast","overall.v.mean",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$LIM_2_0_contrast,hcp.dat$overall.v.mean)
quantile(boot.CI.cor("LIM_2_0_contrast","overall.v.mean",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$DMN_2_0_contrast,hcp.dat$overall.v.mean)
quantile(boot.CI.cor("DMN_2_0_contrast","overall.v.mean",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))


##############################################
## Supplemental single-condition estimates ###
##############################################

cor(abcd.dat$FPN_0bk,abcd.dat$overall.v.mean)
quantile(boot.CI.cor("FPN_0bk","overall.v.mean",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

cor(abcd.dat$FPN_2bk,abcd.dat$overall.v.mean)
quantile(boot.CI.cor("FPN_2bk","overall.v.mean",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))


cor(abcd.dat$DAN_0bk,abcd.dat$overall.v.mean)
quantile(boot.CI.cor("DAN_0bk","overall.v.mean",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

cor(abcd.dat$DAN_2bk,abcd.dat$overall.v.mean)
quantile(boot.CI.cor("DAN_2bk","overall.v.mean",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

cor(abcd.dat$SMN_0bk,abcd.dat$overall.v.mean)
quantile(boot.CI.cor("SMN_0bk","overall.v.mean",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))

cor(abcd.dat$SMN_2bk,abcd.dat$overall.v.mean)
quantile(boot.CI.cor("SMN_2bk","overall.v.mean",abcd.dat,"Site"),
         probs=c(0.025,0.50,0.975))



cor(hcp.dat$FPN_0bk,hcp.dat$overall.v.mean)
quantile(boot.CI.cor("FPN_0bk","overall.v.mean",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$FPN_2bk,hcp.dat$overall.v.mean)
quantile(boot.CI.cor("FPN_2bk","overall.v.mean",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$DAN_0bk,hcp.dat$overall.v.mean)
quantile(boot.CI.cor("DAN_0bk","overall.v.mean",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$DAN_2bk,hcp.dat$overall.v.mean)
quantile(boot.CI.cor("DAN_2bk","overall.v.mean",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$SMN_0bk,hcp.dat$overall.v.mean)
quantile(boot.CI.cor("SMN_0bk","overall.v.mean",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))

cor(hcp.dat$SMN_2bk,hcp.dat$overall.v.mean)
quantile(boot.CI.cor("SMN_2bk","overall.v.mean",hcp.dat,"Family_ID"),
         probs=c(0.025,0.50,0.975))


# r-value difference CIs for overlapping CIs

quantile((boot.CI.cor("overall.v.mean.resid","DAN_2_0_contrast.resid",abcd.dat,"Site")-
            boot.CI.cor("overall.v.mean.resid","DAN_2bk.resid",abcd.dat,"Site")),
         probs=c(0.025,0.50,0.975))


quantile((boot.CI.cor("overall.v.mean.resid","SMN_2_0_contrast.resid",abcd.dat,"Site")-
            boot.CI.cor("overall.v.mean.resid","SMN_2bk.resid",abcd.dat,"Site")),
         probs=c(0.025,0.50,0.975))
