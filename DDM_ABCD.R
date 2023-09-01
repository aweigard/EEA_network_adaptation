rm(list=ls())
source ("dmc/dmc.R")
source("read_abcd.R")

################################################################
####### format n-back data and apply exclusion criteria ########
################################################################

# load trial-level n-back data, merged across all 
# baseline participants (once downloaded from aws)
abcd.nback<-read.csv("nback_revised.csv")

# make DMC-based columns

abcd.nback$S<-as.character(abcd.nback$enback_targettype)
abcd.nback[abcd.nback$S=="target",]$S<-"tar"
abcd.nback[abcd.nback$S=="lure",]$S<-"lur"
abcd.nback[abcd.nback$S=="nonlure",]$S<-"non"
abcd.nback$S<-factor(abcd.nback$S)

abcd.nback$R<-NA
abcd.nback[abcd.nback$S%in%c("tar") & abcd.nback$enback_stim_acc==1 & 
             !is.na(abcd.nback$enback_stim_resp),]$R<-"YES"
abcd.nback[abcd.nback$S%in%c("tar") & abcd.nback$enback_stim_acc==0 & 
             !is.na(abcd.nback$enback_stim_resp),]$R<-"NO"
abcd.nback[abcd.nback$S%in%c("lur","non") & abcd.nback$enback_stim_acc==1 & 
             !is.na(abcd.nback$enback_stim_resp),]$R<-"NO"
abcd.nback[abcd.nback$S%in%c("lur","non") & abcd.nback$enback_stim_acc==0 & 
             !is.na(abcd.nback$enback_stim_resp),]$R<-"YES"
abcd.nback$R<-factor(abcd.nback$R)

abcd.nback$RT<-NA
abcd.nback[!is.na(abcd.nback$R),]$RT<-abcd.nback[!is.na(abcd.nback$R),]$enback_stim_rt/1000

abcd.nback$acc<-NA
abcd.nback[!is.na(abcd.nback$R),]$acc<-abcd.nback[!is.na(abcd.nback$R),]$enback_stim_acc

# separate into baseline and year 2

nback.base<-abcd.nback[abcd.nback$eventname=="baseline_year_1_arm_1",]
nback.y2<-abcd.nback[abcd.nback$eventname=="2_year_follow_up_y_arm_1",]

#compute summary stats

sum.base<-data.frame(unique(nback.base$subject));colnames(sum.base)<-"s"
sum.y2<-data.frame(unique(nback.y2$subject));colnames(sum.y2)<-"s"

sum.base$mrt.0.tar<-NA
sum.base$mrt.0.lur<-NA
sum.base$mrt.0.non<-NA
sum.base$sdrt.0.tar<-NA
sum.base$sdrt.0.lur<-NA
sum.base$sdrt.0.non<-NA
sum.base$acc.0.tar<-NA
sum.base$acc.0.lur<-NA
sum.base$acc.0.non<-NA

sum.base$acc.0.total<-NA
sum.base$p.omit.0<-NA

sum.base$mrt.2.tar<-NA
sum.base$mrt.2.lur<-NA
sum.base$mrt.2.non<-NA
sum.base$sdrt.2.tar<-NA
sum.base$sdrt.2.lur<-NA
sum.base$sdrt.2.non<-NA
sum.base$acc.2.tar<-NA
sum.base$acc.2.lur<-NA
sum.base$acc.2.non<-NA

sum.base$acc.2.total<-NA
sum.base$p.omit.2<-NA

for (r in 1:length(sum.base$s)){
  tmp0<-nback.base[nback.base$subject==sum.base$s[r] & nback.base$enback_loadcon=="0-Back",]
  tmp2<-nback.base[nback.base$subject==sum.base$s[r] & nback.base$enback_loadcon=="2-Back",]
  
  sum.base$mrt.0.tar[r]<-mean(tmp0[tmp0$S=="tar",]$RT,na.rm=TRUE)
  sum.base$mrt.0.lur[r]<-mean(tmp0[tmp0$S=="lur",]$RT,na.rm=TRUE)
  sum.base$mrt.0.non[r]<-mean(tmp0[tmp0$S=="non",]$RT,na.rm=TRUE)
  sum.base$sdrt.0.tar[r]<-sd(tmp0[tmp0$S=="tar",]$RT,na.rm=TRUE)
  sum.base$sdrt.0.lur[r]<-sd(tmp0[tmp0$S=="lur",]$RT,na.rm=TRUE)
  sum.base$sdrt.0.non[r]<-sd(tmp0[tmp0$S=="non",]$RT,na.rm=TRUE)
  sum.base$acc.0.tar[r]<-mean(tmp0[tmp0$S=="tar",]$acc,na.rm=TRUE)
  sum.base$acc.0.lur[r]<-mean(tmp0[tmp0$S=="lur",]$acc,na.rm=TRUE)
  sum.base$acc.0.non[r]<-mean(tmp0[tmp0$S=="non",]$acc,na.rm=TRUE)
  
  sum.base$acc.0.total[r]<-mean(tmp0$acc,na.rm=TRUE)
  sum.base$p.omit.0[r]<-mean(is.na(tmp0$acc))
  
  sum.base$mrt.2.tar[r]<-mean(tmp2[tmp2$S=="tar",]$RT,na.rm=TRUE)
  sum.base$mrt.2.lur[r]<-mean(tmp2[tmp2$S=="lur",]$RT,na.rm=TRUE)
  sum.base$mrt.2.non[r]<-mean(tmp2[tmp2$S=="non",]$RT,na.rm=TRUE)
  sum.base$sdrt.2.tar[r]<-sd(tmp2[tmp2$S=="tar",]$RT,na.rm=TRUE)
  sum.base$sdrt.2.lur[r]<-sd(tmp2[tmp2$S=="lur",]$RT,na.rm=TRUE)
  sum.base$sdrt.2.non[r]<-sd(tmp2[tmp2$S=="non",]$RT,na.rm=TRUE)
  sum.base$acc.2.tar[r]<-mean(tmp2[tmp2$S=="tar",]$acc,na.rm=TRUE)
  sum.base$acc.2.lur[r]<-mean(tmp2[tmp2$S=="lur",]$acc,na.rm=TRUE)
  sum.base$acc.2.non[r]<-mean(tmp2[tmp2$S=="non",]$acc,na.rm=TRUE)
  
  sum.base$acc.2.total[r]<-mean(tmp2$acc,na.rm=TRUE)
  sum.base$p.omit.2[r]<-mean(is.na(tmp2$acc))
  
}

save(sum.base,file="sum_base_fits.RData")
load("sum_base_fits.RData")

# inclusion for each task

inc.0b<-sum.base[sum.base$acc.0.total>=.55 & sum.base$p.omit.0<=.25,"s"]
inc.2b<-sum.base[sum.base$acc.2.total>=.55 & sum.base$p.omit.2<=.25,"s"]

# final groups

base.0.dmc<-nback.base[nback.base$enback_loadcon=="0-Back" & nback.base$subject%in%inc.0b,]
base.0.dmc<-base.0.dmc[,c("subject","S","R","RT")]
colnames(base.0.dmc)<-c("s","S","R","RT")
base.0.dmc$s<-factor(base.0.dmc$s)

base.2.dmc<-nback.base[nback.base$enback_loadcon=="2-Back" & nback.base$subject%in%inc.2b,]
base.2.dmc<-base.2.dmc[,c("subject","S","R","RT")]
colnames(base.2.dmc)<-c("s","S","R","RT")
base.2.dmc$s<-factor(base.2.dmc$s)

# exclude fast guesses

base.0.dmc<-base.0.dmc[base.0.dmc$RT>=0.200 | is.na(base.0.dmc$RT),]
base.2.dmc<-base.2.dmc[base.2.dmc$RT>=0.200 | is.na(base.2.dmc$RT),]

# save out

save(base.0.dmc,file="base_0_dat_DDM.RData")
save(base.2.dmc,file="base_2_dat_DDM.RData")


######################################################
####### select hierarchical subsamples (for priors) ##
######################################################

# load summary data

load("sum_base_fits.RData")

# exclude data from individuals with accuracy <55% and
# omissions > 25% on either task

sum.base[sum.base$acc.0.overall<0.55 | 
           sum.base$acc.0.omit>0.25 | is.na(sum.base$acc.0.overall),
         colnames(sum.base)[grepl(".0.",colnames(sum.base))] ]<-NA

sum.base[sum.base$acc.2.overall<0.55 | 
           sum.base$acc.2.omit>0.25 | is.na(sum.base$acc.2.overall),
         colnames(sum.base)[grepl(".2.",colnames(sum.base))] ]<-NA

# identify subjects with acceptable data from both tasks

sum.base$include<-(is.na(sum.base$acc.0.omit)+is.na(sum.base$acc.2.omit))<1

# load family data

acspsw03<-read.delim("acspsw03.txt")
acspsw03.labs<-acspsw03[1,]
acspsw03<-data.frame(acspsw03[-1,])
row.names(acspsw03)<-1:length(acspsw03[,1])

fam.dat<-acspsw03[acspsw03$eventname=="baseline_year_1_arm_1",]
fam.dat<-data.frame(s=fam.dat$subject,
                    fam=fam.dat$rel_family_id)

# identify singletons 

singles<-table(fam.dat$fam)
singles<-singles[singles==1]
singles<-names(singles)

# load fMRI inclusion data

imaging<-read.csv("ABCD_task_general.csv")
imaging$s<-imaging$subjectkey

# merge fam.dat and imaging with sum.base

sum.fam<-merge(sum.base,fam.dat,by="s",all.x = TRUE)
sum.fam<-merge(sum.fam,imaging,by="s",all.x = TRUE)

# only look at individuals not included in the fMRI analyses

sum.fam<-sum.fam[!sum.fam$Include.nback,]

# select subsample IDs at random

sub.ids<-sum.fam[sum.fam$fam%in%singles & sum.fam$include,]

sub.ids<-sample(sub.ids$s,size = 300,replace = FALSE)

write.csv(sub.ids,
          file="subsample_fmri_project.csv",
          row.names = FALSE)


#####################################################################
###### individual-level estimation: broad priors ####################
#####################################################################

### load data ###

load("base_0_dat_DDM.RData")

load("base_2_dat_DDM.RData")

#### Design and contaminant (gf) ----
load_model ("DDM","ddm_omit.R")

model <- model.dmc(
  p.map = list(a="1",v=c("S"),z="1",d="1",sz="1",
               sv="1",t0="1",st0="1",
               censor="1",gf="1"), 
  responses = c("NO","YES"),
  match.map = list(M=list(tar="YES",lur="NO",non="NO")),
  factors=list(S=c("lur","non","tar")),
  constants = c(sv=0,sz=0,d=0,censor=2.00), 
  type="rd")

# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "a"     "v.lur" "v.non" "v.tar" "z"     "t0"   
# [7] "st0"   "gf"   
# 
# Constants are (see attr(,"constants") ):
#   sz     sv      d censor 
# 0      0      0      2 


p1 <- c(a=1,v.lur=3, v.non = 3, v.tar= 3,
        z=0.5,t0=0.3, st0=0.1,gf=0)

p.prior <- prior.p.dmc(
  dists = rep("tnorm",8),
  p1=p1,                           
  p2=c(.5,1,1,1,.1,.1,.05,1),
  lower=c(0,NA,NA,NA,0,0,0,NA),upper=c(NA,NA,NA,NA,1,2,2,NA)
)
# par(mfcol=c(1,3)); for (i in names(p.prior)) plot.prior(i,p.prior)


zero.base.mdi <- data.model.dmc(base.0.dmc, model)

two.base.mdi <- data.model.dmc(base.2.dmc, model)

sZero.base  <- h.samples.dmc(nmc = 120, p.prior=p.prior,data = zero.base.mdi, thin = 5)
sZero.base.first<-sZero.base[1:4500]
save(sZero.base.first,file="zero_ddm_broad_first.RData")
sZero.base.second<-sZero.base[4501:length(sZero.base)]
save(sZero.base.second,file="zero_ddm_broad_second.RData")

sTwo.base  <- h.samples.dmc(nmc = 120, p.prior=p.prior,data = two.base.mdi, thin = 5)
sTwo.base.first<-sTwo.base[1:4500]
save(sTwo.base.first,file="two_ddm_broad_first.RData")
sTwo.base.second<-sTwo.base[4501:length(sTwo.base)]
save(sTwo.base.second,file="two_ddm_broad_second.RData")

# sampling can be conducted efficiently by using a high-performance 
# computer to estimate parameters for each participant in parallel
# using the RUN.dmc() function (with default parameters) as described
# in the Dynamic Models of Choice tutorial: https://osf.io/pbwx8/


####################################################################
### hierarchical model fits in independent subsample (for priors) ###
#####################################################################

### load data ###

load("base_0_dat_DDM.RData")

load("base_2_dat_DDM.RData")

subsample<-read.csv("subsample_fmri_project.csv")

base.0.dmc<-base.0.dmc[base.0.dmc$s%in%subsample$x,]
base.0.dmc$s<-as.factor(as.character(base.0.dmc$s))

base.2.dmc<-base.2.dmc[base.2.dmc$s%in%subsample$x,]
base.2.dmc$s<-as.factor(as.character(base.2.dmc$s))

#### Design and contaminant (gf) ----
load_model ("DDM","ddm_omit.R")

model <- model.dmc(
  p.map = list(a="1",v=c("S"),z="1",d="1",sz="1",
               sv="1",t0="1",st0="1",
               censor="1",gf="1"), 
  responses = c("NO","YES"),
  match.map = list(M=list(tar="YES",lur="NO",non="NO")),
  factors=list(S=c("lur","non","tar")),
  constants = c(sz = 0, sv = 0, d=0,censor=2.00), 
  type="rd")

p1 <- c(a=1,v.lur=3, v.non = 3, v.tar= 3,
        z=0.5,t0=0.3, st0=.1, gf=0)

p.prior <- prior.p.dmc(
  dists = rep("tnorm",8),
  p1=p1,                           
  p2=c(.5,1,1,1,.1,.1,.05,1),
  lower=c(0,NA,NA,NA,0,0,0,NA),upper=c(NA,NA,NA,NA,1,2,2,NA)
)
# par(mfcol=c(1,3)); for (i in names(p.prior)) plot.prior(i,p.prior)

# same prior for mu
p.prior.mu<-p.prior 

# gamma priors for sigma
p1[1:length(p1)] <- rep(1,length(p1))
p.prior.sigma<-prior.p.dmc(
  dists=rep("gamma",length(p1)),
  p1=p1,p2=p1)
# par(mfcol=c(1,3)); for (i in names(p.prior.sigma)) plot.prior(i,p.prior.sigma)

# combine priors for hyperparameters
pp.prior<-list(p.prior.mu,p.prior.sigma)


# make mdi

zero.hier.mdi <- data.model.dmc(base.0.dmc, model)

two.hier.mdi <- data.model.dmc(base.2.dmc, model)


#start points
load("nback_ddm_broad_samples.RData")

hstart.zero <- make.hstart(sZero.base[subsample$x])
theta1.zero <- make.theta1(sZero.base[subsample$x])

hstart.two <- make.hstart(sTwo.base[subsample$x])
theta1.two <- make.theta1(sTwo.base[subsample$x])


# make sampling objects
sZero.hier <- h.samples.dmc(nmc = 40, p.prior=p.prior, 
                            data=zero.hier.mdi, 
                            pp.prior=pp.prior,
                            thin=10,
                            hstart.prior=hstart.zero,
                            theta1=theta1.zero)

save(sZero.hier ,file="zero_hier_fmri_project.RData")


sTwo.hier <- h.samples.dmc(nmc = 40, p.prior=p.prior, 
                           data=two.hier.mdi, 
                           pp.prior=pp.prior,
                           thin=10,
                           hstart.prior=hstart.two,
                           theta1=theta1.two)

save(sTwo.hier ,file="two_hier_fmri_project.RData")


# run the following code on a high-performance computer 
# to efficiently generate samples with parallel processing:

rm(list=ls())
source ("dmc/dmc.R")
load_model ("DDM","ddm_omit.R") 

load("zero_hier_fmri_project.RData")

cores=32

sZero.hier <- h.run.unstuck.dmc(sZero.hier,cores=cores,report=10,p.migrate=0.025,h.p.migrate=0.025)

save.image("zero_hier_fmri_project.RData")

sZero.hier1 <- h.run.converge.dmc(samples=h.samples.dmc(samples=sZero.hier,nmc=40,thin=10),
                                  nmc=40,cores=cores,report=10,verbose=TRUE)

save.image("zero_hier_fmri_project.RData")

sZero.hier2 <- h.run.converge.dmc(samples=h.samples.dmc(samples=sZero.hier1,nmc=120,thin=10),
                                  nmc=40,cores=cores,report=10,verbose=TRUE)

save.image("zero_hier_fmri_project.RData")

rm(list=ls())
source ("dmc/dmc.R")
load_model ("DDM","ddm_omit.R") 

load("two_hier_fmri_project.RData")

cores=32

sTwo.hier <- h.run.unstuck.dmc(sTwo.hier,cores=cores,report=10,p.migrate=0.025,h.p.migrate=0.025)

save.image("two_hier_fmri_project.RData")

sTwo.hier1 <- h.run.converge.dmc(samples=h.samples.dmc(samples=sTwo.hier,nmc=40,thin=10),
                                 nmc=40,cores=cores,report=10,verbose=TRUE)

save.image("two_hier_fmri_project.RData")

sTwo.hier2 <- h.run.converge.dmc(samples=h.samples.dmc(samples=sTwo.hier1,nmc=120,thin=10),
                                 nmc=40,cores=cores,report=10,verbose=TRUE)

save.image("two_hier_fmri_project.RData")


###### results ######

# read in data

load("zero_hier_fmri_project.RData")
load("two_hier_fmri_project.RData")

plot.dmc(sZero.hier2,hyper=TRUE,layout=c(2,2))
gelman.diag.dmc(sZero.hier2,hyper=TRUE)

plot.dmc(sTwo.hier2,hyper=TRUE,layout=c(2,2))
gelman.diag.dmc(sTwo.hier2,hyper=TRUE)

#### model fit

sim.Zero.hier <- h.post.predict.dmc(sZero.hier2,cores=3)
sim.Two.hier <- h.post.predict.dmc(sTwo.hier2,cores=3)

save(sim.Zero.hier,sim.Two.hier,file="nback_hier_fmri_project_pp.RData")
#

plot.pp.dmc(sim.Zero.hier,layout = c(1,3))

plot.pp.dmc(sim.Two.hier,layout = c(1,3))

##### make informed priors

library(MASS)

# combine all individuals' samples into a single vector for each parameter

zero.samps<-vector("list", length = length(attr(model,"p.vector"))) 
names(zero.samps)<-names(attr(model,"p.vector"))

for (s in 1:length(sZero.hier2)){
  tmp<-sZero.hier2[[s]]
  for (par in names(attr(model,"p.vector"))){
    zero.samps[[par]]<-c(zero.samps[[par]],as.vector(tmp$theta[,par,]))
  }  
}


two.samps<-vector("list", length = length(attr(model,"p.vector"))) 
names(two.samps)<-names(attr(model,"p.vector"))

for (s in 1:length(sTwo.hier2)){
  tmp<-sTwo.hier2[[s]]
  for (par in names(attr(model,"p.vector"))){
    two.samps[[par]]<-c(two.samps[[par]],as.vector(tmp$theta[,par,]))
  }  
}


# upper and lower bounds

lower<-c(0,-Inf,-Inf,-Inf,0,0,0,-Inf)
upper<-c(Inf,Inf,Inf,Inf,1,2,2,Inf) 
names(lower)<-names(p.prior)
names(upper)<-names(p.prior)


h1.gauss.fits.zero<-data.frame(names(attr(model,"p.vector")),NA,NA)
h1.gauss.fits.two<-data.frame(names(attr(model,"p.vector")),NA,NA)
colnames(h1.gauss.fits.zero)<-c("par","m","sd")
colnames(h1.gauss.fits.two)<-c("par","m","sd")

for (par in names(attr(model,"p.vector"))){
  dtnorm0 <- function(X, mean, sd, log = FALSE,low,up) {
    dtnorm(X, mean, sd, low, up,log)}
  
  G.tmp<-fitdistr(zero.samps[[par]], dtnorm0, 
                  start=list(mean=mean(zero.samps[[par]]), 
                             sd=sd(zero.samps[[par]])),
                  low=lower[par],up=upper[par])
  h1.gauss.fits.zero[h1.gauss.fits.zero$par==par,2:3]<-G.tmp$estimate
  
  G.tmp<-fitdistr(two.samps[[par]], dtnorm0, 
                  start=list(mean=mean(two.samps[[par]]), 
                             sd=sd(two.samps[[par]])),
                  low=lower[par],up=upper[par])
  h1.gauss.fits.two[h1.gauss.fits.two$par==par,2:3]<-G.tmp$estimate
}

h1.gauss.fits.zero$m<-as.numeric(h1.gauss.fits.zero$m)
h1.gauss.fits.zero$sd<-as.numeric(h1.gauss.fits.zero$sd)
h1.gauss.fits.two$m<-as.numeric(h1.gauss.fits.two$m)
h1.gauss.fits.two$sd<-as.numeric(h1.gauss.fits.two$sd)


# write out
write.csv(h1.gauss.fits.zero,file = "zero_informed_priors_fmri_project.csv",row.names = FALSE)
write.csv(h1.gauss.fits.two,file = "two_informed_priors_fmri_project.csv",row.names = FALSE)


#plot fits

jpeg(filename = "zero_fits_dists_fmri_project.jpeg",
     units = "in", res = 300,
     width = 12,height = 9)
par(mfrow=c(3,4))
for (par in names(attr(model,"p.vector"))){
  msd<-as.numeric(h1.gauss.fits.zero[h1.gauss.fits.zero$par==par,2:3])
  hist(zero.samps[[par]], prob = TRUE,main = par,xlab= par)
  curve(dtnorm0(x, msd[1], msd[2],low=lower[par],up=upper[par]), 
        col = "red", add = TRUE)
}
dev.off()

jpeg(filename = "two_fits_dists_fmri_project.jpeg",
     units = "in", res = 300,
     width = 12,height = 9)
par(mfrow=c(3,4))
for (par in names(attr(model,"p.vector"))){
  msd<-as.numeric(h1.gauss.fits.two[h1.gauss.fits.two$par==par,2:3])
  hist(two.samps[[par]], prob = TRUE,main = par,xlab= par)
  curve(dtnorm0(x, msd[1], msd[2],low=lower[par],up=upper[par]), 
        col = "red", add = TRUE)
}
dev.off()

#####################################################################
###### individual-level estimation: informed priors #################
#####################################################################

### load data ###

load("base_0_dat_DDM.RData")

load("base_2_dat_DDM.RData")

subsample<-read.csv("nback_subsample_fmri_project.csv")

base.0.dmc<-base.0.dmc[!base.0.dmc$s%in%subsample$x,]
base.0.dmc$s<-as.factor(as.character(base.0.dmc$s))

base.2.dmc<-base.2.dmc[!base.2.dmc$s%in%subsample$x,]
base.2.dmc$s<-as.factor(as.character(base.2.dmc$s))

#### Design and contaminant (gf) ----
load_model ("DDM","ddm_omit.R")

model <- model.dmc(
  p.map = list(a="1",v=c("S"),z="1",d="1",sz="1",
               sv="1",t0="1",st0="1",
               censor="1",gf="1"), 
  responses = c("NO","YES"),
  match.map = list(M=list(tar="YES",lur="NO",non="NO")),
  factors=list(S=c("lur","non","tar")),
  constants = c(sz = 0, sv = 0, d=0,censor=2.00), 
  type="rd")

# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "a"     "v.lur" "v.non" "v.tar" "z"     "t0"   
# [7] "st0"   "gf"   
# 
# Constants are (see attr(,"constants") ):
#   sz     sv      d censor 
# 0      0      0      2 


hp0<-read.csv("zero_informed_priors_fmri_project.csv")
rownames(hp0)<-hp0$par

p1.zero <- c(a=hp0["a","m"],v.lur=hp0["v.lur","m"],v.non=hp0["v.non","m"],
             v.tar=hp0["v.tar","m"],z=hp0["z","m"],t0=hp0["t0","m"],
             st0=hp0["st0","m"],gf=hp0["gf","m"])
p2.zero <- c(a=hp0["a","sd"],v.lur=hp0["v.lur","sd"],v.non=hp0["v.non","sd"],
             v.tar=hp0["v.tar","sd"],z=hp0["z","sd"],t0=hp0["t0","sd"],
             st0=hp0["st0","sd"],gf=hp0["gf","sd"])

p.prior.zero <- prior.p.dmc(
  dists = rep("tnorm",8),
  p1=p1.zero,                           
  p2=p2.zero,
  lower=c(0,NA,NA,NA,0,0,0,NA),upper=c(NA,NA,NA,NA,1,2,2,NA)
)
# par(mfcol=c(1,3)); for (i in names(p.prior.zero)) plot.prior(i,p.prior.zero)


hp2<-read.csv("two_informed_priors_fmri_project.csv")
rownames(hp2)<-hp2$par

p1.two <- c(a=hp2["a","m"],v.lur=hp2["v.lur","m"],v.non=hp2["v.non","m"],
            v.tar=hp2["v.tar","m"],z=hp2["z","m"],t0=hp2["t0","m"],
            st0=hp2["st0","m"],gf=hp2["gf","m"])
p2.two <- c(a=hp2["a","sd"],v.lur=hp2["v.lur","sd"],v.non=hp2["v.non","sd"],
            v.tar=hp2["v.tar","sd"],z=hp2["z","sd"],t0=hp2["t0","sd"],
            st0=hp2["st0","sd"],gf=hp2["gf","sd"])

p.prior.two <- prior.p.dmc(
  dists = rep("tnorm",8),
  p1=p1.two,                           
  p2=p2.two,
  lower=c(0,NA,NA,NA,0,0,0,NA),upper=c(NA,NA,NA,NA,1,2,2,NA)
)
# par(mfcol=c(1,3)); for (i in names(p.prior.two)) plot.prior(i,p.prior.two)

zero.base.mdi <- data.model.dmc(base.0.dmc, model)

two.base.mdi <- data.model.dmc(base.2.dmc, model)

sZero.base  <- h.samples.dmc(nmc = 120, p.prior=p.prior.zero,data = zero.base.mdi, thin = 5)
sTwo.base  <- h.samples.dmc(nmc = 120, p.prior=p.prior.two,data = two.base.mdi, thin = 5)


sZero.base.first<-sZero.base[1:4500]
save(sZero.base.first,file="zero_ddm_inf_fmri_first.RData")
sZero.base.second<-sZero.base[4501:length(sZero.base)]
save(sZero.base.second,file="zero_ddm_inf_fmri_second.RData")

sTwo.base.first<-sTwo.base[1:4500]
save(sTwo.base.first,file="two_ddm_inf_fmri_first.RData")
sTwo.base.second<-sTwo.base[4501:length(sTwo.base)]
save(sTwo.base.second,file="two_ddm_inf_fmri_second.RData")

# sampling can be conducted efficiently by using a high-performance 
# computer to estimate parameters for each participant in parallel
# using the RUN.dmc() function (with default parameters) as described
# in the Dynamic Models of Choice tutorial: https://osf.io/pbwx8/


##############################################################
##### generate point estimates (posterior medians)############
##############################################################

#### load in parts, merge them, save out

source ("dmc/dmc.R")
load_model ("DDM","ddm_omit.R")

load("zero_ddm_inf_fmri_first.RData")
load("zero_ddm_inf_fmri_second.RData")

sZero.base<-sZero.base.first
sZero.base[4501:(4500+length(sZero.base.second))]<-sZero.base.second
names(sZero.base)<-c(names(sZero.base.first),names(sZero.base.second))

load("two_ddm_inf_fmri_first.RData")
load("two_ddm_inf_fmri_second.RData")

sTwo.base<-sTwo.base.first
sTwo.base[4501:(4500+length(sTwo.base.second))]<-sTwo.base.second
names(sTwo.base)<-c(names(sTwo.base.first),names(sTwo.base.second))


save(sZero.base,sTwo.base,file="nback_ddm_inf_fmri_samples.RData")

#### check convergence, remove non-converged if needed

# load("nback_ddm_inf_fmri_samples.RData")

rhat.Zero<-gelman.diag.dmc(sZero.base, hyper = FALSE)
rhat.Two<-gelman.diag.dmc(sTwo.base, hyper = FALSE)
save(rhat.Zero,rhat.Two,file="nback_ddm_inf_fmri_rhat.RData")

#### parameter estimates

# loop to make summary parameters

sZero.base.pars<-sZero.base

for (s in 1:length(sZero.base)){
  t<-sZero.base[[s]]$theta; d<-dim(t)
  Ter<-t[,"t0",]+(t[,"st0",]/2)
  gf.natural<-pnorm(t[,"gf",])
  D <- d + c(0, 2, 0)
  t2<-array(sapply(1:D[3], function(x) cbind(t[,,x], Ter[,x],
                                             gf.natural[,x])), D)
  dimnames(t2)[[2]]<-c(dimnames(t)[[2]],"Ter","gf.natural")
  sZero.base.pars[[s]]$theta<-t2
}

sTwo.base.pars<-sTwo.base

for (s in 1:length(sTwo.base)){
  t<-sTwo.base[[s]]$theta; d<-dim(t)
  Ter<-t[,"t0",]+(t[,"st0",]/2)
  gf.natural<-pnorm(t[,"gf",])
  D <- d + c(0, 2, 0)
  t2<-array(sapply(1:D[3], function(x) cbind(t[,,x], Ter[,x],
                                             gf.natural[,x])), D)
  dimnames(t2)[[2]]<-c(dimnames(t)[[2]],"Ter","gf.natural")
  sTwo.base.pars[[s]]$theta<-t2
}


# save out point estimates 

sZero.base_medians<-lapply(sZero.base.pars,FUN=function(x) apply(x$theta,2,median))
sZero.base_medians<-as.data.frame(t(as.data.frame(sZero.base_medians)))

sTwo.base_medians<-lapply(sTwo.base.pars,FUN=function(x) apply(x$theta,2,median))
sTwo.base_medians<-as.data.frame(t(as.data.frame(sTwo.base_medians)))

colnames(sZero.base_medians)<-paste0("zero.",colnames(sZero.base_medians))
colnames(sTwo.base_medians)<-paste0("two.",colnames(sTwo.base_medians))

# merge
base.medians<-merge(sZero.base_medians, sTwo.base_medians,
                    by="row.names",all = TRUE)

# exclude anyone with data that failed QC
base.medians<-base.medians[!is.na(base.medians$zero.a) & !is.na(base.medians$two.a),]

colnames(base.medians)[1]<-"s"

base.medians$zero.v.mean<-rowMeans(base.medians[,c("zero.v.lur","zero.v.tar","zero.v.non")])
base.medians$two.v.mean<-rowMeans(base.medians[,c("two.v.lur","two.v.tar","two.v.non")])
base.medians$overall.v.mean<-rowMeans(base.medians[,c("zero.v.lur","zero.v.tar","zero.v.non",
                                                      "two.v.lur","two.v.tar","two.v.non")])
row.names(base.medians)<-base.medians$s

base.medians<-base.medians[,-1]

write.csv(base.medians,file="abcd_base_nback_fmri_postmedians.csv")

##############################################################
##### model fit plots ########################################
##############################################################

load("nback_ddm_inf_fmri_samples.RData")

# final sample with included fMRI data
final<-read.csv("abcd_base_nback_postmedians_v4_asw.csv")

sZero.paper<-sZero.base[names(sZero.base)%in%final$subjectkey]
sTwo.paper<-sTwo.base[names(sTwo.base)%in%final$subjectkey]
save(sZero.paper,sTwo.paper,file="nback_fmri_inf_paper_samples.RData")

sim.Zero.paper <- h.post.predict.dmc(sZero.paper,cores=20)
sim.Two.paper <- h.post.predict.dmc(sTwo.paper,cores=20)

save(sim.Zero.paper,sim.Two.paper,file="nback_fmri_inf_paper_pp.RData")


plot.pp.dmc(sim.Zero.paper,layout = c(1,3),
            model.legend = F,show.fits=FALSE,mar=c(4,5,3,1),pos=NA)

plot.pp.dmc(sim.Two.paper,layout = c(1,3),
            model.legend = F,show.fits=FALSE,mar=c(4,5,3,1),pos=NA)

