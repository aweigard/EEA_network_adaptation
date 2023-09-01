rm(list=ls())
source ("dmc/dmc.R")
source("read_abcd.R")

################################################################
####### format n-back data and apply exclusion criteria ########
################################################################

# load data
zero.dat<-read.csv("HCP_0bk_trials.csv")
two.dat<-read.csv("HCP_2bk_trials.csv")

# exclude session 2 data
zero.dat<-zero.dat[zero.dat$session==1,]
two.dat<-two.dat[two.dat$session==1,]

# make and format DMC columns

zero.dat$s<-zero.dat$subject_num
zero.dat$S<-zero.dat$targettype
zero.dat$R<-zero.dat$resp

zero.dat[zero.dat$targettype=="target","S"]<-"tar"
zero.dat[zero.dat$targettype=="lure","S"]<-"lur"
zero.dat[zero.dat$targettype=="nonlure","S"]<-"non"

zero.dat[zero.dat$resp==2 & !is.na(zero.dat$resp),"R"]<-"YES"
zero.dat[zero.dat$resp==3 & !is.na(zero.dat$resp),"R"]<-"NO"

zero.dat[is.na(zero.dat$resp),"RT"]<-NA

two.dat$s<-two.dat$subject_num
two.dat$S<-two.dat$targettype
two.dat$R<-two.dat$resp

two.dat[two.dat$targettype=="target","S"]<-"tar"
two.dat[two.dat$targettype=="lure","S"]<-"lur"
two.dat[two.dat$targettype=="nonlure","S"]<-"non"

two.dat[two.dat$resp==2 & !is.na(two.dat$resp),"R"]<-"YES"
two.dat[two.dat$resp==3 & !is.na(two.dat$resp),"R"]<-"NO"

two.dat[is.na(two.dat$resp),"RT"]<-NA

#trim fast guesses
zero.dat.trim<-zero.dat[zero.dat$RT>=0.200 | is.na(zero.dat$RT),]

two.dat.trim<-two.dat[two.dat$RT>=0.200 | is.na(two.dat$RT),]

#exclude two responses on an extra key
two.dat.trim<-two.dat.trim[two.dat.trim$resp!=4 | is.na(two.dat.trim$resp),]

# summary accuracy stats
sum.acc.dat<-data.frame(unique(zero.dat$subject_num))
colnames(sum.acc.dat)<-"s"


sum.acc.dat$zero.acc<-NA
sum.acc.dat$zero.tar.acc<-NA
sum.acc.dat$zero.lur.acc<-NA
sum.acc.dat$zero.non.acc<-NA
sum.acc.dat$zero.omit<-NA
sum.acc.dat$zero.tar.n<-NA
sum.acc.dat$zero.lur.n<-NA
sum.acc.dat$zero.non.n<-NA
sum.acc.dat$zero.n<-NA

sum.acc.dat$two.acc<-NA
sum.acc.dat$two.tar.acc<-NA
sum.acc.dat$two.lur.acc<-NA
sum.acc.dat$two.non.acc<-NA
sum.acc.dat$two.omit<-NA
sum.acc.dat$two.tar.n<-NA
sum.acc.dat$two.lur.n<-NA
sum.acc.dat$two.non.n<-NA
sum.acc.dat$two.n<-NA


for (s in sum.acc.dat$s){
  
  zero.tmp<-zero.dat.trim[zero.dat.trim$s==s & !is.na(zero.dat.trim$RT),]
  two.tmp<-two.dat.trim[two.dat.trim$s==s & !is.na(two.dat.trim$RT),]
  
  sum.acc.dat[sum.acc.dat$s==s,]$zero.acc<-mean(zero.tmp$accuracy)
  sum.acc.dat[sum.acc.dat$s==s,]$zero.tar.acc<-mean(zero.tmp[zero.tmp$S=="tar",]$accuracy)
  sum.acc.dat[sum.acc.dat$s==s,]$zero.lur.acc<-mean(zero.tmp[zero.tmp$S=="lur",]$accuracy)
  sum.acc.dat[sum.acc.dat$s==s,]$zero.non.acc<-mean(zero.tmp[zero.tmp$S=="non",]$accuracy)
  sum.acc.dat[sum.acc.dat$s==s,]$zero.omit<-mean(is.na(zero.dat.trim[zero.dat.trim$s==s,"RT"]))
  sum.acc.dat[sum.acc.dat$s==s,]$zero.tar.n<-length(zero.tmp[zero.tmp$S=="tar",]$accuracy)
  sum.acc.dat[sum.acc.dat$s==s,]$zero.lur.n<-length(zero.tmp[zero.tmp$S=="lur",]$accuracy)
  sum.acc.dat[sum.acc.dat$s==s,]$zero.non.n<-length(zero.tmp[zero.tmp$S=="non",]$accuracy)
  sum.acc.dat[sum.acc.dat$s==s,]$zero.n<-length(zero.tmp$accuracy)
  
  sum.acc.dat[sum.acc.dat$s==s,]$two.acc<-mean(two.tmp$accuracy)
  sum.acc.dat[sum.acc.dat$s==s,]$two.tar.acc<-mean(two.tmp[two.tmp$S=="tar",]$accuracy)
  sum.acc.dat[sum.acc.dat$s==s,]$two.lur.acc<-mean(two.tmp[two.tmp$S=="lur",]$accuracy)
  sum.acc.dat[sum.acc.dat$s==s,]$two.non.acc<-mean(two.tmp[two.tmp$S=="non",]$accuracy)
  sum.acc.dat[sum.acc.dat$s==s,]$two.omit<-mean(is.na(two.dat.trim[two.dat.trim$s==s,"RT"]))
  sum.acc.dat[sum.acc.dat$s==s,]$two.tar.n<-length(two.tmp[two.tmp$S=="tar",]$accuracy)
  sum.acc.dat[sum.acc.dat$s==s,]$two.lur.n<-length(two.tmp[two.tmp$S=="lur",]$accuracy)
  sum.acc.dat[sum.acc.dat$s==s,]$two.non.n<-length(two.tmp[two.tmp$S=="non",]$accuracy)
  sum.acc.dat[sum.acc.dat$s==s,]$two.n<-length(two.tmp$accuracy)
}

write.csv(sum.acc.dat,file="sum_hcp_nback.csv",row.names=F)

# dmc data

zero.dat.dmc<-zero.dat.trim[,c("s","S","R","RT")]
# zero.dat.dmc[is.na(zero.dat.dmc$R),"R"]<-sample(c("NO","YES"),
#                                             length(zero.dat.dmc[is.na(zero.dat.dmc$R),"R"]),TRUE)
zero.dat.dmc$s<-as.factor(zero.dat.dmc$s)
zero.dat.dmc$S<-as.factor(zero.dat.dmc$S)
zero.dat.dmc$R<-as.factor(zero.dat.dmc$R)

two.dat.dmc<-two.dat.trim[,c("s","S","R","RT")]
# two.dat.dmc[is.na(two.dat.dmc$R),"R"]<-sample(c("NO","YES"),
#                                                 length(two.dat.dmc[is.na(two.dat.dmc$R),"R"]),TRUE)
two.dat.dmc$s<-as.factor(two.dat.dmc$s)
two.dat.dmc$S<-as.factor(two.dat.dmc$S)
two.dat.dmc$R<-as.factor(two.dat.dmc$R)


save(zero.dat.dmc,two.dat.dmc,sum.acc.dat,
     file="HCP_DMC_dat.RData")

write.csv(zero.dat.dmc,file="HCP_zero_DMC.csv",row.names = FALSE)
write.csv(two.dat.dmc,file="HCP_two_DMC.csv",row.names = FALSE)

#######################################################################
###### individual-level estimation: priors informed by ABCD sample ####
#######################################################################

### load data ###

hcp.0.dmc<-read.csv("HCP_zero_DMC.csv")
hcp.2.dmc<-read.csv("HCP_two_DMC.csv")

hcp.0.dmc$s<-as.factor(hcp.0.dmc$s)
hcp.0.dmc$S<-as.factor(hcp.0.dmc$S)
hcp.0.dmc$R<-as.factor(hcp.0.dmc$R)

hcp.2.dmc$s<-as.factor(hcp.2.dmc$s)
hcp.2.dmc$S<-as.factor(hcp.2.dmc$S)
hcp.2.dmc$R<-as.factor(hcp.2.dmc$R)

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

# increase scale 50% because these priors were from 9-10 year olds:
p2.zero<-p2.zero*1.5

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

# increase scale 50% because these priors were from 9-10 year olds:
p2.two<-p2.two*1.5

p.prior.two <- prior.p.dmc(
  dists = rep("tnorm",8),
  p1=p1.two,                           
  p2=p2.two,
  lower=c(0,NA,NA,NA,0,0,0,NA),upper=c(NA,NA,NA,NA,1,2,2,NA)
)
# par(mfcol=c(1,3)); for (i in names(p.prior.two)) plot.prior(i,p.prior.two)

zero.hcp.mdi <- data.model.dmc(hcp.0.dmc, model)

two.hcp.mdi <- data.model.dmc(hcp.2.dmc, model)

sZero.hcp  <- h.samples.dmc(nmc = 120, p.prior=p.prior.zero,data = zero.hcp.mdi, thin = 5)
sTwo.hcp  <- h.samples.dmc(nmc = 120, p.prior=p.prior.two,data = two.hcp.mdi, thin = 5)

save(sZero.hcp,file="zero_DDM_inf_HCP.RData")
save(sTwo.hcp,file="two_DDM_inf_HCP.RData")

# sampling can be conducted efficiently by using a high-performance 
# computer to estimate parameters for each participant in parallel
# using the RUN.dmc() function (with default parameters) as described
# in the Dynamic Models of Choice tutorial: https://osf.io/pbwx8/

##############################################################
##### generate point estimates (posterior medians)############
##############################################################

load("zero_DDM_inf_HCP.RData")
load("two_DDM_inf_HCP.RData")

#### check convergence

# load("nback_ddm_inf_fmri_samples.RData")

rhat.Zero<-gelman.diag.dmc(sZero.hcp, hyper = FALSE)
rhat.Two<-gelman.diag.dmc(sTwo.hcp, hyper = FALSE)
save(rhat.Zero,rhat.Two,file="nback_hcp_inf_rhat.RData")

#### model fit

load("nback_hcp_inf_samples.RData")

# all particiapnts with included fMRI data
final<-read.csv("hcp_nback_postmedians_v4_asw.csv")

sZero.paper<-sZero.hcp[names(sim.Zero.base)%in%final$Subject]
sTwo.paper<-sTwo.hcp[names(sim.Two.base)%in%final$Subject]
save(sZero.paper,sTwo.paper,file="nback_hcp_inf_paper_samples.RData")

sim.Zero.paper <- h.post.predict.dmc(sZero.paper,cores=3)
sim.Two.paper <- h.post.predict.dmc(sTwo.paper,cores=3)

save(sim.Zero.paper,sim.Two.paper,file="nback_ddm_inf_hcp_paper_pp.RData")

plot.pp.dmc(sim.Zero.paper,layout = c(1,3),
            model.legend = F,show.fits=FALSE,mar=c(4,5,3,1),pos=NA)

plot.pp.dmc(sim.Two.paper,layout = c(1,3),
            model.legend = F,show.fits=FALSE,mar=c(4,5,3,1),pos=NA)

#### parameter estimates

# loop to make summary parameters

sZero.hcp.pars<-sZero.hcp

for (s in 1:length(sZero.hcp)){
  t<-sZero.hcp[[s]]$theta; d<-dim(t)
  Ter<-t[,"t0",]+(t[,"st0",]/2)
  gf.natural<-pnorm(t[,"gf",])
  D <- d + c(0, 2, 0)
  t2<-array(sapply(1:D[3], function(x) cbind(t[,,x], Ter[,x],
                                             gf.natural[,x])), D)
  dimnames(t2)[[2]]<-c(dimnames(t)[[2]],"Ter","gf.natural")
  sZero.hcp.pars[[s]]$theta<-t2
}

sTwo.hcp.pars<-sTwo.hcp

for (s in 1:length(sTwo.hcp)){
  t<-sTwo.hcp[[s]]$theta; d<-dim(t)
  Ter<-t[,"t0",]+(t[,"st0",]/2)
  gf.natural<-pnorm(t[,"gf",])
  D <- d + c(0, 2, 0)
  t2<-array(sapply(1:D[3], function(x) cbind(t[,,x], Ter[,x],
                                             gf.natural[,x])), D)
  dimnames(t2)[[2]]<-c(dimnames(t)[[2]],"Ter","gf.natural")
  sTwo.hcp.pars[[s]]$theta<-t2
}


# save out point estimates 

sZero.hcp_medians<-lapply(sZero.hcp.pars,FUN=function(x) apply(x$theta,2,median))
sZero.hcp_medians<-as.data.frame(t(as.data.frame(sZero.hcp_medians)))

sTwo.hcp_medians<-lapply(sTwo.hcp.pars,FUN=function(x) apply(x$theta,2,median))
sTwo.hcp_medians<-as.data.frame(t(as.data.frame(sTwo.hcp_medians)))

colnames(sZero.hcp_medians)<-paste0("zero.",colnames(sZero.hcp_medians))
colnames(sTwo.hcp_medians)<-paste0("two.",colnames(sTwo.hcp_medians))

hcp.medians<-cbind(sZero.hcp_medians,sTwo.hcp_medians)

hcp.medians$zero.v.mean<-rowMeans(hcp.medians[,c("zero.v.lur","zero.v.tar","zero.v.non")])
hcp.medians$two.v.mean<-rowMeans(hcp.medians[,c("two.v.lur","two.v.tar","two.v.non")])
hcp.medians$overall.v.mean<-rowMeans(hcp.medians[,c("zero.v.lur","zero.v.tar","zero.v.non",
                                                    "two.v.lur","two.v.tar","two.v.non")])

row.names(hcp.medians)<-gsub("X","",row.names(hcp.medians))

write.csv(hcp.medians,file="hcp_nback_postmedians.csv")


