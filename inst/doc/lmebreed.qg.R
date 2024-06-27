## -----------------------------------------------------------------------------
library(lme4breeding)
data(DT_example)
DT <- DT_example
A <- A_example

ans1 <- lmebreed(Yield~ (1|Name) + (1|Env) + 
                   (1|Env:Name) + (1|Env:Block),
             data=DT)
vc <- VarCorr(ans1); print(vc,comp=c("Variance"))
ve <- attr(VarCorr(ans1), "sc")^2
n.env <- length(levels(DT$Env))
H2=vc$Name / ( vc$Name + (vc$`Env:Name`/n.env) + (ve/(n.env*2)) )
H2

## -----------------------------------------------------------------------------
data(DT_cpdata)
DT <- DT_cpdata
GT <- GT_cpdata
MP <- MP_cpdata
#### create the variance-covariance matrix
A <- A.mat(GT) # additive relationship matrix
A <- A + diag(1e-4, ncol(A), ncol(A))
#### look at the data and fit the model
head(DT)
mix1 <- lmebreed(Yield~ (1|id) + (1|Rowf) + (1|Colf),
                 relmat=list(id=A),
                 control = lmerControl(
                   check.nobs.vs.nlev = "ignore",
                   check.nobs.vs.rankZ = "ignore",
                   check.nobs.vs.nRE="ignore"
                 ),
                 data=DT)
vc <- VarCorr(mix1); print(vc,comp=c("Variance"))
ve <- attr(VarCorr(mix1), "sc")^2
h2= vc$id / ( vc$id + ve )
as.numeric(h2)

## ----fig.show='hold'----------------------------------------------------------
data(DT_example)
DT <- DT_example
A <- A_example
head(DT)
## Compound simmetry (CS) + Diagonal (DIAG) model
Z <- with(DT, smm(Env))
csdiagFormula <- paste0( "Yield ~ Env + (", paste(colnames(Z), collapse = "+"), "|| Name)")
for(i in 1:ncol(Z)){DT[,colnames(Z)[i]] <- Z[,i]}
ansCSDG <- lmebreed(as.formula(csdiagFormula),
                    relmat = list(Name = A ),
                    data=DT)
vc <- VarCorr(ansCSDG); print(vc,comp=c("Variance"))

## -----------------------------------------------------------------------------
data(DT_cornhybrids)
DT <- DT_cornhybrids
DTi <- DTi_cornhybrids
GT <- GT_cornhybrids

modFD <- lmebreed(Yield~Location + (1|GCA1)+(1|GCA2)+(1|SCA),
              data=DT)

vc <- VarCorr(modFD); print(vc,comp=c("Variance"))
Vgca <- vc$GCA1 + vc$GCA2
Vsca <- vc$SCA
Ve <- attr(vc, "sc")^2
Va = 4*Vgca
Vd = 4*Vsca
Vg <- Va + Vd
(H2 <- Vg / (Vg + (Ve)) )
(h2 <- Va / (Vg + (Ve)) )

## -----------------------------------------------------------------------------
data("DT_halfdiallel")
DT <- DT_halfdiallel
head(DT)
DT$femalef <- as.factor(DT$female)
DT$malef <- as.factor(DT$male)
DT$genof <- as.factor(DT$geno)
# overlay matrix to be added to the addmat argument
Z <- with(DT, overlay(femalef,malef) )
# create inital values for incidence matrix but irrelevant
# since these will be replaced by admat argument
fema <- (rep(colnames(Z), nrow(DT)))[1:nrow(DT)]
#### model using overlay without relationship matrix
modh <- lmebreed(sugar ~ (1|genof) + (1|fema),
                 addmat = list(fema=Z),
             data=DT)
vc <- VarCorr(modh); print(vc,comp=c("Variance"))
ve <- attr(vc, "sc")^2;ve

## -----------------------------------------------------------------------------
# data(DT_wheat)
# DT <- DT_wheat
# GT <- GT_wheat[,1:200]
# colnames(DT) <- paste0("X",1:ncol(DT))
# DT <- as.data.frame(DT);DT$line <- as.factor(rownames(DT))
# # select environment 1
# rownames(GT) <- rownames(DT)
# K <- A.mat(GT) # additive relationship matrix
# colnames(K) <- rownames(K) <- rownames(DT)
# # GBLUP pedigree-based approach
# set.seed(12345)
# y.trn <- DT
# vv <- sample(rownames(DT),round(nrow(DT)/5))
# y.trn[vv,"X1"] <- NA
# head(y.trn)
# ## GBLUP
# K <- K + diag(1e-4, ncol(K), ncol(K) )
# ans <- lmebreed(X1 ~ (1|line), 
#                 relmat = list(line=K),
#                 control = lmerControl(
#                   check.nobs.vs.nlev = "ignore",
#                   check.nobs.vs.rankZ = "ignore",
#                   check.nobs.vs.nRE="ignore"
#                 ),
#                 data=y.trn)
# vc <- VarCorr(ans); print(vc,comp=c("Variance"))
# 
# # take a extended dataset and fit a dummy model 
# # just to get required matrices
# y.tst <- y.trn; y.tst$X1 <- imputev(y.tst$X1)
# ans2 <- update(ans, 
#                start = getME(ans, "theta"),
#                data = y.tst,
#                control = lmerControl(check.nobs.vs.nlev = "ignore",
#                                      check.nobs.vs.rankZ = "ignore",
#                                      check.nobs.vs.nRE="ignore",
#                                      optCtrl = list(maxeval= 1),
#                                      calc.derivs = FALSE))
# # compute predictive ability
# cor(ranef(ans2)$line[vv,],DT[vv,"X1"], use="complete")
# # # other approach
# # mme <- getMME(ans2, vc=vc, recordsToKeep = which(!is.na(y.trn$X1)))
# # cor(mme$bu[vv,],DT[vv,"X1"], use="complete")
# 
# ## rrBLUP
# M <- tcrossprod(GT)
# xx <- with(y.trn, redmm(x=line, M=M, nPC=100, returnLam = TRUE))
# custom <- (rep(colnames(Z), nrow(DT)))[1:nrow(DT)]
# ansRRBLUP <- lmebreed(X1 ~ (1|custom),
#                       addmat = list(custom=Z),
#                       data=y.trn)
# re <- ranef(ansRRBLUP)$custom
# u = tcrossprod(xx$Lam, t(as.matrix( re[colnames(xx$Lam),] ) ))
# cor(u[vv,],DT[vv,"X1"], use="complete")

## -----------------------------------------------------------------------------

# data(DT_ige)
# DT <- DT_ige
# A_ige <- A_ige + diag(1e-4, ncol(A_ige), ncol(A_ige) )
# # Define 2 dummy variables to make a fake covariance
# # for two different random effects
# DT$fn <- DT$nn <- 1
# # Create the incidence matrix for the first random effect
# Zf <- Matrix::sparse.model.matrix( ~ focal-1, data=DT )
# colnames(Zf) <- gsub("focal","", colnames(Zf))
# # Create the incidence matrix for the second random effect
# Zn <- Matrix::sparse.model.matrix( ~ neighbour-1, data=DT )
# colnames(Zn) <- gsub("neighbour","", colnames(Zn))
# # Make inital values for incidence matrix but irrelevant
# # since these will be replaced by the addmat argument
# both <- (rep(colnames(Zf), nrow(DT)))[1:nrow(DT)]
# # Fit the model
# modIGE <- lmebreed(trait ~ block + (0+fn+nn|both),
#                    addmat = list(both=list(Zf,Zn)),
#                    relmat = list(both=A_ige),
#                    data = DT)
# vc <- VarCorr(modIGE); print(vc,comp=c("Variance"))
# blups <- ranef(modIGE)
# pairs(blups$both)
# cov2cor(vc$both)


## -----------------------------------------------------------------------------
# data(DT_technow)
# DT <- DT_technow
# Md <- (Md_technow*2) - 1
# Mf <- (Mf_technow*2) - 1
# Ad <- A.mat(Md)
# Af <- A.mat(Mf)
# Ad <- Ad + diag(1e-4, ncol(Ad), ncol(Ad))
# Af <- Af + diag(1e-4, ncol(Af), ncol(Af))
# # simulate some missing hybrids to predict
# y.trn <- DT
# vv1 <- which(!is.na(DT$GY))
# vv2 <- sample(DT[vv1,"hy"], 100)
# y.trn[which(y.trn$hy %in% vv2),"GY"] <- NA
# ans2 <- lmebreed(GY ~ (1|dent) + (1|flint),
#                  relmat = list(dent=Ad,
#                                flint=Af),
#                  data=y.trn)
# vc <- VarCorr(ans2); print(vc,comp=c("Variance"))
# 
# # take a extended dataset and fit a dummy model 
# # just to get required matrices
# y.tst <- y.trn; y.tst$GY <- imputev(y.tst$GY)
# ans2p <- update(ans2, 
#                 start = getME(ans2, "theta"),
#                 data = y.tst,
#                 control = lmerControl(check.nobs.vs.nlev = "ignore",
#                                       check.nobs.vs.rankZ = "ignore",
#                                       check.nobs.vs.nRE="ignore",
#                                       optCtrl = list(maxeval= 1),
#                                       calc.derivs = FALSE))
# 
# re <- ranef(ans2p)
# 
# Pdent <- as.matrix(re$dent[,1,drop=FALSE]) %*% Matrix(1, ncol=nrow(re$flint), nrow=1)
# Pflint <- as.matrix(re$flint[,1,drop=FALSE]) %*% Matrix(1, ncol=nrow(re$dent), nrow=1)
# P <- Pdent + t(Pflint); colnames(P) <- rownames(re$flint)
# 
# preds <- real <- numeric()
# for(iHyb in vv2){ 
#   parents <- strsplit(iHyb,":")[[1]]
#   preds[iHyb] <- P[which(rownames(P) %in% parents),which(colnames(P) %in% parents)]
#   real[iHyb] <- DT[which(DT$hy == iHyb),"GY"]
# }
# plot(preds, real)
# cor(preds, real)

## -----------------------------------------------------------------------------
data(DT_cpdata)
DT <- DT_cpdata
# add the units column
DT$units <- as.factor(1:nrow(DT))
# get spatial incidence matrix
Zs <- with(DT, tps(Row, Col))$All
rownames(Zs) <- DT$units
# reduce the matrix to its PCs
Z = with(DT, redmm(x=units, M=Zs, nPC=100))
# create dummy variable
spatial <- (rep(colnames(Z), nrow(DT)))[1:nrow(DT)]
# fit model
mix1 <- lmebreed(Yield~ (1|Rowf) + (1|Colf) + (1|spatial),
                 addmat =list(spatial=Z),
                 control = lmerControl(
                   check.nobs.vs.nlev = "ignore",
                   check.nobs.vs.rankZ = "ignore",
                   check.nobs.vs.nRE="ignore"
                 ),
                 data=DT)
vc <- VarCorr(mix1); print(vc,comp=c("Variance"))


## -----------------------------------------------------------------------------
# data(DT_cpdata)
# DT <- DT_cpdata
# GT <- GT_cpdata
# MP <- MP_cpdata
# #### create the variance-covariance matrix
# A <- A.mat(GT) # additive relationship matrix
# A <- A + diag(1e-4, ncol(A), ncol(A))
# #### look at the data and fit the model
# head(DT)
# DT2 <- stackTrait(data=DT, traits = c("Yield","color"))
# head(DT2$long)
# 
# mix1 <- lmebreed(valueS~ (0+trait|id),
#                  relmat=list(id=A),
#                  control = lmerControl(
#                    check.nobs.vs.nlev = "ignore",
#                    check.nobs.vs.rankZ = "ignore",
#                    check.nobs.vs.nRE="ignore"
#                  ),
#                  data=DT2$long)
# vc <- VarCorr(mix1); print(vc,comp=c("Variance"))

## -----------------------------------------------------------------------------
# cov2cor(vc$id)

## -----------------------------------------------------------------------------

# data("DT_cpdata")
# DT <- as.data.frame(DT_cpdata) 
# M <- GT_cpdata
# 
# ################
# # PARTITIONED GBLUP MODEL
# ################
# 
# MMT <-tcrossprod(M) ## MM' = additive relationship matrix 
# MMTinv<-solve(MMT) ## inverse
# MTMMTinv<-t(M)%*%MMTinv # M' %*% (M'M)-
# 
# mix.part <- lmebreed(color ~ (1|id),
#                      relmat = list(id=MMT),
#                      control = lmerControl(
#                        check.nobs.vs.nlev = "ignore",
#                        check.nobs.vs.rankZ = "ignore",
#                        check.nobs.vs.nRE="ignore"
#                      ),
#                      data=DT)
# 
# #convert BLUPs to marker effects me=M'(M'M)- u
# re <- ranef(mix.part)$id
# me.part<-MTMMTinv[,rownames(re)]%*%matrix(re[,1],ncol=1)
# plot(me.part)


## -----------------------------------------------------------------------------

# 
# data("DT_wheat")
# rownames(GT_wheat) <- rownames(DT_wheat)
# G <- A.mat(GT_wheat)
# Y <- data.frame(DT_wheat)
# 
# # make the decomposition
# UD<-eigen(G) # get the decomposition: G = UDU'
# U<-UD$vectors
# D<-diag(UD$values)# This will be our new 'relationship-matrix'
# rownames(D) <- colnames(D) <- rownames(G)
# X<-model.matrix(~1, data=Y) # here: only one fixed effect (intercept)
# UX<-t(U)%*%X # premultiply X and y by U'
# UY <- t(U) %*% as.matrix(Y) # multivariate
# 
# # dataset for decomposed model
# DTd<-data.frame(id = rownames(G) ,UY, UX =UX[,1])
# DTd$id<-as.character(DTd$id)
# head(DTd)
# 
# modeld <- lmebreed(X1~ UX + (1|id),
#                  relmat=list(id=D),
#                  control = lmerControl(
#                    check.nobs.vs.nlev = "ignore",
#                    check.nobs.vs.rankZ = "ignore",
#                    check.nobs.vs.nRE="ignore"
#                  ),
#                  data=DTd)
# vc <- VarCorr(modeld); print(vc,comp=c("Variance"))
# 
# # dataset for normal model
# DTn<-data.frame(id = rownames(G) , DT_wheat)
# DTn$id<-as.character(DTn$id)
# 
# modeln <- lmebreed(X1~ (1|id),
#                    relmat=list(id=G),
#                    control = lmerControl(
#                      check.nobs.vs.nlev = "ignore",
#                      check.nobs.vs.rankZ = "ignore",
#                      check.nobs.vs.nRE="ignore"
#                    ),
#                    data=DTn)
# 
# ## compare regular and transformed blups
# red <- ranef(modeld)$id
# ren <- ranef(modeln)$id
# plot(x=(solve(t(U)))%*%  red[colnames(D),],
#      y=ren[colnames(D),], 
#      xlab="UDU blup", ylab="blup")
# 



## -----------------------------------------------------------------------------


data(DT_expdesigns)
DT <- DT_expdesigns$car1
DT <- aggregate(yield~set+male+female+rep, data=DT, FUN = mean)
DT$setf <- as.factor(DT$set)
DT$repf <- as.factor(DT$rep)
DT$malef <- as.factor(DT$male)
DT$femalef <- as.factor(DT$female)
#levelplot(yield~male*female|set, data=DT, main="NC design I")
##############################
## Expected Mean Square method
##############################
mix1 <- lm(yield~ setf + setf:repf + femalef:malef:setf + malef:setf, data=DT)
MS <- anova(mix1); MS
ms1 <- MS["setf:malef","Mean Sq"]
ms2 <- MS["setf:femalef:malef","Mean Sq"]
mse <- MS["Residuals","Mean Sq"]
nrep=2
nfem=2
Vfm <- (ms2-mse)/nrep
Vm <- (ms1-ms2)/(nrep*nfem)

## Calculate Va and Vd
Va=4*Vm # assuming no inbreeding (4/(1+F))
Vd=4*(Vfm-Vm) # assuming no inbreeding(4/(1+F)^2)
Vg=c(Va,Vd); names(Vg) <- c("Va","Vd"); Vg
##############################
## REML method
##############################
mix2 <- lmebreed(yield~ setf + setf:repf +
                   (1|femalef:malef:setf) + (1|malef:setf), 
             data=DT)
vc <- VarCorr(mix2); print(vc,comp=c("Variance"))
Vfm <- vc$`femalef:malef:setf`
Vm <- vc$`malef:setf`

## Calculate Va and Vd
Va=4*Vm # assuming no inbreeding (4/(1+F))
Vd=4*(Vfm-Vm) # assuming no inbreeding(4/(1+F)^2)
Vg=c(Va,Vd); names(Vg) <- c("Va","Vd"); Vg



## -----------------------------------------------------------------------------


DT <- DT_expdesigns$car2
DT <- aggregate(yield~set+male+female+rep, data=DT, FUN = mean)
DT$setf <- as.factor(DT$set)
DT$repf <- as.factor(DT$rep)
DT$malef <- as.factor(DT$male)
DT$femalef <- as.factor(DT$female)
#levelplot(yield~male*female|set, data=DT, main="NC desing II")
head(DT)

N=with(DT,table(female, male, set))
nmale=length(which(N[1,,1] > 0))
nfemale=length(which(N[,1,1] > 0))
nrep=table(N[,,1])
nrep=as.numeric(names(nrep[which(names(nrep) !=0)]))

##############################
## Expected Mean Square method
##############################

mix1 <- lm(yield~ setf + setf:repf + 
             femalef:malef:setf + malef:setf + 
             femalef:setf, 
           data=DT)
MS <- anova(mix1); MS
ms1 <- MS["setf:malef","Mean Sq"]
ms2 <- MS["setf:femalef","Mean Sq"]
ms3 <- MS["setf:femalef:malef","Mean Sq"]
mse <- MS["Residuals","Mean Sq"]
nrep=length(unique(DT$rep))
nfem=length(unique(DT$female))
nmal=length(unique(DT$male))
Vfm <- (ms3-mse)/nrep; 
Vf <- (ms2-ms3)/(nrep*nmale); 
Vm <- (ms1-ms3)/(nrep*nfemale); 

Va=4*Vm; # assuming no inbreeding (4/(1+F))
Va=4*Vf; # assuming no inbreeding (4/(1+F))
Vd=4*(Vfm); # assuming no inbreeding(4/(1+F)^2)
Vg=c(Va,Vd); names(Vg) <- c("Va","Vd"); Vg

##############################
## REML method
##############################

mix2 <- lmebreed(yield~ setf + setf:repf +
               (1|femalef:malef:setf) + (1|malef:setf) + 
               (1|femalef:setf),
             data=DT)
vc <- VarCorr(mix2); print(vc,comp=c("Variance"))
Vfm <- vc$`femalef:malef:setf`
Vm <- vc$`malef:setf`
Vf <- vc$`femalef:setf`

Va=4*Vm; # assuming no inbreeding (4/(1+F))
Va=4*Vf; # assuming no inbreeding (4/(1+F))
Vd=4*(Vfm); # assuming no inbreeding(4/(1+F)^2)
Vg=c(Va,Vd); names(Vg) <- c("Va","Vd"); Vg



## -----------------------------------------------------------------------------
data(DT_cpdata)
DT <- DT_cpdata
GT <- GT_cpdata#[,1:200]
MP <- MP_cpdata
M<- GT
n <- nrow(DT) # to be used for degrees of freedom
k <- 1 # to be used for degrees of freedom (number of levels in fixed effects)

## -----------------------------------------------------------------------------
# ###########################
# #### GWAS by GBLUP approach
# ###########################
# MMT <-tcrossprod(M) ## MM' = additive relationship matrix 
# MMT <- MMT + diag(1e-4, ncol(MMT), ncol(MMT) )
# MMTinv<-solve( MMT ) ## inverse
# MTMMTinv<-t(M)%*%MMTinv # M' %*% (M'M)-
# 
# mix.part <- lmebreed(color ~ (1|id) + (1|Rowf) + (1|Colf),
#                      relmat = list(id=MMT),
#                      control = lmerControl(
#                        check.nobs.vs.nlev = "ignore",
#                        check.nobs.vs.rankZ = "ignore",
#                        check.nobs.vs.nRE="ignore"
#                      ),
#                      data=DT)
# vc <- VarCorr(mix.part); print(vc,comp=c("Variance"))
# mme <- getMME(object=mix.part)
# #convert BLUPs to marker effects me=M'(M'M)- u
# re <- ranef(mix.part)$id
# a.from.g<-MTMMTinv[,rownames(re)]%*%matrix(re[,1],ncol=1)
# var.g <- kronecker(MMT[rownames(re),rownames(re)],vc$id) - 
#   mme$Ci[rownames(re),rownames(re) ]
# var.a.from.g <- t(M)%*%MMTinv[,rownames(re)]%*% (var.g) %*% t(MMTinv[,rownames(re)])%*%M
# se.a.from.g <- sqrt(diag(var.a.from.g))
# t.stat.from.g <- a.from.g/se.a.from.g # t-statistic
# pvalGBLUP <- dt(t.stat.from.g,df=n-k-1) # -log10(pval)

## -----------------------------------------------------------------------------
# plot(-log(pvalGBLUP), main="GWAS by GBLUP")

