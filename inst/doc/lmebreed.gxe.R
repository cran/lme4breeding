## ----setup, include=FALSE-----------------------------------------------------
# knitr::opts_chunk$set(warning = FALSE, message = FALSE) 
library(lme4breeding)

## -----------------------------------------------------------------------------
data(DT_example)
DT <- DT_example
A <- A_example

ansMain <- lmebreed(Yield ~ Env + (1|Name),
                        relmat = list(Name = A ),
                         verbose = 0L, trace=0L, data=DT)
vc <- VarCorr(ansMain); print(vc,comp=c("Variance"))


## -----------------------------------------------------------------------------

ansDG <- lmebreed(Yield ~ Env + (0+Env || Name),
                      relmat = list(Name = A ),
                       verbose = 0L, trace=0L, data=DT)
vc <- VarCorr(ansDG); print(vc,comp=c("Variance"))
ve <- attr(vc, "sc")^2; ve


## -----------------------------------------------------------------------------

DT$EnvName <- paste(DT$Env, DT$Name, sep = ":")
E <- Matrix::Diagonal(length(unique(DT$Env)));
colnames(E) <- rownames(E) <- unique(DT$Env);E
EA <- Matrix::kronecker(E,A, make.dimnames = TRUE)
ansCS <- lmebreed(Yield ~ Env + (1|Name) + (1|EnvName),
                    relmat = list(Name = A, EnvName=  EA ),
                     verbose = 0L, trace=0L, data=DT)
vc <- VarCorr(ansCS); print(vc,comp=c("Variance"))
ve <- attr(vc, "sc")^2; ve

## -----------------------------------------------------------------------------
ansCSDG <- lmebreed(Yield ~ Env + (Env || Name),
                      relmat = list(Name = A ),
                       verbose = 0L, trace=0L, data=DT)
vc <- VarCorr(ansCSDG); print(vc,comp=c("Variance"))
ve <- attr(vc, "sc")^2; ve

## -----------------------------------------------------------------------------

ansUS <- lmebreed(Yield ~ Env + (0+Env | Name),
                    relmat = list(Name = A ),
                     verbose = 0L, trace=0L, data=DT)
vc <- VarCorr(ansUS); print(vc,comp=c("Variance"))
ve <- attr(vc, "sc")^2; ve


## -----------------------------------------------------------------------------

# library(orthopolynom)
# DT$EnvN <- as.numeric(as.factor(DT$Env))
# # add new columns to dataset
# Z <- with(DT, smm(leg(EnvN,1)) )
# for(i in 1:ncol(Z)){DT[,colnames(Z)[i]] <- Z[,i]}
# ansRR <- lmebreed(Yield ~ Env + (0+leg0+leg1| Name),
#                   relmat = list(Name = A ),
#                   verbose = 0L, trace=0L,
#                   data=DT)
# vc <- VarCorr(ansRR); print(vc,comp=c("Variance"))
# ve <- attr(vc, "sc")^2; ve


## -----------------------------------------------------------------------------
# ansRR2 <- lmebreed(Yield ~ Env + (0+leg0+leg1|| Name),
#                   relmat = list(Name = A ),
#                   verbose = 0L, trace=0L,
#                   data=DT)
# vc <- VarCorr(ansRR2); print(vc,comp=c("Variance"))
# ve <- attr(vc, "sc")^2; ve


## -----------------------------------------------------------------------------

data(DT_h2)
DT <- DT_h2
## build the environmental index
ei <- aggregate(y~Env, data=DT,FUN=mean)
colnames(ei)[2] <- "envIndex"
ei$envIndex <- ei$envIndex - mean(ei$envIndex,na.rm=TRUE) # center the envIndex to have clean VCs
ei <- ei[with(ei, order(envIndex)), ]
## add the environmental index to the original dataset
DT2 <- merge(DT,ei, by="Env")
DT2 <- DT2[with(DT2, order(Name)), ]

ansFW <- lmebreed(y~ Env + (envIndex || Name),  verbose = 0L, trace=0L, data=DT2)
vc <- VarCorr(ansFW); print(vc,comp=c("Variance"))
ve <- attr(vc, "sc")^2; ve


## -----------------------------------------------------------------------------

ansFW2 <- lmebreed(y~ Env + (envIndex | Name),  verbose = 0L, trace=0L, data=DT2)
vc <- VarCorr(ansFW2); print(vc,comp=c("Variance"))
ve <- attr(vc, "sc")^2; ve


## -----------------------------------------------------------------------------

data(DT_h2)
DT <- DT_h2
DT=DT[with(DT, order(Env)), ]
head(DT)
indNames <- na.omit(unique(DT$Name))
A <- diag(length(indNames))
rownames(A) <- colnames(A) <- indNames

# fit diagonal model first to produce a two-way table of genotype by env BLUPs
Z <- with(DT, smm(Env))
diagFormula <- paste0( "y ~ Env + (0+", paste(colnames(Z), collapse = "+"), "|| Name)")
for(i in 1:ncol(Z)){DT[,colnames(Z)[i]] <- Z[,i]}
ans1a <- lmebreed(as.formula(diagFormula),
                  relmat = list(Name = A ),
                  verbose = 0L, trace=0L, 
                  data=DT)
vc <- VarCorr(ans1a); 
H0 <- ranef(ans1a)$Name # GxE table

# Use the GxE table to obtain loadings and build loadings columns in dataset
nPC=3
Z <- with(DT,  rrm(Env, H = H0, nPC = nPC) ) 
for(i in 1:ncol(Z)){DT[,colnames(Z)[i]] <- Z[,i]}
# fit the FA model (rr + diag) to calculate factor scores
ansFA <- lmebreed(y ~ Env + (0+PC1+PC2+PC3|Name) + (0+Env|| Name),
                  relmat = list(Name = A ),
                  verbose = 0L, trace=0L, data=DT)

u <- ranef(ansFA)$Name # all BLUPs
scores <- as.matrix(u[,1:nPC])  # extract factor scores
loadings=with(DT, rrm(Env, nPC = nPC, H = H0, # lodings (latent covars)
                      returnGamma = TRUE) )$Gamma

vc <- VarCorr(ansFA); print(vc,comp=c("Variance")) # extract all varcomps
vcFA <- vc[[1]] # G of PCs only
vcDG <- diag( unlist(lapply(vc[2:16], function(x){x[[1]]})) ) # G of diag model
vcUS <- loadings %*% vcFA %*% t(loadings) # G of unstructured model
G <- vcUS + vcDG # total G as the sum of both Gs
colfunc <- colorRampPalette(c("steelblue4","springgreen","yellow"))
hv <- heatmap(cov2cor(G), col = colfunc(100), symm = TRUE)

uFA <- scores %*% t(loadings) # recover UNS effects
uDG <- as.matrix(u[,(nPC+1):ncol(u)]) # recover DIAG effects 
u <- uFA + uDG # total effects


## -----------------------------------------------------------------------------

##########
## stage 1
##########
data(DT_h2)
DT <- DT_h2
head(DT)
envs <- unique(DT$Env)
vals <- list()
for(i in 1:length(envs)){
  ans1 <- lmebreed(y~Name + (1|Block), trace = 0L, verbose = 0L,
                   data= droplevels(DT[which(DT$Env == envs[i]),]) )
  b <- fixef(ans1) 
  b[2:length(b)] <- b[2:length(b)] + b[1]
  ids <- colnames(model.matrix(~Name-1, data=droplevels(DT[which(DT$Env == envs[i]),]) ))
  ids <- gsub("Name","",ids)
  vals[[i]] <- data.frame(Estimate=b , stdError= diag( vcov(ans1)), Effect= ids, Env= envs[i])
}
DT2 <- do.call(rbind, vals)

##########
## stage 2
##########
DT2$w <- 1/DT2$stdError
ans2 <- lmebreed(Estimate~Env + (1|Effect) + (1|Env:Effect), weights = w,
                 data=DT2,  verbose = 0L, trace=0L
                 )
vc <- VarCorr(ans2); print(vc,comp=c("Variance"))
ve <- attr(vc, "sc")^2; ve


## -----------------------------------------------------------------------------

# load data (phenotypes and genotypes)
data(DT_big)
DT = DT_big
M = apply(M_big,2,as.numeric) # change from bits to numeric
rownames(M) <- rownames(M_big)
DT[,"envf_repf"] = paste(DT[,"envf"],DT[,"repf"],sep="_")

# we will make up missing data
indsNa <- sample(unique(DT$id), 200)
nas <- which(DT$id %in% indsNa) # sample(1:nrow(DT), round(nrow(DT)*.25))
DT$value[nas]=NA
DT$value <- imputev(DT$value, by=DT$envf_repf)
DT$nas <- 0; DT$nas[nas]=1

# compute the relationship matrix
MMT <- tcrossprod(M)
m <- sum(diag(MMT))/nrow(MMT)
MMT <- MMT/m
MMT <- MMT + diag(1e-05, ncol(MMT), ncol(MMT))

## Fit the main effect model
lmod <- lmebreed(formula=value~(1|id), data = DT, 
                  relmat = list(id = MMT), 
                  verbose = 0L, trace=0L, rotation=TRUE
)

# extract GEBVs
ran0 <- ranef(lmod, condVar=FALSE)
xx=as.data.frame(as.matrix(ran0$id) %*% matrix(rep(1,50), nrow=1) )
xx$id <- rownames(xx)

# Compare against
xxl <- reshape(xx, idvar = "id", varying = list(1:(ncol(xx)-1)),
               v.names = "gebv", direction = "long", timevar = "envf_repf",
               times = unique(DT$envf_repf) )
DTc <- merge(DT, xxl, by=c("id","envf_repf"), all.x = TRUE)
with(DTc[which(DTc$nas ==1),], plot(gebv,gv))
with(DTc[which(DTc$nas ==1),], cor(gebv,gv)) # accuracy is 0.39 



## -----------------------------------------------------------------------------
## Fit the diagonal model
lmod <- lmebreed(formula=value~(0+envf_repf||id), data = DT, 
                 relmat = list(id = MMT), nIters = 1000,
                 verbose = 0L, trace=0L, rotation=TRUE
)

# Extract GEBVs
ran0 <- ranef(lmod, condVar=FALSE)
xx=as.data.frame(ran0$id)# u=H0
xx$id <- rownames(xx)

# Compare against
xxl <- reshape(xx, idvar = "id", varying = list(1:(ncol(xx)-1)),
               v.names = "gebv", direction = "long", timevar = "envf_repf",
               times = colnames(xx)[1:(ncol(xx)-1)])
DTc <- merge(DT, xxl, by=c("id","envf_repf"), all.x = TRUE)
# with(DTc[which(DTc$nas ==1),], plot(gebv,gv))
with(DTc[which(DTc$nas ==1),], cor(gebv,gv)) # accuracy is ~0.62


## -----------------------------------------------------------------------------

# lmodRes <- lmebreed(formula=value~(0+envf_repf||id) + (0+envf_repf||unitsR), data = DT,
#                  relmat = list(id = MMT), nIters=1000,
#                  verbose = 0L, trace=0L, rotation=TRUE
# )
# 
# # Extract GEBVs
# ran0 <- ranef(lmodRes, condVar=FALSE)
# xx=as.data.frame(ran0$id)# u=H0
# xx$id <- rownames(xx)
# 
# # Compare against
# xxl <- reshape(xx, idvar = "id", varying = list(1:(ncol(xx)-1)),
#                v.names = "gebv", direction = "long", timevar = "envf_repf",
#                times = colnames(xx)[1:(ncol(xx)-1)])
# DTc <- merge(DT, xxl, by=c("id","envf_repf"), all.x = TRUE)
# with(DTc[which(DTc$nas ==1),], plot(gebv,gv))
# with(DTc[which(DTc$nas ==1),], cor(gebv,gv)) # accuracy is ~0.62


## -----------------------------------------------------------------------------
# # get latent covariates and add them to the dataset
# Gammas <- with(DT,  rrm(envf_repf, H = ran0$id, nPC = 3, returnGamma = TRUE))
# Zrr <- Gammas$Zstar
# for(i in 1:ncol(Zrr)){DT[,colnames(Zrr)[i]] <- Zrr[,i]}
# # fit the model
# ansRR <- lmebreed(value~(0+PC1+PC2+PC3|id) + (0+envf_repf||id), data = DT, 
#                  relmat = list(id = MMT), nIters = 1000,
#                  verbose = 0L, trace=0L, rotation=TRUE
# )
# vc <- VarCorr(ansRR); print(vc,comp=c("Variance"))
# ve <- attr(vc, "sc")^2; ve
# 
# loadings=Gammas$Gamma
# Gint <- loadings %*% vc$id %*% t(loadings)
# Gspec <- diag( unlist(lapply(vc[2:(length(vc))], function(x){x[[1]]})) )
# G <- Gint + Gspec
# u <- ranef(ansRR, condVar=FALSE)$id
# uInter <- as.matrix(u[,1:3]) %*% t(as.matrix(loadings))
# uSpec <- as.matrix(u[,-c(1:3)])
# u <- uSpec + uInter
# 
# xx=as.data.frame(u)# u=H0
# xx$id <- rownames(xx)
# xxl <- reshape(xx, idvar = "id", varying = list(1:(ncol(xx)-1)),
#                v.names = "gebv", direction = "long", timevar = "envf_repf",
#                times = colnames(xx)[1:(ncol(xx)-1)])
# DTc <- merge(DT, xxl, by=c("id","envf_repf"), all.x = TRUE)
# with(DTc[which(DTc$nas ==1),], plot(gebv,gv))
# with(DTc[which(DTc$nas ==1),], cor(gebv,gv)) # accuracy is 0.62


