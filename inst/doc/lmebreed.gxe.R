## ----setup, include=FALSE-----------------------------------------------------
# knitr::opts_chunk$set(warning = FALSE, message = FALSE) 
library(lme4breeding)

## -----------------------------------------------------------------------------
data(DT_example)
DT <- DT_example
A <- A_example

ansMain <- lmebreed(Yield ~ Env + (1|Name),
                        relmat = list(Name = A ),
                        trace = 0L, data=DT)
vc <- VarCorr(ansMain); print(vc,comp=c("Variance"))


## -----------------------------------------------------------------------------

Z <- with(DT, smm(Env))
diagFormula <- paste0( "Yield ~ Env + (0+", paste(colnames(Z), collapse = "+"), "|| Name)")
for(i in 1:ncol(Z)){DT[,colnames(Z)[i]] <- Z[,i]}
print(as.formula(diagFormula))
ansDG <- lmebreed(as.formula(diagFormula),
                      relmat = list(Name = A ),
                      trace = FALSE, data=DT)
vc <- VarCorr(ansDG); print(vc,comp=c("Variance"))
ve <- attr(vc, "sc")^2; ve


## -----------------------------------------------------------------------------

DT$EnvName <- paste(DT$Env, DT$Name, sep = ":")
E <- Matrix::Diagonal(length(unique(DT$Env)));
colnames(E) <- rownames(E) <- unique(DT$Env);E
EA <- Matrix::kronecker(E,A, make.dimnames = TRUE)
ansCS <- lmebreed(Yield ~ Env + (1|Name) + (1|EnvName),
                    relmat = list(Name = A, EnvName=  EA ),
                    trace = FALSE, data=DT)
vc <- VarCorr(ansCS); print(vc,comp=c("Variance"))
ve <- attr(vc, "sc")^2; ve

## -----------------------------------------------------------------------------
Z <- with(DT, smm(Env))
csdiagFormula <- paste0( "Yield ~ Env + (", paste(colnames(Z), collapse = "+"), "|| Name)")
for(i in 1:ncol(Z)){DT[,colnames(Z)[i]] <- Z[,i]}
print(as.formula(csdiagFormula))
ansCSDG <- lmebreed(as.formula(csdiagFormula),
                      relmat = list(Name = A ),
                      trace = FALSE, data=DT)
vc <- VarCorr(ansCSDG); print(vc,comp=c("Variance"))
ve <- attr(vc, "sc")^2; ve

## -----------------------------------------------------------------------------

Z <- with(DT, smm(Env))
usFormula <- paste0( "Yield ~ Env + (0+", paste(colnames(Z), collapse = "+"), "| Name)")
for(i in 1:ncol(Z)){DT[,colnames(Z)[i]] <- Z[,i]}
print(as.formula(usFormula))
ansDG <- lmebreed(as.formula(usFormula),
                    relmat = list(Name = A ),
                    trace = FALSE, data=DT)
vc <- VarCorr(ansDG); print(vc,comp=c("Variance"))
ve <- attr(vc, "sc")^2; ve


## -----------------------------------------------------------------------------

# library(orthopolynom)
# DT$EnvN <- as.numeric(as.factor(DT$Env))
#  
# Z <- with(DT, smm(leg(EnvN,1)) )
# rrFormula <- paste0( "Yield ~ Env + (0+", paste(colnames(Z), collapse = "+"), "| Name)")
# for(i in 1:ncol(Z)){DT[,colnames(Z)[i]] <- Z[,i]}
# ansRR <- lmebreed(as.formula(rrFormula),
#                   relmat = list(Name = A ),
#                   data=DT)
# vc <- VarCorr(ansRR); print(vc,comp=c("Variance"))
# ve <- attr(vc, "sc")^2; ve


## -----------------------------------------------------------------------------
# library(orthopolynom)
# DT$EnvN <- as.numeric(as.factor(DT$Env))
#  Z <- with(DT, smm(leg(EnvN,1)) )
# rrFormula <- paste0( "Yield ~ Env + (0+", paste(colnames(Z), collapse = "+"), "|| Name)")
# for(i in 1:ncol(Z)){DT[,colnames(Z)[i]] <- Z[,i]}
# ansRR <- lmebreed(as.formula(rrFormula),
#                   relmat = list(Name = A ),
#                   data=DT)
# vc <- VarCorr(ansRR); print(vc,comp=c("Variance"))
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

ansFW <- lmebreed(y~ Env + (envIndex || Name), trace = FALSE, data=DT2)
vc <- VarCorr(ansFW); print(vc,comp=c("Variance"))
ve <- attr(vc, "sc")^2; ve


## -----------------------------------------------------------------------------

ansFW2 <- lmebreed(y~ Env + (envIndex | Name), trace = FALSE, data=DT2)
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
                  trace = FALSE,data=DT)
vc <- VarCorr(ans1a); 
H0 <- ranef(ans1a)$Name # GxE table

# Use the GxE table to obtain loadings and build loadings columns
nPC=3
Z <- with(DT,  rrm(Env, H = H0, nPC = nPC) ) 
Zd <- with(DT, smm(Env) )
faFormula <- paste0( "y ~ Env + (0+", paste(colnames(Z), collapse = "+"),
                     "| Name) + (0+",paste(colnames(Zd), collapse = "+"), "|| Name)")
for(i in 1:ncol(Z)){DT[,colnames(Z)[i]] <- Z[,i]}
# fit the FA model (rr + diag) to calculate factor scores
ansFA <- lmebreed(as.formula(faFormula),
                  relmat = list(Name = A ),
                  trace = FALSE, data=DT)

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
  ans1 <- lmebreed(y~Name + (1|Block), trace = FALSE, 
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
                 data=DT2, trace = FALSE
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

# compute the relationship matrix
MMT <- tcrossprod(M)
m <- sum(diag(MMT))/nrow(MMT)
MMT <- MMT/m
MMT <- MMT + diag(1e-05, ncol(MMT), ncol(MMT))

## Fit the main effect model
lmod <- lmebreed(formula=value~(1|id), data = DT, 
                  relmat = list(id = MMT), 
                  rotation=TRUE, 
                  trace=0L
)

# extract GEBVs
ran0 <- ranef(lmod)
xx=as.data.frame(as.matrix(ran0$id) %*% matrix(rep(1,50), nrow=1) )
xx$id <- rownames(xx)

# Compare against
xxl <- reshape(xx, idvar = "id", varying = list(1:(ncol(xx)-1)),
               v.names = "gebv", direction = "long", timevar = "envf_repf",
               times = unique(DT$envf_repf) )
DTc <- merge(DT, xxl, by=c("id","envf_repf"), all.x = TRUE)
with(DTc, plot(gebv,gv))
with(DTc, cor(gebv,gv)) # accuracy is 0.39 



## -----------------------------------------------------------------------------

# Add dummy variables to the dataset to estimate effects per environment
Z <- with(DT, smm(envf_repf))
for(i in 1:ncol(Z)){DT[,colnames(Z)[i]] <- Z[,i]}
csdiagFormula <- paste0( "value ~ (1|envf) + (0+", paste(colnames(Z), collapse = "+"), "|| id)")

## Fit the model
lmod <- lmebreed(formula=as.formula(csdiagFormula), data = DT, 
                  relmat = list(id = MMT), 
                  control = lmerControl( # how to control n iterations
                   calc.derivs = FALSE,
                   restart_edge = FALSE,
                   optCtrl = list(maxfun = 500, maxeval = 500)
                  ),
                  rotation=TRUE, #returnParams = T,
                  trace=0L
)

# Extract GEBVs
ran0 <- ranef(lmod)
xx=as.data.frame(ran0$id)# u=H0
xx$id <- rownames(xx)

# Compare against
xxl <- reshape(xx, idvar = "id", varying = list(1:(ncol(xx)-1)),
               v.names = "gebv", direction = "long", timevar = "envf_repf",
               times = colnames(xx)[1:(ncol(xx)-1)])
DTc <- merge(DT, xxl, by=c("id","envf_repf"), all.x = TRUE)
with(DTc, plot(gebv,gv))
with(DTc, cor(gebv,gv)) # accuracy is 0.62


## -----------------------------------------------------------------------------

# Create a custom matrix for residuals at each environment
DT$units <- NA
for(i in unique(DT$envf_repf)){
  sam <- which(DT$envf_repf == i)
  DT$units[sam] <- 1:length(sam)
}
DT$units <- as.factor(DT$units)
Zu <- sparse.model.matrix(~units-1, data=DT)
DT$res <- (rep(colnames(Zu), nrow(DT)))[1:nrow(DT)]
residualFormula <- paste0( "(0+", paste(colnames(Z), collapse = "+"), "|| res)")
csdiagFormula <- paste(csdiagFormula , residualFormula, sep="+")

# lmodRes <- lmebreed(formula=as.formula(csdiagFormula), data = DT,
#                  relmat = list(id = MMT),
#                  addmat = list(res=Zu),
#                  control = lmerControl( # how to control n iterations
#                    calc.derivs = FALSE,
#                    restart_edge = FALSE,
#                    optCtrl = list(maxfun = 500, maxeval = 500)
#                  ),
#                  rotation=TRUE,
#                  trace=0L
# )


