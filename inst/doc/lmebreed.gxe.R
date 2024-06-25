## -----------------------------------------------------------------------------
library(lme4breeding)
data(DT_example)
DT <- DT_example
A <- A_example

ansMain <- lmebreed(Yield ~ Env + (1|Name),
                        relmat = list(Name = A ),
                        data=DT)
vc <- VarCorr(ansMain); print(vc,comp=c("Variance"))


## -----------------------------------------------------------------------------

Z <- with(DT, dsc(Env))$Z
diagFormula <- paste0( "Yield ~ Env + (0+", paste(colnames(Z), collapse = "+"), "|| Name)")
for(i in 1:ncol(Z)){DT[,colnames(Z)[i]] <- Z[,i]}
print(as.formula(diagFormula))
ansDG <- lmebreed(as.formula(diagFormula),
                      relmat = list(Name = A ),
                      data=DT)
vc <- VarCorr(ansDG); print(vc,comp=c("Variance"))
ve <- attr(vc, "sc")^2; ve


## -----------------------------------------------------------------------------

DT$EnvName <- paste(DT$Env, DT$Name, sep = ":")
E <- Matrix::Diagonal(length(unique(DT$Env)));
colnames(E) <- rownames(E) <- unique(DT$Env);E
EA <- Matrix::kronecker(E,A, make.dimnames = TRUE)
ansCS <- lmebreed(Yield ~ Env + (1|Name) + (1|EnvName),
                    relmat = list(Name = A, EnvName=  EA ),
                    data=DT)
vc <- VarCorr(ansCS); print(vc,comp=c("Variance"))
ve <- attr(vc, "sc")^2; ve

## -----------------------------------------------------------------------------
Z <- with(DT, dsc(Env))$Z
csdiagFormula <- paste0( "Yield ~ Env + (", paste(colnames(Z), collapse = "+"), "|| Name)")
for(i in 1:ncol(Z)){DT[,colnames(Z)[i]] <- Z[,i]}
print(as.formula(csdiagFormula))
ansCSDG <- lmebreed(as.formula(csdiagFormula),
                      relmat = list(Name = A ),
                      data=DT)
vc <- VarCorr(ansCSDG); print(vc,comp=c("Variance"))
ve <- attr(vc, "sc")^2; ve

## -----------------------------------------------------------------------------

Z <- with(DT, dsc(Env))$Z
usFormula <- paste0( "Yield ~ Env + (0+", paste(colnames(Z), collapse = "+"), "| Name)")
for(i in 1:ncol(Z)){DT[,colnames(Z)[i]] <- Z[,i]}
print(as.formula(usFormula))
ansDG <- lmebreed(as.formula(usFormula),
                    relmat = list(Name = A ),
                    data=DT)
vc <- VarCorr(ansDG); print(vc,comp=c("Variance"))
ve <- attr(vc, "sc")^2; ve


## -----------------------------------------------------------------------------

# library(orthopolynom)
# DT$EnvN <- as.numeric(as.factor(DT$Env))
#  
# Z <- with(DT, dsc(leg(EnvN,1)) )$Z
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
#  Z <- with(DT, dsc(leg(EnvN,1)) )$Z
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

ansFW <- lmebreed(y~ Env + (envIndex || Name), data=DT2)
vc <- VarCorr(ansFW); print(vc,comp=c("Variance"))
ve <- attr(vc, "sc")^2; ve


## -----------------------------------------------------------------------------

ansFW2 <- lmebreed(y~ Env + (envIndex | Name), data=DT2)
vc <- VarCorr(ansFW2); print(vc,comp=c("Variance"))
ve <- attr(vc, "sc")^2; ve


## -----------------------------------------------------------------------------

# data(DT_h2)
# DT <- DT_h2
# DT=DT[with(DT, order(Env)), ]
# 
# # fit diagonal model first to produce H matrix
# Z <- with(DT, dsc(Env))$Z
# diagFormula <- paste0( "y ~ Env + (0+", paste(colnames(Z), collapse = "+"), "|| Name)")
# for(i in 1:ncol(Z)){DT[,colnames(Z)[i]] <- Z[,i]}
# print(as.formula(diagFormula))
# ans1a <- lmebreed(as.formula(diagFormula),
#                   data=DT)
# vc <- VarCorr(ans1a); print(vc,comp=c("Variance"))
# H0 <- ranef(ans1a)$Name # GxE table
# 
# # reduced rank model
# Z <- with(DT,  dsc(rrc(Env, H = H0, nPC = 3)) )$Z
# Zd <- with(DT, dsc(Env))$Z
# faFormula <- paste0( "y ~ Env + (0+", paste(colnames(Z), collapse = "+"), "| Name) + (0+",paste(colnames(Zd), collapse = "+"), "|| Name)")
# for(i in 1:ncol(Z)){DT[,colnames(Z)[i]] <- Z[,i]}
# print(as.formula(faFormula))
# ansFA <- lmebreed(as.formula(faFormula),
#                   data=DT)
# vc <- VarCorr(ansFA); print(vc,comp=c("Variance"))
# ve <- attr(vc, "sc")^2; ve
# 
# loadings=with(DT, rrc(Env, nPC = 3, H = H0, returnGamma = T) )$Gamma
# Gint <- loadings %*% vc$Name %*% t(loadings)
# Gspec <- diag( unlist(lapply(vc[2:16], function(x){x[[1]]})) )
# G <- Gint + Gspec
# # lattice::levelplot(cov2cor(G))
# # colfunc <- colorRampPalette(c("steelblue4","springgreen","yellow"))
# # hv <- heatmap(cov2cor(G), col = colfunc(100), symm = TRUE)
# 
# u <- ranef(ansFA)$Name
# uInter <- as.matrix(u[,1:3]) %*% t(as.matrix(loadings))
# uSpec <- as.matrix(u[,-c(1:3)])
# u <- uSpec + uInter


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
  ans1 <- lmebreed(y~Name + (1|Block), data= droplevels(DT[which(DT$Env == envs[i]),]) )
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
ans2 <- lmebreed(Estimate~Env + (1|Effect) + (1|Env:Effect), weights = w,data=DT2, 
                 control = lmerControl(
                   # optimizer="bobyqa",
                   check.nobs.vs.nlev = "ignore",
                   check.nobs.vs.rankZ = "ignore",
                   check.nobs.vs.nRE="ignore"
                 )
                 )
vc <- VarCorr(ans2); print(vc,comp=c("Variance"))
ve <- attr(vc, "sc")^2; ve


