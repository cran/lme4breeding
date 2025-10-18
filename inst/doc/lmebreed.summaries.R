## ----setup, include=FALSE-----------------------------------------------------
# knitr::opts_chunk$set(warning = FALSE, message = FALSE) 
library(lme4breeding)

## -----------------------------------------------------------------------------
data(DT_example)
DT <- DT_example
A <- A_example

ans1 <- lmebreed(Yield~ (1|Name) + (1|Env) + 
                   (1|Env:Name) + (1|Env:Block),
             verbose = FALSE, data=DT)

## -----------------------------------------------------------------------------
vc <- VarCorr(ans1); print(vc,comp=c("Variance"))
ve <- attr(VarCorr(ans1), "sc")^2; ve

## -----------------------------------------------------------------------------
BLUP <- ranef(ans1, condVar=TRUE)
PEV <- lapply(BLUP, function(x){attr(x, which="postVar")}) # take sqrt() for SEs
head(BLUP$Name)
head(PEV$Name)

## -----------------------------------------------------------------------------
Ci <- getCi(ans1)
Matrix::image(Ci)

## -----------------------------------------------------------------------------
pp <- predict(ans1, classify="Name")
head(pp$pvals)

## -----------------------------------------------------------------------------
image(pp$D)

## -----------------------------------------------------------------------------
pp$hyperTable

