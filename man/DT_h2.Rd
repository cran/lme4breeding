\name{DT_h2}
\alias{DT_h2}
\docType{data}
\title{
Broad sense heritability calculation.
}
\description{
This dataset contains phenotpic data for 41 potato lines evaluated in 5 locations across 3 years in an RCBD design. The phenotypic trait is tuber quality and we show how to obtain an estimate of DT_h2 for the trait.

}
\usage{data("DT_h2")}
\format{
  The format is:
 chr "DT_h2"
}
\source{
This data was generated by a potato study.
}
\references{

Giovanny Covarrubias-Pazaran (2024).  lme4breeding: enabling genetic evaluation in the age of genomic data. To be submitted to Bioinformatics.

Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48.

}
\examples{

data(DT_h2)
DT <- DT_h2
head(DT)

\donttest{

DT=DT[with(DT, order(Env)), ]

Z <- with(DT, smm(Env))
for(i in 1:ncol(Z)){DT[,colnames(Z)[i]] <- Z[,i]}
csdiagFormula <- paste0( "y ~ Env + (", paste(colnames(Z), collapse = "+"), "|| Name)")
blockFormula <- paste0( "(0+", paste(colnames(Z), collapse = "+"), "|| Block)")
form <- paste(csdiagFormula , blockFormula, sep="+")

ans1b <- lmebreed(as.formula(form),
                  data=DT)
vc <- VarCorr(ans1b); print(vc,comp=c("Variance"))

}


}
\keyword{datasets}
\seealso{The core function of the package \code{\link{lmebreed}} }
