\name{DT_augment}
\alias{DT_augment}
\docType{data}
\title{
DT_augment design example.
}
\description{
This dataset contains phenotpic data for one trait evaluated in the experimental design known as augmented design. This model allows to obtain BLUPs for genotypes that are unreplicated by dividing the field in blocks and replicating 'check genotypes' in the blocks and unreplicated genotypes randomly within the blocks. The presence of check genotypes (usually cultivars) allows the adjustment of unreplicated genotypes.

}
\usage{data("DT_augment")}
\format{
  The format is:
 chr "DT_augment"
}
\source{
This data was generated by a potato study.
}
\references{

Giovanny Covarrubias-Pazaran (2024).  lme4breeding: enabling genetic evaluation in the age of genomic data. To be submitted to Bioinformatics.

Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48.

}
\examples{

## AUGMENTED DESIGN EXAMPLE
data(DT_augment)
DT <- DT_augment
head(DT)

## fit the mixed model and check summary
mix1 <- lmebreed(TSW ~ Check.Gen + (1|Block) + (1|Genotype:Check),
             data=DT)
vc <- VarCorr(mix1); print(vc,comp=c("Variance"))


}
\keyword{datasets}
\seealso{The core functions of the package \code{\link{lmebreed}} }
