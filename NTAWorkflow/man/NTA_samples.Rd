\name{NTA_samples}
\alias{NTA_samples}
\docType{data}
\title{A subset of data after preprocessing with MSDial}
\description{The subset contains a small number of columns and LC-MS samples.}
\usage{data("NTA_samples")}
\format{
  A data frame with 10 observations on the following 17 variables.
  \describe{
    \item{\code{Alignment.ID}}{Increasing number for features}
    \item{\code{Average.Rt.min.}}{Retention time}
    \item{\code{Average.Mz}}{m/z value}
    \item{\code{Metabolite.name}}{Name found in database used in MSDial}
    \item{\code{MS.MS.assigned}}{logical}
    \item{\code{RT.matched}}{logical. True if retention time matched with database used in MSDial}
    \item{\code{m.z.matched}}{logical. True if m.z. matched with database used in MSDial}
    \item{\code{MS.MS.matched}}{logical. True if MS/MS matched with database used in MSDial}
    \item{\code{Spectrum.reference.file.name}}{After this column the samples come}
    \item{\code{X230828_DLR_neg_extBlank1_001}}{blank sample}
    \item{\code{X230828_DLR_neg_extBlank2_001}}{blank sample}
    \item{\code{X230828_DLR_neg_Pellet_Ground_11_001}}{sample of class "Ground"}
    \item{\code{X230828_DLR_neg_Pellet_Ground_12_001}}{sample of class "Ground"}
    \item{\code{X230828_DLR_neg_Pellet_SMG_1_001}}{sample of class "SMG"}
    \item{\code{X230828_DLR_neg_Pellet_SMG_2_001}}{sample of class "SMG"}
    \item{\code{X230828_DLR_neg_QC1_001}}{quality control}
    \item{\code{X230828_DLR_neg_QC2_001}}{quality control}
  }
}
\examples{
data(NTA_samples)
}
\keyword{datasets}
