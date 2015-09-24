\encoding{utf8}
\docType{data}
\name{saheart}
\alias{saheart}
\title{South African Hearth data set}
\format{A data matrix 462 rows and 10 columns.}
\description{
Data are available at http://sci2s.ugr.es/keel/dataset.php?cod=184. This data set
is a retrospective sample of males in a heart-disease high-risk region of the Western Cape,
South Africa. Many of the coronary heart disease (CHD) positive men have undergone
blood pressure reduction treatment and other programs to reduce their risk factors after
their CHD event. In some cases the measurements were made after these treatments. The
class label indicates if the person has a coronary heart disease (negative or positive) and
is hidden for our analysis. Individuals are described by the following nine variables. The
continuous variables are systolic blood pressure (sbp), cumulative tobacco (tobacco), low
density lipoprotein cholesterol (ldl), adiposity, obesity and current alcohol consumption
(alcohol). The integer variables are type-A behavior (typea) and age at onset (age).
Finally, the binary variable indicates the presence or not of heart disease in the family
history (famhist).
}
\examples{
data("saheart")
summary(saheart)
}
\references{
Rousseauw et al, 1983, South African Medical Journal. 
}
\keyword{datasets}

