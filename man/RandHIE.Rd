\name{RandHIE}
\docType{data}
\alias{RandHIE}

\title{RAND Health Insurance Experiment}

\description{
'The RAND Health Insurance Experiment (RAND HIE) was a comprehensive study
of health care cost, utilization and outcome in the United States.
It is the only randomized study of health insurance,
and the only study which can give definitive evidence as to the causal effects
of different health insurance plans.
[...]
Although the fieldwork of the study was conducted between 1974 and 1982,
the results are still highly relevant,
since RAND HIE is the only study which can make causal statements.'
(Wikipedia, RAND Health Insurance Experiment, \url{http://en.wikipedia.org/w/index.php?title=RAND_Health_Insurance_Experiment&oldid=110166949}, accessed April 8, 2007).
}

\usage{data(RandHIE)}

\format{
   This data frame contains the following columns:
   \describe{
      \item{plan}{HIE plan number.}
      \item{site}{Participant's  place of residence
         when the participant was initially enrolled.}
      \item{coins}{Coinsurance rate.}
      \item{tookphys}{Took baseline physical.}
      \item{year}{Study year.}
      \item{zper}{Person identifier.}
      \item{black}{1 if race of household head is black.}
      \item{income}{Family income.}
      \item{xage}{Age in years.}
      \item{female}{1 if person is female.}
      \item{educdec}{Education of household head in years.}
      \item{time}{Time eligible during the year.}
      \item{outpdol}{Outpatient expenses:
         all covered outpatient medical services
         excluding dental care, outpatient psychotherapy,
         outpatient drugs or supplies.}
      \item{drugdol}{Drug expenses:
         all covered outpatient and dental drugs.}
      \item{suppdol}{Supply expenses:
         all covered outpatient supplies including dental.}
      \item{mentdol}{Psychotherapy expenses:
         all covered outpatient psychotherapy services including injections
         excluding charges for visits in excess of 52 per year,
         prescription drugs, and inpatient care.}
      \item{inpdol}{Inpatient expenses:
         all covered inpatient expenses in a hospital, mental hospital,
         or nursing home,
         excluding outpatient care and renal dialysis.}
      \item{meddol}{Medical expenses:
         all covered inpatient and outpatient services,
         including drugs, supplies, and inpatient costs of newborns
         excluding dental care and outpatient psychotherapy.}
      \item{totadm}{Hospital admissions:
         annual number of covered hospitalizations.}
      \item{inpmis}{Incomplete Hospital Records:
         missing inpatient records.}
      \item{mentvis}{Psychotherapy visits:
         indicates the annual number of outpatient visits for psychotherapy.
         It includes billed visits only.
         The limit was 52 covered visits per person per year.
         The count includes an initial visit to a psychiatrist or psychologist.}
      \item{mdvis}{Face-to-Face visits to physicians:
         annual covered outpatient visits with physician providers
         (excludes dental, psychotherapy, and
         radiology/anesthesiology/pathology-only visits).}
      \item{notmdvis}{Face-to-Face visits to nonphysicians:
         annual covered outpatient visits with nonphysician providers
         such as speech and physical therapists, chiropractors,
         podiatrists, acupuncturists, Christian Science etc.
         (excludes dental, healers, psychotherapy,
         and radiology/anesthesiology/pathology-only visits).}
      \item{num}{Family size.}
      \item{mhi}{Mental health index.}
      \item{disea}{Number of chronic diseases.}
      \item{physlm}{Physical limitations.}
      \item{ghindx}{General health index.}
      \item{mdeoff}{Maximum expenditure offer.}
      \item{pioff}{Participation incentive payment.}
      \item{child}{1 if age is less than 18 years.}
      \item{fchild}{\code{female * child}.}
      \item{lfam}{log of \code{num} (family size).}
      \item{lpi}{log of \code{pioff} (participation incentive payment).}
      \item{idp}{1 if individual deductible plan.}
      \item{logc}{\code{log(coins+1)}.}
      \item{fmde}{0 if \code{idp=1},
         \code{ln(max(1,mdeoff/(0.01*coins)))} otherwise.}
      \item{hlthg}{1 if self-rated health is good
         -- baseline is excellent self-rated health.}
      \item{hlthf}{1 if self-rated health is fair
         -- baseline is excellent self-rated health.}
      \item{hlthp}{1 if self-rated health is poor
         -- baseline is excellent self-rated health.}
      \item{xghindx}{\code{ghindx} (general healt index)
         with imputations of missing values.}
      \item{linc}{log of \code{income} (family income).}
      \item{lnum}{log of \code{num} (family size).}
      \item{lnmeddol}{log of \code{meddol} (medical expenses).}
      \item{binexp}{1 if \code{meddol} > 0.}
   }
}

\source{
   Data sets of Cameron and Trivedi (2005),
   \url{http://cameron.econ.ucdavis.edu/mmabook/mmadata.html}.

   Additional information of variables
   from Table 20.4 of Cameron and Trivedi (2005)
   and from Newhouse (1999).
}

\references{
   Cameron, A. C. and Trivedi, P. K. (2005)
   \emph{Microeconometrics: Methods and Applications},
   Cambridge University Press.

   Newhouse, J. P. (1999)
   \emph{RAND Health Insurance Experiment [in Metropolitan and Non-Metropolitan
   Areas of the {United} {States}], 1974--1982},
   ICPSR Inter-university Consortium for Political and Social Research,
   Aggregated Claims Series,
   Volume 1: Codebook for Fee-for-Service Annual Expenditures and Visit Counts,
   ICPSR 6439.

   Wikipedia, \emph{RAND Health Insurance Experiment},
   \url{ http://en.wikipedia.org/wiki/RAND_Health_Insurance_Experiment}.
}

\keyword{datasets}
