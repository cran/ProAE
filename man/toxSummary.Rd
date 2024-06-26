% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toxSummary.R
\name{toxSummary}
\alias{toxSummary}
\title{Create patient-level and group-level summary statistics.}
\usage{
toxSummary(
  dsn,
  id_var,
  cycle_var,
  summary_measure,
  baseline_val = NA,
  arm_var = NA
)
}
\arguments{
\item{dsn}{A data.frame object with PRO-CTCAE data.}

\item{id_var}{A character string. Name of ID variable differentiating each
PRO-CTCAE survey/participant entered as a quoted string.}

\item{cycle_var}{A character string. Name of variable differentiating one
longitudinal/repeated. PRO-CTCAE survey from another, within an individual
ID.}

\item{summary_measure}{A character string. Type of summary statistic to be
used. Please consult current literature for appropriate interpretations of
the summary measure selected and suitable analysis procedures for comparing
groups. Options include: \code{"max"} = Use subjects'
maximum score. \code{"max_post_bl"} = Use subjects' maximum score
post-baseline visit. \code{"bl_adjusted"} = Use subjects' baseline adjusted
score over the study period. The baseline adjusted score is derived by the
following: If the maximum score post-baseline is more severe than the
baseline score, then the use maximum score post-baseline is used as the
adjusted score. Otherwise, if the maximum score post-baseline is the same
or less serve than the baseline score, then zero (0) is used as the
adjusted score. \code{"toxicity_index"} = Construct patient-level toxicity
index. \code{"AUC_worsening"} = Calculate group-level AUC describing}

\item{baseline_val}{A number indicating the expected baseline cycle/time
point.}

\item{arm_var}{A character string. Name of arm variable differentiating
treatment arms or other grouping factor. Required for group-level
summary measures.}
}
\value{
A data.frame with only the id and PRO-CTCAE variables being summarized.
  Each subject will now only have 1 observation (PRO-CTCAE variables are now the summary measure value).
}
\description{
Data format should be in 'long' format, where each PRO-CTCAE item is a
	  variable/column.
}
\examples{
toxSummary(dsn=ProAE::tox_acute,
id_var="id",
cycle_var="Cycle",
baseline_val=1,
summary_measure = "max")
}
