% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CVXSynth.R
\name{sc_solve}
\alias{sc_solve}
\title{Solve for synthetic control weights in CVXR}
\usage{
sc_solve(y_t, y_c, solv = "MOSEK")
}
\arguments{
\item{y_t}{t_0 X 1   matrix of pre-treatment outcomes for treatment units}

\item{y_c}{t_0 X n_0 matrix of pre-treatment outcomes for donor units}

\item{solv}{what solver to use. default is mosek}
}
\value{
vector of weights
}
\description{
Solve for synthetic control weights in CVXR
}
