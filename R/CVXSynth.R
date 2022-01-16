#' @title \pkg{CVXSynth} provides a way to run synthetic control using CVXR.
#'
#' @description \pkg{CVXSynth} provides a way to run synthetic control using CVXR.
#'
#' @name CVXSynth
#' @importFrom utils head tail
#' @importFrom CVXR  Minimize sum_squares Problem solve Variable
#' @importFrom devtools install_github
#' @docType package
NULL

# %%
#' Solve for synthetic control weights in CVXR
#' @param y_t t_0 X 1   matrix of pre-treatment outcomes for treatment units
#' @param y_c t_0 X n_0 matrix of pre-treatment outcomes for donor units
#' @param solv what solver to use. default is mosek
#' @return vector of weights
#' @import CVXR
#' @export
sc_solve = function(y_t, y_c, solv = 'MOSEK'){
  o = Variable(ncol(y_c))
  objective = Minimize(sum_squares(y_t - y_c %*% o))
  constraints = list( # no intercept
    o >= 0,            # positive weights
    sum(o) == 1        # weights sum to 1
    )
  problem = Problem(objective, constraints)
  result = solve(problem, solver = solv)
  o_hat = result$getValue(o)
  return(o_hat)
}

# %%
#' Solve for Imbens-Doudchenko elastic net synthetic control weights in CVXR
#' @param y_t t_0 X 1   matrix of pre-treatment outcomes for treatment units
#' @param y_c t_0 X n_0 matrix of pre-treatment outcomes for donor units
#' @param lambdas       vector of penalty values
#' @param alpha         scalar mixing value between L1 and L2 regularisation
#' @param t             number of lambdas to try (when range not manually specified)
#' @param solv what solver to use. default is mosek
#' @return vector of weights
#' @import CVXR
#' @export
en_sc_solve = function(y_t, y_c, lambdas = NULL, alpha = 0.5, t = 10,
    solv = 'MOSEK'){
  # sequence of lambdas
  if (is.null(lambdas)) lambdas = 10^seq(-2, log10(max(y_t)), length.out = t)
  # penalty term
  elastic_penalty = function(o, lam = 0, alph = 0) {
      lasso =  cvxr_norm(o, 1) * alph                   # l1 norm
      ridge =  cvxr_norm(o, 2) * (1 - alph) / 2         # l2 norm
      lam * (lasso + ridge)
  }
  # targets : intercept and weights
  mu = Variable(1) ; o = Variable(ncol(y_c))
  en_solve = function(l) {
    obj = (sum_squares(y_t - mu - y_c %*% o)/(2 * nrow(y_c)) +
            elastic_penalty(o, l, alpha) )
    # unconstrained problem, just apply L1 and L2 regularisation to weights
    prob = Problem(Minimize(obj))
    sol = solve(prob, solver = solv)
    # concatenate intercept and weights
    sols = c(sol$getValue(mu), sol$getValue(o))
  }
  # solve for each value of lambda
  solution_mat = lapply(lambdas, en_solve)
  # convert to (n_0 + 1) × |lambdas| matrix
  solution_mat = do.call(cbind, solution_mat)
  return(solution_mat)
}

# %% ####################################################
#' Choose regularisation parameter that minimises pseudo-treatment prediction error for
#' Doudchenko Imbens (2016) estimator
#' @param j  index of pseudo-treatment unit
#' @param Y  matrix of control outcomes (wide, where each column is a series)
#' @param lambdas sequence of lambda parameters
#' @param T0 number of pre-treatment periods
#' @return value of lambda that minimises prediction error for pseudo-treatment
#' @export

pick_lambda = function(j, Y, lambdas, T0){
  # pre period data
  ypre = head(Y, T0)
  y_j  = ypre[, j]; y_nj = ypre[, -j]
  # fit mu, o for pseudo-treatment unit
  o_tilde = en_sc_solve(y_j, y_nj, lambdas)
  # compute prediction error on post-period
  ypost = tail(Y,  (nrow(Y) - T0))
  mse_j = function(w) mean(ypost[, j] - cbind(1, ypost[, -j])  %*% w)
  # each column in o_tilde is a set of weights and intercepts for a given lam, so
  # summarise across columns
  mses = apply(o_tilde, 2, mse_j)
  # lambda that minimises error
  lambdas[which(mses == min(mses))]
}

# %%
