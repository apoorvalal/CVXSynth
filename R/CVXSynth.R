# %%
#' Solve for synthetic control weights in CVXR
#' @param y_t t_0 X 1   matrix of pre-treatment outcomes for treatment units
#' @param y_c t_0 X n_0 matrix of pre-treatment outcomes for donor units
#' @return vector of weights
#' @import CVXR
#' @export
sc_solve = function(y_t, y_c){
  require(CVXR)
  ω = Variable(ncol(y_c))
  objective = Minimize(sum_squares(y_t - y_c %*% ω))
  constraints = list( # no intercept
    ω >= 0,            # positive weights
    sum(ω) == 1        # weights sum to 1
    )
  problem = Problem(objective, constraints)
  result = solve(problem, solver = 'MOSEK')
  ω_hat = result$getValue(ω)
  return(ω_hat)
}

# %%
#' Solve for Imbens-Doudchenko elastic net synthetic control weights in CVXR
#' @param y_t t_0 X 1   matrix of pre-treatment outcomes for treatment units
#' @param y_c t_0 X n_0 matrix of pre-treatment outcomes for donor units
#' @param y_c t_0 X n_0 matrix of pre-treatment outcomes for donor units
#' @param lambdas       vector of penalty values
#' @param alpha         scalar mixing value between L1 and L2 regularisation
#' @param t             number of lambdas to try (when range not manually specified)
#' @return vector of weights
#' @import CVXR
#' @export
en_sc_solve = function(y_t, y_c, lambdas = NULL, alpha = 0.5, t = 10){
  # sequence of lambdas
  if (is.null(lambdas)) lambdas = 10^seq(-2, log10(max(y_t)), length.out = t)
  # penalty term
  elastic_penalty = function(ω, λ = 0, α = 0) {
      lasso =  cvxr_norm(ω, 1) * α
      ridge =  cvxr_norm(ω, 2) * (1 - α) / 2
      λ * (lasso + ridge)
  }
  # targets : intercept and weights
  require(CVXR)
  μ = Variable(1) ; ω = Variable(ncol(y_c))
  en_solve = function(l) {
    obj = (sum_squares(y_t - μ - y_c %*% ω)/(2 * nrow(y_c)) +
            elastic_penalty(ω, l, alpha) )
    # unconstrained problem, just apply L1 and L2 regularisation to weights
    prob = Problem(Minimize(obj))
    sol = solve(prob, solver = 'MOSEK')
    # concatenate intercept and weights
    sols = c(sol$getValue(μ), sol$getValue(ω))
  }
  # solve for each value of lambda
  solution_mat = lapply(lambdas, en_solve)
  # convert to (n_0 + 1) × |lambdas| matrix
  solution_mat = do.call(cbind, solution_mat)
  return(solution_mat)
}

# %% ####################################################
#' Choose regularisation parameter that minimises pseudo-treatment prediction error
#' @param j  index of pseudo-treatment unit
#' @param Y  matrix of control outcomes (wide, where each column is a series)
#' @param lambdas sequence of lambda parameters
#' @param T0 number of pre-treatment periods
#' @return value of lambda that minimises prediction error
#' @export
pick_lambda = function(j, Y, lambdas, T0){
  # pre period data
  ypre = head(Y, T0)
  y_j  = ypre[, j]; y_nj = ypre[, -j]
  # fit μ, ω for pseudo-treatment unit
  ω_tilde = en_sc_solve(y_j, y_nj, lambdas)
  # compute prediction error on post-period
  ypost = tail(Y, T1)
  mse_j = function(w) mean(ypost[, j] - cbind(1, ypost[, -j])  %*% w)
  # each column in ω_tilde is a set of weights and intercepts for a given λ, so
  # summarise across columns
  mses = apply(ω_tilde, 2, mse_j)
  # lambda that minimises error
  lambdas[which(mses == min(mses))]
}

# %%
