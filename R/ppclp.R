#' 2D prbabilistic principal curve with length penalty
#'
#' This function applies the probabilistic principal curve algorithm with length penalty
#'
#' @param x The coordinate of the x axis
#' @param y The coordinate of the y axis
#' @param x_fix The cordinate of the starting point of the curve
#' @param y_fix The coordinate of the ending point of the curve
#' @param K the number of knots of the curve
#' @param degree_free the degree of freedom of the B spline
#' @param lambda The magnitude value added to the length penalty
#' @param T The total number of iterations to do the EM algorithm
#' @import tidyverse splines splines2
#' @export ppclp2D
#' @examples
#' data("threeExample)
#' tmpCurve = ppclp2D(threeExample$x, threeExample$y, threeExample$x_fix, threeExample$y_fix)
#' plot(threeExample$x, threeExample$y, xlim = c(0,1), ylim = c(0,1), pch = 16, cex = 0.8)
#' lines(tmpCurve,$xFit tmpCurve$yFit, type = "l", col = "red", lwd = 5)



ppclp2D <- function(x, y, x_fix, y_fix, K = 50, degree_free = 10, lambda = 0.2, T = 100){
  ## Total number of points in the curve
  N = length(x)

  ## Basis for the spline
  w <- seq(0, 2 * pi, by = 2 * pi / K)[1 : K]
  B = cbind(1, bs(w, df = degree_free))
  B_der = cbind(0, dbs(w, df = degree_free))
  B_der_der = cbind(0, dbs(w, df = degree_free, derivs = 2))
  B_tilde = B[c(1, nrow(B)), ]


  ## PI is the matrix of p_ik
  PI = matrix(1 / K, nrow = N, ncol = K)
  PI_sum_old = rep(1 / K, K)

  ## sigma is the variance of the noise in the gaussian distribution
  sigma_old = 1

  ## beta_x and beta_y are the coefficients for the splines to the x-axis and y-axis
  beta_x_old = runif(degree_free + 1, -5, 5)
  beta_y_old = runif(degree_free + 1, -5, 5)

  likelihood_store = c()

  length_penalty = t(B_der) %*% B_der / K * 2 * pi
  smooth_penalty = t(B_der_der) %*% B_der_der / K * 2 * pi

  ## The procedure of the EM-Algorithm
  for(t in 1 : T){

    ## items used during the EM procedure
    x.i.matrix = matrix(x,nrow=length(x),ncol=K,byrow=FALSE)
    x.k.matrix = matrix(B %*% beta_x_old, nrow = N, ncol = length(B %*% beta_x_old), byrow = TRUE)
    y.i.matrix = matrix(y,nrow = length(y), ncol = K, byrow = FALSE)
    y.k.matrix = matrix(B %*% beta_y_old, nrow = N, ncol = length(B %*% beta_y_old), byrow = TRUE)

    ## E-step
    PI = exp(-1 / as.numeric((2 * sigma_old)) * ((x.i.matrix - x.k.matrix) ^ 2 + (y.i.matrix - y.k.matrix)^ 2)) %*% diag(PI_sum_old)
    PI = PI / apply(PI, 1, sum)

    ## M-step
    ## Update PI_sum
    PI_sum_new = 1 / N * apply(PI, 2, sum)

    ## Update sigma
    sigma_temp = 0
    sigma_temp = sum(((x.i.matrix - x.k.matrix)^2 + (y.i.matrix - y.k.matrix)^2 ) * PI)
    sigma_new =  sigma_temp / (2 * N)

    ## Update beta_x and beta_y
    B_XX = 0
    B_YY = 0
    B_ZZ = 0
    B_XY = 0

    for(i in 1 : N){
      B_XX = B_XX + t(B) %*% as.matrix(PI[i, ]) * x[i]
      B_YY = B_YY + t(B) %*% as.matrix(PI[i, ]) * y[i]
    }

    diag_B = apply(PI, 2, sum) %>% diag
    B_XY = t(B) %*% diag_B %*% B

    ## Inverse matrix for the estimation
    Inverse_M = solve(B_XY + lambda1 * length_penalty )

    beta_x_new  = Inverse_M %*% B_XX
    beta_y_new  = Inverse_M %*% B_YY

    ## Psu-inverse matrix for the estimation of coefficient

    Inverse_P = Inverse_M %*% t(B_tilde) %*%
      solve(B_tilde %*% Inverse_M %*% t(B_tilde))

    beta_x_new = beta_x_new -
      Inverse_P %*%
      (B_tilde %*% beta_x_new - x_fix)

    beta_y_new = beta_y_new -
      Inverse_P %*%
      (B_tilde %*% beta_y_new - y_fix)

    ## Computation of the log likelihood
    likelihood = 0
    for(i in 1 : N){
      likelihood_temp = 0
      for(k in 1 : K){
        likelihood_temp = likelihood_temp + PI_sum_new[k] * 1 / sigma_new * exp(-1/(2 * sigma_new) * ((x[i] - B[k, ] %*% beta_x_new )^2 + (y[i] - B[k, ] %*% beta_y_new)^2))
      }
      likelihood = likelihood +  log(likelihood_temp)
    }

    PI_sum_old = PI_sum_new
    sigma_old = sigma_new
    beta_x_old = beta_x_new
    beta_y_old = beta_y_new

    #print(likelihood)
    likelihood_store = c(likelihood_store, likelihood)

  }

  plot(x,y, xlim = c(0,1), ylim = c(0,1), pch = 16, cex = 0.8)
  lines(B %*% beta_x_new, B %*% beta_y_new, type = "l", col = "red", lwd = 5)
  dev.off()

  return(list(xFit = B %*% beta_x_new, yFit = B %*% beta_y_new))
}
