x = 1:25
y = 1:25

coord_mat = cbind(rep(x, 25), rep(y, each = 25))

correlation_func <- function(h){return(exp(-0.15*h))}

distance_mat = as.matrix(dist(coord_mat))

correlation_mat = correlation_func(distance_mat)

covariance_mat = correlation_mat

obs_cov_mat = 25*diag(25)

mean_var = rep(NA, 25)
for (j in 1:25){
  F_j = matrix(0, ncol = 625, nrow = 25)
  for (i in 1:25){
    F_j[i,j+25*(i-1)] = 1
  }
  Sigma_j = covariance_mat - covariance_mat%*%t(F_j)%*%solve(F_j%*%covariance_mat%*%t(F_j)+obs_cov_mat)%*%F_j%*%covariance_mat
  mean_var[j] = mean(Sigma_j)
}
