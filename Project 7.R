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


plot(mean_var, xlab = "j", ylab = "Mean Variance")


# a) second part


covariance_mat = correlation_mat
selected_lines = rep(NA, 10)
for (iter in 1:10){
  mean_var = rep(NA, 25)
  Sigma_j = list()
  for (j in 1:25){
    F_j = matrix(0, ncol = 625, nrow = 25)
    for (i in 1:25){
      F_j[i,j+25*(i-1)] = 1
    }
    Sigma_j[[j]] = covariance_mat - covariance_mat%*%t(F_j)%*%solve(F_j%*%covariance_mat%*%t(F_j)+obs_cov_mat)%*%F_j%*%covariance_mat
    mean_var[j] = mean(Sigma_j[[j]])
  }
  plot(mean_var)
  min_index = which(mean_var == min(mean_var, na.rm = T))
  if(length(min_index)>1){min_index = min_index[2]}
  covariance_mat = Sigma_j[[min_index]]
  selected_lines[iter] = min_index
  print(iter)
}

# b) first part
covariance_mat = correlation_mat
r_wj = rep(NA, 25)
R_j = list()
voi_j = rep(NA, 25)
for (j in 1:25){
  F_j = matrix(0, ncol = 625, nrow = 25)
  for (i in 1:25){
    F_j[i,j+25*(i-1)] = 1
  }
  R_j[[j]] = covariance_mat%*%t(F_j)%*%solve(F_j%*%covariance_mat%*%t(F_j)+obs_cov_mat)%*%F_j%*%covariance_mat
  r_wj2 = sum(R_j[[j]])
  r_wj[j] = sqrt(r_wj2)
  observed_indexes = seq(j, 600 + j, l = 25)
  cov_mat_temp = covariance_mat[observed_indexes, observed_indexes]
  voi = 0
  
  y = mvrnorm(200000, mu = rep(0,25), Sigma = cov_mat_temp) + matrix(rnorm(25*200000, 0, 5), ncol = 25)
  mu_j = covariance_mat%*%t(F_j)%*%solve(F_j%*%covariance_mat%*%t(F_j)+obs_cov_mat)%*%t(y)
  mu_w = apply(mu_j, 2, sum)
  voi = mu_w*pnorm(mu_w/r_wj[j]) + r_wj[j]*dnorm(mu_w/r_wj[j])
  voi_j[j] = mean(voi)
  print(j)
}
plot(voi_j, xlab = "j", ylab = "VoI")

# b) second part