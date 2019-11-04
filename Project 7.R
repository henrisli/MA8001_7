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
  voi_j[j] = r_wj[j]*dnorm(0)
  #observed_indexes = seq(j, 600 + j, l = 25)
  #cov_mat_temp = covariance_mat[observed_indexes, observed_indexes]
  #voi = 0
  
  # y = mvrnorm(200000, mu = rep(0,25), Sigma = cov_mat_temp) + matrix(rnorm(25*200000, 0, 5), ncol = 25)
  # mu_j = covariance_mat%*%t(F_j)%*%solve(F_j%*%covariance_mat%*%t(F_j)+obs_cov_mat)%*%t(y)
  # mu_w = apply(mu_j, 2, sum)
  # voi = mu_w*pnorm(mu_w/r_wj[j]) + r_wj[j]*dnorm(mu_w/r_wj[j])
  # voi_j[j] = mean(voi)
  print(j)
}
plot(voi_j, xlab = "j", ylab = "VoI")

# b) second part
simulation_run <- function(){
  covariance_mat = correlation_mat
  x = mvrnorm(1, mu = rep(0,625), Sigma = covariance_mat)
  r_wj = rep(NA, 25)
  voi_j = rep(NA, 25)
  for (j in 1:25){
    F_j = matrix(0, ncol = 625, nrow = 25)
    for (i in 1:25){
      F_j[i,j+25*(i-1)] = 1
    }
    R_j[[j]] = covariance_mat%*%t(F_j)%*%solve(F_j%*%covariance_mat%*%t(F_j)+obs_cov_mat)%*%F_j%*%covariance_mat
    r_wj2 = sum(R_j[[j]])
    r_wj[j] = sqrt(r_wj2)
    voi_j[j] = r_wj[j]*dnorm(0)
  }
  survey_line = list()
  survey_line[[1]] = which.max(voi_j)
  r_wj = r_wj[survey_line[[1]]]
  F_j = matrix(0, ncol = 625, nrow = 25)
  for (i in 1:25){
    F_j[i,survey_line[[1]]+25*(i-1)] = 1
  }
  y = F_j%*%x + rnorm(25, 0, 5)
  mu_new = covariance_mat%*%t(F_j)%*%solve(obs_cov_mat + F_j%*%covariance_mat%*%t(F_j))%*%y
  covariance_mat = covariance_mat - R_j[[survey_line[[1]]]]
  mu_w = sum(mu_new)
  iter = 1
  while(max(0,mu_w) < mu_w*pnorm(mu_w/r_wj) + r_wj*pnorm(mu_w/r_wj) - 0.5){
    iter = iter+1
    r_wj = rep(NA, 25)
    voi_j = rep(NA, 25)
    for (j in 1:25){
      F_j = matrix(0, ncol = 625, nrow = 25)
      for (i in 1:25){
        F_j[i,j+25*(i-1)] = 1
      }
      R_j[[j]] = covariance_mat%*%t(F_j)%*%solve(F_j%*%covariance_mat%*%t(F_j)+obs_cov_mat)%*%F_j%*%covariance_mat
      r_wj2 = sum(R_j[[j]])
      r_wj[j] = sqrt(r_wj2)
      voi_j[j] = r_wj[j]*dnorm(0)
    }
    survey_line[[iter]] = which.max(voi_j)
    r_wj = r_wj[survey_line[[iter]]]
    F_j = matrix(0, ncol = 625, nrow = 25)
    for (i in 1:25){
      F_j[i,survey_line[[iter]]+25*(i-1)] = 1
    }
    y = F_j%*%x + rnorm(25, 0, 5)
    mu_new = covariance_mat%*%t(F_j)%*%solve(obs_cov_mat + F_j%*%covariance_mat%*%t(F_j))%*%y
    covariance_mat = covariance_mat - R_j[[survey_line[[iter]]]]
    mu_w = sum(mu_new)
    #print(mu_w)
    #print(iter)
  }
  return(unlist(survey_line))
}
survey_lines = list()
for (i in 1:10){
  print(i)
  survey_lines[[i]] = simulation_run()
}
#write.table(survey_lines, "C://Users//henri//Documents//GitHub//MA8001_7//survey_lines.csv")
num_iter = sapply(survey_lines, length)
df = data.frame(x = num_iter)
ggplot(data = df, aes(x=x)) + geom_histogram() + xlab("Number of survey lines") + ylab("Counts")
#hist(num_iter)



max_iterations = max(num_iter)

line_number = rep(NA, sum(num_iter))
run_number = rep(NA, sum(num_iter))
plot_matrix = matrix(0, nrow = 100, ncol = 25)
iter = 1
for(i in 1:204){
  for (j in 1:num_iter[i]){
    line_number[iter] = survey_lines[[i]][j]
    run_number[iter] = i
    plot_matrix[i, survey_lines[[i]][j]] = plot_matrix[i, survey_lines[[i]][j]] + 1
    iter = iter + 1
  }
}

plot(run_number, line_number, xlab = "Run number", ylab = "Survey Line")
image.plot(matrix(apply(plot_matrix, 2, sum), ncol = 1), yaxt = 'n', x = 1:25, y = 1, ylab = "", xlab = "Survey Line")
plot(apply(plot_matrix, 2, sum), xlab = "Survey Line", ylab = "Number of surveys")
