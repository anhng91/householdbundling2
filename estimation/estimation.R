args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) { 
  if (Sys.info()[['sysname']] == 'Windows') {
    numcores = 10;
  } else {
    numcores = 4;
  }
  job_index = 2;  
} else {
  job_index = as.numeric(args[1]);
  numcores = as.numeric(args[2]); 
}
if (dir.exists('work/teckyongtan/tecktan/Oil/Data/R_lib')) {
  .libPaths('work/teckyongtan/tecktan/Oil/Data/R_lib')
  remote = TRUE
} else {
  remote = FALSE
}
options("install.lock"=FALSE)
library(knitr)
library(tidyverse)
library(lfe)
library(sandwich)
library(plm)
library(stargazer)
library(parallel)
library(randomForest)
library(randtoolbox)
library(fastDummies)

# setwd('./familyenrollment')
devtools::install(upgrade='never')
library(familyenrollment)

message('constructing the list of Compulsory households')

Com_HH_list_index = lapply(1:length(data_hh_list), function(hh_index) {
  data = data_hh_list[[hh_index]]; 
  if (0 == nrow(data %>% filter(Bef_sts + Com_sts + Std_w_ins == 0))) {
    return(hh_index);
  }
  else {
    return(NA); 
  }
}) %>% unlist(); 

Com_HH_list_index = Com_HH_list_index[!(is.na(Com_HH_list_index))]

message('constructing the list of Voluntary households')
Vol_HH_list_index = lapply(1:length(data_hh_list), function(hh_index) {
  data = data_hh_list[[hh_index]]; 
  if (nrow(data) > nrow(data %>% filter(Bef_sts + Com_sts + Std_w_ins == 1)) & !(data$Year[1] == 2006 & data$HHsize_s[1] > 1)) {
    return(hh_index);
  }
  else {
    return(NA); 
  }
}) %>% unlist(); 

Vol_HH_list_index = Vol_HH_list_index[!(is.na(Vol_HH_list_index))] 

full_index = 1:length(data_hh_list); 

message('Identify the sample for theta')
sample_identify_theta = full_index[which((lapply(full_index, function(sample_index_i) {
      data_mini = data_hh_list[[sample_index_i]]
      if ((data_mini$Year[1] == 2008) & (data_mini$Income[1] < 0) & (data_mini$HHsize_s[1] == 0)) {
        return(1)
      } else {
        return(0)
      }
    }) %>% unlist()) == 1)]

sample_no_sick = full_index[which((lapply(full_index, function(sample_index_i) {
      data_mini = data_hh_list[[sample_index_i]]
      if (((data_mini$HHsize_s[1] == 0) & (sum(data_mini$sick_dummy == 0)))) {
        return(1)
      } else {
        return(0)
      }
    }) %>% unlist()) == 1)]

message('Identify the sample for preference')
sample_identify_pref = lapply(Com_HH_list_index, function(x) ifelse(x %in% c(sample_identify_theta, sample_no_sick), NA, x)) %>% unlist()
sample_identify_pref = sample_identify_pref[!(is.na(sample_identify_pref))]

# Bootstrapping indices 
message('bootstrapping indices')
set.seed(job_index);
sample_index = sample(1:length(data_hh_list), length(data_hh_list), replace=TRUE)
sample_r_theta = Vol_HH_list_index
if (remote) {
  sample_r_theta = sample(sample_r_theta, length(sample_r_theta), replace=TRUE)
  sample_identify_pref = sample(sample_identify_pref, length(sample_identify_pref), replace=TRUE)
  sample_identify_theta = sample(sample_identify_theta, length(sample_identify_theta), replace=TRUE)
} else {
  sample_r_theta = sample(sample_r_theta, 1000, replace=TRUE)
  sample_identify_pref = sample(sample_identify_pref, length(sample_identify_pref), replace=TRUE)
  sample_identify_theta = sample(sample_identify_theta, length(sample_identify_theta), replace=TRUE)
}


message('merging household frames')
data = do.call('rbind', data_hh_list[sample_index])

# Estimate the probability of getting sick 
message('computing the sick parameters')
sick_parameters = optim(rep(0, ncol(var_ind(data_hh_list[[1]]))), fn = function(x) llh_sick(x, data), gr = function(x) grad_llh_sick(x, data), control=list(maxit = 1e4), method='BFGS')

message('computing the coverage parameters')
xi_parameters = optim(rep(0, 2 * ncol(var_ind(data_hh_list[[1]]))), fn = function(x) llh_xi(x, data), gr = function(x) grad_llh_xi(x, data), control=list(maxit = 1e4), method='BFGS')

param_trial = transform_param(return_index=TRUE, init=TRUE);
transform_param_trial = transform_param(param_trial, return_index=TRUE)
x_transform = transform_param_trial

var_list = c('sigma_thetabar', 'beta_theta', 'beta_theta_ind');
active_index = unlist(transform_param_trial[[2]][var_list]);
active_index = setdiff(unlist(transform_param_trial[[2]][var_list]),tail(transform_param_trial[[2]][['beta_theta']], n=4));

message('Preparing data for estimation')
data_hh_list_theta = list()
for (index in 1:length(sample_identify_theta)) {
  message(paste0('constructing data for index', index))
  data_hh_list_theta[[index]] = household_draw_theta_kappa_Rdraw(hh_index=sample_identify_theta[index], param=transform_param_trial[[1]], n_draw_halton = 1000, n_draw_gauss = 10, sick_parameters, xi_parameters, short=FALSE)
}


message('Estimating theta')

if (Sys.info()[['sysname']] == 'Windows') {
  cl = makeCluster(numcores);
  clusterEvalQ(cl, library('tidyverse'))
  clusterEvalQ(cl, library('familyenrollment'))
  clusterExport(cl, c('active_index', 'param_trial', 'data_hh_list_theta'))
}

n_draw_halton = 100; 

n_halton_at_r = 100; 


if (Sys.info()[['sysname']] == 'Windows') {
  clusterExport(cl, c('active_index', 'param_trial', 'transform_param_trial', 'sample_identify_pref', 'xi_parameters', 'sick_parameters'))
  clusterExport(cl,'n_draw_halton')
}


message('compute covariates before estimation')

if (Sys.info()[['sysname']] == 'Windows') {
  data_hh_list_pref = parLapply(cl, sample_identify_pref,function(index) {
    output = tryCatch(household_draw_theta_kappa_Rdraw(hh_index=index, param=transform_param_trial[[1]], n_draw_halton = n_draw_halton, n_draw_gauss = 10, sick_parameters, xi_parameters, short=FALSE),error=function(e) e)
    return(output)
  })
} else {
  data_hh_list_pref = mclapply(sample_identify_pref, function(index) tryCatch(household_draw_theta_kappa_Rdraw(hh_index=index, param=transform_param_trial[[1]], n_draw_halton = n_draw_halton, n_draw_gauss = 10, sick_parameters, xi_parameters, short=FALSE), error=function(e) e), mc.cores=numcores)
}


mat_M = do.call('rbind', lapply(data_hh_list_pref, function(x) {
      return(cbind(x$data$M_expense, x$data$M_expense^2))
    }))

mat_M = cbind(mat_M, mat_M[,1] > quantile(mat_M[,1], 0.025) & mat_M[,1] < quantile(mat_M[,1], 0.975))

mat_Year = do.call('rbind', lapply(data_hh_list_pref, function(x) {
        return(cbind(x$data$Year == 2004, x$data$Year == 2006, x$data$Year == 2010, x$data$Year == 2012))
      }))
X_ind_pref = do.call('rbind', lapply(data_hh_list_pref, function(x) x$X_ind))
X_hh_pref = do.call('rbind', lapply(data_hh_list_pref, function(x) x$X_hh))
X_ind_pref_with_year = cbind(X_ind_pref, mat_Year)
mat_M_rtheta = do.call('rbind', lapply(data_hh_list[sample_r_theta], function(x) {
    return(cbind(x$M_expense, x$M_expense^2))
  }))

mat_Y_rtheta = do.call('c', lapply(data_hh_list[sample_r_theta], function(x) x$Income[1]))

mat_Y_rtheta = pnorm(mat_Y_rtheta, mean = mean(mat_Y_rtheta), sd = sd(mat_Y_rtheta))

full_insurance_indicator = do.call('c', lapply(data_hh_list[sample_r_theta], function(x) {
    return(x$HHsize_s[1] == x$N_vol[1])
  }))

full_insurance_indicator_ind_level = do.call('c', lapply(data_hh_list[sample_r_theta], function(x) {
    return(rep(x$HHsize_s[1] == x$N_vol[1], x$HHsize[1]))
  }))

no_insurance_indicator = do.call('c', lapply(data_hh_list[sample_r_theta], function(x) {
    return(x$N_vol[1] == 0)
  }))

no_insurance_indicator_ind_level = do.call('c', lapply(data_hh_list[sample_r_theta], function(x) {
    return(rep(x$N_vol[1] == 0, x$HHsize[1]))
  }))

X_r_all = do.call('rbind', lapply(data_hh_list[sample_r_theta], function(x) var_hh(x)))

X_hh_theta_r = do.call('rbind',lapply(sample_r_theta, function(output_hh_index) var_hh(data_hh_list[[output_hh_index]])))

n_involuntary = do.call('c', lapply(sample_r_theta, function(output_hh_index) data_hh_list[[output_hh_index]]$N_com[1] + data_hh_list[[output_hh_index]]$N_bef[1] + data_hh_list[[output_hh_index]]$N_std_w_ins[1]))

# initial_param_trial = init_param
initial_param_trial = rep(0, length(init_param))
# initial_param_trial[[x_transform[[2]]$beta_delta[1]]] = 0.1
# initial_param_trial[[x_transform[[2]]$beta_theta_ind[1]]] = 1
# initial_param_trial[[x_transform[[2]]$beta_theta[1]]] = -3
# initial_param_trial[[x_transform[[2]]$beta_omega[1]]] = 1
# initial_param_trial[[x_transform[[2]]$beta_theta_ind[1]]] = 0.5
# initial_param_trial[[x_transform[[2]]$beta_omega[1]]] = 0.5
# initial_param_trial[[x_transform[[2]]$sigma_thetabar[1]]] = 1

compute_inner_loop = function(x_stheta, return_result=FALSE, estimate_theta=TRUE, estimate_pref=TRUE) {
  print(paste0('sigma_theta value is', x_stheta))
  param_trial_here = initial_param_trial
  param_trial_here[transform_param_trial[[2]]$sigma_theta] = x_stheta; 
  transform_param_trial = transform_param(param_trial_here, return_index=TRUE);
  var_list = c('beta_omega','beta_delta', 'beta_gamma', 'sigma_delta', 'sigma_gamma', 'sigma_omega')
  x_transform = transform_param(param_trial_here, return_index=TRUE)

  aggregate_moment_pref = function(x_transform, silent=TRUE) {
    if (Sys.info()[['sysname']] == 'Windows') {
      clusterExport(cl, c('x_transform', 'n_halton_at_r'),envir=environment())
      data_hh_list_pref = parLapply(cl, sample_identify_pref,function(index) {
        output = tryCatch(household_draw_theta_kappa_Rdraw(hh_index=index, param=x_transform[[1]], n_draw_halton = n_draw_halton, n_draw_gauss = 10, sick_parameters, xi_parameters, short=FALSE),error=function(e) e)
        return(output)
      })
      mat_YK = do.call('cbind', parLapply(cl, data_hh_list_pref, function(x) {
            output = rbind(colMeans(matrix(x$kappa_draw[[1]], nrow=n_draw_halton)), colMeans(matrix(x$kappa_draw[[1]]^2, nrow=n_draw_halton)), x$income[1], x$income[1]^2, colMeans(matrix(x$kappa_draw[[1]], nrow=n_draw_halton))*x$income[1])
            return(output)
          }))
    } else {
      data_hh_list_pref = mclapply(sample_identify_pref, function(index) tryCatch(household_draw_theta_kappa_Rdraw(hh_index=index, param=x_transform[[1]], n_draw_halton = n_draw_halton, n_draw_gauss = 10, sick_parameters, xi_parameters, short=FALSE), error=function(e) e), mc.cores=numcores)
      mat_YK = do.call('cbind', mclapply(data_hh_list_pref, function(x) {
            output = rbind(colMeans(matrix(x$kappa_draw[[1]], nrow=n_draw_halton)), colMeans(matrix(x$kappa_draw[[1]]^2, nrow=n_draw_halton)), x$income[1], x$income[1]^2, colMeans(matrix(x$kappa_draw[[1]], nrow=n_draw_halton))*x$income[1])
            return(output)
          }, mc.cores=numcores))
    }

    mat_YK = rbind(mat_YK, do.call('c', lapply(data_hh_list[sample_identify_pref], function(x) x$Income)))

    if (Sys.info()[['sysname']] == 'Windows') {
      clusterExport(cl, 'x_transform',envir=environment())
      moment_ineligible_hh_output = parLapply(cl, data_hh_list_pref,function(mini_data) {
        output = tryCatch(moment_ineligible_hh(mini_data, x_transform[[1]]),error=function(e) e)
        return(output)
      })
    } else {
      moment_ineligible_hh_output = mclapply(data_hh_list_pref, function(mini_data) tryCatch(moment_ineligible_hh(mini_data, x_transform[[1]]), error=function(e) e), mc.cores=numcores)
    }

    output_1 = do.call('c', lapply(moment_ineligible_hh_output, function(x) x[[1]]))
    # if (!(is.na(sum(output_1)))) {
    #   if (max(output_1) == 0) {
    #     return(list(NA, rep(NA, length(active_index_pref))))
    #   }
    # }
    output_2 = do.call('c', lapply(moment_ineligible_hh_output, function(x) x[[2]]))
    d_output_1 = list();
    d_output_2 = list();
    for (varname in var_list) {
      if (grepl('sigma', varname)) {
        d_output_1[[varname]] = do.call('c', lapply(moment_ineligible_hh_output, function(x) x[[3]][[varname]]))
        d_output_2[[varname]] = do.call('c', lapply(moment_ineligible_hh_output, function(x) x[[4]][[varname]]))
      } else {
        d_output_1[[varname]] = do.call('cbind', lapply(moment_ineligible_hh_output, function(x) x[[3]][[paste0('mean_',varname)]]))
        d_output_2[[varname]] = do.call('cbind', lapply(moment_ineligible_hh_output, function(x) x[[4]][[paste0('mean_',varname)]]))
      }
    }

    moment = NULL
    d_moment = list()
    output = list(); 
    output[[3]] = list();
    # if (!silent) {
      print('fit of ineligible HHs')
      print(summary(output_1))
      print(summary(mat_M[,1]))
    # }
    for (moment_index in c(1:nrow(mat_YK))) {
      output[[3]][[moment_index]] = ((output_1 - mat_M[,1]) * mat_YK[moment_index,])^2
      moment[moment_index] = mean(output[[3]][[moment_index]]); 
      d_moment[[moment_index]] = 2*((output_1 - mat_M[,1]) * mat_YK[moment_index,]);
    }
    for (moment_index in c(1:nrow(mat_YK))) {
      output[[3]][[moment_index + nrow(mat_YK)]] = ((output_2 - mat_M[,2]) * mat_YK[moment_index,])^2
      moment[moment_index + nrow(mat_YK)] = mean(output[[3]][[moment_index + nrow(mat_YK)]]); 
      d_moment[[moment_index + nrow(mat_YK)]] = 2*((output_2 - mat_M[,2]) * mat_YK[moment_index,]);
    }

    d_moment = do.call('cbind', d_moment)

    output[[1]] = sum(moment); 
    active_index_pref = x_transform[[2]][var_list] %>% unlist()
    output[[2]] = rep(0, length(active_index_pref)); 

    output[[2]][x_transform[[2]][['beta_delta']]] = apply(d_output_1[['beta_delta']], 1, function(x) sum((x %*% d_moment[,1:nrow(mat_YK)])/length(x))) + apply(d_output_2[['beta_delta']], 1, function(x) sum((x %*% d_moment[,(nrow(mat_YK) + 1):(2 * nrow(mat_YK))])/length(x)))
    output[[2]][x_transform[[2]][['beta_omega']]] = apply(d_output_1[['beta_omega']], 1, function(x) sum((x %*% d_moment[,1:nrow(mat_YK)])/length(x))) + apply(d_output_2[['beta_omega']], 1, function(x) sum((x %*% d_moment[,(nrow(mat_YK) + 1):(2 * nrow(mat_YK))])/length(x)))
    output[[2]][x_transform[[2]][['beta_gamma']]] = apply(d_output_1[['beta_gamma']], 1, function(x) sum((x %*% d_moment[,1:nrow(mat_YK)])/length(x))) + apply(d_output_2[['beta_gamma']], 1, function(x) sum((x %*% d_moment[,(nrow(mat_YK) + 1):(2 * nrow(mat_YK))])/length(x)))

    output[[2]][x_transform[[2]][['sigma_delta']]] =  sum((d_output_1[['sigma_delta']] %*% d_moment[,1:nrow(mat_YK)])/length(active_index)) + sum((d_output_2[['sigma_delta']] %*% d_moment[,(nrow(mat_YK) + 1):(2 * nrow(mat_YK))])/length(active_index))
    output[[2]][x_transform[[2]][['sigma_gamma']]] =  sum((d_output_1[['sigma_gamma']] %*% d_moment[,1:nrow(mat_YK)])/length(active_index)) + sum((d_output_2[['sigma_gamma']] %*% d_moment[,(nrow(mat_YK) + 1):(2 * nrow(mat_YK))])/length(active_index))
    output[[2]][x_transform[[2]][['sigma_omega']]] =  sum((d_output_1[['sigma_omega']] %*% d_moment[,1:nrow(mat_YK)])/length(active_index)) + sum((d_output_2[['sigma_omega']] %*% d_moment[,(nrow(mat_YK) + 1):(2 * nrow(mat_YK))])/length(active_index))

    output[[2]] = output[[2]][active_index_pref]
    output[[3]] = do.call('cbind', output[[3]])
    output[[4]] = abs(mean(output_1) - mean(mat_M[,1]))
    return(output)
  }

  aggregate_moment_theta = function(x_transform) {
    # computing moment using realized expenses
    if (Sys.info()[['sysname']] == 'Windows') {
      clusterExport(cl, 'x_transform',envir=environment())
      moment_realized_expense = parLapply(cl, data_hh_list_theta,function(mini_data) {
        output = tryCatch(identify_theta(mini_data, x_transform[[1]], n_draw_halton = n_draw_halton),error=function(e) e)
        return(output)
      })
    } else {
      moment_realized_expense = mclapply(data_hh_list_theta, function(mini_data) tryCatch(identify_theta(mini_data, x_transform[[1]], n_draw_halton = n_draw_halton), error=function(e) e), mc.cores=numcores)
    }
    
    moment_realized_expense_val = do.call('c', lapply(moment_realized_expense, function(x) x[[1]])) %>% mean
    # print(paste0('moment from theta = ', moment_realized_expense_val))
    moment_realized_expense_deriv = rep(0, length(param_trial));
    moment_realized_expense_deriv[x_transform[[2]]$beta_theta] = do.call('rbind', lapply(moment_realized_expense, function(x) x[[2]]$beta_theta)) %>% colMeans
    moment_realized_expense_deriv[x_transform[[2]]$beta_theta_ind] = do.call('rbind', lapply(moment_realized_expense, function(x) x[[2]]$beta_theta_ind)) %>% colMeans
    moment_realized_expense_deriv[x_transform[[2]]$sigma_thetabar] = do.call('c', lapply(moment_realized_expense, function(x) x[[2]]$sigma_thetabar)) %>% mean

    return(list(moment_realized_expense_val, moment_realized_expense_deriv))
  }

  active_index = x_transform[[2]][c('beta_theta','sigma_thetabar', 'beta_theta_ind', var_list)] %>% unlist()
  tol = 1e-4

  iteration = 1; 

  init_pref = aggregate_moment_pref(transform_param(param_trial_here, return_index=TRUE))
  init_theta = aggregate_moment_theta(transform_param(param_trial_here, return_index=TRUE))

  if (is.na(init_pref[[1]]) | is.nan(init_pref[[1]]) | is.na(init_theta[[1]]) | is.nan(init_theta[[1]])) {
    print('stop here')
    return(NA)
  }

  if (estimate_theta) {
    optim_pref_theta = splitfngr::optim_share(param_trial_here[active_index], function(x_pref_theta) {
      param_trial_inner_theta = param_trial_here
      param_trial_inner_theta[active_index] = x_pref_theta 
      x_transform = transform_param(param_trial_inner_theta, return_index = TRUE)
      if (max(abs(x_pref_theta)) > 4) {
        return(list(NA, rep(NA, length(active_index))))
      }
      output_theta = aggregate_moment_theta(x_transform)
      if (is.nan(output_theta[[1]])) {
        return(list(NA, rep(NA, length(active_index))))
      }
      pref_moment = aggregate_moment_pref(transform_param(param_trial_inner_theta, return_index=TRUE))
      deriv_theta_index = c(x_transform[[2]]$beta_theta[1], x_transform[[2]]$beta_theta_ind[1], x_transform[[2]]$sigma_thetabar)
      if (is.na(pref_moment[[1]]) | any(is.nan(pref_moment[[2]]))) {
        return(list(NA, rep(NA, length(x_pref_theta))))
      }

      pref_derivative = rep(0, length(param_trial_here));
      pref_derivative[x_transform[[2]][var_list] %>% unlist()] = pref_moment[[2]]
      for (i in deriv_theta_index) {
        param_trial_i = param_trial_inner_theta; param_trial_i[i] = param_trial_inner_theta[i] + tol
        fi = aggregate_moment_pref(transform_param(param_trial_i, return_index=TRUE), silent=TRUE)
        if (i == deriv_theta_index[1]) {
          pref_derivative[i] = (rowMeans((fi[[3]] - pref_moment[[3]])/tol) %*% X_ind_pref_with_year)/nrow(X_ind_pref)
        } else if (i == deriv_theta_index[2]) {
          pref_derivative[i] = (rowMeans((fi[[3]] - pref_moment[[3]])/tol) %*% X_ind_pref)/nrow(X_ind_pref)
        } else {
          pref_derivative[i] = (fi[[1]] - pref_moment[[1]])/tol
        }
      }
      output = list()
      output[[1]] = pref_moment[[1]] * length(sample_identify_pref) + output_theta[[1]] 
      output[[2]] = pref_derivative[active_index] * length(sample_identify_pref) + output_theta[[2]][active_index] 
      return(output)
    }, control=list(maxit=1e3), method='BFGS')

    param_trial_here[active_index] = optim_pref_theta$par
  }
  
  active_index = c(x_transform[[2]]$sigma_r,  x_transform[[2]]$sigma_theta, x_transform[[2]]$beta_r); 

  x_transform = transform_param(param_trial_here, return_index=TRUE)


  min_theta_R = do.call('c', lapply(sample_r_theta, function(output_hh_index) min(cbind(var_ind(data_hh_list[[output_hh_index]]), data_hh_list[[output_hh_index]]$Year == 2004, data_hh_list[[output_hh_index]]$Year == 2006, data_hh_list[[output_hh_index]]$Year == 2010, data_hh_list[[output_hh_index]]$Year == 2012) %*% x_transform[[1]]$beta_theta)))

  max_theta_R = do.call('c', lapply(sample_r_theta, function(output_hh_index) max(cbind(var_ind(data_hh_list[[output_hh_index]]), data_hh_list[[output_hh_index]]$Year == 2004, data_hh_list[[output_hh_index]]$Year == 2006, data_hh_list[[output_hh_index]]$Year == 2010, data_hh_list[[output_hh_index]]$Year == 2012) %*% x_transform[[1]]$beta_theta)))

  mean_theta_R = do.call('c', lapply(sample_r_theta, function(output_hh_index) mean(cbind(var_ind(data_hh_list[[output_hh_index]]), data_hh_list[[output_hh_index]]$Year == 2004, data_hh_list[[output_hh_index]]$Year == 2006, data_hh_list[[output_hh_index]]$Year == 2010, data_hh_list[[output_hh_index]]$Year == 2012) %*% x_transform[[1]]$beta_theta)))

  message('start estimation of r')

  if (Sys.info()[['sysname']] == 'Windows') {
    clusterExport(cl, c('x_transform', 'sick_parameters', 'xi_parameters'),envir=environment())
    moment_eligible_hh_output = parLapply(cl, sample_r_theta, function(mini_data_index) tryCatch(household_draw_theta_kappa_Rdraw(mini_data_index, x_transform[[1]], n_halton_at_r, 10, sick_parameters, xi_parameters, u_lowerbar = -1, derivative_r_threshold = TRUE), error=function(e) e))
  } else {
    moment_eligible_hh_output = mclapply(sample_r_theta, function(mini_data_index) tryCatch(household_draw_theta_kappa_Rdraw(mini_data_index, x_transform[[1]], n_halton_at_r, 10, sick_parameters, xi_parameters, u_lowerbar = -1, derivative_r_threshold = TRUE), error=function(e) e), mc.cores=numcores)
  }

  root_r = do.call('c',lapply(moment_eligible_hh_output, function(output_hh) output_hh$root_r))
  hh_theta = do.call('c',lapply(moment_eligible_hh_output, function(output_hh) output_hh$hh_theta))
  relevant_index = which(full_insurance_indicator + no_insurance_indicator == 1)
  fx_r = function(x_transform, derivative=FALSE, silent=TRUE) {
    root_r = do.call('c',lapply(moment_eligible_hh_output, function(output_hh) output_hh$root_r))
    hh_theta = do.call('c',lapply(moment_eligible_hh_output, function(output_hh) output_hh$hh_theta))
    sd_r = exp(x_transform[[1]]$sigma_r);
    correlation = x_transform[[1]]$correlation
    mean_vec = rep(X_hh_theta_r %*% x_transform[[1]]$beta_r, each = n_halton_at_r) + correlation * hh_theta

    denominator = (pnorm((5 - mean_vec)/sd_r) - pnorm((0 - mean_vec)/sd_r))

    full_insurance_prob = (pnorm((5 - mean_vec)/sd_r) - pnorm((root_r - mean_vec)/sd_r))/(denominator + 1e-20)

    no_insurance_prob = (pnorm((root_r - mean_vec)/sd_r) - pnorm((0 - mean_vec)/sd_r))/(denominator + 1e-20)

    output = -sum(colMeans(matrix(full_insurance_prob, nrow=n_halton_at_r))[full_insurance_indicator]^2) - sum(colMeans(matrix(no_insurance_prob, nrow=n_halton_at_r))[no_insurance_indicator]^2)

    if (derivative) {
      derivative_root_r = list(); 
      derivative_root_r$beta_theta = do.call('rbind', lapply(moment_eligible_hh_output, function(output_hh) (t(matrix(output_hh$derivative_root_r$theta_bar, nrow = output_hh$HHsize)) %*% output_hh$X_ind_year))); 
      derivative_root_r$beta_theta_ind = do.call('rbind', lapply(moment_eligible_hh_output, function(output_hh) (t(matrix(output_hh$derivative_root_r$theta_bar, nrow = output_hh$HHsize)) %*% output_hh$X_ind)))
      derivative_root_r$sigma_thetabar = do.call('c', lapply(moment_eligible_hh_output, function(output_hh) (rowSums(t(matrix(output_hh$derivative_root_r$theta_bar, nrow = output_hh$HHsize)))) * exp(x_transform[[1]][['sigma_thetabar']])))
      derivative_root_r$correlation = do.call('c', lapply(moment_eligible_hh_output, function(output_hh) (rowSums(t(matrix(output_hh$derivative_root_r$theta_bar, nrow = output_hh$HHsize))) * output_hh$hh_theta) * (x_transform[[1]][['correlation']])));
      derivative_root_r$beta_gamma = do.call('rbind', lapply(moment_eligible_hh_output, function(output_hh) (t(matrix(output_hh$derivative_root_r$gamma, nrow = output_hh$HHsize)) %*% output_hh$X_ind)))
      derivative_root_r$sigma_gamma = do.call('c', lapply(moment_eligible_hh_output, function(output_hh) (rowSums(t(matrix(output_hh$derivative_root_r$gamma, nrow = output_hh$HHsize)) * output_hh$gamma_draw)) * exp(x_transform[[1]][['sigma_gamma']])))
      derivative_root_r$beta_delta = do.call('rbind', lapply(moment_eligible_hh_output, function(output_hh) (t(matrix(output_hh$derivative_root_r$delta, nrow = output_hh$HHsize)) %*% output_hh$X_ind)))
      derivative_root_r$sigma_delta = do.call('c', lapply(moment_eligible_hh_output, function(output_hh) (rowSums(t(matrix(output_hh$derivative_root_r$delta, nrow = output_hh$HHsize)) * output_hh$delta_draw)) * exp(x_transform[[1]][['sigma_delta']])))
      derivative_root_r$beta_omega = do.call('rbind', lapply(moment_eligible_hh_output, function(output_hh) (apply(output_hh$X_hh, 2, function(x) x * output_hh$derivative_root_r$omega))))
      derivative_root_r$sigma_omega = do.call('c', lapply(moment_eligible_hh_output, function(output_hh) (output_hh$derivative_root_r$omega * output_hh$omega_draw) * exp(x_transform[[1]][['sigma_omega']])))

      denominator_derivative = list(); 
      denominator_derivative$mean = (-dnorm((5 - mean_vec)/sd_r)/sd_r + dnorm((0 - mean_vec)/sd_r)/sd_r)
      denominator_derivative$sd =  (-dnorm((5 - mean_vec)/sd_r)/sd_r^2 + dnorm((0 - mean_vec)/sd_r)/sd_r^2)

      full_insurance_prob_derivative = list(); 
      full_insurance_prob_derivative$mean = (-dnorm((5 - mean_vec)/sd_r)/sd_r + dnorm((root_r - mean_vec)/sd_r)/sd_r)/(denominator + 1e-20) - (pnorm((5 - mean_vec)/sd_r) - pnorm((root_r - mean_vec)/sd_r))/(denominator + 1e-20)^2 * denominator_derivative$mean; 
      full_insurance_prob_derivative$sd = (-dnorm((5 - mean_vec)/sd_r)/sd_r^2 + dnorm((root_r - mean_vec)/sd_r)/sd_r^2)/(denominator + 1e-20) - (pnorm((5 - mean_vec)/sd_r) - pnorm((root_r - mean_vec)/sd_r))/(denominator + 1e-20)^2 * denominator_derivative$sd;
      full_insurance_prob_derivative$root_r = list(); 
      for (name_i in c('beta_theta','beta_theta_ind', 'beta_gamma', 'beta_delta', 'beta_omega', 'sigma_thetabar', 'sigma_delta', 'sigma_omega', 'sigma_gamma', 'correlation')) {
        if (grepl('beta',name_i)) {
          full_insurance_prob_derivative$root_r[[name_i]] = matrix(0, nrow = length(sample_r_theta) * n_halton_at_r, ncol = ncol(derivative_root_r[[name_i]]));
          for (col_index in 1:ncol(full_insurance_prob_derivative$root_r[[name_i]])) {
            full_insurance_prob_derivative$root_r[[name_i]][,col_index] = - dnorm((root_r - mean_vec)/sd_r)/sd_r * derivative_root_r[[name_i]][,col_index]/(denominator + 1e-20)
          }
        } else {
          full_insurance_prob_derivative$root_r[[name_i]] =- dnorm((root_r - mean_vec)/sd_r)/sd_r * derivative_root_r[[name_i]]/(denominator + 1e-20);
        }
      }  

      no_insurance_prob_derivative = list(); 
      no_insurance_prob_derivative$mean = (- dnorm((root_r - mean_vec)/sd_r)/sd_r + dnorm((0 - mean_vec)/sd_r)/sd_r)/(denominator + 1e-20) - (pnorm((root_r - mean_vec)/sd_r) - pnorm((0 - mean_vec)/sd_r))/(denominator + 1e-20)^2 * denominator_derivative$mean; 
      no_insurance_prob_derivative$sd = (- dnorm((root_r - mean_vec)/sd_r)/sd_r^2 + dnorm((0 - mean_vec)/sd_r)/sd_r^2)/(denominator + 1e-20) - (pnorm((root_r - mean_vec)/sd_r) - pnorm((0 - mean_vec)/sd_r))/(denominator + 1e-20)^2 * denominator_derivative$sd; 

      no_insurance_prob_derivative$root_r = list()

      for (name_i in c('beta_theta','beta_theta_ind', 'beta_gamma', 'beta_delta', 'beta_omega', 'sigma_thetabar', 'sigma_delta', 'sigma_omega', 'sigma_gamma', 'correlation')) {
        if (grepl('beta',name_i)) {
          no_insurance_prob_derivative$root_r[[name_i]] = matrix(0, nrow = length(sample_r_theta) * n_halton_at_r, ncol = ncol(derivative_root_r[[name_i]]));
          for (col_index in 1:ncol(no_insurance_prob_derivative$root_r[[name_i]])) {
            no_insurance_prob_derivative$root_r[[name_i]][,col_index] = - dnorm((root_r - mean_vec)/sd_r)/sd_r * derivative_root_r[[name_i]][,col_index]/(denominator + 1e-20)
          }
        } else {
          no_insurance_prob_derivative$root_r[[name_i]] =- dnorm((root_r - mean_vec)/sd_r)/sd_r * derivative_root_r[[name_i]]/(denominator + 1e-20);
        }
      }  

      d_output = list(); 

      d_output$beta_r = -colSums(apply(X_hh_theta_r[full_insurance_indicator,], 2, function(x) x * colMeans(matrix(full_insurance_prob * 2 * full_insurance_prob_derivative$mean, nrow = n_halton_at_r))[full_insurance_indicator])) - colSums(apply(X_hh_theta_r[no_insurance_indicator,], 2, function(x) x * colMeans(matrix(no_insurance_prob * 2 * no_insurance_prob_derivative$mean, nrow = n_halton_at_r))[no_insurance_indicator]))
      d_output$sigma_r = -sum(colMeans(matrix(full_insurance_prob * 2 * full_insurance_prob_derivative$sd, nrow = n_halton_at_r))[full_insurance_indicator] * exp(x_transform[[1]]$sigma_r)) -sum(colMeans(matrix(no_insurance_prob * 2 * no_insurance_prob_derivative$sd, nrow = n_halton_at_r))[no_insurance_indicator] * exp(x_transform[[1]]$sigma_r)) 

      for (name_i in c('beta_theta','beta_theta_ind', 'beta_gamma', 'beta_delta', 'beta_omega', 'sigma_thetabar', 'sigma_delta', 'sigma_omega', 'sigma_gamma', 'correlation')) {
        if (grepl('beta',name_i)) {
          d_output[[name_i]] = -colSums(apply(full_insurance_prob_derivative$root_r[[name_i]], 2, function(x) colMeans(matrix(full_insurance_prob * 2 * x, nrow = n_halton_at_r))[full_insurance_indicator])) -colSums(apply(no_insurance_prob_derivative$root_r[[name_i]], 2, function(x) colMeans(matrix(no_insurance_prob * 2 * x, nrow = n_halton_at_r))[no_insurance_indicator])) 
        } else {
          d_output[[name_i]] = -sum(colMeans(matrix(full_insurance_prob * 2 * full_insurance_prob_derivative$root_r[[name_i]], nrow = n_halton_at_r))[full_insurance_indicator]) - sum(colMeans(matrix(no_insurance_prob * 2 * no_insurance_prob_derivative$root_r[[name_i]], nrow = n_halton_at_r))[no_insurance_indicator])
        }
      }  
    }
    
    if (!(silent)) {
      print(summary(colMeans(matrix(full_insurance_prob, nrow=n_halton_at_r))[relevant_index]))
      print(summary(full_insurance_indicator[relevant_index]))
      print(summary(colMeans(matrix(no_insurance_prob, nrow=n_halton_at_r))[relevant_index]))
      print(summary(no_insurance_indicator[relevant_index]))
    }

    if (!(derivative)) {
      if (is.nan(output) | is.infinite(output)) {
        return(NA)
      }
      else {
        return(output)
      }
    } else {
      d_output_vec = rep(0, length(x_transform[[1]]))
      for (name_i in names(d_output)) {
        d_output_vec[x_transform[[2]][[name_i]]] = d_output[[name_i]]
      }
      d_output_vec[which(is.na(d_output_vec))] = 0 
      if (is.nan(output) | is.infinite(output)) {
        return(list(NA, d_output_vec))
      }
      else {
        return(list(output, d_output_vec))
      }
    }

  }

  index_r = x_transform[[2]][c('beta_r', 'sigma_r', 'correlation')] %>% unlist()

  optim_r = splitfngr::optim_share(param_trial_here[index_r], function(x) {
    x_with_new_r = param_trial_here; 
    x_with_new_r[index_r] = x; 
    output = fx_r(transform_param(x_with_new_r, return_index=TRUE), silent=TRUE, derivative=TRUE)
    return(list(output[[1]], output[[2]][index_r]))
  }, control=list(maxit=1e4), method='BFGS') 

  print('-------fit--------')
  param_trial_here[index_r] = optim_r$par;
  print('output of r'); fx_r(transform_param(param_trial_here, return_index =TRUE), silent=FALSE)

  param_trial_here[c(x_transform[[2]]$beta_r, x_transform[[2]]$sigma_r, x_transform[[2]]$correlation)] = optim_r$par
  x_transform = transform_param(param_trial_here,return_index=TRUE); 
  print('output of preference moments'); aggregate_moment_pref(x_transform, silent=FALSE) 

  print(x_transform[[1]])
  output_2 = lapply(moment_eligible_hh_output, function(output_hh) {
      sd = exp(optim_r$par[length(optim_r$par)-1])
      correlation = optim_r$par[length(optim_r$par)]
      mean_vec = rep(X_hh_theta_r %*% optim_r$par[1:(length(optim_r$par)-2)], each = n_halton_at_r) + correlation * hh_theta
      full_insurance_prob = (pnorm((5 - mean_vec)/sd) - pnorm((output_hh$root_r - mean_vec)/sd))/(pnorm((5 - mean_vec)/sd) - pnorm((0 - mean_vec)/sd))
      full_insurance_prob[is.nan(full_insurance_prob)] = 1; 
      no_insurance_prob = 1 - full_insurance_prob 
      Em = list()
      Em$full_insurance = colSums(apply(matrix(output_hh$m, nrow = n_draw_halton), 2, function(x) x * full_insurance_prob/(sum(full_insurance_prob) + 1e-20)))
      Em$no_insurance = colSums(apply(matrix(output_hh$m, nrow = n_draw_halton), 2, function(x) x * no_insurance_prob/(sum(no_insurance_prob) + 1e-20)))
      return(Em)
    })

  Em_full_insurance = do.call('c', lapply(output_2, function(x) x$full_insurance))
  Em_no_insurance = do.call('c', lapply(output_2, function(x) x$no_insurance))

  output_2 =  sum(((Em_full_insurance - mat_M_rtheta[,1]) * full_insurance_indicator_ind_level)^2 + ((Em_no_insurance - mat_M_rtheta[,1]) * no_insurance_indicator_ind_level)^2, na.rm=TRUE)

  summary(Em_full_insurance * full_insurance_indicator_ind_level + Em_no_insurance * no_insurance_indicator_ind_level) %>% print
  summary(mat_M_rtheta[,1]) %>% print

  if (!(estimate_theta)) {
    optim_pref_theta = list()
    optim_pref_theta$value = 0 
  }

  if (return_result) {
    return(param_trial_here)
  } else {
    print(paste0('output_2 = ', output_2));
    print(paste0('optim_r$value = ', optim_r$value));
    print(paste0('optim_pref_theta$value = ', optim_pref_theta$value))
    initial_param_trial <<- param_trial_here; 
    return(output_2 + optim_r$value/length(relevant_index) + optim_pref_theta$value) 
  }
  
}

param_trial = compute_inner_loop(-1, return_result=TRUE, estimate_theta=TRUE, estimate_pref = TRUE)

stop()

estimate_r_thetabar = optimize(function(xs) {
  output = compute_inner_loop(log(xs))
  return(output)
}, c(1e-3,1)) 


param_trial = compute_inner_loop(log(estimate_r_thetabar$minimum), return_result=TRUE, estimate_theta=TRUE, estimate_pref = TRUE)

message('computing final param_trial')

param_final = list(); 
param_final$other = param_trial; 
param_final$xi = xi_parameters;
param_final$sick = sick_parameters

param = param_final 
transform_param_final = transform_param(param_final$other)

fit_sample = sample(Vol_HH_list_index, 1000)

for (seed_number in c(1:10)) {
  if (Sys.info()[['sysname']] == 'Windows') {
    clusterExport(cl, c('transform_param_final', 'param','counterfactual_household_draw_theta_kappa_Rdraw'))
    mini_fit_values = parLapply(cl, c(fit_sample), function(id) {
      output = counterfactual_household_draw_theta_kappa_Rdraw(id, transform_param_final, n_draw_halton, n_draw_gauss, param$sick, param$xi, u_lowerbar = -1, policy_mat_hh = policy_mat[[id]], seed_number = seed_number, constraint_function = function(x) x)
      output = as.data.frame(output)
      output$Y = data_hh_list[[id]]$Income; 
      output$m_observed = data_hh_list[[id]]$M_expense; 
      output$id = id;
      output$HHsize_s = data_hh_list[[id]]$HHsize_s; 
      return(output)
    })
  } else {
    mini_fit_values = mclapply(c(fit_sample), function(id) {
    output = tryCatch(counterfactual_household_draw_theta_kappa_Rdraw(id, transform_param_final, n_draw_halton, n_draw_gauss, param$sick, param$xi, u_lowerbar = -1, policy_mat_hh = policy_mat[[id]], seed_number = seed_number, constraint_function = function(x) x), error=function(e) e)
    output = as.data.frame(output)
    output$Y = data_hh_list[[id]]$Income; 
    output$m_observed = data_hh_list[[id]]$M_expense; 
    output$id = id; 
    output$HHsize_s = data_hh_list[[id]]$HHsize_s; 
    return(output)}, mc.cores=numcores)
  }
  if (seed_number == 1) {
    fit_values = do.call('rbind', mini_fit_values)
  } else {
    fit_values = rbind(fit_values, do.call('rbind', mini_fit_values))
  }
}


observed_data_voluntary = do.call('rbind', data_hh_list[c(fit_sample)])

fit_values = as.data.frame(fit_values)

observed_data_voluntary = as.data.frame(observed_data_voluntary)

predicted_data_summary = fit_values %>% filter(Com_sts + Bef_sts + Std_w_ins == 0) %>% mutate(Y2 = as.numeric(Hmisc::cut2(Y, g=5))) %>% group_by(Y2) %>% summarise(mean_Vol_sts = mean(vol_sts_counterfactual), mean_m = mean(m, na.rm=TRUE))
predicted_data_summary$type = 'predicted'
actual_data_summary = observed_data_voluntary %>% filter(Com_sts + Bef_sts + Std_w_ins == 0)%>% mutate(Y2 = as.numeric(Hmisc::cut2(Income, g=5))) %>% group_by(Y2) %>% summarise(mean_Vol_sts = mean(Vol_sts), mean_m = mean(M_expense, na.rm=TRUE))
actual_data_summary$type = 'actual'

plot_1 = ggplot(data = rbind(predicted_data_summary, actual_data_summary), aes(x = Y2, y = mean_Vol_sts, color=type)) + geom_line() 
plot_2 = ggplot(data = rbind(predicted_data_summary, actual_data_summary), aes(x = Y2, y = mean_m, color=type)) + geom_line() 
plot = gridExtra::grid.arrange(plot_1, plot_2, nrow=1)

if (dir.exists('../../householdbundling_estimate')) {
  saveRDS(param_final, file=paste0('../../householdbundling_estimate/estimate_',job_index,'.rds'))
  saveRDS(plot, file=paste0('../../householdbundling_estimate/estimate_fit_',job_index,'.rds'))
} else {
  dir.create('../../householdbundling_estimate') 
  saveRDS(param_final, file=paste0('../../householdbundling_estimate/estimate_',job_index,'.rds'))
  saveRDS(plot, file=paste0('../../householdbundling_estimate/estimate_fit_',job_index,'.rds'))
}