args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) { 
  if (Sys.info()[['sysname']] == 'Windows') {
    numcores = 10;
  } else {
    numcores = 8;
  }
  job_index = as.integer(Sys.time());  
} else {
  job_index = as.numeric(args[1]);
  numcores = as.numeric(args[2]); 
}

if (numcores < 6) {
  mini=TRUE
} else {
  mini=FALSE
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
      if (((data_mini$HHsize_s[1] == 0) & (sum(data_mini$sick_dummy) == 0))) {
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
if (mini) {
  message('estimating in mini mode')
  sample_r_theta = sample(sample_r_theta, 500, replace=TRUE)
  sample_identify_pref = sample(sample_identify_pref, length(sample_identify_pref), replace=TRUE)
  sample_identify_theta = sample(sample_identify_theta, length(sample_identify_theta), replace=TRUE)

  n_draw_halton = 20;

  n_halton_at_r = 20;

  n_draw_gauss = 10;
} else {
  sample_r_theta = sample(sample_r_theta, 3000, replace=TRUE)
  sample_identify_pref = sample(sample_identify_pref, length(sample_identify_pref), replace=TRUE)
  sample_identify_theta = sample(sample_identify_theta, length(sample_identify_theta), replace=TRUE)

  n_draw_halton = 100;

  n_halton_at_r = 100;

  n_draw_gauss = 10;
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

message('Preparing data for estimation')
data_hh_list_theta = list()
for (index in 1:length(sample_identify_theta)) {
  message(paste0('constructing data for index', index))
  data_hh_list_theta[[index]] = household_draw_theta_kappa_Rdraw(hh_index=sample_identify_theta[index], param=transform_param_trial[[1]], n_draw_halton = 1000, n_draw_gauss = n_draw_gauss, sick_parameters, xi_parameters, short=FALSE)
}


message('Estimating theta')

if (Sys.info()[['sysname']] == 'Windows') {
  cl = makeCluster(numcores);
  clusterEvalQ(cl, library('tidyverse'))
  clusterEvalQ(cl, library('familyenrollment'))
  clusterExport(cl, c('param_trial', 'data_hh_list_theta'))
}
 


if (Sys.info()[['sysname']] == 'Windows') {
  clusterExport(cl, c('param_trial', 'transform_param_trial', 'sample_identify_pref', 'xi_parameters', 'sick_parameters'))
  clusterExport(cl,'n_draw_halton')
}


message('compute covariates before estimation')

if (Sys.info()[['sysname']] == 'Windows') {
  data_hh_list_pref = parLapply(cl, sample_identify_pref,function(index) {
    output = tryCatch(household_draw_theta_kappa_Rdraw(hh_index=index, param=transform_param_trial[[1]], n_draw_halton = n_draw_halton, n_draw_gauss = n_draw_gauss, sick_parameters, xi_parameters, short=FALSE),error=function(e) e)
    return(output)
  })
} else {
  data_hh_list_pref = mclapply(sample_identify_pref, function(index) tryCatch(household_draw_theta_kappa_Rdraw(hh_index=index, param=transform_param_trial[[1]], n_draw_halton = n_draw_halton, n_draw_gauss = n_draw_gauss, sick_parameters, xi_parameters, short=FALSE), error=function(e) e), mc.cores=numcores)
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
initial_param_trial[x_transform[[2]]$beta_theta[1]] = -2;
initial_param_trial[x_transform[[2]]$sigma_theta] = -1.4;


iteration = 1;
save_output = list();
tol = 1e-4

pref_list = c('beta_omega', 'beta_delta', 'beta_gamma', 'sigma_omega', 'sigma_delta', 'sigma_gamma') 
list_theta_var = c('beta_theta', 'beta_theta_ind', 'sigma_thetabar', 'sigma_theta')

dir.create('../../householdbundling_estimate') 

aggregate_moment_pref = function(x_transform, silent=TRUE, recompute_pref=FALSE) {
  param_trial_here = as.vector(unlist(x_transform[[1]]))
  print(param_trial_here)
  if (Sys.info()[['sysname']] == 'Windows') {
    clusterExport(cl, c('x_transform', 'n_halton_at_r'),envir=environment())
    data_hh_list_pref = parLapply(cl, sample_identify_pref,function(index) {
      output = tryCatch(household_draw_theta_kappa_Rdraw(hh_index=index, param=x_transform[[1]], n_draw_halton = n_draw_halton, n_draw_gauss = n_draw_gauss, sick_parameters, xi_parameters, short=FALSE),error=function(e) e)
      return(output)
    })
    mat_YK = do.call('cbind', parLapply(cl, data_hh_list_pref, function(x) {
          output = rbind(colMeans(matrix(x$kappa_draw[[1]], nrow=n_draw_halton)), colMeans(matrix(x$kappa_draw[[1]]^2, nrow=n_draw_halton)), x$income[1], colMeans(matrix(x$kappa_draw[[1]], nrow=n_draw_halton))*x$income[1])
          return(output)
        }))
  } else {
    data_hh_list_pref = mclapply(sample_identify_pref, function(index) tryCatch(household_draw_theta_kappa_Rdraw(hh_index=index, param=x_transform[[1]], n_draw_halton = n_draw_halton, n_draw_gauss = n_draw_gauss, sick_parameters, xi_parameters, short=FALSE), error=function(e) e), mc.cores=numcores)
    mat_YK = do.call('cbind', mclapply(data_hh_list_pref, function(x) {
          output = rbind(colMeans(matrix(x$kappa_draw[[1]], nrow=n_draw_halton)), colMeans(matrix(x$kappa_draw[[1]]^2, nrow=n_draw_halton)), x$income[1], colMeans(matrix(x$kappa_draw[[1]], nrow=n_draw_halton))*x$income[1])
          return(output)
        }, mc.cores=numcores))
  }


  mat_YK = rbind(mat_YK, do.call('c', lapply(data_hh_list[sample_identify_pref], function(x) x$Income)))
  mat_YK = rbind(1, mat_YK)

  mini_f = function(x_transform) {
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
    output_2 = do.call('c', lapply(moment_ineligible_hh_output, function(x) x[[2]]))

    print('fit of ineligible hh')
    print(summary(output_1))
    print(summary(mat_M[,1]))
    d_output_1 = list();
    d_output_2 = list();
    
    for (varname in pref_list) {
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
    for (moment_index in c(1:nrow(mat_YK))) {
      output[[3]][[moment_index]] = ((output_1 - mat_M[,1]) * mat_YK[moment_index,])^2
      moment[moment_index] = mean(output[[3]][[moment_index]]); 
      d_moment[[moment_index]] = 2*((output_1 - mat_M[,1]) * mat_YK[moment_index,]^2);
    }
    for (moment_index in c(1:nrow(mat_YK))) {
      output[[3]][[moment_index + nrow(mat_YK)]] = ((output_2 - mat_M[,2]) * mat_YK[moment_index,])^2
      moment[moment_index + nrow(mat_YK)] = mean(output[[3]][[moment_index + nrow(mat_YK)]]); 
      d_moment[[moment_index + nrow(mat_YK)]] = 2*((output_2 - mat_M[,2]) * mat_YK[moment_index,]^2);
    }

    d_moment = do.call('cbind', d_moment)

    output[[1]] = sum(moment); 
    output[[2]] = rep(0, length(initial_param_trial)); 

    output[[2]][x_transform[[2]][['beta_delta']]] = apply(d_output_1[['beta_delta']], 1, function(x) sum((x %*% d_moment[,1:nrow(mat_YK)])/length(x))) + apply(d_output_2[['beta_delta']], 1, function(x) sum((x %*% d_moment[,(nrow(mat_YK) + 1):(2 * nrow(mat_YK))])/length(x)))
    output[[2]][x_transform[[2]][['beta_omega']]] = apply(d_output_1[['beta_omega']], 1, function(x) sum((x %*% d_moment[,1:nrow(mat_YK)])/length(x))) + apply(d_output_2[['beta_omega']], 1, function(x) sum((x %*% d_moment[,(nrow(mat_YK) + 1):(2 * nrow(mat_YK))])/length(x)))
    output[[2]][x_transform[[2]][['beta_gamma']]] = apply(d_output_1[['beta_gamma']], 1, function(x) sum((x %*% d_moment[,1:nrow(mat_YK)])/length(x))) + apply(d_output_2[['beta_gamma']], 1, function(x) sum((x %*% d_moment[,(nrow(mat_YK) + 1):(2 * nrow(mat_YK))])/length(x)))

    output[[2]][x_transform[[2]][['sigma_delta']]] =  sum((d_output_1[['sigma_delta']] %*% d_moment[,1:nrow(mat_YK)])/nrow(mat_YK) ) + sum((d_output_2[['sigma_delta']] %*% d_moment[,(nrow(mat_YK) + 1):(2 * nrow(mat_YK))])/nrow(mat_YK) )
    output[[2]][x_transform[[2]][['sigma_gamma']]] =  sum((d_output_1[['sigma_gamma']] %*% d_moment[,1:nrow(mat_YK)])/nrow(mat_YK) ) + sum((d_output_2[['sigma_gamma']] %*% d_moment[,(nrow(mat_YK) + 1):(2 * nrow(mat_YK))])/nrow(mat_YK) )
    output[[2]][x_transform[[2]][['sigma_omega']]] =  sum((d_output_1[['sigma_omega']] %*% d_moment[,1:nrow(mat_YK)])/nrow(mat_YK) ) + sum((d_output_2[['sigma_omega']] %*% d_moment[,(nrow(mat_YK) + 1):(2 * nrow(mat_YK))])/nrow(mat_YK))

    output[[3]] = do.call('cbind', output[[3]])
    output[[4]] = output_1
    return(output)
  }

  if (recompute_pref) {
    output_initial = mini_f(x_transform);
    if (mean(output_initial[[4]]) > 5 * mean(mat_M[,1])) {
      return(list(par = initial_param_trial[x_transform[[2]][pref_list] %>% unlist()]))
    }
    mini_param = splitfngr::optim_share(initial_param_trial[x_transform[[2]][pref_list] %>% unlist()], function(x) {
      xmini_new = x;
      xmini_new[(length(x)-2):length(x)] = 1 - exp(x[(length(x)-2):length(x)])
      x_new = param_trial_here; 
      x_new[x_transform[[2]][pref_list] %>% unlist()] = xmini_new; 
      output = mini_f(transform_param(x_new, return_index = TRUE)); 
      print(paste0('output of aggregate_moment_pref  = ', output[[1]])); 
      print('x = '); print(xmini_new)
      deriv = output[[2]][x_transform[[2]][pref_list] %>% unlist()];
      deriv[(length(x)-2):length(x)] = deriv[(length(x)-2):length(x)] * exp(x[(length(x)-2):length(x)]) * (-1)
      print('derivative = '); print(deriv);
      return(list(output[[1]], deriv))}, method='BFGS')
    mini_param$par[(length(mini_param$par)-2):length(mini_param$par)] = 1 - exp(mini_param$par[(length(mini_param$par)-2):length(mini_param$par)])
    return(mini_param)
  } else {
    print(x_transform)
    return(mini_f(x_transform))
  }
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
   moment_realized_expense_deriv[x_transform[[2]]$sigma_theta] = do.call('c', lapply(moment_realized_expense, function(x) x[[2]]$sigma_theta)) %>% mean

  return(list(moment_realized_expense_val, moment_realized_expense_deriv))
}

index_theta_only = x_transform[[2]][list_theta_var] %>% unlist(); 
index_pref_only = x_transform[[2]][pref_list] %>% unlist();

save_output = list()
iter = 1; 
optim_f =  function(x_pref_theta) {
    print('x_pref_theta'); print(x_pref_theta)
    # if (max(abs(x_pref_theta)) > 10) {
    #   return(list(NA, rep(NA, length(x_pref_theta))))
    # }
    param_trial_here = initial_param_trial; 
    param_trial_here[c(index_theta_only)] = x_pref_theta 
    optim_pref = aggregate_moment_pref(transform_param(param_trial_here, return_index = TRUE), recompute_pref=TRUE); 
    param_trial_here[index_pref_only] = optim_pref$par; 

    x_transform = transform_param(param_trial_here, return_index = TRUE)

    output_theta = aggregate_moment_theta(x_transform)

    if (is.nan(output_theta[[1]])) {
      return(list(NA, rep(NA, length(c(index_theta_only, index_pref_only)))))
    }
    pref_moment = aggregate_moment_pref(x_transform)
      
    if (Sys.info()[['sysname']] == 'Windows') {
      clusterExport(cl, c('x_transform', 'sick_parameters', 'xi_parameters'),envir=environment())
      moment_eligible_hh_output = parLapply(cl, sample_r_theta, function(mini_data_index) tryCatch(household_draw_theta_kappa_Rdraw(mini_data_index, x_transform[[1]], n_halton_at_r, 10, sick_parameters, xi_parameters, u_lowerbar = -1, derivative_r_threshold = TRUE), error=function(e) e))
    } else {
      moment_eligible_hh_output = mclapply(sample_r_theta, function(mini_data_index) tryCatch(household_draw_theta_kappa_Rdraw(mini_data_index, x_transform[[1]], n_halton_at_r, 10, sick_parameters, xi_parameters, u_lowerbar = -1, derivative_r_threshold = TRUE), error=function(e) e), mc.cores=numcores)
    }

    X_ind_pref_with_year_r = do.call('rbind', lapply(moment_eligible_hh_output, function(x) x$X_ind_year));

    X_ind_pref_r = do.call('rbind', lapply(moment_eligible_hh_output, function(x) x$X_ind));

    output_wrt_r = function(x_transform, silent = TRUE) {
      param_trial_inner_r = as.vector(unlist(x_transform[[1]]))
      if (Sys.info()[['sysname']] == 'Windows') {
        clusterExport(cl, c('x_transform', 'sick_parameters', 'xi_parameters'),envir=environment())
        moment_eligible_hh_output = parLapply(cl, sample_r_theta, function(mini_data_index) tryCatch(household_draw_theta_kappa_Rdraw(mini_data_index, x_transform[[1]], n_halton_at_r, 10, sick_parameters, xi_parameters, u_lowerbar = -1, derivative_r_threshold = FALSE), error=function(e) e))
      } else {
        moment_eligible_hh_output = mclapply(sample_r_theta, function(mini_data_index) tryCatch(household_draw_theta_kappa_Rdraw(mini_data_index, x_transform[[1]], n_halton_at_r, 10, sick_parameters, xi_parameters, u_lowerbar = -1, derivative_r_threshold = FALSE), error=function(e) e), mc.cores=numcores)
      }

      root_r = do.call('rbind',lapply(moment_eligible_hh_output, function(output_hh) output_hh$root_r))
      hh_theta = do.call('c',lapply(moment_eligible_hh_output, function(output_hh) output_hh$hh_theta))
      relevant_index = which(full_insurance_indicator + no_insurance_indicator == 1)
      fx_r = function(x_transform, derivative=FALSE, silent=TRUE) {
        sd_r = exp(x_transform[[1]]$sigma_r);
        correlation = x_transform[[1]]$correlation
        mean_vec = rep(X_hh_theta_r %*% x_transform[[1]]$beta_r, each = n_halton_at_r) + correlation * hh_theta

        denominator = (pnorm((5 - mean_vec)/sd_r) - pnorm((0 - mean_vec)/sd_r))

        output = -sum(colMeans(matrix((pnorm((root_r[,2] - mean_vec)/sd_r) - pnorm((root_r[,1] - mean_vec)/sd_r))/denominator, nrow = n_halton_at_r))) 
        
        output[which(is.nan(output))] = 0;
        if (!(silent)) {
          print(summary(colMeans(matrix((pnorm((root_r[,2] - mean_vec)/sd_r) - pnorm((root_r[,1] - mean_vec)/sd_r))/denominator, nrow = n_halton_at_r))))
        }

        if (!(derivative)) {
          if (is.nan(output) | is.infinite(output)) {
            return(NA)
          }
          else {
            return(output)
          }
        }
      }

      index_r = x_transform[[2]][c('beta_r', 'sigma_r', 'correlation')] %>% unlist()

      optim_r = optim(param_trial_inner_r[index_r], function(x) {
        x_with_new_r = param_trial_inner_r; 
        x_with_new_r[index_r] = x; 
        output = fx_r(transform_param(x_with_new_r, return_index=TRUE), silent=silent, derivative=FALSE)
        return(output[[1]])
      }, control=list(maxit=1e4), method='BFGS') 

      param_trial_inner_r[c(x_transform[[2]]$beta_r, x_transform[[2]]$sigma_r, x_transform[[2]]$correlation)] = optim_r$par
      x_transform = transform_param(param_trial_inner_r,return_index=TRUE); 


      output_2 = lapply(moment_eligible_hh_output, function(output_hh) {
          sd_r = exp(x_transform[[1]]$sigma_r)
          correlation = x_transform[[1]]$correlation
          mean_vec = rep(output_hh$X_hh %*% x_transform[[1]]$beta_r, each = n_halton_at_r) + correlation * output_hh$hh_theta
          denominator = (pnorm((5 - mean_vec)/sd_r) - pnorm((0 - mean_vec)/sd_r)) + 1e-20;

          prob_optimal = (matrix((pnorm((output_hh$root_r[,2] - mean_vec)/sd_r) - pnorm((output_hh$root_r[,1] - mean_vec)/sd_r))/denominator, nrow = n_halton_at_r))

          Em = apply(matrix(output_hh$m, nrow = n_halton_at_r), 2, function(x) sum(x * prob_optimal/sum(prob_optimal)))
          return(list(Em, rep(mean(prob_optimal), length(Em))))
        })

      Em = do.call('c', lapply(output_2, function(x) x[[1]]))
      print('summary(Em)');print(summary(Em));

      prob_optimal = do.call('c', lapply(output_2, function(x) x[[2]]))
      m_observed = do.call('c', lapply(data_hh_list[sample_r_theta], function(x) x$M_expense))
      print('m_observed');print(summary(m_observed));

      return(list(Em = Em, m_observed = m_observed, insurance_match = optim_r$value, optim = optim_r, prob_optimal = prob_optimal))
    }

    f0 = output_wrt_r(x_transform, silent = FALSE); 

    param_trial_here[c(x_transform[[2]]$beta_r, x_transform[[2]]$sigma_r, x_transform[[2]]$correlation)] = f0$optim$par
    output_r =  sum((f0$Em - f0$m_observed)^2, na.rm=TRUE)  - sum(f0$prob_optimal) 

    if (is.na(output_r) | max(f0$Em[which(!(is.na(f0$Em)))]) == 0) {
      return(list(NA, rep(NA, length(index_theta_only))))
    }
    r_derivative = rep(0, length(param_trial_here))

    for (name_i in list_theta_var) {
      i = x_transform[[2]][[name_i]][1]
      param_trial_i = param_trial_here; param_trial_i[i] = param_trial_here[i] + tol
      fi = output_wrt_r(transform_param(param_trial_i, return_index=TRUE))
      output_2i =  sum((fi$Em - fi$m_observed)^2)  - sum(f0$prob_optimal) 
      if (name_i == 'beta_theta') {
        r_derivative[x_transform[[2]][[name_i]]] = colSums(apply(X_ind_pref_with_year_r, 2, function(x) x * 2 * ((f0$Em - f0$m_observed) * (fi$Em - f0$Em)/tol)), na.rm=TRUE) -  colSums(apply(X_ind_pref_with_year_r, 2, function(x) x * (fi$prob_optimal - f0$prob_optimal)/tol), na.rm=TRUE)
      } else if (name_i == 'beta_theta_ind') {
        r_derivative[x_transform[[2]][[name_i]]] = colSums(apply(X_ind_pref_r, 2, function(x) x * 2 * ((f0$Em - f0$m_observed) * (fi$Em - f0$Em)/tol)),na.rm=TRUE) -  colSums(apply(X_ind_pref_r, 2, function(x) x * (fi$prob_optimal - f0$prob_optimal)/tol),na.rm=TRUE)
      } else {
        r_derivative[x_transform[[2]][[name_i]]] = sum(2 * (f0$Em - f0$m_observed) * (fi$Em - f0$Em)/tol, na.rm=TRUE) -  sum(fi$prob_optimal - f0$prob_optimal)/tol
      }
    }

    print('----------FINAL OUTPUT---------');
    print(pref_moment[[1]] + output_theta[[1]] + output_r/length(f0$Em))
    print('--------------------')
    save_output[[iter]] <<- param_trial_here
    iter <<- iter + 1; 
    saveRDS(save_output, file='../../householdbundling_estimate/save_output_',job_index,'.rds')
    return(list(pref_moment[[1]] * 100 + output_theta[[1]] + output_r/length(f0$Em), pref_moment[[2]][index_theta_only] * 100 + output_theta[[2]][index_theta_only] + r_derivative[index_theta_only]/length(f0$Em), param_trial_here))
}


optim_pref_theta = splitfngr::optim_share(initial_param_trial[index_theta_only], function(x) {
  output = try(optim_f(x))
  if("try-error" %in% class(output)) {
    return(list(NA, rep(NA, length(index_theta_only))))
  } else {
    return(output[1:2])
  }
}, control=list(maxit=1e3), method='BFGS')

param_trial = optim_f(optim_pref_theta$par)[[3]]; 



message('computing final param_trial')
 
param_final = list(); 
param_final$other = param_trial; 
param_final$xi = xi_parameters;
param_final$sick = sick_parameters

param = param_final 
transform_param_final = transform_param(param_final$other)

fit_sample = Com_HH_list_index

for (seed_number in c(1:1)) {
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
  if (seed_number == min(seed_number)) {
    fit_values = do.call('rbind', mini_fit_values)
  } else {
    fit_values = rbind(fit_values, do.call('rbind', mini_fit_values))
  }
}


observed_data_voluntary = do.call('rbind', data_hh_list[c(fit_sample)])

fit_values = as.data.frame(fit_values)

observed_data_voluntary = as.data.frame(observed_data_voluntary)

predicted_data_summary = fit_values  %>% mutate(Y2 = as.numeric(Hmisc::cut2(Y, g=5))) %>% group_by(Y2) %>% summarise(mean_Vol_sts = mean(vol_sts_counterfactual), mean_m = mean(m, na.rm=TRUE))
predicted_data_summary$type = 'predicted'
actual_data_summary = observed_data_voluntary %>% mutate(Y2 = as.numeric(Hmisc::cut2(Income, g=5))) %>% group_by(Y2) %>% summarise(mean_Vol_sts = mean(Vol_sts), mean_m = mean(M_expense, na.rm=TRUE))
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