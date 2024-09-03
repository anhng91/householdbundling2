args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) { 
  if (Sys.info()[['sysname']] == 'Windows') {
    numcores = 20;
  } else if (Sys.info()[['sysname']] == 'Linux') {
    numcores = 12;
  } else {
    numcores = 6;
  }
} else {
  numcores = as.numeric(args[1]); 
}

if (Sys.info()[['nodename']] == 'Anh-Macbook-3.local' | grepl("vpn", Sys.info()[['nodename']])  | grepl("macbook", Sys.info()[['nodename']]) ) {
  mini=TRUE
  numcores = 4; 
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
sample_identify_theta_f = full_index[which((lapply(full_index, function(sample_index_i) {
      data_mini = data_hh_list[[sample_index_i]]
      if ((data_mini$Year[1] == 2008) & (data_mini$Income[1] < 0) & (data_mini$HHsize_s[1] == 0)) {
        return(1)
      } else {
        return(0)
      }
    }) %>% unlist()) == 1)]

sample_no_sick_f = full_index[which((lapply(full_index, function(sample_index_i) {
      data_mini = data_hh_list[[sample_index_i]]
      if (((data_mini$HHsize_s[1] == 0) & (sum(data_mini$sick_dummy) == 0))) {
        return(1)
      } else {
        return(0)
      }
    }) %>% unlist()) == 1)]

message('Identify the sample for preference')
sample_identify_pref_f = lapply(Com_HH_list_index, function(x) ifelse(x %in% c(sample_identify_theta_f, sample_no_sick_f), NA, x)) %>% unlist()
sample_identify_pref_f = sample_identify_pref_f[!(is.na(sample_identify_pref_f))]

message('Change pref sample')
sample_identify_pref_f = Com_HH_list_index

out_sample_index_f = lapply(1:length(data_hh_list), function(hh_index) {
  data = data_hh_list[[hh_index]]; 
  if (nrow(data) > nrow(data %>% filter(Bef_sts + Com_sts + Std_w_ins == 1)) & (data$Year[1] == 2006 & data$HHsize_s[1] > 1)) {
    return(hh_index);
  }
  else {
    return(NA); 
  }
}) %>% unlist(); 

out_sample_index_f = out_sample_index_f[!(is.na(out_sample_index_f))] 
message('iterating over job indices') 

sample_r_theta_f = Vol_HH_list_index

for (job_index_iter in c(1:100)) {
  job_index = as.integer(Sys.time());  

  # Bootstrapping indices 
  message('bootstrapping indices')
  set.seed(job_index);
  sample_index = sample(1:length(data_hh_list), length(data_hh_list), replace=TRUE)
  if (mini) {
    message('estimating in mini mode')
    sample_r_theta = sample(sample_r_theta_f, round(length(sample_r_theta_f)/50), replace=TRUE)
    sample_identify_pref = sample(sample_identify_pref_f, round(length(sample_identify_pref_f)/50), replace=TRUE)
    out_sample_index = sample(out_sample_index_f, round(length(out_sample_index_f)/50), replace=TRUE)
    sample_identify_theta = sample(sample_identify_theta_f, length(sample_identify_theta_f), replace=TRUE)

    n_draw_halton = 10;

    n_halton_at_r = 10;

    n_draw_gauss = 10;
  } else {
    sample_r_theta = sample(sample_r_theta_f, round(length(sample_r_theta_f)/10), replace=TRUE)
    sample_identify_pref = sample(sample_identify_pref_f, round(length(sample_identify_pref_f)/10), replace=TRUE)
    out_sample_index = sample(out_sample_index_f, round(length(out_sample_index_f)/10), replace=TRUE)
    sample_identify_theta = sample(sample_identify_theta_f, length(sample_identify_theta_f), replace=TRUE)
    n_draw_halton = 10;

    n_halton_at_r = 10;

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
    data_hh_list_theta[[index]] = household_draw_theta_kappa_Rdraw(hh_index=sample_identify_theta[index], param=transform_param_trial[[1]], n_draw_halton = 10, n_draw_gauss = n_draw_gauss, sick_parameters, xi_parameters, short=FALSE)
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

  mat_M2 = do.call('c', lapply(data_hh_list_pref, function(mini_data) c(mini_data$data$M_expense %*% t(mini_data$data$M_expense))))

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

  initial_param_trial = init_param
  initial_param_trial[x_transform[[2]]$beta_theta_ind[1]] = log(initial_param_trial[x_transform[[2]]$beta_theta_ind[1]])
  # initial_param_trial = rep(0, length(init_param))
  # initial_param_trial[x_transform[[2]]$beta_theta[1]] = -0.2;
  # initial_param_trial[x_transform[[2]]$sigma_theta] = log(0.2);
  # initial_param_trial[x_transform[[2]]$beta_delta[1]] = 0;
  # initial_param_trial[x_transform[[2]]$beta_theta_ind[1]] = log(0.1);
  # initial_param_trial[x_transform[[2]]$sigma_thetabar] = log(0.2);
  # initial_param_trial[x_transform[[2]]$beta_omega[1]] = 1;
  # initial_param_trial[x_transform[[2]]$beta_gamma[1]] = 0;
  # initial_param_trial[x_transform[[2]]$sigma_gamma[1]] = -1;
  # initial_param_trial[x_transform[[2]]$sigma_omega[1]] = -1;

  if (Sys.info()[['sysname']] == 'Windows') {
    clusterExport(cl, c('initial_param_trial'))
  }

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
        output = tryCatch(household_draw_theta_kappa_Rdraw(hh_index=index, param=x_transform[[1]], n_draw_halton = n_draw_halton, n_draw_gauss = n_draw_gauss, sick_parameters, xi_parameters, short=FALSE, realized_sick = TRUE),error=function(e) e)
        return(output)
      })
      n_draw_here = data_hh_list_pref[[1]]$theta_draw %>% nrow()
      mat_YK = do.call('cbind', parLapply(cl, data_hh_list_pref, function(x) {
            output = rbind(colMeans(matrix(x$kappa_draw[[1]], nrow=n_draw_here)), colMeans(matrix(x$kappa_draw[[1]]^2, nrow=n_draw_here)), (x$income[1] + 2), colMeans(matrix(x$kappa_draw[[1]], nrow=n_draw_here))*(x$income[1] + 2))
            return(output)
          }))
    } else {
      data_hh_list_pref = mclapply(sample_identify_pref, function(index) tryCatch(household_draw_theta_kappa_Rdraw(hh_index=index, param=x_transform[[1]], n_draw_halton = n_draw_halton, n_draw_gauss = n_draw_gauss, sick_parameters, xi_parameters, short=FALSE, realized_sick = TRUE), error=function(e) e), mc.cores=numcores)
      n_draw_here = data_hh_list_pref[[1]]$theta_draw %>% nrow()
      mat_YK = do.call('cbind', mclapply(data_hh_list_pref, function(x) {
            output = rbind(colMeans(matrix(x$kappa_draw[[1]], nrow=n_draw_here)), colMeans(matrix(x$kappa_draw[[1]]^2, nrow=n_draw_here)), (x$income[1] + 2), colMeans(matrix(x$kappa_draw[[1]], nrow=n_draw_here))*(x$income[1] + 2))
            return(output)
          }, mc.cores=numcores))
    }


    mat_YK = rbind(1, mat_YK)
    mat_YK = apply(mat_YK, 2, function(x) {x[which(is.nan(x))] = 0; return(x)})
    mat_YK = mat_YK[c(1,2,4),];

    mini_f = function(x_transform, return_long = FALSE) {
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
      print(summary(mat_M))

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

      output[[2]][x_transform[[2]][['sigma_delta']]] =  sum((d_output_1[['sigma_delta']] %*% d_moment[,1:nrow(mat_YK)])/ncol(mat_YK) ) + sum((d_output_2[['sigma_delta']] %*% d_moment[,(nrow(mat_YK) + 1):(2 * nrow(mat_YK))])/ncol(mat_YK) )
      output[[2]][x_transform[[2]][['sigma_gamma']]] =  sum((d_output_1[['sigma_gamma']] %*% d_moment[,1:nrow(mat_YK)])/ncol(mat_YK) ) + sum((d_output_2[['sigma_gamma']] %*% d_moment[,(nrow(mat_YK) + 1):(2 * nrow(mat_YK))])/ncol(mat_YK) )
      output[[2]][x_transform[[2]][['sigma_omega']]] =  sum((d_output_1[['sigma_omega']] %*% d_moment[,1:nrow(mat_YK)])/ncol(mat_YK) ) + sum((d_output_2[['sigma_omega']] %*% d_moment[,(nrow(mat_YK) + 1):(2 * nrow(mat_YK))])/ncol(mat_YK))

      output[[3]] = do.call('cbind', output[[3]])
      output[[4]] = output_1

      if (return_long) { 
        output_mat = list();
        output_mat[[1]] = output_1
        output_mat[[2]] = output_2;
        output_mat[[3]] = do.call('c', lapply(moment_ineligible_hh_output, function(x) {output_mini = x[[1]] %*% t(x[[1]]); diag(output_mini) = x[[2]]; return(c(output_mini))}));
        return(output_mat)
      }
      return(output)
    }

    if (recompute_pref) {
      output_initial = mini_f(x_transform);
      mini_param = splitfngr::optim_share(initial_param_trial[x_transform[[2]][pref_list] %>% unlist()], function(x) {
        x_new = param_trial_here; 
        x_new[x_transform[[2]][pref_list] %>% unlist()] = x; 
        if (max(abs(x)) > 3) {
          return(list(NA, rep(NA, length(x))))
        }
        output = mini_f(transform_param(x_new, return_index = TRUE)); 
        print(paste0('output of aggregate_moment_pref  = ', output[[1]])); 
        print('x = '); print(x)
        deriv = output[[2]][x_transform[[2]][pref_list] %>% unlist()];
        print('derivative = '); print(deriv);
        # deriv[(length(x)-2):length(x)] = 0
        return(list(output[[1]], deriv))}, method='BFGS', control=list(reltol=1e-4))
      return(mini_param)
    } else {
        output_short = list()
        fid = function(index) {
            output_hh = tryCatch(household_draw_theta_kappa_Rdraw(hh_index=index, param=x_transform[[1]], n_draw_halton = n_draw_halton, n_draw_gauss = n_draw_gauss, sick_parameters, xi_parameters, short=FALSE, realized_sick = FALSE),error=function(e) e)
            f0 = moment_ineligible_hh(output_hh, x_transform[[1]])
            realized_sick = data_hh_list[[index]]$sick_dummy
            off_diag = function(x) {
              # diag(x) = 0; 
              x = x  
              return(x)
            }

            f_on_sqr = function(x) x^2;


            mat_mimj = f0[[1]] %*% t(f0[[1]]); diag(mat_mimj) = f0[[2]];
            output_0 =  sum((f0[[1]] - data_hh_list[[index]]$M_expense)^2) + sum(f_on_sqr(c(off_diag(mat_mimj) - off_diag(data_hh_list[[index]]$M_expense %*% t(data_hh_list[[index]]$M_expense))))) 

            deriv = rep(0, length(initial_param_trial));

            for (name_i in c('beta_theta', 'beta_theta_ind')) {
              for (i in 1:length(f0[[1]])) {
                numerical_derivative = rep(0, length(f0[[1]])); numerical_derivative[i] = tol;
                
                if (name_i == 'beta_theta') {
                  output_hh_i = household_draw_theta_kappa_Rdraw(index, x_transform[[1]], n_draw_halton, 10, sick_parameters, xi_parameters, u_lowerbar = -1, derivative_r_threshold = FALSE, derivative=TRUE, numerical_derivative = numerical_derivative, short=FALSE, option_derivative = name_i, realized_sick = FALSE);
                  f1 = moment_ineligible_hh(output_hh_i, x_transform[[1]])
                  mat_mimj_1 = f1[[1]] %*% t(f1[[1]]); diag(mat_mimj_1) = f1[[2]];
                  output_i =  sum((f1[[1]] - data_hh_list[[index]]$M_expense)^2) + sum(f_on_sqr(c(off_diag(mat_mimj_1) - off_diag(data_hh_list[[index]]$M_expense %*% t(data_hh_list[[index]]$M_expense))))) 
                  deriv[x_transform[[2]][[name_i]]] = deriv[x_transform[[2]][[name_i]]] + (output_i - output_0)/tol * output_hh$X_ind_year[i,];
                } else {
                  
                  output_hh_i = household_draw_theta_kappa_Rdraw(index, x_transform[[1]], n_draw_halton, 10, sick_parameters, xi_parameters, u_lowerbar = -1, derivative_r_threshold = FALSE, derivative=TRUE, numerical_derivative = numerical_derivative, short=FALSE, option_derivative = name_i, realized_sick = FALSE);
                  f1 = moment_ineligible_hh(output_hh_i, x_transform[[1]])
                  mat_mimj_1 = f1[[1]] %*% t(f1[[1]]); diag(mat_mimj_1) = f1[[2]];
                  output_i =  sum((f1[[1]] - data_hh_list[[index]]$M_expense)^2) + sum(f_on_sqr(c(off_diag(mat_mimj_1) - off_diag(data_hh_list[[index]]$M_expense %*% t(data_hh_list[[index]]$M_expense)))))  
                  deriv[x_transform[[2]][[name_i]]] = deriv[x_transform[[2]][[name_i]]] + (output_i - output_0)/tol * output_hh$X_ind[i,];
                }
              }
            }


            for (name_i in c('sigma_theta', 'sigma_thetabar')) { 
              x_transform_i = x_transform; x_transform_i[[1]][[name_i]] = x_transform_i[[1]][[name_i]] + tol
              output_hh_i = household_draw_theta_kappa_Rdraw(index, x_transform_i[[1]], n_draw_halton, 10, sick_parameters, xi_parameters, u_lowerbar = -1, derivative_r_threshold = FALSE, derivative=FALSE, numerical_derivative = NA, short=FALSE, realized_sick = FALSE);
                f1 = moment_ineligible_hh(output_hh_i, x_transform_i[[1]])
                mat_mimj_1 = f1[[1]] %*% t(f1[[1]]); diag(mat_mimj_1) = f1[[2]];
                output_i =  sum((f1[[1]] - data_hh_list[[index]]$M_expense)^2) + sum(f_on_sqr(c(off_diag(mat_mimj_1) - off_diag(data_hh_list[[index]]$M_expense %*% t(data_hh_list[[index]]$M_expense))))) 
                deriv[x_transform[[2]][[name_i]]] = (output_i - output_0)/tol;
            }
            return(list(deriv, output_0, f0[[1]]))
          }
        if (Sys.info()[['sysname']] == 'Windows') {
          clusterExport(cl, c('x_transform', 'n_draw_halton','fid'),envir=environment())
          deriv_list = parLapply(cl, sample_identify_pref,fid)
        } else {
          deriv_list = mclapply(sample_identify_pref,fid, mc.cores = numcores)
        }

        output_short[[2]] = colSums(do.call('rbind', lapply(deriv_list, function(x) x[[1]])))
        output_short[[1]] = sum(do.call('c', lapply(deriv_list, function(x) x[[2]])))
        print(summary(do.call('c', lapply(deriv_list, function(x) x[[3]]))))
        
        return(list(output_short[[1]], output_short[[2]]))
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
    
    moment_realized_expense_val = do.call('c', lapply(moment_realized_expense, function(x) x[[1]])) %>% sum
    # print(paste0('moment from theta = ', moment_realized_expense_val))
    moment_realized_expense_deriv = rep(0, length(param_trial));
    moment_realized_expense_deriv[x_transform[[2]]$beta_theta] = do.call('rbind', lapply(moment_realized_expense, function(x) x[[2]]$beta_theta)) %>% colSums
    moment_realized_expense_deriv[x_transform[[2]]$beta_theta_ind] = do.call('rbind', lapply(moment_realized_expense, function(x) x[[2]]$beta_theta_ind)) %>% colSums
    moment_realized_expense_deriv[x_transform[[2]]$sigma_thetabar] = do.call('c', lapply(moment_realized_expense, function(x) x[[2]]$sigma_thetabar)) %>% sum
     moment_realized_expense_deriv[x_transform[[2]]$sigma_theta] = do.call('c', lapply(moment_realized_expense, function(x) x[[2]]$sigma_theta)) %>% sum

    return(list(moment_realized_expense_val, moment_realized_expense_deriv))
  }


  # index_theta_only = c(x_transform[[2]]$beta_theta[1], x_transform[[2]]$beta_theta_ind[1], x_transform[[2]]$sigma_thetabar, x_transform[[2]]$sigma_theta)
  index_pref_only = x_transform[[2]][pref_list] %>% unlist();

  save_output = list()
  iter = 1; 
  current_output = 0; 

  optim_f =  function(x_pref_theta, include_r=TRUE, include_pref=TRUE, include_theta_raw=FALSE) {
      print('START OF A NEW ITERATION'); 

      param_trial_here = initial_param_trial; 
      param_trial_here[c(index_theta_only)] = x_pref_theta
      param_trial_here[x_transform[[2]]$beta_theta_ind[1]] = exp(param_trial_here[x_transform[[2]]$beta_theta_ind[1]]);

      x_transform = transform_param(param_trial_here, return_index = TRUE)

      output_theta = list(0, rep(0, length(initial_param_trial)))
      if (include_theta_raw) {
        output_theta = aggregate_moment_theta(x_transform)
      }
      if (is.infinite(output_theta[[1]])) {
        return(list(NA, rep(NA, length(x_pref_theta))))
      }
      
      if (include_pref) {
        optim_pref = aggregate_moment_pref(transform_param(param_trial_here, return_index = TRUE), recompute_pref=TRUE); 

        param_trial_here[index_pref_only] = optim_pref$par; 

        x_transform = transform_param(param_trial_here, return_index = TRUE)

        # output_theta = aggregate_moment_theta(x_transform)
        output_theta = list(0, rep(0, length(initial_param_trial)))
        if (is.nan(output_theta[[1]])) {
          return(list(NA, rep(NA, length(c(index_theta_only, index_pref_only)))))
        }
        pref_moment = aggregate_moment_pref(x_transform)
      } else {
        pref_moment = list(0, rep(0, length(initial_param_trial)))
      } 
        
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
        fx_r = function(x_transform_fxr, derivative=FALSE, silent=TRUE) {
          sd_r = exp(x_transform_fxr[[1]]$sigma_r);
          correlation = x_transform_fxr[[1]]$correlation
          mean_vec = rep(X_hh_theta_r %*% x_transform_fxr[[1]]$beta_r, each = n_halton_at_r) + correlation * hh_theta

          denominator = (pnorm(-(5 - mean_vec)/sd_r, lower.tail=FALSE) - pnorm(-(0 - mean_vec)/sd_r, lower.tail=FALSE)) + 1e-20

          output = -sum(colMeans(matrix((pnorm(-(root_r[,2] - mean_vec)/sd_r, lower.tail=FALSE) - pnorm(-(root_r[,1] - mean_vec)/sd_r, lower.tail=FALSE))/denominator, nrow = n_halton_at_r))) 
          
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
          output = try(fx_r(transform_param(x_with_new_r, return_index=TRUE), silent=silent, derivative=FALSE))
          if("try-error" %in% class(output)) {
            return(NA)
          }
          return(output[[1]])
        }, control=list(maxit=1e4), method='BFGS') 

        param_trial_inner_r[c(x_transform[[2]]$beta_r, x_transform[[2]]$sigma_r, x_transform[[2]]$correlation)] = optim_r$par

        x_transform = transform_param(param_trial_inner_r,return_index=TRUE); 

        # compute derivative 
        compute_deriv = function(mini_data_index) {
          data_hh_i = data_hh_list[[mini_data_index]]

          if (data_hh_i$Year[1] == 2006) {
            constraint_function = function(x) {x_new = x; x_new[-c(1, length(x))] = -Inf; return(x_new) }
          } else {
            constraint_function = function(x) x
          }
          output_hh = counterfactual_household_draw_theta_kappa_Rdraw(mini_data_index, param = x_transform[[1]], n_draw_halton, n_draw_gauss, sick_parameters, xi_parameters, u_lowerbar = -1, policy_mat_hh = policy_mat[[mini_data_index]], seed_number = 1, constraint_function = constraint_function, compute_WTP = FALSE)
          
          realized_vol_sts = data_hh_i$Vol_sts
          m_observed = data_hh_i$M_expense
          HHsize = length(m_observed)
          X_ind = var_ind(data_hh_i)
          X_hh = var_hh(data_hh_i)
          X_ind_year = cbind(var_ind(data_hh_i), data_hh_i$Year == 2004, data_hh_i$Year == 2006, data_hh_i$Year == 2010, data_hh_i$Year == 2012)
          f_output = function(output_hh) {
            return(list(sum((output_hh$m - m_observed)^2) + sum((c(output_hh$m %*% t(output_hh$m)) - c(m_observed %*% t(m_observed)))^2), sum(output_hh$vol_sts_counterfactual),output_hh$m))
          }

          output_0 = f_output(output_hh)
          deriv = list(rep(0, length(initial_param_trial)), rep(0, length(initial_param_trial)), rep(0, length(initial_param_trial)))

          for (name_i in c('beta_theta', 'beta_theta_ind')) {
            for (i in 1:HHsize) {
              numerical_derivative = rep(0, HHsize); numerical_derivative[i] = tol;
              if (name_i == 'beta_theta') {
                output_hh_i = counterfactual_household_draw_theta_kappa_Rdraw(mini_data_index, param = x_transform[[1]], n_draw_halton, n_draw_gauss, sick_parameters, xi_parameters, u_lowerbar = -1, policy_mat_hh = policy_mat[[mini_data_index]], seed_number = 1, constraint_function = constraint_function, compute_WTP = FALSE, derivative=TRUE, numerical_derivative = numerical_derivative, option_derivative = name_i)
                output_i = f_output(output_hh_i);
                deriv[[1]][x_transform[[2]][[name_i]]] = deriv[[1]][x_transform[[2]][[name_i]]] + (output_i[[1]] - output_0[[1]])/tol * X_ind_year[i,];
                deriv[[2]][x_transform[[2]][[name_i]]] = deriv[[2]][x_transform[[2]][[name_i]]] + (output_i[[2]] - output_0[[2]])/tol * X_ind_year[i,];
                deriv[[3]][x_transform[[2]][[name_i]]] = deriv[[3]][x_transform[[2]][[name_i]]] + (output_i[[3]] - output_0[[3]])/tol * X_ind_year[i,];
              } else {
                output_hh_i = counterfactual_household_draw_theta_kappa_Rdraw(mini_data_index, param = x_transform[[1]], n_draw_halton, n_draw_gauss, sick_parameters, xi_parameters, u_lowerbar = -1, policy_mat_hh = policy_mat[[mini_data_index]], seed_number = 1, constraint_function = constraint_function, compute_WTP = FALSE, derivative=TRUE, numerical_derivative = numerical_derivative, option_derivative = name_i);
                output_i = f_output(output_hh_i);
                deriv[[1]][x_transform[[2]][[name_i]]] = deriv[[1]][x_transform[[2]][[name_i]]] + (output_i[[1]] - output_0[[1]])/tol * X_ind[i,];
                deriv[[2]][x_transform[[2]][[name_i]]] = deriv[[2]][x_transform[[2]][[name_i]]] + (output_i[[2]] - output_0[[2]])/tol * X_ind[i,];
                deriv[[3]][x_transform[[2]][[name_i]]] = deriv[[3]][x_transform[[2]][[name_i]]] + (output_i[[3]] - output_0[[3]])/tol * X_ind[i,];
              }
            }
          }


          for (name_i in c('sigma_theta', 'sigma_thetabar')) { 
            x_transform_i = x_transform; x_transform_i[[1]][[name_i]]=x_transform_i[[1]][[name_i]]+ tol
            output_hh_i = counterfactual_household_draw_theta_kappa_Rdraw(mini_data_index, x_transform_i[[1]], n_draw_halton, n_draw_gauss, sick_parameters, xi_parameters, u_lowerbar = -1, policy_mat_hh = policy_mat[[mini_data_index]], seed_number = 1, constraint_function = constraint_function, compute_WTP = FALSE, derivative=FALSE, numerical_derivative = NA, option_derivative = NA);
                output_i = f_output(output_hh_i);
            deriv[[1]][x_transform[[2]][[name_i]]] = (output_i[[1]] - output_0[[1]])/tol;
            deriv[[2]][x_transform[[2]][[name_i]]] = (output_i[[2]] - output_0[[2]])/tol;
            deriv[[3]][x_transform[[2]][[name_i]]] = (output_i[[3]] - output_0[[3]])/tol;
          }
          return(list(output_0[[1]], deriv, cbind(output_hh$m, data_hh_i$M_expense, output_hh$vol_sts_counterfactual, data_hh_i$Vol_sts)))
        }

        if (Sys.info()[['sysname']] == 'Windows') {
          clusterExport(cl, c('x_transform', 'sick_parameters', 'xi_parameters'),envir=environment())
          deriv_list = parLapply(cl, sample_r_theta, compute_deriv)
        } else {
          deriv_list = mclapply(sample_r_theta, compute_deriv, mc.cores=numcores)
        }
        output_1 = do.call('c', lapply(deriv_list, function(x) x[[1]])) %>% sum
        deriv_1 = do.call('rbind', lapply(deriv_list, function(x) x[[2]][[1]])) %>% colSums()
        output_2_sqrt = (do.call('c', lapply(deriv_list, function(x) x[[3]][,3])) %>% mean - do.call('c', lapply(deriv_list, function(x) x[[3]][,4])) %>% mean)
        output_3_sqrt = (do.call('c', lapply(deriv_list, function(x) x[[3]][,1])) %>% mean - do.call('c', lapply(deriv_list, function(x) x[[3]][,2])) %>% mean)
        deriv_2 = 2 * output_2_sqrt * (do.call('rbind', lapply(deriv_list, function(x) x[[2]][[2]])) %>% colMeans())
        deriv_3 = 2 * output_2_sqrt * (do.call('rbind', lapply(deriv_list, function(x) x[[2]][[3]])) %>% colMeans())
        # output = output_1 + output_2_sqrt^2*1000; 
        # deriv = deriv_1 + deriv_2*1000;
        output = (output_3_sqrt)^2 * length(sample_r_theta) + (output_2_sqrt)^2 * length(sample_r_theta);
        deriv = output_3_sqrt * 2 * deriv_3 * length(sample_r_theta) + output_2_sqrt * 2 * deriv_2 * length(sample_r_theta);
        print('------VOLUNTARY HH UNDER BD---------')
        print(summary(do.call('rbind', lapply(deriv_list, function(x) x[[3]]))))

        if (Sys.info()[['sysname']] == 'Windows') {
          clusterExport(cl, c('x_transform', 'sick_parameters', 'xi_parameters'),envir=environment())
          deriv_list = parLapply(cl, out_sample_index, compute_deriv)
        } else {
          deriv_list = mclapply(out_sample_index, compute_deriv, mc.cores=numcores)
        }
        output_1 = do.call('c', lapply(deriv_list, function(x) x[[1]])) %>% sum
        deriv_1 = do.call('rbind', lapply(deriv_list, function(x) x[[2]][[1]])) %>% colSums()
        output_2_sqrt = (do.call('c', lapply(deriv_list, function(x) x[[3]][,3])) %>% mean - do.call('c', lapply(deriv_list, function(x) x[[3]][,4])) %>% mean)
        output_3_sqrt = (do.call('c', lapply(deriv_list, function(x) x[[3]][,1])) %>% mean - do.call('c', lapply(deriv_list, function(x) x[[3]][,2])) %>% mean)
        deriv_2 = 2 * output_2_sqrt * (do.call('rbind', lapply(deriv_list, function(x) x[[2]][[2]])) %>% colMeans())
        deriv_3 = 2 * output_2_sqrt * (do.call('rbind', lapply(deriv_list, function(x) x[[2]][[3]])) %>% colMeans())
        # output = output_1 + output_2_sqrt^2*1000; 
        # deriv = deriv_1 + deriv_2*1000;
        output = output + (output_3_sqrt)^2 * length(out_sample_index) + (output_2_sqrt)^2 * length(out_sample_index);
        deriv = deriv + output_3_sqrt * 2 * deriv_3 * length(out_sample_index) + output_2_sqrt * 2 * deriv_2 * length(out_sample_index);
        print('------VOLUNTARY HH UNDER PURE BUNDLING---------')
        print(summary(do.call('rbind', lapply(deriv_list, function(x) x[[3]]))))
        # if (max(do.call('rbind', lapply(deriv_list, function(x) x[[3]]))[,3]) == 0) {
        #   message('0 insured --- eliminate this value')
        #   output = NA 
        #   deriv = rep(NA, length(initial_param_trial))
        # }
        return(list(output, deriv, param_trial_inner_r))
      }

      if (include_r) {
        output_r = output_wrt_r(x_transform, silent = TRUE);   
      } else {
        output_r = list(0, rep(0, length(initial_param_trial)), param_trial_here)
      }
      
      print('----------FINAL OUTPUT---------');
      print('at x value = '); print(output_r[[3]])
      print(pref_moment[[1]] + output_theta[[1]] + output_r[[1]])
      print(pref_moment[[1]])
      print(output_r[[1]])
      print('--------------------')
      param_trial_here = output_r[[3]];
      save_output[[iter]] <<- param_trial_here
      current_output <<- pref_moment[[1]] + output_theta[[1]] + output_r[[1]]
      iter <<- iter + 1; 

      pref_moment[[2]][x_transform[[2]]$beta_theta_ind[1]] = pref_moment[[2]][x_transform[[2]]$beta_theta_ind[1]]  * exp(param_trial_here[x_transform[[2]]$beta_theta_ind[1]])
      output_theta[[2]][x_transform[[2]]$beta_theta_ind[1]] = output_theta[[2]][x_transform[[2]]$beta_theta_ind[1]]  * exp(param_trial_here[x_transform[[2]]$beta_theta_ind[1]])
      output_r[[2]][x_transform[[2]]$beta_theta_ind[1]] = output_r[[2]][x_transform[[2]]$beta_theta_ind[1]]  * exp(param_trial_here[x_transform[[2]]$beta_theta_ind[1]])

      total_output = pref_moment[[1]] + output_r[[1]]; 
      total_derivative = pref_moment[[2]]  + output_r[[2]];

      total_derivative[x_transform[[2]]$beta_theta_ind] = total_derivative[x_transform[[2]]$beta_theta_ind] * exp(x_transform[[1]]$beta_theta_ind[1])
      print('pref_moment[[2]]'); print(pref_moment[[2]][index_theta_only])
      print('output_theta[[2]]'); print(output_theta[[2]][index_theta_only])
      print('output_r[[2]]'); print( output_r[[2]][index_theta_only])
      return(list(total_output, total_derivative[index_theta_only], param_trial_here))
  }


  index_theta_only = x_transform[[2]][list_theta_var] %>% unlist(); 

  optim_pref_theta = splitfngr::optim_share(initial_param_trial[index_theta_only], function(x) {
    print('x at optim_pref_theta = '); print(x)
    if (max(abs(x)) > 10) {
      return(list(NA, rep(NA, length(index_theta_only))))
    }
    output = try(optim_f(x, include_r = TRUE, include_pref=TRUE))
    if("try-error" %in% class(output)) {
      print(output)
      return(list(NA, rep(NA, length(index_theta_only))))
    } else {
      return(output[1:2])
    }
  }, control=list(maxit=1e3,reltol=1e-4), method='BFGS')

  param_trial = optim_f(optim_pref_theta$par)[[3]]; 



  message('computing final param_trial')
   
  param_final = list(); 
  param_final$other = param_trial; 
  param_final$xi = xi_parameters;
  param_final$sick = sick_parameters

  param = param_final 
  transform_param_final = transform_param(param_final$other)

  fit_sample = sample(Vol_HH_list_index, 1000)

  for (seed_number in c(1:1)) {
    if (Sys.info()[['sysname']] == 'Windows') {
      clusterExport(cl, c('transform_param_final', 'param','counterfactual_household_draw_theta_kappa_Rdraw'))
      mini_fit_values = parLapply(cl, c(fit_sample), function(id) {
        output = counterfactual_household_draw_theta_kappa_Rdraw(id, transform_param_final, n_draw_halton, n_draw_gauss, param$sick, param$xi, u_lowerbar = -1, policy_mat_hh = policy_mat[[id]], seed_number = seed_number + job_index, constraint_function = function(x) x)
        output = as.data.frame(output)
        output$Y = data_hh_list[[id]]$Income; 
        output$m_observed = data_hh_list[[id]]$M_expense; 
        output$id = id;
        output$HHsize_s = data_hh_list[[id]]$HHsize_s; 
        return(output)
      })
    } else {
      mini_fit_values = mclapply(c(fit_sample), function(id) {
      output = tryCatch(counterfactual_household_draw_theta_kappa_Rdraw(id, transform_param_final, n_draw_halton, n_draw_gauss, param$sick, param$xi, u_lowerbar = -1, policy_mat_hh = policy_mat[[id]], seed_number = seed_number + job_index, constraint_function = function(x) x), error=function(e) e)
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
  # plot = gridExtra::grid.arrange(plot_1, plot_2, nrow=1)

  if (dir.exists('../../householdbundling_estimate')) {
    saveRDS(param_final, file=paste0('../../householdbundling_estimate/estimate_',job_index,'.rds'))
  } else {
    dir.create('../../householdbundling_estimate') 
    saveRDS(param_final, file=paste0('../../householdbundling_estimate/estimate_',job_index,'.rds'))
  }

  if (Sys.info()[['sysname']] == 'Windows') {
    stopCluster(cl)
  }
}