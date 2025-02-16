args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) { 
  numcores = 4;  
  job_index_iter = 1:100; 
} else {
  numcores = as.numeric(args[2]); 
  job_index_iter = as.numeric(args[1])
}

library(knitr)
library(tidyverse)
library(lfe)
library(sandwich)
library(plm)
library(stargazer)
library(parallel)
library(randomForest)
library(randtoolbox)
library(Hmisc)

# setwd('./familyenrollment')
# devtools::install(upgrade='never')
library(familyenrollment)

Vol_HH_list_index = lapply(1:length(data_hh_list), function(hh_index) {
  data = data_hh_list[[hh_index]]; 
  if (nrow(data) > nrow(data %>% filter(Bef_sts + Com_sts + Std_w_ins == 1)) & data$Year[1] == 2012 & !(data$Year[1] == 2006 & data$HHsize_s[1] > 1)) {
    return(hh_index);
  }
  else {
    return(NA); 
  }
}) %>% unlist(); 

Vol_HH_list_index = Vol_HH_list_index[!(is.na(Vol_HH_list_index))] 

Com_HH_list_index = lapply(1:length(data_hh_list), function(hh_index) {
  data = data_hh_list[[hh_index]]; 
  if (0 == nrow(data %>% filter(Bef_sts + Com_sts + Std_w_ins == 0)) & data$Year[1] == 2012) {
    return(hh_index);
  }
  else {
    return(NA); 
  }
}) %>% unlist(); 

out_sample_index = lapply(1:length(data_hh_list), function(hh_index) {
  data = data_hh_list[[hh_index]]; 
  if (nrow(data) > nrow(data %>% filter(Bef_sts + Com_sts + Std_w_ins == 1)) & (data$Year[1] == 2006 & data$HHsize_s[1] > 1)) {
    return(hh_index);
  }
  else {
    return(NA); 
  }
}) %>% unlist(); 

out_sample_index = out_sample_index[!(is.na(out_sample_index))] 

Com_HH_list_index = Com_HH_list_index[!(is.na(Com_HH_list_index))] 

if (Sys.info()[['sysname']] == 'Windows') {
  numcores = 10; 
  cl = makeCluster(numcores);
  clusterEvalQ(cl, library('tidyverse'))
  clusterEvalQ(cl, library('familyenrollment'))
  clusterExport(cl,c('Vol_HH_list_index', 'Com_HH_list_index', 'out_sample_index'))
}
	
job_index_list = as.numeric(gsub("\\D", "", list.files('../../householdbundling_estimate/'))) %>% unique()
# job_index_list = c('1724060072')

iter_list = c(1:1)

if (Sys.info()[['sysname']] == 'Windows') {
  clusterExport(cl,c('iter_list'))
}

if (!(dir.exists('../../Obj_for_manuscript/'))) {
	dir.create('../../Obj_for_manuscript/')
}

# if (!file.exists('../../Obj_for_manuscript/all_param_final_2.rds')) {
	all_param_final = list()
	job_index_normalized = 1;
	for (job_index in job_index_list) {
		if (file.exists(paste0('../../householdbundling_estimate/estimate_',job_index,'.rds'))) {
			param_final <- readRDS(paste0('../../householdbundling_estimate/estimate_',job_index,'.rds'))
			all_param_final[[job_index_normalized]] = transform_param(param_final$other)
			all_param_final[[job_index_normalized]]$beta_sick = param_final$sick$par
			job_index_normalized = job_index_normalized + 1;
		}
	}

	name_summary_param = names(all_param_final[[1]]); 
	all_param_final_2 = lapply(name_summary_param, function(name_x) do.call('rbind', lapply(all_param_final, function(x) x[[name_x]])))
	names(all_param_final_2) = name_summary_param
	all_param_final_2$beta_theta = as.data.frame(all_param_final_2$beta_theta); names(all_param_final_2$beta_theta) = c(names(var_ind(data_hh_list[[1]]) %>% as.data.frame), 'year_dummy_2004', 'year_dummy_2006', 'year_dummy_2010', 'year_dummy_2012')
	all_param_final_2$beta_theta_ind = as.data.frame(all_param_final_2$beta_theta_ind); names(all_param_final_2$beta_theta_ind) = c(names(var_ind(data_hh_list[[1]]) %>% as.data.frame))
	all_param_final_2$beta_gamma = as.data.frame(all_param_final_2$beta_gamma); names(all_param_final_2$beta_gamma) = c(names(var_ind(data_hh_list[[1]]) %>% as.data.frame))
	all_param_final_2$beta_sick = as.data.frame(all_param_final_2$beta_sick); names(all_param_final_2$beta_sick) = c(names(var_ind(data_hh_list[[1]]) %>% as.data.frame))
	all_param_final_2$beta_delta = as.data.frame(all_param_final_2$beta_delta); names(all_param_final_2$beta_delta) = c(names(var_ind(data_hh_list[[1]]) %>% as.data.frame))
	all_param_final_2$beta_r = as.data.frame(all_param_final_2$beta_r); names(all_param_final_2$beta_r) = c(names(var_hh(data_hh_list[[1]]) %>% as.data.frame))
	all_param_final_2$beta_omega = as.data.frame(all_param_final_2$beta_omega); names(all_param_final_2$beta_omega) = c(names(var_hh(data_hh_list[[1]]) %>% as.data.frame))
	for (name_i in c('sigma_theta', 'sigma_thetabar', 'sigma_r', 'sigma_delta', 'sigma_gamma', 'sigma_omega')) {
		all_param_final_2[[name_i]] = exp(all_param_final_2[[name_i]])
	}

	saveRDS(all_param_final_2, file='../../Obj_for_manuscript/all_param_final_2.rds')
# }



# if (!file.exists('../../Obj_for_manuscript/correlation_values.rds')) {
	message('computing correlation values')
	# Compute covariance matrices 
	correlation_values = lapply(all_param_final, function(x) {
		output = do.call('rbind',lapply(data_hh_list[c(Vol_HH_list_index, Com_HH_list_index)], 
			function(data_id) {
				x_W = var_ind(data_id) %*% x$beta_theta_ind; 
				if (nrow(data_id) > 1) {
					cov_matrix = x_W %*% t(x_W) + diag(rep(exp(x$sigma_thetabar)^2, nrow(data_id)))
				} else {
					return(NA)
				}

				# diag(cor_matrix) = NA; 
				cov_matrix_vec = c(cov_matrix)
				cov_relationship = matrix(NA, nrow=5, ncol=5)
				for (i in 1:5) {
					for (j in 1:i) {
						cov_relationship[i,j] = mean(c(cov_matrix[which(data_id$relationship == i), which(data_id$relationship == j)]), na.rm=TRUE)
					}
				}
				cov_relationship[which(is.nan(cov_relationship))] = NA; 
				return(c(cov_relationship)); 
			}))
		return(matrix(apply(output, 2, function(y) mean(y, na.rm=TRUE)), nrow=5))
	})

	saveRDS(correlation_values, file='../../Obj_for_manuscript/correlation_values.rds')
# }



fit_values_all = NULL; 
no_heterogeneity_values_all = NULL; 
n_draw_halton = 10;
iter=1

for (job_index in job_index_list[job_index_iter]) {
	print(paste0('computing at index = ', job_index))
	if (file.exists(paste0('../../Obj_for_manuscript/no_heterogeneity_values',job_index,'.rds'))) {
		next 
	}
	if (file.exists(paste0('../../householdbundling_estimate/estimate_',job_index,'.rds'))) {
		param_final <- readRDS(paste0('../../householdbundling_estimate/estimate_',job_index,'.rds'))
		transform_param_final = transform_param(param_final$other)
		fid = function(id) {
				output = counterfactual_household_draw_theta_kappa_Rdraw(id, transform_param_final, n_draw_halton, 10, param_final$sick, param_final$xi, u_lowerbar = -1, policy_mat_hh = policy_mat[[id]], seed_number = iter + job_index, constraint_function = function(x) x, within_hh_heterogeneity = list(omega=FALSE, gamma=FALSE, delta=FALSE, theta_bar=FALSE))
				output = as.data.frame(output)
				output$Y = data_hh_list[[id]]$Income; 
				output$m_observed = data_hh_list[[id]]$M_expense; 
				output$fit_type = ifelse(id %in% Vol_HH_list_index, 2, ifelse(id %in% Com_HH_list_index, 1, 3))
				output$id = id; 
				output$relationship = data_hh_list[[id]]$relationship
				output$iter = iter;
				return(output)
				}
		if (Sys.info()[['sysname']] == 'Windows') {
		  clusterExport(cl, c('transform_param_final', 'param_final','counterfactual_household_draw_theta_kappa_Rdraw', 'n_draw_halton'))
			no_heterogeneity_values = parLapply(cl, c(Vol_HH_list_index), fid)
		} else {
			no_heterogeneity_values = mclapply(c(Vol_HH_list_index), fid, mc.cores = numcores)
		}
		
		no_heterogeneity_values = do.call('rbind', no_heterogeneity_values)
		no_heterogeneity_values = as.data.frame(no_heterogeneity_values)
		no_heterogeneity_values$job_index = job_index; 

		saveRDS(no_heterogeneity_values, file=paste0('../../Obj_for_manuscript/no_heterogeneity_values', job_index, '.rds'))
	}
}