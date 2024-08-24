args = commandArgs(trailingOnly=TRUE)
if (length(args)<3) { 
  numcores = 4;  
  contraction_variance = 1;
  within_hh_heterogeneity = list(omega=TRUE, gamma=TRUE, delta=TRUE, theta_bar=TRUE);
  file_name_label = paste0('contraction_variance_', contraction_variance, '_full_heterogeneity')
  job_index_iter = 1:100; 
} else {
  numcores = as.numeric(args[1]); 
  contraction_variance = as.numeric(args[2]);
  if (as.numeric(args[3]) == 1) {
  	within_hh_heterogeneity = list(omega=TRUE, gamma=TRUE, delta=TRUE, theta_bar=TRUE);
  	file_name_label = paste0('contraction_variance_', contraction_variance, '_full_heterogeneity')
  } else {
  	within_hh_heterogeneity = list(omega=FALSE, gamma=FALSE, delta=FALSE, theta_bar=TRUE);
  	file_name_label = paste0('contraction_variance_', contraction_variance, 'no_pref_heterogeneity')
  }
  job_index_iter = as.numeric(args[3])
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
  if (nrow(data) > nrow(data %>% filter(Bef_sts + Com_sts + Std_w_ins == 1)) & !(data$Year[1] == 2006 & data$HHsize_s[1] > 1)) {
    return(hh_index);
  }
  else {
    return(NA); 
  }
}) %>% unlist(); 

Vol_HH_list_index = Vol_HH_list_index[!(is.na(Vol_HH_list_index))] 

Com_HH_list_index = lapply(1:length(data_hh_list), function(hh_index) {
  data = data_hh_list[[hh_index]]; 
  if (0 == nrow(data %>% filter(Bef_sts + Com_sts + Std_w_ins == 0))) {
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


list_hh_2012 = unlist(lapply(1:length(data_hh_list), function(hh_index) ifelse(data_hh_list[[hh_index]]$Year[1] == 2012 & data_hh_list[[hh_index]]$HHsize_s[1] == 2, hh_index, NA)))
list_hh_2012 = list_hh_2012[which(!(is.na(list_hh_2012)))] %>% sample(1000)

list_p1 = seq(0, 0.06, by = 0.005) 
list_p2 = seq(0, 1, by = 0.05)


if (Sys.info()[['sysname']] == 'Windows') {
  numcores = numcores; 
  cl = makeCluster(numcores);
  clusterEvalQ(cl, library('tidyverse'))
  clusterEvalQ(cl, library('familyenrollment'))
  clusterExport(cl,c('Vol_HH_list_index', 'Com_HH_list_index', 'out_sample_index'))
}

job_index_list = as.numeric(gsub("\\D", "", list.files('../../householdbundling_estimate/'))) %>% unique()

iter_list = c(1:1);

if (Sys.info()[['sysname']] == 'Windows') {
  clusterExport(cl,c('iter_list'))
}

bd_output = list()
pb_output = list()

counterfactual_values_bd = list();
counterfactual_values_pb = list();

for (job_index in job_index_list[job_index_iter]) {
	if (file.exists(paste0('../../Obj_for_manuscript/bd_output', file_name_label,'_',job_index,'.rds'))) {
		next;
	}
	if (file.exists(paste0('../../householdbundling_estimate/estimate_',job_index,'.rds'))) {
		param_final <- readRDS(paste0('../../householdbundling_estimate/estimate_',job_index,'.rds'))
		transform_param_final = transform_param(param_final$other)

		counterfactual_premium = function(premium_param, type, id) {
			MW = (Income_net_premium[[id]][1] - Income_net_premium[[id]][2]) * 4/3 / 0.06 # 2012 premium is at 4.5\% MW for the first member
			if (type == 'bundle discount') {
				hhsize_s = length(Income_net_premium[[id]]) - 1; 
				income_vec = Income_net_premium[[id]][1] - c(0, premium_param[[hhsize_s]][1:hhsize_s]) * MW;
			} else if (type == 'pure bundling') {
				income_vec = rep(Income_net_premium[[id]][1], hhsize_s); 
				income_vec[hhsize_s] = Income_net_premium[[id]][1] - premium_param[[hhsize_s]] * MW; 
			} else {
				income_vec =  Income_net_premium[[id]][1] - premium_param[[hhsize_s]] * c(0:hhsize_s) * MW;
			}
			return(income_vec)
		}

		if (Sys.info()[['sysname']] == 'Windows') {
		  clusterExport(cl,c('counterfactual_premium'))
		}

		bd_prem = list()
		pb_prem = list()
		prem_id = 0; 
		for (i1 in list_p1) {
			for (i2 in list_p2) {
				prem_id = prem_id + 1; 
				bd_prem[[prem_id]] = list(); 
				bd_prem[[prem_id]][[2]] = c(i1, i1 * (1 + i2))
				pb_prem[[prem_id]] = list(); 
				pb_prem[[prem_id]][[2]] = c(0, i1 * (1 + i2))
				if (i1 == 0) {
					break
				}
			}
		} 

		bd_output[[job_index]] = list(); 
		bd_output[[job_index]]$budget = NULL; 
		bd_output[[job_index]]$surplus = NULL

		pb_output[[job_index]] = list(); 
		pb_output[[job_index]]$budget = NULL; 
		pb_output[[job_index]]$surplus = NULL

		print_index = 0; 


		counterfactual_values_bd[[job_index]] = list();
		counterfactual_values_pb[[job_index]] = list();

		f_prem = function(prem, constraint_function, iter_list, job_index, return_long = FALSE) {
			f_id = function(id) {
				income_vec = counterfactual_premium(prem, 'bundle discount', id)
				output = do.call('rbind', lapply(iter_list, function(iter) {
					output = as.data.frame(counterfactual_household_draw_theta_kappa_Rdraw(id, transform_param_final, 100, 10, param_final$sick, param_final$xi, u_lowerbar = -1, policy_mat_hh = policy_mat[[id]], seed_number = iter, constraint_function = constraint_function, income_vec = income_vec, contraction_variance = contraction_variance))
					output$iter = iter; 
					output$Y = data_hh_list[[id]]$Income; 
					output$m_observed = data_hh_list[[id]]$M_expense; 
					output$fit_type = ifelse(id %in% Vol_HH_list_index, 2, ifelse(id %in% Com_HH_list_index, 1, 3))
					output$id = id
					output$premium_optimal = income_vec[1] - income_vec[sum(output$vol_sts_counterfactual) + 1];
					output$HHsize_s = data_hh_list[[id]]$HHsize_s; 
					output$job_index = job_index
					return(output) 
					}))
				return(output)}
			if (Sys.info()[['sysname']] == 'Windows') {
				clusterExport(cl, c('prem', 'param_final', 'transform_param_final', 'job_index', 'constraint_function', 'contraction_variance'))
				counterfactual_values_bd = parLapply(cl, c(list_hh_2012), f_id)
			} else {
				counterfactual_values_bd = mclapply(c(list_hh_2012), f_id, mc.cores=numcores)
			}
			counterfactual_values_bd = do.call('rbind', counterfactual_values_bd)
			bd_output = list()
			bd_output$surplus = counterfactual_values_bd  %>% filter(Com_sts + Bef_sts + Std_w_ins == 0) %>% group_by(id, iter) %>% slice(1) %>% pull(wtp) %>% sum(na.rm = TRUE) - sum(counterfactual_values_bd %>% group_by(id, iter) %>% slice(1) %>% pull(premium_optimal)) 
			bd_output$budget = sum(counterfactual_values_bd  %>% filter(Com_sts + Bef_sts + Std_w_ins == 0) %>% group_by(id, iter) %>% slice(1) %>% pull(premium_optimal)) - (counterfactual_values_bd %>% filter(Com_sts + Bef_sts + Std_w_ins == 0) %>% pull(cost_to_insurance) %>% sum())
			bd_output$cost = (counterfactual_values_bd %>% filter(Com_sts + Bef_sts + Std_w_ins == 0) %>% pull(cost_to_insurance) %>% sum())
			bd_output$demand = sum(counterfactual_values_bd  %>% filter(Com_sts + Bef_sts + Std_w_ins == 0) %>% group_by(id, iter) %>% mutate(N_vol = sum(vol_sts_counterfactual)) %>% slice(1) %>% ungroup() %>% pull(N_vol))/sum(counterfactual_values_bd  %>% filter(Com_sts + Bef_sts + Std_w_ins == 0) %>% group_by(id, iter) %>% mutate(N_vol = sum(vol_sts_counterfactual)) %>% slice(1) %>% ungroup() %>% pull(HHsize_s))
			bd_output$p1 = prem[[2]][1];
			bd_output$p2 = prem[[2]][1];
			if (return_long) {
				return(counterfactual_values_bd)
			} else {
				return(bd_output)
			}
		}

		constraint_function = function(x) x;
		bd_output[[job_index]] = do.call('rbind', lapply(bd_prem, function(prem) f_prem(prem, constraint_function, iter_list, job_index, return_long = FALSE)));

		constraint_function = function(x) {x_new = x; x_new[-c(1, length(x))] = -Inf; return(x_new)}
		pb_output[[job_index]] = do.call('rbind', lapply(bd_prem, function(prem) f_prem(prem, constraint_function, iter_list, job_index, return_long = FALSE)));

		bd_output[[job_index]] = as.data.frame(bd_output[[job_index]])
		pb_output[[job_index]] = as.data.frame(pb_output[[job_index]])

		names(bd_output[[job_index]]) = c('surplus', 'budget', 'cost', 'demand', 'p1', 'p2')
		names(pb_output[[job_index]]) = c('surplus', 'budget', 'cost', 'demand', 'p1', 'p2')

		bd_output[[job_index]] = bd_output[[job_index]] %>% mutate(p_ratio = p2/p1); 
		bd_output[[job_index]]$p_ratio[is.nan(bd_output[[job_index]]$p_ratio)] = 1; 

		saveRDS(bd_output[[job_index]], file = paste0('../../Obj_for_manuscript/bd_output', file_name_label,'_',job_index,'.rds'))
		saveRDS(pb_output[[job_index]], file = paste0('../../Obj_for_manuscript/pb_output', file_name_label,'_',job_index,'.rds'))
	}
}




