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
devtools::install(upgrade='never')
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

n_draw_halton = 10;

if (Sys.info()[['sysname']] == 'Windows') {
  numcores = numcores; 
  cl = makeCluster(numcores);
  clusterEvalQ(cl, library('tidyverse'))
  clusterEvalQ(cl, library('familyenrollment'))
  clusterExport(cl,c('Vol_HH_list_index', 'Com_HH_list_index', 'out_sample_index', 'n_draw_halton'))
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

		f_prem = function(prem, constraint_function, iter_list, job_index, return_long = FALSE, penalty = 0) {
			f_id = function(id) {
				income_vec = counterfactual_premium(prem, 'bundle discount', id)
				output = do.call('rbind', lapply(iter_list, function(iter) {
					output = as.data.frame(counterfactual_household_draw_theta_kappa_Rdraw(id, transform_param_final, n_draw_halton, 10, param_final$sick, param_final$xi, u_lowerbar = -1, policy_mat_hh = policy_mat[[id]], seed_number = iter + job_index, constraint_function = constraint_function, income_vec = income_vec, contraction_variance = contraction_variance, penalty = penalty))
					output$iter = iter; 
					output$Y = data_hh_list[[id]]$Income; 
					output$m_observed = data_hh_list[[id]]$M_expense; 
					output$fit_type = ifelse(id %in% Vol_HH_list_index, 2, ifelse(id %in% Com_HH_list_index, 1, 3))
					output$id = id
					output$premium_optimal = income_vec[1] - income_vec[sum(output$vol_sts_counterfactual) + 1];
					output$HHsize_s = data_hh_list[[id]]$HHsize_s; 
					output$eligible = (data_hh_list[[id]]$Com_sts + data_hh_list[[id]]$Bef_sts + data_hh_list[[id]]$Std_w_ins) == 0  
					output$job_index = job_index
					return(as.data.frame(output))
					}))
				return(output)}
			if (Sys.info()[['sysname']] == 'Windows') {
				clusterExport(cl, c('prem', 'param_final', 'transform_param_final', 'job_index', 'constraint_function', 'contraction_variance'))
				counterfactual_values_bd = parLapply(cl, c(list_hh_2012), f_id)
			} else {
				counterfactual_values_bd = mclapply(c(list_hh_2012), f_id, mc.cores=numcores) 
			}
			counterfactual_values_bd = do.call('rbind', counterfactual_values_bd) %>% as.data.frame()
			bd_output = list()
			bd_output$surplus = counterfactual_values_bd  %>% filter(eligible) %>% group_by(id, iter) %>% slice(1) %>% pull(wtp) %>% sum(na.rm = TRUE) - sum(counterfactual_values_bd %>% group_by(id, iter) %>% slice(1) %>% pull(premium_optimal)) 
			bd_output$budget = sum(counterfactual_values_bd  %>% filter(eligible) %>% group_by(id, iter) %>% slice(1) %>% pull(premium_optimal)) - (counterfactual_values_bd %>% filter(eligible) %>% pull(cost_to_insurance) %>% sum())
			bd_output$cost = (counterfactual_values_bd %>% filter(eligible) %>% pull(cost_to_insurance) %>% sum())
			bd_output$demand = sum(counterfactual_values_bd  %>% filter(eligible) %>% group_by(id, iter) %>% mutate(N_vol = sum(vol_sts_counterfactual)) %>% slice(1) %>% ungroup() %>% pull(N_vol))/sum(counterfactual_values_bd  %>% filter(eligible) %>% group_by(id, iter) %>% mutate(N_vol = sum(vol_sts_counterfactual)) %>% slice(1) %>% ungroup() %>% pull(HHsize_s))
			bd_output$p1 = prem[[2]][1];
			bd_output$p2 = prem[[2]][2];
			if (return_long) {
				return(counterfactual_values_bd)
			} else {
				return(bd_output)
			}
		}

		constraint_function = function(x) x;
		prem_2012 = list(); prem_2012[[2]] = c(0.045, 1.9 * 0.045);
		benchmark_2012 = f_prem(prem_2012, constraint_function, iter_list, job_index, return_long = FALSE)

		compute_best_ps = function(prem_normalized, benchmark = benchmark_2012, constraint_function = function(x) x, type = c('ip','pb','bd'), return_output = FALSE, penalty = 0) {
			if (type == 'ip') {
				if (penalty == 0) {
					prem = list(); 
					prem[[2]] = c(NA, NA); 
					prem[[2]][1] = exp(prem_normalized)/(1 + exp(prem_normalized)) * 0.06;
					prem[[2]][2] = 2 * prem[[2]][1]; 
				} else {
					prem = list(); 
					prem[[2]] = c(NA, NA); 
					prem[[2]][1] = exp(prem_normalized[1])/(1 + exp(prem_normalized[1])) * 0.06;
					prem[[2]][2] = 2 * prem[[2]][1]; 
					penalty = exp(prem_normalized[2])/(1 + exp(prem_normalized[2])) * 0.06;
				}
					
			} else if (type == 'bd') {
				prem = list(); 
				prem[[2]] = c(NA, NA); 
				prem[[2]][1] = exp(prem_normalized[1])/(1 + exp(prem_normalized[1])) * 0.06;
				prem[[2]][2] = (1 + exp(prem_normalized[2])/(1 + exp(prem_normalized[2]))) * prem[[2]][1]; 
			} else {
				prem = list(); 
				prem[[2]] = c(NA, NA); 
				prem[[2]][1] = 0 ;
				prem[[2]][2] = exp(prem_normalized)/(1 + exp(prem_normalized)) * 0.12; 
			}
			
			if (return_output) {
				output = f_prem(prem, constraint_function, iter_list, job_index, return_long = TRUE, penalty = penalty); 
			} else {
				output = f_prem(prem, constraint_function, iter_list, job_index, return_long = FALSE, penalty = penalty); 
				if (output$budget > benchmark$budget) {
					print('at'); print(prem[[2]])
					print(output)
				}
				output_summary = - output$surplus + 100 * (output$budget < benchmark$budget);
				return(output_summary)
			}
		}

		optimal_bd = optim(c(log(1/(1 - 0.75) - 1), log(1/(1 - 0.9) - 1)), function(prem_normalized) compute_best_ps(prem_normalized, benchmark_2012, function(x) x, type='bd') , control = list(reltol = 1e-2))
		output_optimal_bd = compute_best_ps(optimal_bd$par, benchmark_2012, function(x) x, type='bd', return_output=TRUE)

		constraint_function = function(x) {x_new = x; x_new[-c(1, length(x))] = -Inf; return(x_new)}
		optimal_pb = optimize(function(prem_normalized) compute_best_ps(prem_normalized, benchmark_2012, constraint_function, type='pb') , c(-5,5), tol = 1e-3)

		output_optimal_pb = compute_best_ps(optimal_pb$minimum, benchmark_2012, constraint_function, type='pb', return_output=TRUE)

		constraint_function = function(x) x;
		optimal_ip = optimize(function(prem_normalized) compute_best_ps(prem_normalized, benchmark_2012, constraint_function, type='ip'), c(0, 4), tol=1e-3)
		output_optimal_ip = compute_best_ps(optimal_ip$minimum, benchmark_2012, constraint_function, type='ip', return_output=TRUE)

		constraint_function = function(x) {x_new = x; x_new[-length(x)] = -Inf; return(x_new)}
		optimal_mandate = optimize(function(prem_normalized) compute_best_ps(prem_normalized, benchmark_2012, constraint_function, type='ip'), c(-4, 4), tol=1e-3)
		output_optimal_mandate = compute_best_ps(optimal_mandate$minimum, benchmark_2012, constraint_function, type='ip', return_output=TRUE)

		constraint_function = function(x) x;
		optimal_penalty = optim(c(log(1/(1 - 0.75) - 1), log(1/(1 - 0.9) - 1)), function(prem_normalized) compute_best_ps(prem_normalized, benchmark_2012, constraint_function, type='ip',penalty = 1), control=list(reltol = 1e-2))
		output_optimal_penalty = compute_best_ps(optimal_penalty$par, benchmark_2012, constraint_function, type='ip', return_output=TRUE, penalty = 1)

		constraint_function = function(x) x;
		optimal_risk_rating = optim(c(log(1/(1 - 0.75) - 1), 1), function(prem_normalized)  compute_best_ps(prem_normalized, benchmark_2012, constraint_function, type='ip',penalty = 1), control=list(reltol = 1e-2))

		# compute across ip prices 
		ip_price = seq(0,1,0.005) * 0.06 
		optimal_bd_all = list();
		benchmark_ip_all = list();
		index = 0
		for (p1 in ip_price) {
			index = index + 1; 
			prem = list(); 
			prem[[2]] = c(p1, p1 * 2)
			benchmark_ip_all[[index]] = f_prem(prem, constraint_function, iter_list, job_index, return_long = FALSE)
			optimal_bd_all[[index]] = optim(c(log(1/(1 - p1/0.06) - 1), 4.5), function(prem_normalized) compute_best_ps(prem_normalized, benchmark_ip_all[[index]], function(x) x, type='bd') , control = list(reltol = 1e-2))
			optimal_pb_all[[index]] = optim(c(log(1/(1 - p1/0.06) - 1), 4.5), function(prem_normalized) compute_best_ps(prem_normalized, benchmark_ip_all[[index]], function(x) x, type='pb') , control = list(reltol = 1e-2))
		}

		save(c('optimal_bd', 'output_optimal_bd', 'optimal_pb', 'output_optimal_pb', 'optimal_ip', 'output_optimal_ip', 'optimal_bd_all', 'optimal_pb_all', 'output_optimal_mandate', 'output_optimal_penalty', 'optimal_penalty', 'optimal_mandate'), file = paste0('../../Obj_for_manuscript/counterfactual', file_name_label,'_',job_index,'.rdata'))

		stop()

		next 

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




