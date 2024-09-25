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

if (Sys.info()[['sysname']] == 'Windows') {
    numcores = 20;
} else if (Sys.info()[['sysname']] == 'Linux') {
  numcores = 12;
} else {
  numcores = 6;
}

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
list_hh_2012 = list_hh_2012[which(!(is.na(list_hh_2012)))] 

list_p1 = seq(0, 0.06, by = 0.005) 
list_p2 = seq(0, 1, by = 0.05)

n_draw_halton = 100;

if (Sys.info()[['sysname']] == 'Windows') {
  numcores = numcores; 
  cl = makeCluster(numcores);
  clusterEvalQ(cl, library('tidyverse'))
  clusterEvalQ(cl, library('familyenrollment'))
  clusterExport(cl,c('Vol_HH_list_index', 'Com_HH_list_index', 'out_sample_index', 'n_draw_halton'))
}

job_index_list = as.numeric(gsub("\\D", "", list.files('../../householdbundling_estimate/'))) %>% unique()
job_index_list = job_index_list[which(!(is.na(job_index_list)))]

iter_list = c(1:1);

if (Sys.info()[['sysname']] == 'Windows') {
  clusterExport(cl,c('iter_list'))
}

bd_output = list()
pb_output = list()

counterfactual_values_bd = list();
counterfactual_values_pb = list();


if (Sys.info()[['nodename']] == 'Anh-Macbook-3.local' | grepl("vpn", Sys.info()[['nodename']])  | grepl("macbook", tolower(Sys.info()[['nodename']])) ) {
  mini=FALSE
  numcores = 4; 
} else {
  mini=FALSE
}

if (mini) {
	list_hh_2012 = sample(list_hh_2012, 100)
}

f_counterfactual = function(file_name_label, within_hh_heterogeneity, job_index_iter, contraction_variance, compute_mandate = FALSE, compute_penalty = FALSE, compute_risk = FALSE, compute_frontier = FALSE, max_obj = 'ss', recompute_basic = TRUE) {
	for (job_index in job_index_list[job_index_iter]) {
		print(paste0('at job_index = ', job_index))
		if (file.exists(paste0('../../Obj_for_manuscript/counterfactual', file_name_label,'_',job_index,'.rdata'))) {
			next;
		}
		if (file.exists(paste0('../../householdbundling_estimate/estimate_',job_index,'.rds'))) {
			filename = paste0('../../Obj_for_manuscript/',file_name_label,'_',job_index,'.rdata')
			if (file.exists(filename)) {
				next;
			}
			if (job_index == min(job_index_list[job_index_iter], na.rm=TRUE)) {
				grid_search = FALSE
			} else {
				grid_search = FALSE
			}
			param_final <- readRDS(paste0('../../householdbundling_estimate/estimate_',job_index,'.rds'))
			transform_param_final = transform_param(param_final$other)

			counterfactual_premium = function(premium_param, id) {
				MW = (Income_net_premium[[id]][1] - Income_net_premium[[id]][2]) * 4/3 / 0.06 # 2012 premium is at 4.5\% MW for the first member

				hhsize_s = length(Income_net_premium[[id]]) - 1; 
				income_vec = Income_net_premium[[id]][1] - c(0, premium_param[1:hhsize_s]) * MW;
				
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

			f_prem = function(prem, constraint_function, iter_list, job_index, return_long = FALSE, penalty = 0, compute_WTP = FALSE, wtp_mat = NA, cost_mat = NA, price_var = NA, risk_discrimination = function(x) 0, risk_discrimination_dummy = FALSE) {
				f_id = function(id, prem, penalty = 0, compute_counterfactual_cost = TRUE) {
					income_vec = counterfactual_premium(prem, id)
					output = do.call('rbind', lapply(iter_list, function(iter) {
						output = as.data.frame(counterfactual_household_draw_theta_kappa_Rdraw(id, transform_param_final, n_draw_halton, 10, param_final$sick, param_final$xi, u_lowerbar = -1, policy_mat_hh = policy_mat[[id]], seed_number = iter + job_index, constraint_function = constraint_function, income_vec = income_vec, contraction_variance = contraction_variance, penalty = penalty, always_covered = FALSE, compute_wtp2 = TRUE, within_hh_heterogeneity = within_hh_heterogeneity, compute_WTP = TRUE, compute_counterfactual_cost = compute_counterfactual_cost, afford = FALSE, risk_discrimination = risk_discrimination, risk_discrimination_dummy = risk_discrimination_dummy))
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
				if (compute_WTP) {
					if (Sys.info()[['sysname']] == 'Windows') {
						clusterExport(cl, c('prem', 'param_final', 'transform_param_final', 'job_index', 'constraint_function', 'contraction_variance'))
						counterfactual_values_bd = parLapply(cl, c(list_hh_2012), function(id) f_id(id, prem = c(0.045, 0.045 * 1.9)))
					} else {
						counterfactual_values_bd = mclapply(c(list_hh_2012), function(id) f_id(id, prem = c(0.045, 0.045 * 1.9)), mc.cores=numcores) 
					}
					wtp_frame = lapply(counterfactual_values_bd, function(mini_data) list(wtp = c(mini_data$wtp_uni[which(mini_data$eligible)], mini_data$wtp_2[1]), cost_to_insurance = mini_data$cost_to_insurance[which(mini_data$eligible)]))
					wtp_mat = (do.call('rbind', lapply(wtp_frame, function(mini_data) mini_data$wtp)))
					cost_mat = (do.call('rbind', lapply(wtp_frame, function(mini_data) mini_data$cost)))
					counterfactual_values_bd = do.call('rbind', counterfactual_values_bd) %>% as.data.frame()
					return(list(counterfactual_values_bd = counterfactual_values_bd, wtp_mat = wtp_mat, cost_mat = cost_mat))
				} else {
					outcome_at_premium_f = function(prem, constraint_function, price_var = NA) {
						output = do.call('rbind', lapply(1:nrow(wtp_mat), function(row_index) {
							if (is.na(price_var)) {
								prem_hh = counterfactual_premium(prem, id = list_hh_2012[row_index])
								optim_index = which.max(constraint_function(c(0 - 2 * penalty, wtp_mat[row_index,1] - penalty, wtp_mat[row_index,2] - penalty, wtp_mat[row_index,3])) - c(0, prem_hh[1] - prem_hh[2], prem_hh[1] - prem_hh[2], prem_hh[1] - prem_hh[3]))
								optim_index = max(optim_index)
								ss = c(0, wtp_mat[row_index,1], wtp_mat[row_index,2], wtp_mat[row_index, 3])[optim_index] - c(0,cost_mat[row_index,], sum(cost_mat[row_index,]))[optim_index]
								cs = c(0 - 2 * penalty, wtp_mat[row_index,1] - penalty, wtp_mat[row_index,2] - penalty, wtp_mat[row_index, 3])[optim_index] - c(0, prem_hh[1] - prem_hh[2], prem_hh[1] - prem_hh[2], prem_hh[1] - prem_hh[3])[optim_index]
								budget =  c(2 * penalty, prem_hh[1] - prem_hh[2] + penalty, prem_hh[1] - prem_hh[2] + penalty, prem_hh[1] - prem_hh[3])[optim_index] - c(0,cost_mat[row_index,],sum(cost_mat[row_index,]))[optim_index]
								vol_sts = (optim_index == 1) * c(0,0) + (optim_index == 2) * c(1,0) + (optim_index == 3) * c(0,1) + (optim_index == 4) * c(1,1)
							} else {
								# must be a binary variable
								data_hh_i_eligible = output_base$counterfactual_values_bd %>% filter(id == list_hh_2012[row_index]) %>% filter(Bef_sts + Com_sts + Std_w_ins == 0)
								prem_hh = NULL; 
								type_mem = data_hh_i_eligible[[price_var]]
								prem_hh[1] = type_mem[1] * prem[1] + (1 - type_mem[1]) * prem[2];
								prem_hh[2] = type_mem[2] * prem[1] + (1 - type_mem[2]) * prem[2];
								y = counterfactual_premium(c(prem_hh[1], 2 * prem_hh[1]), id = list_hh_2012[row_index]); 
								prem_hh[1] = y[1] - y[2]
								y = counterfactual_premium(c(prem_hh[2], 2 * prem_hh[2]), id = list_hh_2012[row_index]); 
								prem_hh[2] = y[1] - y[2]

								optim_index = which.max(constraint_function(c(0 - 2 * penalty,wtp_mat[row_index,1] - penalty, wtp_mat[row_index,2] - penalty, wtp_mat[row_index,3])) - c(0, prem_hh[1], prem_hh[2], prem_hh[1] + prem_hh[2]))
								optim_index = max(optim_index)
								ss = c(0, wtp_mat[row_index,1], wtp_mat[row_index,2], wtp_mat[row_index, 3])[optim_index] - c(0,cost_mat[row_index,],sum(cost_mat[row_index,]))[optim_index]
								cs = c(0 - 2 * penalty, wtp_mat[row_index,1] - penalty, wtp_mat[row_index,2] - penalty, wtp_mat[row_index, 3])[optim_index] -  c(0, prem_hh[1], prem_hh[2], prem_hh[2] + prem_hh[1])[optim_index]
								budget =  c(2 * penalty, prem_hh[1] + penalty, prem_hh[2] + penalty, prem_hh[2] + prem_hh[1])[optim_index] - c(0,cost_mat[row_index,],sum(cost_mat[row_index,]))[optim_index]
								vol_sts = (optim_index == 1) * c(0,0) + (optim_index == 2) * c(1,0) + (optim_index == 3) * c(0,1) + (optim_index == 4) * c(1,1)
							}
							return(c(ss, cs, budget, vol_sts))
						}))
						names(output) = c('ss','cs','budget','vol_sts_1','vol_sts_2')
						return(output)
					}

					if (return_long) {
						if (is.na(price_var)) {
							if (Sys.info()[['sysname']] == 'Windows') {
								clusterExport(cl, c('prem', 'param_final', 'transform_param_final', 'job_index', 'constraint_function', 'contraction_variance'))
								counterfactual_values_bd = parLapply(cl, c(list_hh_2012), function(id) f_id(id = id, prem = prem, penalty = penalty, compute_counterfactual_cost = FALSE))
							} else {
								counterfactual_values_bd = mclapply(c(list_hh_2012), function(id) f_id(id = id, prem = prem, penalty = penalty, compute_counterfactual_cost = FALSE), mc.cores=numcores) 
							}
							counterfactual_values_bd = do.call('rbind', counterfactual_values_bd) %>% as.data.frame()
							return(counterfactual_values_bd)
						} else {
							if (Sys.info()[['sysname']] == 'Windows') {
								clusterExport(cl, c('prem', 'param_final', 'transform_param_final', 'job_index', 'constraint_function', 'contraction_variance'))
								counterfactual_values_bd = parLapply(cl, c(list_hh_2012), function(id) f_id(id = id, prem = prem, penalty = penalty, risk_discrimination = risk_discrimination, compute_counterfactual_cost = FALSE))
							} else {
								counterfactual_values_bd = mclapply(c(list_hh_2012), function(id) f_id(id = id, prem = prem, penalty = penalty, risk_discrimination = risk_discrimination, compute_counterfactual_cost = FALSE), mc.cores=numcores) 
							}
							counterfactual_values_bd = do.call('rbind', counterfactual_values_bd) %>% as.data.frame()
							return(counterfactual_values_bd)
						}
					} else {
						output = colSums(outcome_at_premium_f(prem, constraint_function, price_var));
						names(output) = c('ss','cs', 'budget','vol_sts_1','vol_sts_2')
						return(output)
					}
				}		
			}

			prem_2012 = c(0.045, 1.9 * 0.045);
			output_base = f_prem(prem_2012, constraint_function = function(x) x, iter_list = 1, job_index, return_long = TRUE, penalty = 0, compute_WTP = TRUE, wtp_mat = NA, cost_mat = NA); 

			benchmark_2012 = f_prem(prem_2012, function(x) x, 1, job_index = job_index, wtp_mat = output_base$wtp_mat, cost_mat = output_base$cost_mat)
				benchmark_2012[['job_index']] = job_index; 

			if (recompute_basic) {
				optimal_bd = optim(c(log(1/(1 - 0.75) - 1), log(1/(1 - 0.9) - 1)), function(prem_normalized) {
					prem = NULL
					prem[1] = exp(prem_normalized[1])/(1 + exp(prem_normalized[1])) * 0.06;
					prem[2] = prem[1] * (1 + exp(prem_normalized[2])/(1 + exp(prem_normalized[2])));
					print('prem = '); print(prem)
					output = f_prem(prem, function(x) x, 1, job_index = job_index, wtp_mat = output_base$wtp_mat, cost_mat = output_base$cost_mat)
					final_output = -output[['ss']] + 1e4 * (output[['budget']] < benchmark_2012[['budget']])
					print(output)
					return(final_output)
				})

				optimal_bd_summary = f_prem(c(exp(optimal_bd$par[1])/(1 + exp(optimal_bd$par[1])) * 0.06,  exp(optimal_bd$par[1])/(1 + exp(optimal_bd$par[1])) * 0.06 * (1 + exp(optimal_bd$par[2])/(1 + exp(optimal_bd$par[2])))), function(x) x, 1, job_index = job_index, wtp_mat = output_base$wtp_mat, cost_mat = output_base$cost_mat)

				output_optimal_bd = f_prem(prem = c(exp(optimal_bd$par[1])/(1 + exp(optimal_bd$par[1])) * 0.06,  exp(optimal_bd$par[1])/(1 + exp(optimal_bd$par[1])) * 0.06 * (1 + exp(optimal_bd$par[2])/(1 + exp(optimal_bd$par[2])))), constraint_function = function(x) x, 1, job_index = job_index, wtp_mat = output_base$wtp_mat, cost_mat = output_base$cost_mat, return_long=TRUE)

				optimal_pb = optimize(function(prem_normalized){
					prem = c(0, exp(prem_normalized)/(1 + exp(prem_normalized)) * 0.12)
					output = f_prem(prem, function(x) c(x[1], -Inf, -Inf, x[4]), 1, job_index = job_index, wtp_mat = output_base$wtp_mat, cost_mat = output_base$cost_mat)
					print(output)
					print(-output[[max_obj]] + 1e2 * (output[['budget']] < benchmark_2012[['budget']]))
					return(-output[[max_obj]] + 1e2 * (output[['budget']] < benchmark_2012[['budget']]))
				}, c(0, 5), tol = 1e-4)

				optimal_pb_summary = f_prem(c(0, exp(optimal_pb$minimum)/(1 + exp(optimal_pb$minimum)) * 0.12), function(x) c(x[1], -Inf, -Inf, x[4]), 1, job_index = job_index, wtp_mat = output_base$wtp_mat, cost_mat = output_base$cost_mat)

				output_optimal_pb = f_prem(c(0, exp(optimal_pb$minimum)/(1 + exp(optimal_pb$minimum)) * 0.12), function(x) c(x[1], -Inf, x[3]), 1, job_index = job_index, wtp_mat = output_base$wtp_mat, cost_mat = output_base$cost_mat, return_long = TRUE)

				optimal_ip = optimize(function(prem_normalized){
					prem = c(exp(prem_normalized)/(1 + exp(prem_normalized)) * 0.06, exp(prem_normalized)/(1 + exp(prem_normalized)) * 0.06 * 2)
					output = f_prem(prem, function(x) x, 1, job_index = job_index, wtp_mat = output_base$wtp_mat, cost_mat = output_base$cost_mat)
					print(output)
					return(-output[[max_obj]] + 1e4 * (output[['budget']] < benchmark_2012[['budget']]))
				}, c(-5, 5), tol = 1e-3)

				optimal_ip_summary = f_prem(c(exp(optimal_ip$minimum)/(1 + exp(optimal_ip$minimum)) * 0.06, exp(optimal_ip$minimum)/(1 + exp(optimal_ip$minimum)) * 0.12), function(x) x, 1, job_index = job_index, wtp_mat = output_base$wtp_mat, cost_mat = output_base$cost_mat)

				output_optimal_ip = f_prem(c(exp(optimal_ip$minimum)/(1 + exp(optimal_ip$minimum)) * 0.06, exp(optimal_ip$minimum)/(1 + exp(optimal_ip$minimum)) * 0.12), function(x) x, 1, job_index = job_index, wtp_mat = output_base$wtp_mat, cost_mat = output_base$cost_mat, return_long = TRUE)
			}

			if (compute_mandate) {
				optimal_mandate = optimize(function(prem_normalized){
					prem = c(0, exp(prem_normalized)/(1 + exp(prem_normalized)) * 0.12)
					output = f_prem(prem, function(x) c(-Inf, -Inf, -Inf, x[4]), 1, job_index = job_index, wtp_mat = output_base$wtp_mat, cost_mat = output_base$cost_mat)
					return(-output[['cs']] + 1e4 * (output[['budget']] < benchmark_2012[['budget']]))
				}, c(-5, 5), tol = 1e-3)

				optimal_mandate_summary = f_prem(c(0, exp(optimal_mandate$minimum)/(1 + exp(optimal_mandate$minimum)) * 0.12), function(x) c(-Inf, -Inf, -Inf, x[4]), 1, job_index = job_index, wtp_mat = output_base$wtp_mat, cost_mat = output_base$cost_mat)
				output_optimal_mandate = f_prem(c(0, exp(optimal_mandate$minimum)/(1 + exp(optimal_mandate$minimum)) * 0.12), function(x) c(-Inf, -Inf, x[3]), 1, job_index = job_index, wtp_mat = output_base$wtp_mat, cost_mat = output_base$cost_mat, return_long = TRUE)
			} else {
				optimal_mandate = list()
				output_optimal_mandate = list()
				optimal_mandate_summary = list() 
			}


			if (compute_penalty) {
				optimal_penalty = optim(c(log(1/(1 - 0.75) - 1), log(1/(1 - 0.9) - 1)), function(prem_normalized) {
					prem = NULL
					prem[1] = exp(prem_normalized[1])/(1 + exp(prem_normalized[1])) * 0.06;
					prem[2] = prem[1] * 2;
					output = f_prem(prem, function(x) x, 1, job_index = job_index, wtp_mat = output_base$wtp_mat, cost_mat = output_base$cost_mat, penalty = exp(prem_normalized[2])/(1 + exp(prem_normalized[2])) * 0.06)
					return(-output[[max_obj]] + 1e4 * (output[['budget']] < benchmark_2012[['budget']]))
				})

				output_optimal_penalty = f_prem(c(exp(optimal_penalty$par[1])/(1 + exp(optimal_penalty$par[1])) * 0.06, exp(optimal_penalty$par[1])/(1 + exp(optimal_penalty$par[1])) * 0.12), function(x) x, 1, job_index = job_index, wtp_mat = output_base$wtp_mat, cost_mat = output_base$cost_mat, return_long=TRUE, penalty = exp(optimal_penalty$par[2])/(1 + exp(optimal_penalty$par[2])) * 0.06)

				optimal_penalty_summary = f_prem(c(exp(optimal_penalty$par[1])/(1 + exp(optimal_penalty$par[1])) * 0.06, exp(optimal_penalty$par[1])/(1 + exp(optimal_penalty$par[1])) * 0.12), function(x) x, 1, job_index = job_index, wtp_mat = output_base$wtp_mat, cost_mat = output_base$cost_mat, penalty = exp(optimal_penalty$par[2])/(1 + exp(optimal_penalty$par[2])) * 0.06)
			} else {
				optimal_penalty = list()
				output_optimal_penalty = list()
				optimal_penalty_summary = list()
			}
				

			if (compute_risk) {
				threshold = 0.5;
				output_base$counterfactual_values_bd$high_theta = output_base$counterfactual_values_bd$theta_bar > quantile(output_base$counterfactual_values_bd$theta_bar, threshold)

				optimal_risk = optim(c(0, optimal_ip$minimum), function(prem_normalized) {
					prem = NULL
					prem[2] = exp(prem_normalized[2])/(1 + exp(prem_normalized[2])) * 0.06;
					prem[1] = prem[2] * exp(prem_normalized[1])
					output = f_prem(prem, function(x) x, 1, job_index = job_index, wtp_mat = output_base$wtp_mat, cost_mat = output_base$cost_mat, price_var = 'high_theta')
					print(output)
					print(-output[[max_obj]] + 1e4 * (output[['budget']] < benchmark_2012[['budget']]))
					return(-output[[max_obj]] + 1e4 * (output[['budget']] < benchmark_2012[['budget']]))
				})

				prem = NULL
				prem[2] = exp(optimal_risk$par[2])/(1 + exp(optimal_risk$par[2])) * 0.06;
				prem[1] = prem[2] * exp(optimal_risk$par[1])

				optimal_risk_summary = f_prem(prem, function(x) x, 1, job_index = job_index, wtp_mat = output_base$wtp_mat, cost_mat = output_base$cost_mat, return_long=FALSE)

				risk_discrimination = function(x) (x > quantile(output_base$counterfactual_values_bd$theta_bar, threshold)) + (x <= quantile(output_base$counterfactual_values_bd$theta_bar, threshold)) * (exp(optimal_risk$par[2])/(1 + exp(optimal_risk$par[2])))/(exp(optimal_risk$par[1])/(1 + exp(optimal_risk$par[1])))

				output_optimal_risk = f_prem(prem, function(x) x, 1, job_index = job_index, wtp_mat = output_base$wtp_mat, cost_mat = output_base$cost_mat, return_long=TRUE, risk_discrimination = risk_discrimination, risk_discrimination_dummy = TRUE)


				optimal_risk_pb = optim(c(0, optimal_ip$minimum), function(prem_normalized) {
					prem = NULL
					prem[2] = exp(prem_normalized[2])/(1 + exp(prem_normalized[2])) * 0.06;
					prem[1] = prem[2] * exp(prem_normalized[1])
					output = f_prem(prem, function(x) c(x[1], -Inf, -Inf, x[4]), 1, job_index = job_index, wtp_mat = output_base$wtp_mat, cost_mat = output_base$cost_mat, price_var = 'high_theta')
					print(output)
					print(-output[[max_obj]] + 1e4 * (output[['budget']] < benchmark_2012[['budget']]))
					return(-output[[max_obj]] + 1e4 * (output[['budget']] < benchmark_2012[['budget']]))
				})

				prem = NULL
				prem[2] = exp(optimal_risk_pb$par[2])/(1 + exp(optimal_risk_pb$par[2])) * 0.06;
				prem[1] = prem[2] * exp(optimal_risk_pb$par[1])

				optimal_risk_pb_summary = f_prem(prem, function(x) c(x[1], -Inf, -Inf, x[4]), 1, job_index = job_index, wtp_mat = output_base$wtp_mat, cost_mat = output_base$cost_mat, return_long=FALSE)

				risk_discrimination = function(x) (x > quantile(output_base$counterfactual_values_bd$theta_bar, threshold)) + (x <= quantile(output_base$counterfactual_values_bd$theta_bar, threshold)) * (exp(optimal_risk_pb$par[2])/(1 + exp(optimal_risk_pb$par[2])))/(exp(optimal_risk_pb$par[1])/(1 + exp(optimal_risk_pb$par[1])))

				output_optimal_risk_pb = f_prem(prem, function(x) c(x[1], -Inf, x[3]), 1, job_index = job_index, wtp_mat = output_base$wtp_mat, cost_mat = output_base$cost_mat, return_long=TRUE, risk_discrimination = risk_discrimination, risk_discrimination_dummy = TRUE)

			} else {
				optimal_risk = list()
				optimal_risk_pb = list()
				output_optimal_risk = list()
				output_optimal_risk_pb = list()
				optimal_risk_summary = list()
				optimal_risk_summary_pb = list()
			}
			
			
			if (compute_frontier) {
				ip_price = seq(0.005,0.055,0.01)
				optimal_bd_all = list();
				optimal_pb_all = list()
				benchmark_ip_all = list();
				optimal_bd_all_output = list();
				optimal_pb_all_output = list()
				benchmark_ip_all_output = list();
				index = 0
				for (p1 in ip_price) {
					index = index + 1; 
					prem = list(); 
					prem[[2]] = c(p1, p1 * 2)
					optimal_pb_all_output[[index]] = f_prem(c(p1, p1 * 2), 1, job_index = job_index, wtp_mat = output_base$wtp_mat, cost_mat = output_base$cost_mat, return_long = TRUE, constraint_function = function(x) c(x[1], -Inf, x[3]))
					benchmark_budget = f_prem(c(p1, p1 *2), function(x) c(x[1], -Inf, -Inf, x[4]), 1, job_index = job_index, wtp_mat = output_base$wtp_mat, cost_mat = output_base$cost_mat)
					print('benchmark_budget = '); print(benchmark_budget)
					optimal_bd_all[[index]] = optim(c(log(p1 * 2), 0), function(prem_normalized) {
						prem = NULL
						prem[1] = exp(prem_normalized[1])
						prem[2] = prem[1] * (1 + exp(prem_normalized[2])/(1 + exp(prem_normalized[2])));
						print(prem)
						output = f_prem(prem, function(x) x, 1, job_index = job_index, wtp_mat = output_base$wtp_mat, cost_mat = output_base$cost_mat)
						print(output)
						return(-output[[max_obj]] + 1e2 * (output[['budget']] < benchmark_budget[['budget']]))
					})
					
					# benchmark_ip_all_output[[index]] = f_prem(c(p1, p1 * 2), function(x) x, 1, job_index = job_index, wtp_mat = output_base$wtp_mat, cost_mat = output_base$cost_mat, return_long = TRUE)
					benchmark_ip_all_output[[index]]  = list()
					optimal_bd_all_output[[index]] = f_prem(c(optimal_bd_all[[index]]$par[1] |> exp(), exp(optimal_bd_all[[index]]$par[1]) * (1 + exp(optimal_bd_all[[index]]$par[2])/(1 + exp(optimal_bd_all[[index]]$par[2])))), function(x) x, 1, job_index = job_index, wtp_mat = output_base$wtp_mat, cost_mat = output_base$cost_mat, return_long = TRUE)
					
				}				
			} else {
				ip_price = seq(0.005,0.055,0.01)
				optimal_bd_all = list();
				optimal_pb_all = list()
				benchmark_ip_all = list();
				optimal_bd_all_output = list();
				optimal_pb_all_output = list()
				benchmark_ip_all_output = list();
			}
			# compute across ip prices 
			save(list = c('output_optimal_risk_pb', 'output_optimal_risk', 'optimal_risk','optimal_risk_pb', 'optimal_risk_summary', 'optimal_risk_summary_pb'), file=filename)
			# save(list = c('benchmark_ip_all', 'optimal_bd_all', 'benchmark_ip_all_output', 'optimal_bd_all_output', 'optimal_pb_all_output', 'output_optimal_risk_pb', 'output_optimal_risk', 'optimal_risk','optimal_risk_pb', 'optimal_risk_summary', 'optimal_risk_summary_pb', 'optimal_mandate','output_optimal_mandate','optimal_mandate_summary','optimal_ip_summary','optimal_ip', 'output_optimal_ip','optimal_bd_summary','optimal_bd', 'output_optimal_bd','optimal_pb_summary','optimal_pb', 'output_optimal_pb', 'output_base', 'benchmark_2012'), file = filename)
		}
	}
}

for (i in 1:100) {
f_counterfactual(file_name_label = 'cf_cs_variance_1_full_heterogeneity',within_hh_heterogeneity = list(omega=TRUE, gamma=TRUE, delta=TRUE, theta_bar=TRUE), job_index_iter = c(i), contraction_variance = 1, compute_mandate = TRUE, compute_penalty = TRUE, compute_risk = TRUE, compute_frontier = TRUE, max_obj = 'cs', recompute_basic = TRUE);

f_counterfactual('cf_cs_variance_1_nomh',within_hh_heterogeneity = list(omega=TRUE, gamma=FALSE, delta=TRUE, theta_bar=TRUE), job_index_iter = c(i), contraction_variance = 1, max_obj = 'cs', compute_mandate = TRUE, recompute_basic = TRUE);

f_counterfactual('cf_cs_variance_1_norisk',within_hh_heterogeneity = list(omega=TRUE, gamma=FALSE, delta=TRUE, theta_bar=FALSE), job_index_iter = c(i), contraction_variance = 1, max_obj = 'cs', compute_mandate = TRUE, recompute_basic = TRUE);

# f_counterfactual('cf_cs_variance_1_risk',within_hh_heterogeneity = list(omega=TRUE, gamma=TRUE, delta=TRUE, theta_bar=TRUE), job_index_iter = c(i), contraction_variance = 1, max_obj = 'cs', compute_risk = TRUE, recompute_basic = FALSE);

f_counterfactual('cf_cs_variance_0_full_heterogeneity',within_hh_heterogeneity = list(omega=TRUE, gamma=TRUE, delta=TRUE, theta_bar=TRUE), job_index_iter = c(i), contraction_variance = 0, max_obj = 'cs', compute_mandate = TRUE);
}


