args = commandArgs(trailingOnly=TRUE)
if (length(args)<3) { 
  numcores = 4;  
  contraction_variance = 1;
  within_hh_heterogeneity = list(omega=TRUE, gamma=TRUE, delta=TRUE, theta_bar=TRUE);
  file_name_label = paste0('contraction_variance_', contraction_variance, '_full_heterogeneity')
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

benchmark = readRDS('../../Obj_for_manuscript/fit_values.rds')

list_hh_2012 = unlist(lapply(1:length(data_hh_list), function(hh_index) ifelse(data_hh_list[[hh_index]]$Year[1] == 2012 & data_hh_list[[hh_index]]$HHsize_s[1] == 2, hh_index, NA)))
list_hh_2012 = list_hh_2012[which(!(is.na(list_hh_2012)))]

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

for (job_index in job_index_list) {
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

		for (prem in bd_prem) {
			print_index = print_index + 1; print(paste0('computing bundle discount premium index ', print_index));
			f_id = function(id) {
				income_vec = counterfactual_premium(prem, 'bundle discount', id)
				output = do.call('rbind', lapply(iter_list, function(iter) {
					output = as.data.frame(counterfactual_household_draw_theta_kappa_Rdraw(id, transform_param_final, 100, 10, param_final$sick, param_final$xi, u_lowerbar = -1, policy_mat_hh = policy_mat[[id]], seed_number = iter, constraint_function = function(x) x, income_vec = income_vec))
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
				clusterExport(cl, c('prem', 'param_final', 'transform_param_final'))
				counterfactual_values_bd[[job_index]][[print_index]] = parLapply(cl, c(list_hh_2012), f_id)
			} else {
				counterfactual_values_bd[[job_index]][[print_index]] = mclapply(c(list_hh_2012), f_id, mc.cores=numcores)
			}
			
			counterfactual_values_bd[[job_index]][[print_index]] = do.call('rbind', counterfactual_values_bd[[job_index]][[print_index]])
			bd_output[[job_index]]$surplus = c(bd_output[[job_index]]$surplus, counterfactual_values_bd[[job_index]][[print_index]]  %>% filter(Com_sts + Bef_sts + Std_w_ins == 0) %>% group_by(id, iter) %>% slice(1) %>% pull(wtp) %>% sum(na.rm = TRUE) - sum(counterfactual_values_bd[[job_index]][[print_index]] %>% group_by(id, iter) %>% slice(1) %>% pull(premium_optimal))) 
			bd_output[[job_index]]$budget = c(bd_output[[job_index]]$budget, sum(counterfactual_values_bd[[job_index]][[print_index]]  %>% filter(Com_sts + Bef_sts + Std_w_ins == 0) %>% group_by(id, iter) %>% slice(1) %>% pull(premium_optimal)) - (counterfactual_values_bd[[job_index]][[print_index]] %>% filter(Com_sts + Bef_sts + Std_w_ins == 0) %>% pull(cost_to_insurance) %>% sum()))
			bd_output[[job_index]]$demand = c(bd_output[[job_index]]$demand, sum(counterfactual_values_bd[[job_index]][[print_index]]  %>% filter(Com_sts + Bef_sts + Std_w_ins == 0) %>% group_by(id, iter) %>% mutate(N_vol = sum(vol_sts_counterfactual)) %>% slice(1) %>% ungroup() %>% pull(N_vol))/sum(counterfactual_values_bd[[job_index]][[print_index]]  %>% filter(Com_sts + Bef_sts + Std_w_ins == 0) %>% group_by(id, iter) %>% mutate(N_vol = sum(vol_sts_counterfactual)) %>% slice(1) %>% ungroup() %>% pull(HHsize_s)))
		}

		print_index = 0; 
		for (prem in pb_prem) {
			print_index = print_index + 1; print(paste0('computing pure bundling premium index ', print_index));
			fid = function(id) {
				income_vec = counterfactual_premium(prem, 'bundle discount', id)
				output = do.call('rbind', lapply(iter_list, function(iter) {
				output = counterfactual_household_draw_theta_kappa_Rdraw(id, transform_param_final, 100, 10, param_final$sick, param_final$xi, u_lowerbar = -1, policy_mat_hh = policy_mat[[id]], seed_number = iter, constraint_function = function(x) {x_new = x; x_new[-c(1, length(x))] = -Inf; return(x_new)}, income_vec = income_vec)
				output = as.data.frame(output)
				output$Y = data_hh_list[[id]]$Income; 
				output$m_observed = data_hh_list[[id]]$M_expense; 
				output$fit_type = ifelse(id %in% Vol_HH_list_index, 2, ifelse(id %in% Com_HH_list_index, 1, 3))
				output$id = id
				output$premium_optimal = income_vec[1] - income_vec[sum(output$vol_sts_counterfactual) + 1]
				output$iter = iter; 
				output$HHsize_s = data_hh_list[[id]]$HHsize_s; 
				output$job_index = job_index
				return(output)
				}))	
				return(output)}
			if (Sys.info()[['sysname']] == 'Windows') {
				clusterExport(cl, c('prem', 'param_final','transform_param_final'))
				counterfactual_values_pb[[job_index]][[print_index]] = parLapply(cl, c(list_hh_2012), f_id)
			} else {
				counterfactual_values_pb[[job_index]][[print_index]] = mclapply(c(list_hh_2012), f_id, mc.cores=numcores)
			}
			counterfactual_values_pb[[job_index]][[print_index]] = do.call('rbind', counterfactual_values_pb[[job_index]][[print_index]]) 
			pb_output[[job_index]]$surplus = c(pb_output[[job_index]]$surplus, counterfactual_values_pb[[job_index]][[print_index]] %>% filter(Com_sts + Bef_sts + Std_w_ins == 0) %>% group_by(id, iter) %>% slice(1) %>% ungroup() %>% pull(wtp) %>% sum(na.rm=TRUE) - sum(counterfactual_values_pb[[job_index]][[print_index]] %>% group_by(id, iter) %>% slice(1) %>% ungroup() %>% pull(premium_optimal)))
			pb_output[[job_index]]$budget = c(pb_output[[job_index]]$budget, sum(counterfactual_values_pb[[job_index]][[print_index]] %>% filter(Com_sts + Bef_sts + Std_w_ins == 0) %>% group_by(id, iter) %>% slice(1) %>% ungroup() %>% pull(premium_optimal)) - (counterfactual_values_pb[[job_index]][[print_index]] %>% filter(Com_sts + Bef_sts + Std_w_ins == 0) %>% pull(cost_to_insurance) %>% sum()))
			pb_output[[job_index]]$demand = c(pb_output[[job_index]]$demand, sum(counterfactual_values_pb[[job_index]][[print_index]] %>% filter(Com_sts + Bef_sts + Std_w_ins == 0) %>% group_by(id, iter) %>% mutate(N_vol = sum(vol_sts_counterfactual)) %>% slice(1) %>% ungroup() %>% pull(N_vol))/sum(counterfactual_values_pb[[job_index]][[print_index]]  %>% filter(Com_sts + Bef_sts + Std_w_ins == 0) %>% group_by(id, iter) %>% mutate(N_vol = sum(vol_sts_counterfactual)) %>% slice(1) %>% ungroup() %>% pull(HHsize_s)))
		}	
	}
}


bd_output = do.call('rbind', lapply(bd_output[job_index_list], function(x) as.data.frame(x)))
pb_output = do.call('rbind', lapply(pb_output[job_index_list], function(x) as.data.frame(x)))
pb_premium_vector = as.data.frame(do.call('rbind',lapply(pb_prem, function(x) x[[2]]))); names(pb_premium_vector) = c('p1','p2')
bd_premium_vector = as.data.frame(do.call('rbind',lapply(bd_prem, function(x) x[[2]]))); names(bd_premium_vector) = c('p1','p2')
bd_output = cbind(bd_output, bd_premium_vector); names(bd_output)[c(4,5)] = c('p1','p2')
bd_output = bd_output %>% mutate(p_ratio = p2/p1); 
bd_output$p_ratio[is.nan(bd_output$p_ratio)] = 1; 
pb_output = cbind(pb_output, pb_premium_vector); 

saveRDS(bd_output, file = paste0('../../Obj_for_manuscript/bd_output', file_name_label,'.rds'))
saveRDS(pb_output, file = paste0('../../Obj_for_manuscript/pb_output', file_name_label,'.rds'))

budget_2012_new = (bd_output %>% filter(p1 == 0.045 & p2 == 1.9 * 0.045) %>% pull(budget) %>% sum())

plot_list = list()
plot_list[[1]] = ggplot(data = bd_output , aes(x = p1, y = p_ratio, fill = surplus, color = budget >= (budget_2012_new))) + geom_tile(linewidth=0.3) + geom_text(aes(label = round(surplus,3)), color = 'white',size=2) + theme_bw() + theme(legend.position = "none") + xlab(expression(p[1])) + ylab(expression(p[2]/p[1])) + scale_y_continuous(breaks=seq(1,2,length.out=6))

plot_list[[2]] = ggplot(data = bd_output , aes(x = p1, y = p_ratio, fill = demand, color = budget >= (budget_2012_new))) + geom_tile(linewidth=0.3) + geom_text(aes(label = round(demand,3)), color = 'white',size=2) + theme_bw() + theme(legend.position = "none") + xlab(expression(p[1])) + ylab(expression(p[2]/p[1])) + scale_y_continuous(breaks=seq(1,2,length.out=6))

graph_data = rbind(pb_output, bd_output %>% select(-p_ratio) %>% filter(p2 == 2 * p1), bd_output %>% select(-p_ratio) %>% filter(p2 == 1.9 * p1))

graph_data$type = c(rep('PB', nrow(pb_output)), rep('IP', nrow(bd_output %>% filter(p2 == 2 * p1))), rep('BD', nrow(bd_output %>% filter(p2 == 1.9 * p1))))

plot_list[[3]] = ggplot(data = graph_data, aes(x = p2, y = surplus, linetype= budget >= budget_2012_new, color=type)) + geom_line()

plot_list[[4]] = ggplot(data = graph_data, aes(x = p2, y = demand, linetype= budget >= budget_2012_new, color=type)) + geom_line()


plot_list[[5]] = ggplot(data = data.frame(x = c(bd_output$surplus, pb_output$surplus), y = c(bd_output$budget, pb_output$budget), color = c(rep('BD', nrow(bd_output)), rep('PB', nrow(pb_output)))), aes(x = y, y = x, color = color)) + geom_point()


