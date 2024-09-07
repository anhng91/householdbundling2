deriv_2 = do.call('rbind', lapply(1:ncol(X_ind_sample_r_theta), function(index) 
	2 * output_2_sqrt[index] / ((do.call('c', lapply(deriv_list, function(x) (x[[3]][,3])))[index_type_2] * X_ind_sample_r_theta[index_type_2,index]) %>% mean) * (apply(do.call('rbind', lapply(deriv_list, function(x) x[[2]][[2]]))[index_type_2,], 2, function(xx) xx * X_ind_sample_r_theta[index_type_2,index]) %>% colMeans())
	))

f0 = optim_f(initial_param_trial[index_theta_only], include_r = FALSE, recompute_pref = FALSE)[[1]]
deriv_manual = rep(0, length(initial_param_trial))
for (i in index_theta_only) {
	x_new = initial_param_trial;
	x_new[i] = x_new[i] + 1e-3
	f1 = optim_f(x_new[index_theta_only], include_r = FALSE, recompute_pref = FALSE)[[1]]
	deriv_manual[i] = (f1 - f0)/1e-3
}

f0 = aggregate_moment_pref(x_transform, silent=TRUE, recompute_pref=FALSE)

x_transform_i = x_transform; x_transform_i[[1]]$beta_theta[1] = x_transform_i[[1]]$beta_theta[1] + 1e-3

f1 = aggregate_moment_pref(x_transform_i, silent=TRUE, recompute_pref=FALSE)




f0 = fid(Com_HH_list_index[1], x_transform)

x_transform_i = x_transform; x_transform_i[[1]]$beta_theta[1] = x_transform_i[[1]]$beta_theta[1] + 1e-3

f1 = fid(Com_HH_list_index[1], x_transform_i)

(f1[[2]] - f0[[2]])/1e-3; 