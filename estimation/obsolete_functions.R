#' Compute the expected utility, expected squared utility, and expected medical spending with respect to the preference parameters.
#'
#' @param data_set This is a list that contains all pre-computed elements of each household .
#' @param param This is a list, each element is a coefficient (i.e., a vector) of a preference parameter (beta_gamma, beta_omega, beta_delta, beta_theta, beta_theta_ind) and the standard deviations of the unobserved heterogeneity (sigma_gamma, sigma_omega, sigma_delta, sigma_theta, sigma_thetabar)
#' @param n_draw this is number of Gauss Laguerre points.
#' @param default_sigma is between 0 and 1. Default to 0.5. represents the value of sigma_theta that is used for Taylor second-order approximation. sigma_theta is computed as default_sigma * exp(param$sigma_thetabar)
#' @param taylor_order an integer, equal to the order of Taylor approximation used to approximate exp(-x^2/sigma^2) around (sigma^2 = sigma_0^2).
#' @param u_lowerbar value of utility when the residual income is nonpositive. Default value is -10
#'
#' @return a list that includes the draws of U, U squared, and m. These draws are list with length equal to the order of the Taylor approximation. Each component of U and U2 is a matrix of dimension the number of bundles x length of thetabar. Each component of m is a matrix with dimension equal to length of thetabar x hhsize.
#' 
#' @export
#'
#' @examples
#' mini_data = household_draw_theta_kappa_Rdraw(1, sample_data_and_parameter$param, 1000, 10, sick_parameters_sample, xi_parameters_sample);
#' compute_expected_U_m(mini_data, sample_data_and_parameter$param, n_draw, default_sigma = 0.5, taylor_order=2)



compute_expected_U_m = function(data_set, param, n_draw, default_sigma = 0.5, taylor_order, u_lowerbar = -10) {
	stop('This function is obsolete')
	stopifnot(default_sigma <= 1 & default_sigma >= 0)
	data_hh_i = data_set$data;
	policy_mat_hh_index = data_set$policy_mat;
	hh_index = data_set$index; 
	theta_draw = data_set$theta_draw;
	income_vec = data_set$income;
	X_ind = data_set$X_ind; 
	X_hh = data_set$X_hh;
	sick_p = data_set$sick_p;
	R_draw = data_set$R_draw;
	kappa_draw = data_set$kappa_draw; 
	
	sigma_theta_taylor = default_sigma * exp(param$sigma_thetabar);


	HHsize = data_hh_i$HHsize[1]

	mean_beta_omega = X_hh %*% param$beta_omega; 
	mean_beta_delta = X_ind %*% param$beta_delta;
	mean_beta_gamma = X_ind %*% param$beta_gamma; 
 
	n_draw = 50;

	U = list(); 
	U2 = list(); 
	m_draw = theta_draw; 

	for (insurance_bundle in 1:length(R_draw)) {
		U[[insurance_bundle]] = lapply(R_draw[[insurance_bundle]], function(x) ifelse(x > 0, CRRA_fun_Gauss(log(x), mu = mean_beta_omega, sigma = exp(param$sigma_omega), n_draw = n_draw), 0)) %>% unlist()

		U2[[insurance_bundle]] = lapply(R_draw[[insurance_bundle]], function(x) ifelse(x > 0, CRRA_fun_Gauss(2 * log(x), mu = mean_beta_omega, sigma = exp(param$sigma_omega), n_draw = n_draw), 0)) %>% unlist()  
	}
	
	for (mem_index in 1:HHsize) {
		for (insurance_bundle in 1:length(R_draw)) {
			if (insurance_bundle == 1) {
				obj1 =  mean_E_XW(mean_beta_gamma[mem_index], exp(param$sigma_gamma), 1/(1 + kappa_draw[[1]][,mem_index]))
				obj2 = truncated_normal_mean(mean_beta_delta[mem_index], exp(param$sigma_delta))
				obj3 = mean_E_XW(mean_beta_omega, exp(param$sigma_omega), R_draw[[insurance_bundle]])
				obj4 =  mean_E_XW(mean_beta_omega, exp(param$sigma_omega), R_draw[[insurance_bundle]]);

				obj8 = truncated_normal_variance(mean_beta_delta[mem_index], exp(param$sigma_delta))
				obj9 = mean_E_XW(mean_beta_omega, exp(param$sigma_omega), R_draw[[insurance_bundle]]^2)
				obj11 =  mean_E_XW(mean_beta_gamma[mem_index], exp(param$sigma_gamma), 1/(1 + kappa_draw[[1]][,mem_index])^2)


				m_draw[,mem_index] = theta_draw[,mem_index] * kappa_draw[[1]][,mem_index] + obj2 * theta_draw[, mem_index] * obj4 * obj1 * kappa_draw[[1]][,mem_index]	

			}


			U[[insurance_bundle]] = U[[insurance_bundle]] - unlist(lapply(kappa_draw[[insurance_bundle]][,mem_index],function(x)  CRRA_fun_Gauss(log(1 + x), mu = mean_beta_gamma[mem_index], sigma = exp(param$sigma_gamma), n_draw = n_draw))) * (R_draw[[insurance_bundle]] > 0) * theta_draw[,mem_index] * truncated_normal_mean(mu = mean_beta_delta[mem_index], sigma = exp(param$sigma_delta))


			for (mem_index_j in 1:HHsize) {
				if (mem_index_j == mem_index) {
					U2[[insurance_bundle]] = U2[[insurance_bundle]] + unlist(lapply(kappa_draw[[insurance_bundle]][,mem_index],function(x)  CRRA_fun_Gauss(2*log(1 + x), mu = mean_beta_gamma[mem_index], sigma = exp(param$sigma_gamma), n_draw = n_draw))) * (R_draw[[insurance_bundle]] > 0) * theta_draw[,mem_index]^2 * truncated_normal_variance(mu = mean_beta_delta[mem_index], sigma = exp(param$sigma_delta))

					
				} else {
					U2[[insurance_bundle]] = U2[[insurance_bundle]] + unlist(lapply(kappa_draw[[insurance_bundle]][,mem_index],function(x)  CRRA_fun_Gauss(log(1 + x), mu = mean_beta_gamma[mem_index], sigma = exp(param$sigma_gamma), n_draw = n_draw))) * (R_draw[[insurance_bundle]] > 0) * theta_draw[,mem_index] * truncated_normal_mean(mu = mean_beta_delta[mem_index], sigma = exp(param$sigma_delta)) * unlist(lapply(kappa_draw[[insurance_bundle]][,mem_index_j],function(x)  CRRA_fun_Gauss(log(1 + x), mu = mean_beta_gamma[mem_index_j], sigma = exp(param$sigma_gamma), n_draw = n_draw))) * (R_draw[[insurance_bundle]] > 0) * theta_draw[,mem_index_j] * truncated_normal_mean(mu = mean_beta_delta[mem_index_j], sigma = exp(param$sigma_delta))
				}
			}

			U2[[insurance_bundle]] = U2[[insurance_bundle]] - 2 * unlist(lapply(kappa_draw[[insurance_bundle]][,mem_index],function(x)  CRRA_fun_Gauss(log(1 + x), mu = mean_beta_gamma[mem_index], sigma = exp(param$sigma_gamma), n_draw = n_draw))) * (R_draw[[insurance_bundle]] > 0) * theta_draw[,mem_index] * truncated_normal_mean(mu = mean_beta_delta[mem_index], sigma = exp(param$sigma_delta)) * (lapply(R_draw[[insurance_bundle]], function(x) ifelse(x > 0, CRRA_fun_Gauss(log(x), mu = mean_beta_omega, sigma = exp(param$sigma_omega), n_draw = n_draw), 0)) %>% unlist())

		}
	}

	for (insurance_bundle in 1:length(R_draw)) {
		U[[insurance_bundle]] = U[[insurance_bundle]] * (R_draw[[insurance_bundle]] > 0) + u_lowerbar * (R_draw[[insurance_bundle]] <= 0)
		U2[[insurance_bundle]] = U2[[insurance_bundle]] * (R_draw[[insurance_bundle]] > 0) + u_lowerbar^2 * (R_draw[[insurance_bundle]] <= 0)
	}

	U = do.call('cbind', U); 
	U2 = do.call('cbind', U2);
	data_set$expected_U = list();
	data_set$expected_U2 = list();

	for (n_taylor in c(0:taylor_order)) {
		data_set$expected_U[[n_taylor + 1]] = matrix(NA, nrow = nrow(data_set$Hermite_draw_mat$x), ncol = ncol(U))
		data_set$expected_U2[[n_taylor + 1]] = matrix(NA, nrow = nrow(data_set$Hermite_draw_mat$x), ncol = ncol(U))
		data_set$m[[n_taylor + 1]] = matrix(NA, nrow = nrow(data_set$Hermite_draw_mat$x), ncol = HHsize)
	}

	for (theta_bar_index in c(1:nrow(data_set$Hermite_draw_mat$x))) {
		theta_bar = data_set$Hermite_draw_mat$x[theta_bar_index, ];
		theta_draw_minus_mean = apply(theta_draw, 1, function(x) x - theta_bar) %>% matrix(ncol = HHsize)
		weight = exp(-rowSums(theta_draw_minus_mean)^2/2/sigma_theta_taylor^2)
		for (n_taylor in c(0:taylor_order)) {
			data_set$expected_U[[n_taylor + 1]][theta_bar_index,] = t(U) %*% ((-1/2 * rowSums(theta_draw_minus_mean^2))^n_taylor/factorial(n_taylor) * weight)
			data_set$expected_U2[[n_taylor + 1]][theta_bar_index,] = t(U2) %*% ((-1/2 * rowSums(theta_draw_minus_mean^2))^n_taylor/factorial(n_taylor) * weight)
			data_set$m[[n_taylor + 1]][theta_bar_index,] = t(m_draw) %*% ((-1/2 * rowSums(theta_draw_minus_mean^2))^n_taylor/factorial(n_taylor) * weight)
		}
	} 

	return(data_set);
}	


#' Compute the draws (preference parameters independent).
#'
#' @param hh_index index (from object data_hh_list) of a particular household
#' @param param This is a list, each element is a coefficient (i.e., a vector) of a preference parameter (beta_gamma, beta_omega, beta_delta, beta_theta, beta_theta_ind) and the standard deviations of the unobserved heterogeneity (sigma_gamma, sigma_omega, sigma_delta, sigma_theta, sigma_thetabar)
#' @param n_draw_gauss this is number of Gauss  points.
#' @param n_draw_halton this is the number of halton draws, default to 1000
#'
#' @return a list that includes the draws of household-related objects,taking into account the sick parameters and the distribution of coverage. This does not take into account estimated preference parameters or unconditional distribution of health shocks. See `compute_expected_U_m` for the draws post-estimation of preference parameters and health shocks distribution.
#' 
#' @export
#'
#' @examples
#' household_draw_theta_kappa_Rdraw_obsolete(1, sample_data_and_parameter$param, 1000, 10, sick_parameters_sample, xi_parameters_sample)
#' 
household_draw_theta_kappa_Rdraw_obsolete = function(hh_index, param, n_draw_halton = 1000, n_draw_gauss = 10, sick_parameters, xi_parameters) {
	stop('this function is deprecated. Please use household_draw_theta_kappa_Rdraw')
	set.seed(1);
	data_hh_i = data_hh_list[[hh_index]]; 
	HHsize = nrow(data_hh_i);
	X_ind = var_ind(data_hh_i)
	X_hh = var_hh(data_hh_i)
	draw_p_xi_0 = X_ind %*% xi_parameters$par[1:ncol(X_ind)]
	draw_p_xi_1 = X_ind %*% xi_parameters$par[(ncol(X_ind) + 1):(2*ncol(X_ind))]
	p_0 = exp(draw_p_xi_0)/(1 + exp(draw_p_xi_0) + exp(draw_p_xi_1))
	p_1 = exp(draw_p_xi_1)/(1 + exp(draw_p_xi_1) + exp(draw_p_xi_0))
	sick_p = c(1/(1 + exp(X_ind %*% sick_parameters$par)))
	halton_mat = randtoolbox::halton(n_draw_halton, HHsize * 3 + 1)
	household_random_factor = qnorm(halton_mat[,HHsize * 3 + 1]); 

	
	kappa_draw = list(); 

	theta_draw = matrix(NA, nrow=n_draw_halton, ncol = HHsize) # observed insurance 
	kappa_draw[[1]] = matrix(NA, nrow=n_draw_halton, ncol = HHsize) # observed coinsurance 

	policy_mat_hh_index = list(); 
	policy_mat_hh_index[[1]] = policy_mat[[hh_index]];

	elig_member_index = NULL; 

	if (data_hh_i$HHsize_s[1] != 0) {
		elig_member_index = which(data_hh_i$Bef_sts + data_hh_i$Com_sts + data_hh_i$Std_w_ins == 0)
		for (i in 1:length(elig_member_index)) {
			if (data_hh_i$Vol_sts[elig_member_index[i]] == 0) {
				policy_mat_hh_index[[1]][[1]][,elig_member_index[i]]  = c(1, rep(NA, 9));
				policy_mat_hh_index[[1]][[2]][,elig_member_index[i]]  = c(0, rep(NA, 9)); 
				policy_mat_hh_index[[1]][[3]][,elig_member_index[i]]  = c(0, rep(NA, 9));  
			}
		}
		for (i in 1:length(elig_member_index)) {
			policy_mat_hh_index[[i + 1]] = policy_mat_hh_index[[1]]; 

			if (data_hh_i$Vol_sts[elig_member_index[i]] == 1) {
				policy_mat_hh_index[[i + 1]][[1]][,elig_member_index[i]] = c(1, rep(NA, 9));
				policy_mat_hh_index[[i + 1]][[2]][,elig_member_index[i]] = c(0, rep(NA, 9)); 
				policy_mat_hh_index[[i + 1]][[3]][,elig_member_index[i]] = c(0, rep(NA, 9));  
			} else {
				policy_mat_hh_index[[i + 1]][[1]][,elig_member_index[i]]  = policy_mat[[hh_index]][[1]][,elig_member_index[i]]
				policy_mat_hh_index[[i + 1]][[2]][,elig_member_index[i]]  = policy_mat[[hh_index]][[2]][,elig_member_index[i]]
				policy_mat_hh_index[[i + 1]][[3]][,elig_member_index[i]]  = policy_mat[[hh_index]][[3]][,elig_member_index[i]]
			}

			kappa_draw[[i + 1]] = matrix(NA, nrow=n_draw_halton, ncol = HHsize)
		} 
	}

	

	common_household_factor = matrix(NA, nrow=n_draw_halton, ncol = HHsize);
	individual_factor = matrix(NA, nrow=n_draw_halton, ncol = HHsize); 


	Hermite_draw = pracma::gaussHermite(n_draw_gauss)

	# normalize Hermite draws to take into account the conditional distribution of theta

	Hermite_draw_mat = list()

	if (data_hh_i$HHsize_s[1] != 0) {
		combn_index = all_combn(1:n_draw_gauss, HHsize);
		Hermite_draw_mat$x = do.call('cbind', lapply(1:HHsize, function(member_index) Hermite_draw$x[combn_index[,member_index]] * (exp(param$sigma_thetabar) + (X_ind[member_index,] %*% param$beta_theta_ind)^2) + t(c(X_ind[member_index,], data_hh_i$Year[member_index] == 2004, data_hh_i$Year[member_index] == 2006 , data_hh_i$Year[member_index] == 2010, data_hh_i$Year[member_index] == 2012)) %*% param$beta_theta))
		Hermite_draw_mat$w = apply(do.call('cbind', lapply(1:HHsize, function(member_index) Hermite_draw$w[combn_index[,member_index]])), 1, prod)
	}
	
	for (i in 1:HHsize) {
		common_household_factor[,i] = household_random_factor * (X_ind[i,] %*% param$beta_theta_ind) + t(c(X_ind[i,], data_hh_i$Year[i] == 2006, data_hh_i$Year[i] == 2008, data_hh_i$Year[i] == 2010, data_hh_i$Year[i] == 2012)) %*% param$beta_theta
		individual_factor[,i] = qnorm(halton_mat[,i]);
		theta_draw[,i] = exp(common_household_factor[,i] + individual_factor[,i] * exp(param$sigma_thetabar)) * (halton_mat[,2 * HHsize + i] >= sick_p[i])
		random_xi_draws_i = halton_mat[,HHsize + i]; 
		random_xi_draws = lapply(random_xi_draws_i, function(x) ifelse(x <= p_0[i], 0, ifelse(x <= p_0[i] + p_1[i], 1, x))) %>% unlist()

		for (insurance_bundle in 1:length(policy_mat_hh_index)) {
			kappa_draw[[insurance_bundle]][,i] = (lapply(1:nrow(theta_draw), function(j) policy_mat_hh_index[[insurance_bundle]][[1]][max(which((theta_draw[j,i] * random_xi_draws[j]) >= policy_mat_hh_index[[insurance_bundle]][[2]][,i])),i]) %>% unlist()) * random_xi_draws + 1 - random_xi_draws
		}
	}
	
	R_draw = list()

	for (insurance_bundle in 1:length(policy_mat_hh_index)) {
		R_draw[[insurance_bundle]] =  lapply(Income_net_premium[[hh_index]][ifelse(insurance_bundle == 1, sum(data_hh_i$Vol_sts) + 1, ifelse(data_hh_i$Vol_sts[elig_member_index[insurance_bundle - 1]] == 1, sum(data_hh_i$Vol_sts), sum(data_hh_i$Vol_sts) + 2))] - rowSums(theta_draw * kappa_draw[[insurance_bundle]]), function(x) max(x,0)) %>% unlist()
	}	
	

	output = list(); 
	output$R_draw = R_draw; 
	output$kappa_draw = kappa_draw; 
	output$theta_draw = theta_draw;
	output$common_household_factor = common_household_factor; 
	output$individual_factor = individual_factor;
	output$index = hh_index
	output$data = data_hh_i
	output$policy_mat = policy_mat[[hh_index]]
	output$income = Income_net_premium[[hh_index]]
	output$X_ind = var_ind(data_hh_i)
	output$X_hh = var_hh(data_hh_i)
	output$sick_p = c(1/(1 + exp(output$X_ind %*% sick_parameters$par)))
	output$Hermite_draw_mat = Hermite_draw_mat 
	return(output)	
}

#' Compute the moments of eligible households
#'
#' @param data_set This is a list that contains all pre-computed elements of each household .
#' @param param This is a list, each element is a coefficient (i.e., a vector) of a preference parameter (beta_gamma, beta_omega, beta_delta, beta_theta, beta_theta_ind) and the standard deviations of the unobserved heterogeneity (sigma_gamma, sigma_omega, sigma_delta, sigma_theta, sigma_thetabar)
#' @param deriv_tol Default set to 1e-3. The tolerance level for the numerical computation of derivative. 
#' @param default_sigma Is a number between 0 and 1 to use in the Taylor approximation
#' @param taylor_order is the order of Taylor approximation 
#' @param xi_parameters estimated parameters for coverage
#' @param sick_parameters estimated parameters for probability of being sick 
#'
#' @return a list. The first element is expected value of oop conditional on the observed insurance status being optimal, the second element is the expected value of oop squared conditional on the observed insurance status being optimal, the third element is the probability of observed insurance being optimal. The next 3 elements are the derivative w.r.t the coefficients being optimized of the first 3 elements. 
#' 
#' @export
#'
#' @examples
#' mini_data = household_draw_theta_kappa_Rdraw(1, sample_data_and_parameter$param, 1000, 10, sick_parameters_sample, xi_parameters_sample);
#' mini_data = compute_expected_U_m(mini_data, sample_data_and_parameter$param, 
#' n_draw, default_sigma = 0.5, taylor_order=4)
#' moment_eligible_hh(mini_data, sample_data_and_parameter$param, deriv_tol=1e-3, default_sigma = 0.5, taylor_order=4)
#' 

moment_eligible_hh = function(data_set, param, deriv_tol = 1e-3, default_sigma, taylor_order=3,  sick_parameters, xi_parameters) {
	if (!('expected_U' %in% names(data_set))) {
		stop('data_set must be an object produced by compute_expected_U_m');
	}
	tol = 1e-3;
	data_hh_i = data_set$data;
	policy_mat_hh_index = data_set$policy_mat;
	hh_index = data_set$index; 
	theta_draw = data_set$theta_draw;
	income_vec = data_set$income;
	X_ind = data_set$X_ind; 
	X_hh = data_set$X_hh;
	sick_p = data_set$sick_p;
	R_draw = data_set$R_draw;
	kappa_draw = data_set$kappa_draw; 
	common_household_factor = data_set$common_household_factor; 
	individual_factor = data_set$individual_factor;
	expected_U = data_set$expected_U
	expected_U2 = data_set$expected_U2; 
	draws_m = data_set$m 

	HHsize = data_hh_i$HHsize[1]
	sigma_theta_taylor = default_sigma * exp(param$sigma_thetabar);

	
	inner_f = function(mean_beta_r, sigma_theta, sigma_r) {
		U_draw = 0; 
		U2_draw = 0; 
		m = 0;

		for (n_taylor in c(1:length(expected_U))) {
			U_draw = U_draw + expected_U[[n_taylor]] * (sigma_theta^2 - sigma_theta_taylor^2)^(n_taylor - 1)/sigma_theta/sqrt(2 * pi)
			U2_draw =  U2_draw + expected_U2[[n_taylor]] * (sigma_theta^2 - sigma_theta_taylor^2)^(n_taylor - 1)/sigma_theta/sqrt(2 * pi)
			m = m + draws_m[[n_taylor]] * (sigma_theta^2 - sigma_theta_taylor^2)^(n_taylor - 1)/sigma_theta/sqrt(2 * pi)
		}


		CARA_U = do.call('cbind', lapply(1:ncol(U_draw), function(insurance_bundle) second_taylor_CARA(U_draw[,insurance_bundle], U2_draw[,insurance_bundle], mean_beta_r, sigma_r)));

		prob_insured = apply(CARA_U, 1, function(x) x[1] >= max(x[-1]))

		m_conditional_insurance = colMeans(matrix(apply(m, 2, function(x) prob_insured * x * data_set$Hermite_draw_mat$w), ncol=ncol(m)), na.rm=TRUE)/mean(prob_insured, na.rm=TRUE)

		return(c(mean(prob_insured, na.rm=TRUE) , m_conditional_insurance))
	}

	mean_beta_r = matrix(X_hh, nrow=1) %*% param$beta_r 

	f0 = inner_f(mean_beta_r, exp(param$sigma_theta)/(1 + exp(param$sigma_theta)) * exp(param$sigma_thetabar), exp(param$sigma_r));

	f1 = list(inner_f(mean_beta_r + deriv_tol, exp(param$sigma_theta)/(1 + exp(param$sigma_theta)) * exp(param$sigma_thetabar), exp(param$sigma_r)), inner_f(mean_beta_r, exp(param$sigma_theta + deriv_tol)/(1 + exp(param$sigma_theta + deriv_tol)) * exp(param$sigma_thetabar), exp(param$sigma_r)), inner_f(mean_beta_r, exp(param$sigma_theta)/(1 + exp(param$sigma_theta)) * exp(param$sigma_thetabar), exp(param$sigma_r + deriv_tol))); 

	df = lapply(f1, function(x) (x - f0)/1e-3); 
	names(df) = c('beta_r', 'sigma_theta', 'sigma_r')
	df[['beta_r']] = df[['beta_r']] %*% matrix(X_hh, nrow=1); 

	return(list(f0, df))
} 
