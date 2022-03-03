functions{
	// flatten_lower_tri: function that returns the lower-tri of a matrix, flattened to a vector
	vector flatten_lower_tri(matrix mat) {
		int n_cols = cols(mat);
		int n_uniq = (n_cols * (n_cols - 1)) %/% 2;
		vector[n_uniq] out ;
		int i = 1;
		for(c in 1:(n_cols-1)){
			for(r in (c+1):n_cols){
				out[i] = mat[r,c];
				i += 1;
			}
		}
		return(out) ;
	}

}

data{

	// nI: number of individuals
	int<lower=3> nI ;

	// nXc: number of condition-level predictors
	int<lower=2> nXc ;

	// rXc: number of rows in the condition-level predictor matrix
	int<lower=nXc> gauss_rXc ;
	int<lower=nXc> binom_rXc ;

	// Xc: condition-level predictor matrix
	matrix[gauss_rXc,nXc] gauss_Xc ;
	matrix[binom_rXc,nXc] binom_Xc ;

	// iXc: which individual is associated with each row in Xc
	array[gauss_rXc] int<lower=1,upper=nI> gauss_iXc ;
	array[binom_rXc] int<lower=1,upper=nI> binom_iXc ;

	// nY: num entries in the observation vector
	int<lower=nI> gauss_nY ;
	int<lower=nI> binom_nY ;

	// Y: observations
	vector[gauss_nY] gauss_Y  ;
	array[binom_nY] int<lower=0> binom_zeroes_count  ;
	array[binom_nY] int<lower=0> binom_ones_count  ;

	// yXc: which row in Xc is associated with each observation in Y
	array[gauss_nY] int<lower=1,upper=gauss_rXc> gauss_yXc ;
	array[binom_nY] int<lower=1,upper=binom_rXc> binom_yXc ;

}

transformed data{

	// computing binom_trials_count as a sum makes the automatic debugging data generator work better
	array[binom_nY] int binom_trials_count ;
	array[binom_nY] real binom_p ; // for computing the proportion ones
	for(i in 1:binom_nY){
		binom_trials_count[i] = binom_ones_count[i] + binom_zeroes_count[i] ;
		binom_p[i] = (binom_ones_count[i]*1.0) / (binom_trials_count[i]*1.0) ;
	}

	// logit_intercept_est: empirical mean proportion, log-odds scale
	real logit_intercept_est = logit( mean(binom_p) ) ;

	// tXc: transposed copy of Xc
	matrix[nXc,gauss_rXc] t_gauss_Xc = transpose(gauss_Xc) ;
	matrix[nXc,binom_rXc] t_binom_Xc = transpose(binom_Xc) ;


}

parameters{

	// coef_mean: group-level mean for each coefficient
	matrix[nXc,3] coef_mean_ ;

	// coef_variability_: group-level variability (on log-variance scale) for each coefficient
	matrix[nXc,3] coef_variability_ ; //helper-variable (note underscore suffix) for non-centered parameterization

	// coef_variability_mean_*: parameters for partial-pooling structure on coef_variability_
	vector[3] coef_variability_mean ;
	vector<lower=0>[3] coef_variability_sd ;

	// cors_chol: cholesky-factor of correlation structure associated with variability among individuals on influence of within-individual predictors
	cholesky_factor_corr[nXc] cors_chol ;

	// mvn_ncp_helper: helper-variable for non-centered parameterization multivariate normal
	matrix[nXc,nI] mvn_ncp_helper ;

	// variables encoding the SEM from gauss_locat to both gauss_scale & binom_logit
	matrix<lower=-1,upper=1>[nXc,2] locat_star_cor ;
	array[2] matrix[nXc,nI] star_unique ;

}

model{

	////////
	// Priors
	////////

	// flat prior on correlations
	target += lkj_corr_cholesky_lpdf(cors_chol|1);

	// standard-normal priors on all group-level coefficients
	target += std_normal_lupdf( to_vector(coef_mean_) ) ;

	// priors for partial pooling of coef_variability
	target += std_normal_lupdf( to_vector(coef_variability_mean) ) ;
	target += weibull_lupdf( to_vector(coef_variability_sd) | 2, 1) ; // n.b. non-zero-peaked

	// standard-normal prior on coef_variability_ necessary for non-centered partial-pooling of coef_variability
	target += std_normal_lupdf( to_vector(coef_variability_) ) ;

	// standard-normal prior on mvn_ncp_helper necessary for non-centered partial-pooling of coef_variability
	target += std_normal_lupdf( to_vector(mvn_ncp_helper) ) ;

	// standard-normal prior on star_unique necessary for SEM
	for(i in 1:2){
		target += std_normal_lupdf( to_vector(star_unique[i]) ) ;
	}

	////////
	// Compute derived quantities
	////////

	// coef_sd as partially-pooled
	matrix[nXc,3] coef_sd ;
	for(i in 1:3){
		// coef_sd[,i] = sqrt(exp( coef_variability_[,i] )); // unpooled
		coef_sd[,i] = sqrt(exp(coef_variability_[,i]*coef_variability_sd[i] + coef_variability_mean[i] ));
	}

	// centering the the binomial intercept on the empirical mean
	matrix[nXc,3] coef_mean = coef_mean_ ;
	coef_mean[1,3] += logit_intercept_est ;

	// locat_devs: zero-mean/unit-scale/correlated variates for locat coefficients
	matrix[nXc,nI] locat_devs = cors_chol * mvn_ncp_helper ;

	// Square of SEM weights
	matrix[nXc,2] locat_star_var_common = locat_star_cor.^2 ;
	matrix[nXc,2] locat_star_var_unique = (1-locat_star_var_common) ;

	// apply SEM to get *_devs
	array[2] matrix[nXc,nI] star_devs ;
	for(i in 1:2){
		star_devs[i] = (
			(
				(
					locat_devs
					.* rep_matrix(
						locat_star_cor[,i]
						, nI
					)
				)
				+ (
					star_unique[i]
					.* rep_matrix(
						sqrt(locat_star_var_unique[,i])
						, nI
					)
				)
			)
			// divide by the square-root of the sum of the squared weights to yield unit-scale variates (since the component variates have unit-scale too)
			./ rep_matrix(
				sqrt( locat_star_var_common[,i] + locat_star_var_unique[,i] )
				, nI
			)
		) ;
	}

	// compute condition values implied by coefficients & predictors
	row_vector[gauss_rXc] locat_dot_Xc = columns_dot_product(
		(
			// shift & scale devs
			locat_devs
			.* rep_matrix(coef_sd[,1],nI)
			+ rep_matrix(coef_mean[,1],nI)
		)[,gauss_iXc]
		, t_gauss_Xc
	) ;

	row_vector[gauss_rXc] scale_dot_Xc = columns_dot_product(
		(
			// shift & scale devs
			star_devs[1]
			.* rep_matrix(coef_sd[,2],nI)
			+ rep_matrix(coef_mean[,2],nI)
		)[,gauss_iXc]
		, t_gauss_Xc
	) ;

	row_vector[binom_rXc] logit_dot_Xc = columns_dot_product(
		(
			// shift & scale devs
			star_devs[2]
			.* rep_matrix(coef_sd[,3],nI)
			+ rep_matrix(coef_mean[,3],nI)
		)[,binom_iXc]
		, t_binom_Xc
	) ;

	////////
	// Observation-level structure ("Likelihood")
	////////

	// gauss observations as location-sccale
	target += normal_lupdf(
		gauss_Y
		|
		locat_dot_Xc[gauss_yXc]
		, sqrt(exp( scale_dot_Xc ))[gauss_yXc]
	) ;

	// binom observations
	target += binomial_logit_lpmf(
		binom_ones_count
		|
		binom_trials_count
		, logit_dot_Xc[binom_yXc]
	) ;

}

generated quantities{

	// cors: lower-tri of correlation matrix flattened to a vector
	vector[(nXc*(nXc-1))%/%2] cors = flatten_lower_tri(multiply_lower_tri_self_transpose(cors_chol)) ;

}
