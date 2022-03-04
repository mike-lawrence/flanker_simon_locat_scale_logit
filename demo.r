# preamble (options, installs, imports & custom functions) ----
library(tidyverse)
library(cmdstanr)

for(file in list.files(	path = 'r', full.names = T)){
	source(file)
}

# read in the data ----

dat = readRDS('rds/dat.rds')
print(dat)


# Gaussian inputs ----

(
	dat
	%>% filter(!error)
) ->
	gauss_dat

(
	gauss_dat
	# select individual & any condition-defining columns
	%>% select(individual,flankers)#,simon)
	# collapse down to distinct rows (1 per individual/conditions combo)
	%>% distinct()
	# expand (in case there's any missing data; doesn't hurt if not)
	# %>% exec(expand_grid,!!!.) #should work, but doesn't, so need the next three instead
	%>% as.list()
	%>% map(unique)
	%>% cross_df()
	# arrange (not really necessary, but why not)
	%>% arrange()
	# add the contrast matrix columns
	%>% mutate(
		contrasts = get_contrast_matrix_rows_as_list(
			data = .
			, formula = ~ flankers*simon
			, contrast_kind = halfhelmert_contrasts
		)
	)
) ->
	complete_gauss_Xc_with_vars

# show the unique contrasts
(
	complete_gauss_Xc_with_vars
	%>% select(-individual)
	%>% distinct()
	%>% unnest(contrasts)
)

# subset down to just those individual-condition combos actually present in the data
#  it's ok if there's no missing data and nrow(complete_Xc_with_vars)==nrow(Xc_with_vars)
(
	complete_gauss_Xc_with_vars
	%>% semi_join(gauss_dat)
) ->
	gauss_Xc_with_vars


# join Xc with dat to label observations with corresponding row from Xc
(
	gauss_Xc_with_vars
	# add row identifier
	%>% mutate(Xc_row=1:n())
	# right-join with dat to preserve dat's row order
	%>% right_join(
		mutate(gauss_dat,dat_row=1:n())
	)
	%>% arrange(dat_row)
	# pull the Xc row identifier
	%>% pull(Xc_row)
) ->
	gauss_yXc


# Binomial inputs ----

(
	dat
	%>% group_by(individual,flankers)#,simon)
	%>% summarise(
		zeroes_count = sum(error==0)
		, ones_count = sum(error==1)
		, .groups = 'drop'
	)
) ->
	binom_dat

(
	binom_dat
	# select individual & any condition-defining columns
	%>% select(individual,flankers)#,simon)
	# collapse down to distinct rows (1 per individual/conditions combo)
	%>% distinct()
	# expand (in case there's any missing data; doesn't hurt if not)
	# %>% exec(expand_grid,!!!.) #should work, but doesn't, so need the next three instead
	%>% as.list()
	%>% map(unique)
	%>% cross_df()
	# arrange (not really necessary, but why not)
	%>% arrange()
	# add the contrast matrix columns
	%>% mutate(
		contrasts = get_contrast_matrix_rows_as_list(
			data = .
			, formula = ~ flankers*simon
			, contrast_kind = halfhelmert_contrasts
		)
	)
) ->
	complete_binom_Xc_with_vars

# show the unique contrasts
(
	complete_binom_Xc_with_vars
	%>% select(-individual)
	%>% distinct()
	%>% unnest(contrasts)
)

# subset down to just those individual-condition combos actually present in the data
#  it's ok if there's no missing data and nrow(complete_Xc_with_vars)==nrow(Xc_with_vars)
(
	complete_binom_Xc_with_vars
	%>% semi_join(binom_dat)
) ->
	binom_Xc_with_vars


# join Xc with dat to label observations with corresponding row from Xc
(
	binom_Xc_with_vars
	# add row identifier
	%>% mutate(Xc_row=1:n())
	# right-join with dat to preserve dat's row order
	%>% right_join(
		mutate(binom_dat,dat_row=1:n())
	)
	%>% arrange(dat_row)
	# pull the Xc row identifier
	%>% pull(Xc_row)
) ->
	binom_yXc


# package for stan & sample ----

data_for_stan = lst( # lst permits later entries to refer to earlier entries

	########
	# Entries we need to specify ourselves
	########

	# Y: observations
	gauss_Y = scale(gauss_dat$log_rt)[,1] # note scaling! adjust interpretation of priors accordingly
	, binom_zeroes_count = binom_dat$zeroes_count
	, binom_ones_count = binom_dat$ones_count

	# Xc: condition-level predictor matrix
	, gauss_Xc = (
		gauss_Xc_with_vars
		%>% select(contrasts)
		%>% unnest(contrasts)
		%>% as.matrix()
	)
	, binom_Xc = (
		binom_Xc_with_vars
		%>% select(contrasts)
		%>% unnest(contrasts)
		%>% as.matrix()
	)

	# yXc: which row in Xc is associated with each observation in Y
	, gauss_yXc = gauss_yXc
	, binom_yXc = binom_yXc

	# iXc: which individual is associated with each row in Xc
	, gauss_iXc = as.numeric(factor(gauss_Xc_with_vars$individual))
	, binom_iXc = as.numeric(factor(binom_Xc_with_vars$individual))

	########
	# Entries computable from the above
	########

	# nI: number of individuals
	, nI = max(gauss_iXc)

	# nXc: number of cols in the condition-level predictor matrix
	, nXc = ncol(gauss_Xc)

	# rXc: number of rows in the condition-level predictor matrix
	, gauss_rXc = nrow(gauss_Xc)
	, binom_rXc = nrow(binom_Xc)

	# nY: num entries in the observation vectors
	, gauss_nY = length(gauss_Y)
	, binom_nY = length(binom_ones_count)

)

glimpse(data_for_stan)

mod = cmdstan_model('stan/hierarchical_within_loc_scale_err_sem.stan')
post = mod$sample(
	data = data_for_stan
	, chains = parallel::detectCores()/2
	, parallel_chains = parallel::detectCores()/2
	, refresh = 1
	, init = 1 # default value of 2 tends to yield initialization failure
	# , max_treedepth = 11
	# , metric = 'dense_e' #slow, but might help with divergences?
)

#save outputs
post$save_object('r/post.rds')

#check time
post$time()

#check diagnostics
post$cmdstan_diagnose()

#check diagnostics by hand
(
	post$sampler_diagnostics()
	%>% posterior::as_draws_df()
	%>% group_by(.chain)
	%>% summarise(
		max_treedepth = max(treedepth__)
		, num_divergent = sum(divergent__)
		, rebfmi = var(energy__)/(sum(diff(energy__)^2)/n()) # n.b. reciprocal of typical EBFMI, so bigger=bad, like rhat
	)
)

#check that treedepth distribution doesn't seem bunched at max
(
	post$sampler_diagnostics()
	%>% posterior::as_draws_df()
	%>% select(.chain,.iteration,treedepth__)
	%>% mutate(treedepth__=factor(treedepth__,levels=c(1:max(treedepth__))))
	%>% group_by(.chain)
	%>% count(treedepth__,.drop=F)
	%>% pivot_wider(names_from=.chain,values_from=n)
)


# gather summary for core parameters (inc. r̂ & ess)
(
	post$draws()
	%>% posterior::as_draws_df()
	%>% select(-contains('chol')) #remove cholesky, as we have the lower tri too
	%>% posterior::summarise_draws(
		~ posterior::quantile2(.x, probs = c(0,.065,.25,.5,.75,.935,1))
		, posterior::default_convergence_measures()
		, .cores = parallel::detectCores()
	)
) ->
	par_summary

# show the ranges of r̂/ess's
(
	par_summary
	%>% select(rhat,contains('ess'))
	%T>% {function(x){
		print(summary(x))
	}}()
	# %>% summary()
	%>% mutate(
		log10_ess_bulk = log10(ess_bulk)
		, log10_ess_tail = log10(ess_bulk)
	)
	%>% select(-ess_bulk,-ess_tail)
	%>% pivot_longer(everything())
	%>% ggplot()
	+ facet_wrap(~name,scales='free')
	+ geom_histogram(aes(x=value))
)

# View those with suspect r̂
(
	par_summary
	%>% filter(rhat>1.01)
	%>% (function(suspects){
		if(nrow(suspects)>1){
			View(suspects)
		}
		return(paste('# suspect parameters:',nrow(suspects)))
	})()
)

# Viz group-level coefficient means ----
(
	post$draws('coef_mean_')
	%>% posterior::as_draws_df()
	%>% select(-.draw)
	# %>% select(
	# 	.chain
	# 	, .iteration
	# 	, ends_with(',1]')
	# )
	%>% pivot_longer(
		cols = -c(.chain,.iteration)
		, names_to = 'variable'
	)
	%>% group_nest(variable)
	%>% left_join(par_summary,by='variable')
	%>% separate(
		variable
		, into=c('variable','index')
		, sep='\\['
	)
	%>% separate(
		index
		, into=c('index1','index2')
		, sep=','
	)
	%>% mutate(
		index1 = as.numeric(index1)
		, index2 = as.numeric(str_remove(index2,fixed(']')))
		, variable = (
			complete_binom_Xc_with_vars
			%>% select(contrasts)
			%>% unnest(contrasts)
			%>% names()
		)[index1]
		, type = c('logRT Location','logRT Scale','logit Error')[index2]
	)
	%>% mutate(
		rhat_is_bad = 1.01<rhat
		, ess_bulk_is_bad = 100>ess_bulk
		, ess_tail_is_bad = 100>ess_tail
	)
) ->
	to_plot

(
	to_plot
	%>% ggplot()
	+ facet_wrap(~type,ncol=3)
	+ geom_hline(yintercept = 0,linetype=3,colour='white')
	+ geom_violin(
		data = (. %>% select(type,variable,data) %>% unnest(data))
		, mapping = aes(
			x = variable
			, y = value
		)
		, colour = 'transparent'
		, fill = 'black'
		, scale = 'width'
		, alpha = .5
	)
	+ geom_linerange(
		mapping = aes(
			x = variable
			, ymin = q0
			, ymax = q100
		)
		, colour = 'black'
		# , alpha = .5
		, size = .5
	)
	+ geom_linerange(
		mapping = aes(
			x = variable
			, ymin = q6.5
			, ymax = q93.5
			, colour = ess_tail_is_bad
		)
	)
	+ geom_linerange(
		mapping = aes(
			x = variable
			, ymin = q25
			, ymax = q75
			, colour = ess_bulk_is_bad
		)
		, size = 3
	)
	+ geom_point(
		mapping = aes(
			x = variable
			, y = q50
			, fill = rhat_is_bad
		)
		, shape = 21
		, size = 2
		, colour = 'black'
	)
	+ geom_hline(yintercept=0,linetype=3)
	+ coord_flip()
	+ scale_color_manual(
		values = lst(`TRUE`='red',`FALSE`='white')
		, labels = lst(`TRUE`='<100',`FALSE`='>=100')
	)
	+ scale_fill_manual(
		values = lst(`TRUE`='red',`FALSE`='black')
		, labels = lst(`TRUE`='>1.01',`FALSE`='<=1.01')
	)
	+ labs(
		y = "Posterior Value of Predictor's coefficient"
		, x = 'Predictor'
		, colour = 'ESS'
		, fill = 'Rhat'
	)
	+ theme(
		aspect.ratio = 2/(1+sqrt(5))
	)
)

# Viz logRT Location correlations ----
(
	post$draws('cors')
	%>% posterior::as_draws_df()
	%>% select(-.draw)
	%>% pivot_longer(
		cols = -c(.chain,.iteration)
		, names_to = 'variable'
	)
	%>% group_nest(variable)
	%>% left_join(par_summary,by='variable')
	%>% mutate(
		rhat_is_bad = 1.01<rhat
		, ess_bulk_is_bad = 100>ess_bulk
		, ess_tail_is_bad = 100>ess_tail
	)
) ->
	to_plot

(
	to_plot
	%>% ggplot()
	+ geom_hline(yintercept = 0,linetype=3,colour='white')
	+ geom_violin(
		data = (. %>% select(variable,data) %>% unnest(data))
		, mapping = aes(
			x = variable
			, y = value
		)
		, colour = 'transparent'
		, fill = 'black'
		, scale = 'width'
		, alpha = .5
	)
	+ geom_linerange(
		mapping = aes(
			x = variable
			, ymin = q0
			, ymax = q100
		)
		, colour = 'black'
		# , alpha = .5
		, size = .5
	)
	+ geom_linerange(
		mapping = aes(
			x = variable
			, ymin = q6.5
			, ymax = q93.5
			, colour = ess_tail_is_bad
		)
	)
	+ geom_linerange(
		mapping = aes(
			x = variable
			, ymin = q25
			, ymax = q75
			, colour = ess_bulk_is_bad
		)
		, size = 3
	)
	+ geom_point(
		mapping = aes(
			x = variable
			, y = q50
			, fill = rhat_is_bad
		)
		, shape = 21
		, size = 2
		, colour = 'black'
	)
	+ geom_hline(yintercept=0,linetype=3)
	+ coord_flip()
	+ scale_color_manual(
		values = lst(`TRUE`='red',`FALSE`='white')
		, labels = lst(`TRUE`='<100',`FALSE`='>=100')
	)
	+ scale_fill_manual(
		values = lst(`TRUE`='red',`FALSE`='black')
		, labels = lst(`TRUE`='>1.01',`FALSE`='<=1.01')
	)
	+ scale_y_continuous(limits=c(-1,1),expand=c(0,0),breaks=c(-.5,0.,.5))
	+ labs(
		y = "Posterior Value"
		, x = 'Correlation'
		, colour = 'ESS'
		, fill = 'Rhat'
	)
	+ theme(
		aspect.ratio = 2/(1+sqrt(5))
	)
)


# Viz cross-type correlations ----
(
	post$draws('locat_star_cor')
	%>% posterior::as_draws_df()
	%>% select(-.draw)
	%>% pivot_longer(
		cols = -c(.chain,.iteration)
		, names_to = 'variable'
	)
	%>% group_nest(variable)
	%>% left_join(par_summary,by='variable')
	%>% separate(
		variable
		, into=c('variable','index')
		, sep='\\['
	)
	%>% separate(
		index
		, into=c('index1','index2')
		, sep=','
	)
	%>% mutate(
		index1 = as.numeric(index1)
		, index2 = as.numeric(str_remove(index2,fixed(']')))
		, variable = (
			complete_binom_Xc_with_vars
			%>% select(contrasts)
			%>% unnest(contrasts)
			%>% names()
		)[index1]
		, type = c('Location-Scale','Location-Error')[index2]
	)
	%>% mutate(
		rhat_is_bad = 1.01<rhat
		, ess_bulk_is_bad = 100>ess_bulk
		, ess_tail_is_bad = 100>ess_tail
	)
) ->
	to_plot

(
	to_plot
	%>% ggplot()
	+ facet_wrap(~type,ncol=3)
	+ geom_hline(yintercept = 0,linetype=3,colour='white')
	+ geom_violin(
		data = (. %>% select(type,variable,data) %>% unnest(data))
		, mapping = aes(
			x = variable
			, y = value
		)
		, colour = 'transparent'
		, fill = 'black'
		, scale = 'width'
		, alpha = .5
	)
	+ geom_linerange(
		mapping = aes(
			x = variable
			, ymin = q0
			, ymax = q100
		)
		, colour = 'black'
		# , alpha = .5
		, size = .5
	)
	+ geom_linerange(
		mapping = aes(
			x = variable
			, ymin = q6.5
			, ymax = q93.5
			, colour = ess_tail_is_bad
		)
	)
	+ geom_linerange(
		mapping = aes(
			x = variable
			, ymin = q25
			, ymax = q75
			, colour = ess_bulk_is_bad
		)
		, size = 3
	)
	+ geom_point(
		mapping = aes(
			x = variable
			, y = q50
			, fill = rhat_is_bad
		)
		, shape = 21
		, size = 2
		, colour = 'black'
	)
	+ geom_hline(yintercept=0,linetype=3)
	+ coord_flip()
	+ scale_color_manual(
		values = lst(`TRUE`='red',`FALSE`='white')
		, labels = lst(`TRUE`='<100',`FALSE`='>=100')
	)
	+ scale_fill_manual(
		values = lst(`TRUE`='red',`FALSE`='black')
		, labels = lst(`TRUE`='>1.01',`FALSE`='<=1.01')
	)
	+ scale_y_continuous(limits=c(-1,1),expand=c(0,0),breaks=c(-.5,0.,.5))
	+ labs(
		y = "Posterior Value of Correlation"
		, x = 'Predictor'
		, colour = 'ESS'
		, fill = 'Rhat'
	)
	+ theme(
		aspect.ratio = 2/(1+sqrt(5))
	)
)


#TODO: use generate quantities to get coef_mean, coef_sd, *_dot_Xc, etc

