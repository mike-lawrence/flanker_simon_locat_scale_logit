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
			, formula = ~ flankers#*simon
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
			, formula = ~ flankers#*simon
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
mod$sample(
	data = data_for_stan
	, chains = parallel::detectCores()/2
	, parallel_chains = parallel::detectCores()/2
	, refresh = 1
	, init = 1
	, max_treedepth = 11
	, metric = 'dense_e'
)
