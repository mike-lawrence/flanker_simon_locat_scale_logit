cat('fix_gridlines() now available as a function\n')
fix_gridlines = function(p){
	p<<-p
	pp <<- ggplot_build(p)$layout$panel_params[[1]]
	grid_coords = lst(
		x_major = b$layout$panel_params[[1]]$y$breaks
		, x_minor = b$layout$panel_params[[1]]$y$minor_breaks
		, y_major = b$layout$panel_params[[1]]$x$breaks
		, y_minor = b$layout$panel_params[[1]]$x$minor_breaks
	)
	(
		p
		# + theme(
		# 	panel.grid = element_blank()
		# )
		# + geom_vline(
		# 	xintercept = grid_coords$x_major
		# 	, size = 1
		# 	, colour = 'white'
		# 	# , alpha = .5
		# )
		# + geom_vline(
		# 	xintercept = grid_coords$x_minor
		# 	, size = .5
		# 	, colour = 'white'
		# 	, alpha = .5
		# )
		+ geom_hline(
			yintercept = grid_coords$y_major
			, size = 1
			, colour = 'white'
			, alpha = .5
		)
		+ geom_hline(
			yintercept = grid_coords$y_minor
			, size = .5
			, colour = 'white'
			, alpha = .5
		)
	) ->
		p
	return(p)
}
