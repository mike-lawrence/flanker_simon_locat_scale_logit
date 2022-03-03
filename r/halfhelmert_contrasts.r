cat('halfhelmert_contrasts() now available as a function\n')
halfhelmert_contrasts = function(...){
	contr.helmert(...)*.5
}
