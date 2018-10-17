.onLoad <- function(libname, pkgname) {
	op <- options()
	op.CAGEfightR <- list(
		CAGEfightR.round = 9,
		CAGEfightR.plus = "cornflowerblue",
		CAGEfightR.minus = "tomato",
		CAGEfightR.unstranded = "hotpink"
	)
	toset <- !(names(op.CAGEfightR) %in% names(op))
	if(any(toset)) options(op.CAGEfightR[toset])

	invisible()
}
