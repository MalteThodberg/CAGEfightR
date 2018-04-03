setGeneric("rowsum2", function(x, group, ...){
	standardGeneric("rowsum2")
})

setMethod("rowsum2", signature(x="matrix", group="factor"),
					function(x, group, drop=TRUE, sparse=FALSE){
	# Pre-checks
	stopifnot(nrow(x) == length(group),
						length(sparse) == 1,
						length(drop) == 1)

	if(!drop | sparse){
		warning("Aggregating base::matrix always drops unused levels",
						" and return a base::matrix!")
	}

	# Aggregate
	o <- suppressWarnings(rowsum(x, group))

	# Discard NAs
	if(anyNA(group)){
		o <- o[-nrow(o),]
	}

	# Check output dimensions and names match
	group <- factor(group)
	stopifnot(nrow(o) == length(levels(group)),
						rownames(o) == levels(group))

	# Return
	o
})

setMethod("rowsum2", signature(x="dgCMatrix", group="factor"),
					function(x, group, drop=FALSE, sparse=FALSE){
	# Pre-checks
	stopifnot(nrow(x) == length(group),
						length(sparse) == 1,
						length(drop) == 1)


	# Aggregate
	o <- aggregate.Matrix2(x=x, groupings=group, drop.unused.levels=drop)

	# Discard NAs
	if(anyNA(group)){
		o <- o[-nrow(o),]
	}

	# To matrix
	if(!sparse){
		# To matrix
		o <- Matrix::as.matrix(o)

		# Compress to integer if possible
		if(all(o == floor(o))){
			storage.mode(o) <- "integer"
		}
	}

	# Post-checks
	if(drop){
		stopifnot(nrow(o) <= length(levels(group)),
							all(rownames(o) %in% levels(group)))
	}else{
		stopifnot(nrow(o) == length(levels(group)),
							setequal(rownames(o), levels(group)),
							all(rownames(o) == levels(group)))
	}

	# Return
	o
})

aggregate.Matrix2 <- function (x, groupings = NULL, form = NULL, fun = "sum",
															 drop.unused.levels=TRUE, ...){
	if (!methods::is(x, "Matrix"))
		x <- Matrix::Matrix(Matrix::as.matrix(x), sparse = TRUE)
	if (fun == "count")
		x <- x != 0
	groupings2 <- groupings
	if (!methods::is(groupings2, "data.frame"))
		#groupings2 <- methods::as(groupings2, "data.frame")
		groupings2 <- as.data.frame(groupings2)
	groupings2 <- data.frame(lapply(groupings2, as.factor))
	groupings2 <- data.frame(interaction(groupings2, sep = "_"))
	colnames(groupings2) <- "A"
	if (is.null(form))
		form <- stats::as.formula("~0+.")
	form <- stats::as.formula(form)
	mapping <- Matrix.utils::dMcast(groupings2, form,
																	drop.unused.levels=drop.unused.levels)
	colnames(mapping) <- substring(colnames(mapping), 2)
	result <- t(mapping) %*% x
	if (fun == "mean")
		result@x <- result@x/(Matrix.utils::aggregate.Matrix(x, groupings2,
																					 fun = "count"))@x
	attr(result, "crosswalk") <- grr::extract(groupings, match(rownames(result),
																														 groupings2$A))
	return(result)
}
