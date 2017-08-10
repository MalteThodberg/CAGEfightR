# ### General
#
# has_names <- function(x){
# 	!is.null(names(x))
# }
# on_failure(has_names) <- function(call, env) {
# 	paste0(deparse(call$x), " does not have names (see ?names)")
# }
#
# has_rownames <- function(x){
# 	!is.null(rownames(x))
# }
# on_failure(has_rownames) <- function(call, env) {
# 	paste0(deparse(call$x), " does not have row names (see ?rownames)")
# }
#
# has_colnames <- function(x){
# 	!is.null(rownames(x))
# }
# on_failure(has_colnames) <- function(call, env) {
# 	paste0(deparse(call$x), " does not have column names (see ?colnames)")
# }
#
# has_scores <- function(x){
# 	!is.null(score(x))
# }
# on_failure(has_scores) <- function(call, env) {
# 	paste0(deparse(call$x), " does not have scores (see ?score)")
# }
#
# has_seqlengths <- function(x){
# 	length(seqlengths(x)) != 0
# }
# on_failure(has_seqlengths) <- function(call, env) {
# 	paste0(deparse(call$x), " does not have seqlengths (see ?seqinfo)")
# }
#
# ### GRanges
#
# is_GRanges <- function(x){
# 	isS4(x) & class(x) == "GRanges"
# }
# on_failure(is_GRanges) <- function(call, env) {
# 	paste0(deparse(call$x), " is not a GRanges")
# }
#
# is_disjoint <- function(x){
# 	assert_that(is_GRanges(x),
# 							isDisjoint((x)))
# }
# on_failure(is_disjoint) <- function(call, env) {
# 	paste0(deparse(call$x), " is not disjoint (see ?isDisjoint)")
# }
#
#
# ### GRangesList
#
# is_GRangesList <- function(x){
# 	isS4(x) & class(x) == "GRangesList"
# }
# on_failure(is_GRangesList) <- function(call, env) {
# 	paste0(deparse(call$x), " is not a GRangesList")
# }
#
# ### Deep: CTSS
#
# # is_ctss <- function(x){
# # 	assert_that(is_disjoint(x))
# # 	assert_that(has_scores(x))
# #
# # }
# # on_failure(is_ctss) <- function(call, env) {
# # 	paste0(deparse(call$x), " is not a CTSS.")
# # }
#
# is_ctss <- function(x){
# 	# assert_that(is_GRangesList(x),
# 	# 						has_seqlengths(x))
# 	i <- isDisjoint(x) | vapply(x, has_scores, logical(1))
# 	all(i)
# }
# on_failure(is_ctss) <- function(call, env) {
# 	res <- eval(call[[2]], env)
# 	i <- isDisjoint(res) | vapply(res, has_scores, logical(1))
# 	paste0("Elements ", paste(i, collapse = ", "), " of ",
# 				 deparse(call[[2]]), " are not valid CTSS")
# }
#
# is_ctss(x)
#
# ### BigWigFile
#
# is_BigWigFile <- function(x){
# 	isS4(x) & class(x) == "BigWigFile"
# }
# on_failure(is_BigWigFile) <- function(call, env) {
# 	paste0(deparse(call$x), " is not a BigWigFile")
# }
#
# ### BigWigList
#
# is_BigWigList <- function(x){
# 	isS4(x) & class(x) == "BigWigFileList"
# }
# on_failure(is_BigWigList) <- function(call, env) {
# 	paste0(deparse(call$x), " is not a BigWigFileList")
# }
#
# ### Other
#
# #
# # is_IRanges <- function(x){
# # 	isS4(x) & class(x) == "IRanges"
# # }
# # on_failure(is_IRanges) <- function(call, env) {
# # 	paste0(deparse(call$x), " is not an IRanges")
# # }
# #
# # has_seqlengths <- function(x){
# # 	assert_that(is_GRanges(x))
# # 	!is.null(seqlengths(x) )
# # }
# # on_failure(has_seqlengths) <- function(call, env) {
# # 	paste0(deparse(call$x), " does not have seqlengths")
# # }
# #
# # has_column <- function(x, y){
# # 	assert_that(is_GRanges(x))
# # 	assert_that(is.string(y))
# # 	y %in% colnames(mcols(x))
# #
# # }
# # on_failure(has_column) <- function(call, env) {
# # 	paste0(deparse(call$x), " does not have column: ", deparse(call$y))
# # }
# #
# # has_peaks <- function(x, y="thick"){
# # 	assert_that(has_column(x, y))
# # 	column <- mcols(x)[,y]
# # 	assert_that(is_IRanges(column))
# # 	all(start(x) <= start(column)) & all(end(column) >= end(column))
# # }
# # on_failure(has_peaks) <- function(call, env) {
# # 	paste0(deparse(call$x), " has peaks stored in ", deparse(call$y), "outside of ranges!")
# # }
# #
# # z <-  rowRanges(SE)[1:1000]
