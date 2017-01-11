# library(data.table)
# library(GenomicRanges)
#
# chunkify <- function(l, n){
# 	split(l, factor(rep_len(seq_len(n), length.out=length(l))))
# }
#
# # doubleReduce <- function(X, innerFunction, outerFunction=innerFunction, innerInit=NULL, outerInit=NULL, workers=1){
# # 	if(workers == 1){
# # 		# Only first reduce
# # 		if(!is.null(innerInit)){
# # 			o <- Reduce(f=innerFunction, x=X, init=innerInit)
# # 		}else{
# # 			o <- Reduce(f=innerFunction, x=X)
# # 		}
# # 	}else{
# # 		# Make chunks
# # 		chunks <- chunkify(l=X, n=workers)
# #
# # 		# First Reduce
# # 		if(!is.null(innerInit)){
# # 			o <- parallel::mclapply(X=X, FUN=Reduce, f=innerFunction, init=innerInit, mc.cores=workers, mc.preschedule=TRUE)
# # 		}else{
# # 			o <- parallel::mclapply(X=X, FUN=Reduce, f=innerFunction, mc.cores=workers, mc.preschedule=TRUE)
# # 		}
# #
# # 		# Second reduce
# # 		if(!is.null(outerInit)){
# # 			o <- Reduce(f=outerFunction, x=o, init=outerInit)
# # 		}else{
# # 			o <- Reduce(f=outerFunction, x=o)
# # 		}
# # 	}
# #
# # 	# Return
# # 	o
# # }
#
# # l <- 1:100000
# # f1 <- function(x1, x2){ x1+x2}
# # f2 <- f1
# #
# # doubleReduce(X=l, innerFunction=f1, worker=2)
#
# doubleReduce <- function(X, stackingFunction, preprocessingFunction, init=NULL, workers=1){
# 	# Define new inner function
# 	innerFunction <- function(x1,x2){
# 		stackingFunction(x1, preprocessingFunction(x2))
# 	}
#
# 	if(workers == 1){
# 		# Only first reduce
# 		if(!is.null(init)){
# 			o <- Reduce(f=innerFunction, x=X, init=init)
# 		}else{
# 			o <- Reduce(f=innerFunction, x=X)
# 		}
# 	}else{
# 		# Make chunks
# 		chunks <- chunkify(l=X, n=workers)
#
# 		# First Reduce
# 		if(!is.null(init)){
# 			o <- parallel::mclapply(X=X, FUN=Reduce, f=innerFunction, init=init, mc.cores=workers, mc.preschedule=TRUE)
# 		}else{
# 			o <- parallel::mclapply(X=X, FUN=Reduce, f=innerFunction, mc.cores=workers, mc.preschedule=TRUE)
# 		}
#
# 		# Second reduce
# 		o <- Reduce(f=stackingFunction, x=o)
#
# 	}
#
# 	# Return
# 	o
# }
#
# # Tester files
# bedBoth <- list.files("~/Desktop/to_be_zipped/", full.names=TRUE)
# bedBoth <- rep(bedBoth, 20)
#
# er <- GRanges(score=numeric())
#
# stackGR <- function(gr1, gr2){
# 	c(gr1, gr2)
# }
#
# fastReadBED <- function(fname){
# 	fread(input=fname,
# 				data.table=TRUE,
# 				sep="\t",
# 				col.names=c("seqnames", "start", "end", "score", "strand"),
# 				select=c(1,2,3,5,6),
# 				verbose=FALSE,
# 				showProgress=FALSE,
# 				na.strings=NULL)
# }
#
# DT2GR <- function(dt){
# 	GRanges(seqnames=Rle(dt$seqnames),
# 					ranges=IRanges(dt$start+1, width=1),
# 					strand=Rle(dt$strand),
# 					score=dt$score)
# }
#
# BED2GR <- function(fname){
# 	DT2GR(fastReadBED(fname))
# }
#
# gr <- doubleReduce(X=bedBoth, stackingFunction=stackGR, preprocessingFunction=BED2GR, init=er, workers=3)
