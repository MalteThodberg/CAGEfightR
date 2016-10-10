# Generic for calculating TPM
calcTPM <- function(x) UseMethod("calcTPM")

# For single vector
calcTPM.default <- function(x){
	x / (sum(x) / 1e6)
}

# For GRanges
calcTPM.GRangesList <- function(x){
	# Calculate TPM of score column
	score(x) <- calcTPM(score(x))

	# Return x
	x
}

# For GRanges
calcTPM.GRangesList <- function(x){
	# Turn into data.table
	dts <- as.data.table(x)

	# Calculate TPM for each level


}

