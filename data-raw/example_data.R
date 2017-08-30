library(CAGEfightR)
library(rtracklayer)

# Load data
data("exampleDesign")

bw_plus <- BigWigFileList(system.file("extdata", exampleDesign$BigWigPlus, package = "CAGEfightR"))
bw_minus <- BigWigFileList(system.file("extdata", exampleDesign$BigWigMinus, package = "CAGEfightR"))
names(bw_plus) <- exampleDesign$Name
names(bw_minus) <- exampleDesign$Name

### Wrapper runs

# Read in
CTSS <- readCTSS(bwPlus=bw_plus, bwMinus=bw_minus, genome="mm9")

# Coverage
CTSScoverage <- coverageOfCTSS(CTSS)

# Tuned TCs
tuningResults <- tuneTagClustering(CTSScoverage)
optimalTC <- min(subset(tuningResults, nTCs == max(nTCs), select=threshold))
TCs <- tagClustering(ctssCoverage=CTSScoverage, tpmCutoff=optimalTC)

# BCs
BCs <- bidirectionalClustering(ctssCoverage=CTSScoverage)

# Rename and save
exampleCTSS <- CTSS
exampleCoverage <- CTSScoverage
exampleTagClusters <- TCs
exampleBidirectionalClusters <- BCs

devtools::use_data(exampleCTSS, overwrite = TRUE)
devtools::use_data(exampleCoverage, overwrite = TRUE)
devtools::use_data(exampleTagClusters, overwrite = TRUE)
devtools::use_data(exampleBidirectionalClusters, overwrite = TRUE)

