# Number of 384 well plates
plates <- 5
# Supply data frame with raw CTG values.  
df <- bmn
#Indicate column numbers (starting at 1) for DMSO and drug replicates:
dmsocols <- c(5,6,7)
drugcols <- c(9,10,11)
# Column with plate numbers
platecol <- 1


# Make sure df is sorted by plate
dfo <- df[order(df[,platecol]),]
# Log transform
dfl <- dfo
for (i in c(dmsocols,drugcols)) {
	dfl[,i] <- log(dfl[,i],2)
 
## Get plate medians of log transformed data
# Make a matrix: plates as rows, reps as columns
dmso_plate_medians <- matrix(nrow=plates,ncol=length(dmsocols))
drug_plate_medians <- matrix(nrow=plates,ncol=length(drugcols))
# Populate the matrix with plate medians.  This allows for different numbers of drug/dmso reps.
for (i in 1:plates) {
	for (j in 1:length(dmsocols)) {
		dmso_plate_medians[i,j] <- median(dfl[dfl[,platecol] == i,dmsocols[j]],na.rm=T)
	}
	for (j in 1:length(drugcols)) {
		drug_plate_medians[i,j] <- median(dfl[dfl[,platecol] == i,drugcols[j]],na.rm=T)
	}
}
## Normalise data by plate centering
dfn <- dfl
# counter for replicates
n = 1
for (i in dmsocols) {
        dfn[,i] <- dfn[,i] - dmso_plate_medians[dfn[,platecol],n]
	n = n + 1
}
n = 1
for (i in drugcols) {
        dfn[,i] <- dfn[,i] - drug_plate_medians[dfn[,platecol],n]
	n = n + 1
}
## Take medians
dfn$drugmed <- 0
dfn$dmsomed <- 0
for (i in 1:nrow(dfn)) {
						# Median returns a list for some reason?
	dfn$drugmed[i] <- median(dfn[i,drugcols])[[1]]
	dfn$dmsomed[i] <- median(dfn[i,dmsocols])[[1]]
}
## Calculate drug effect

dfn$ad_dmso <- abs(dfn$dmsomed - median(dfn$dmsomed))
dfn$ad_drug <- abs(dfn$drugmed - median(dfn$drugmed))
dfn$z_dmso <- (dfn$dmsomed - median(dfn$dmsomed))/median(dfn$ad_dmso)
dfn$z_drug <- (dfn$drugmed - median(dfn$drugmed))/median(dfn$ad_drug)
dfn$de <- dfn$z_drug - dfn$z_dmso

