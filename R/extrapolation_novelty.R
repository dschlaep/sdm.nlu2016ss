# Functions developed by DRS for Bell DM, Schlaepfer DR (2016) On the dangers of model complexity without ecological justification in species distribution modeling. Ecological Modelling, 330, 50-59.
# see https://github.com/bellland/SDM.Virtual.Species_Bell.Schlaepfer

#' @export
calc_NT1 <- function(refdat, prodat) {
	stopifnot(identical(colnames(refdat), colnames(prodat)))
	
	#——————————————————————–#
	# NT1 – UNIVARIATE EXTRAPOLATION
	#——————————————————————–#
	# Mesgaran, M. B., R. D. Cousens, B. L. Webber, and J. Franklin. 2014. Here be dragons: a tool for quantifying novelty due to covariate range and correlation change when projecting species distribution models. Diversity and Distributions 20:1147-1159.
	# code based on comment by Matthew Bayly to https://pvanb.wordpress.com/2014/05/13/a-new-method-and-tool-exdet-to-evaluate-novelty-environmental-conditions/

	range_ref <- t(matrixStats::colRanges(refdat))
	#dimnames(range_ref) <- list(c("min", "max"), colnames(refdat))
	range_ref_arr <- array(range_ref, dim = c(dim(range_ref), nrow(prodat)), dimnames = list(c("min", "max"), colnames(refdat), NULL))

	diffs_ref <- matrixStats::colDiffs(range_ref)
	#colnames(diffs_ref) <- colnames(refdat)
	diffs_ref_arr <- matrix(diffs_ref, nrow = nrow(prodat), ncol = ncol(prodat), byrow = TRUE)

	iud <- array(0, dim = c(dim(prodat), 3))
	iud[, , 2] <- prodat - t(range_ref_arr["min", ,])
	iud[, , 3] <- t(range_ref_arr["max", ,]) - prodat

	UDs <- apply(iud, 1:2, min) / diffs_ref_arr
	NT1 <- rowSums(UDs)
}

#' @export
calc_NT2 <- function(refdat, prodat) {
	stopifnot(identical(colnames(refdat), colnames(prodat)))

	#——————————————————————–#
	# NT2 – MULTIVARIATE EXTRAPOLATION
	#——————————————————————–#
	# Mesgaran, M. B., R. D. Cousens, B. L. Webber, and J. Franklin. 2014. Here be dragons: a tool for quantifying novelty due to covariate range and correlation change when projecting species distribution models. Diversity and Distributions 20:1147-1159.
	# code modified from on https://pvanb.wordpress.com/2014/05/13/a-new-method-and-tool-exdet-to-evaluate-novelty-environmental-conditions/

	# Calculate the average and covariance matrix of the variables 
	# in the reference set
	ref.av  <- colMeans(refdat, na.rm = TRUE)
	ref.cov <- var(refdat, na.rm = TRUE)
 
	# Calculate the mahalanobis distance of each raster 
	# cell to the environmental center of the reference 
	# set for both the reference and the projection data 
	# set and calculate the ratio between the two.
	mah.ref <- mahalanobis(x = refdat, center = ref.av, cov = ref.cov)
	mah.pro <- mahalanobis(x = prodat, center = ref.av, cov = ref.cov)
	mah.max <- max(mah.ref[is.finite(mah.ref)])
	NT2 <- mah.pro / mah.max
}
