#' @export
assess_model_performance <- function(obs, pred) {
	res <- list(auc = NA, sensitivity = NA, specificity = NA, tss = NA, tss_max = NA, threshold_tss_max = NA)
	
	if (requireNamespace("ROCR", quietly = TRUE)) {
		# AUC
		temp <- ROCR::prediction(pred, obs)
		res[["auc"]] <- slot(ROCR::performance(temp, "auc"), "y.values")[[1]][1]

		# Sensitivity, specificity, and TSS
		ss <- ROCR::performance(temp, "spec", "sens")
		thresholds <- slot(ss, "alpha.values")[[1]]
		res[["sens"]] <- slot(ss, "x.values")[[1]]
		res[["spec"]] <- slot(ss, "y.values")[[1]]

		# Maximize TSS = sensitivity + specificity - 1
		res[["tss"]] <- res[["sens"]] + res[["spec"]] - 1
		i_max_tss <- which(res[["tss"]] == max(res[["tss"]]))
		res[["tss_max"]] <- mean(tss[i_max_tss])
		res[["tss_max"]] <- mean(thresholds[i_max_tss])
		
	} else {
		warning("The package 'ROCR' is not available; model performance will not be calculated without this package")
	}
	
	res		
}


#' @export
plot_ROC <- function(sens, spec, filename = NULL, ...) {
	if (!is.null(filename)) pdf(file = filename)
	
	plot(1 - spec, sens, type = "l",
			xlim = c(0, 1), ylim = c(0, 1),
			xlab = "1 - specificity", ylab = "sensitivity",
			main = "Receiver Operating Characteristic",
			...)

	if (!is.null(filename)) dev.off()
	
	invisible(TRUE)
}
