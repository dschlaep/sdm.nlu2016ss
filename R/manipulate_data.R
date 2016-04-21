validate_sdm_data <- function(data, obs_coln) {
	# check that observations are binary
	temp <- factor(data[, obs_coln])
	if (!(nlevels(temp) == 2L))
		stop("Observations are not binary!")
		
	if (anyNA(data))
		warning("Observations and/or environmental variables contain NAs!")
		
	invisible(TRUE)
}

#' Import data
#'
#' @param filename A character string representing a path to a csv file.
#' @param obs_coln The name of the column with the species observations.
#' @param var_coln A vector with column names that will be used as predictors.
#' @param na.rm A logical value.
#' @seealso \code{\link{read.csv}}
#' @return A list with the elements 'data', 'subset', ...
#' @export
import_sdm_data <- function(filename, obs_coln, var_coln, na.rm = FALSE) {
	# Load data
	temp_dat <- read.csv(file = filename, header = TRUE)
	
	# Check that requested columns are available
	temp_novar <- is.na(sapply(var_coln, function(x) grep(x, colnames(temp_dat))[1]))
	if (any(temp_novar)) {
		warning("Some of the requested predictors (", paste(var_coln[temp_novar], collapse = ", "), 
				") cannot be found among the fields of ", basename(filename), ": ", paste(colnames(temp_dat), collapse = ", "))
		var_coln2 <- var_coln[!temp_novar]
	} else {
		var_coln2 <- var_coln
	}
	
	# Subset to complete cases
	if (na.rm) {
		temp_subset <- complete.cases(temp_dat[, c(obs_coln, var_coln2)])
	} else {
		temp_subset <- rep(TRUE, nrow(temp_dat))
	}
	
	# Validate observations
	validate_sdm_data(temp_dat[temp_subset, c(obs_coln, var_coln2)], obs_coln)
	
	# Return
	list(data = temp_dat, subset = temp_subset, obs_coln = obs_coln, var_coln = var_coln2)
}


#' @export
download_obs_from_gbif <- function(species) {
	species_records <- as.data.frame(matrix(NA, nrow = 0, ncol = 0))

	if (requireNamespace("spocc", quietly = TRUE)) {
		# ask how many records available
		temp <- spocc::occ(query = species, from = c("gbif"), limit = 5, gbifopts = list(hasCoordinate = TRUE, fields = "all", spatialIssues = FALSE))
		cat(format(Sys.time(), format = ""), species, " get ", temp$gbif$meta$found, " records from gbif\n")
		
		# get all data
		if(temp$gbif$meta$found  ==  0) return(species_records)
		species_occ <- spocc::occ(query = species, from = c("gbif"), limit = temp$gbif$meta$found, gbifopts = list(hasCoordinate = TRUE, fields = "all", spatialIssues = FALSE))
		
		# organize data
		temp <- lapply(species_occ$gbif$data, FUN = function(x) with(x, data.frame(
								ID = gbifID,
								FAMILY = family,
								NAME = scientificName,
								LONGITUDE = as.numeric(longitude),
								LATITUDE = as.numeric(latitude),
								COUNTRY = country,
								YEAR = as.integer(year),
								issues = issues,
								geodeticDatum = geodeticDatum,
								protocol = protocol,
								catalogNumber = catalogNumber,
								taxonRank = taxonRank,
								genericName = genericName,
								specificEpithet = specificEpithet,
								stringsAsFactors = FALSE)))
		species_records <- do.call(rbind, temp)

	} else {
		warning("The package 'spocc' is not available; data from 'gbif' cannot be accessed without this package")
		species_records <- as.data.frame(matrix(NA, nrow = 0, ncol = 0))
	}
	
	species_records
}


#' @export
extract_from_BioClim <- function(xy, dir_data, dir_bc = NULL, ids_bio = NULL, cur_fut = c("current", "future"), rcp = NULL, GCM = NULL){
	# load 19 (ids_bio subset) bioclim variables from the current WorldClim dataset
	#	- variable definition: 'Bioclim | WorldClim - Global Climate Data.pdf'
	env_data <- brick_bioclim <- list()
	
	if (requireNamespace("raster", quietly = TRUE)) {
		if (is.null(ids_bio)) ids_bio <- 1:19
		cur_fut <- match.arg(cur_fut)
		
		# Find directory and files
		if (is.null(dir_bc)) {
			if (cur_fut == "current") {
				dir_bc <- file.path(dir_data, "Grids_Current_1960-1990_10min", "bio_10m_bil")
			} else {
				dir_bc_temp <- file.path(dir_data, "Grids_Future_2061-2080_10min", rcp)
			}		
		} else {
			if (cur_fut == "future") dir_bc_temp <- file.path(dir_bc, "..")
		}
		
		if (cur_fut == "current") {
			file_ext <- ".bil"
			file_tag <- "bio"
		} else {
			file_ext <- ".tif"
			file_tag <- strsplit(tag_gcm, split="_", fixed = TRUE)[[1]][2]
			tag_gcm <- (temp <- list.files(temp))[grepl(GCM, temp)]
			dir_bc <- file.path(dir_bc_temp, tag_gcm)
		}
		
				
		# Load BioClim data
		bio_rasters <- list()
		for (i in ids_bio) {
			ftemp <- file.path(dir_bc, paste0(file_tag, i, ".bil"))
			if (file.exists(ftemp)) bio_rasters[[basename(ftemp)]] <- raster::raster(ftemp)
		}
		
		n <- length(bio_rasters)
		
		if (n < length(ids_bio))
			warning("Not all requested BioClim data layers were available")
		
		if (n > 0) {
			brick_args <- c(list(x = bio_rasters[[1]]), unlist(bio_rasters[-1]), list(nl = n))
			brick_bioclim <- do.call(getFromNamespace("brick", "raster"), args = brick_args)
			env_data <- raster::extract(x = brick_bioclim, y = xy, method = "bilinear", df = TRUE)
		}
	} else {
		warning("The package 'spocc' is not available; data from 'gbif' cannot be accessed without this package")
	}
	
	list(env_data = env_data, brick_bioclim = brick_bioclim)
}
