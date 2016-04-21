



library("sdm.nlu2016ss")


fname_exp <- system.file("extdata", "SOGI_data1.csv", package = "sdm.nlu2016ss")


dat_SOGI <- import_sdm_data(filename = fname_exp, 
				obs_coln = "obs",
				var_coln = c("bio1", "bio4", "bio11"),
				na.rm = TRUE)

obs_OPEL <- download_obs_from_gbif(species = "Ophrys elatior")
env_OPEL <- extract_from_BioClim(xy = my_obs2[, c("LONGITUDE", "LATITUDE")], dir_data = "WorldClim")


my_glm <- fit_sdm_glm(my_data$obs)
my_data <- import_sdm_data()

my_glm <- fit_sdm_glm(my_data$obs)
my_rf <- fit_sdm_rf(my_data$obs)

my_glm_preds <- predict_sdm(my_glm, my_data$newdata)
my_rf_preds <- predict_sdm(my_rf, my_data$newdata)

plot_response_function(my_glm, var = "bio1")
plot_response_function(my_rf, var = "bio1")
