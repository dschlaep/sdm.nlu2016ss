library(sdm.nlu2016ss)

my_data <- prepare_installed_data()
my_glm <- fit_sdm_glm(my_data$obs)
my_rf <- fit_sdm_rf(my_data$obs)
my_glm_preds <- predict_sdm(my_glm, my_data$newdata)
my_rf_preds <- predict_sdm(my_rf, my_data$newdata)
plot_response_function(my_glm, var = “bio1”)
plot_response_function(my_rf, var = “bio1”)
