




library(tidyverse)
packages <- c("lubridate", "survival", "terra", "raster", "tictoc", "TwoStepCLogit")
walk(packages, require, character.only = T)

# import buffalo data
buffalo_data <- read_csv(file = "outputs/buffalo_parametric_popn_covs_20230208.csv")
buffalo_data <- buffalo_data %>% mutate(t1_ = lubridate::with_tz(buffalo_data$t1_, tzone = "Australia/Darwin"),
                                        t2_ = lubridate::with_tz(buffalo_data$t2_, tzone = "Australia/Darwin"))

# create relevant covariates for modelling and predictions
buffalo_CLR <- buffalo_data %>% 
  mutate(id_num = as.numeric(factor(id)), 
         step_id = step_id_, 
         x1 = x1_, x2 = x2_, 
         y1 = y1_, y2 = y2_, 
         t1 = t1_, 
         t1_rounded = round_date(buffalo_data$t1_, "hour"), 
         t2_rounded = round_date(buffalo_data$t2_, "hour"), 
         t2 = t2_, 
         yday = yday(t1_),
         sl = sl_, sl_scaled = scale(sl),
         log_sl = log_sl_, 
         ta = ta_, cos_ta = cos_ta_,
         year = year(t1_), 
         month = month(t1_), 
         hour = hour(t2_rounded),
         ndvi_scaled = scale(ndvi_temporal),
         ndvi_2 = ndvi_scaled[,1]^2,
         canopy_01 = canopy_cover/100,
         herby_scaled = scale(veg_herby),
         canopy_scaled = scale(canopy_01),
         canopy_2 = canopy_scaled[,1]^2,
         elev_scaled = scale(DEM_H_end),
         elev_delta_scaled = scale(elev_delta),
         elev_log_scaled = scale(elev_log),
         slope_scaled = scale(slope_end),
         yday_s1 = sin(2*pi*yday/365),
         yday_s2 = sin(4*pi*yday/365),
         yday_s3 = sin(6*pi*yday/365),
         yday_s4 = sin(8*pi*yday/365),
         yday_c1 = cos(2*pi*yday/365),
         yday_c2 = cos(4*pi*yday/365),
         yday_c3 = cos(6*pi*yday/365),
         yday_c4 = cos(8*pi*yday/365),
         hour_s1 = sin(2*pi*hour/24),
         hour_s2 = sin(4*pi*hour/24),
         hour_s3 = sin(6*pi*hour/24),
         hour_s4 = sin(8*pi*hour/24),
         hour_c1 = cos(2*pi*hour/24),
         hour_c2 = cos(4*pi*hour/24),
         hour_c3 = cos(6*pi*hour/24),
         hour_c4 = cos(8*pi*hour/24)) %>%
  drop_na(c(ndvi_temporal, veg_herby, canopy_01, sl_)) %>% 
  dplyr::select(!(burst_:case_))

# Filtering by animals that have more than a year of data
buffalo_year_ids <- c(2005, 2014, 2018, 2022, 2024, 2154, 2158, 2327, 2354, 2387)
buffalo_CLR_year <- buffalo_CLR %>% filter(id %in% buffalo_year_ids)
unique(buffalo_CLR_year$id)

buffalo_CLR_year_harmonics <- buffalo_CLR_year %>% filter(t1 < "2019-07-25 09:32:42 ACST")

# buffalo_CLR_year %>% ggplot(aes(x = t1, y = factor(id), colour = factor(id))) +
#   geom_point(alpha = 0.1) +
#   scale_y_discrete("Buffalo ID") +
#   scale_x_datetime("Date") +
#   scale_colour_viridis_d() +
#   theme_bw() +
#   theme(legend.position = "none")

# Movement parameters from random sampling
# population
# gamma$params$shape # 0.438167
# gamma$params$scale # 534.3507
# vonmises$params$kappa # 0.1848126
# vonmises$params$mu # 0

movement_parameters <- read_csv("outputs/buffalo_parametric_indv_movement_params_20230208.csv")
movement_parameters_year <- movement_parameters %>% filter(id %in% buffalo_year_ids)

# Import spatial layers (currently in the import spatial layers R file)

ndvi_stack_scaled <- rast("mapping/ndvi_projected_watermask_scaled_20230612.tif")
canopy_scaled <- rast("mapping/cropped rasters/canopy_scaled_by_buffalo_data_20230210.tif")
herby_scaled <- rast("mapping/cropped rasters/herby_scaled_by_buffalo_data_20230210.tif")

xmin <- round(min(buffalo_CLR_year_harmonics$x2), digits = -2)
xmax <- round(max(buffalo_CLR_year_harmonics$x2), digits = -2)
ymin <- round(min(buffalo_CLR_year_harmonics$y2), digits = -2)
ymax <- round(max(buffalo_CLR_year_harmonics$y2), digits = -2)

crop_extent <- ext(xmin, xmax, ymin, ymax)
# subset_extent <- terra::ext(4200, 54300, -1458800, -1409700)

ndvi_stack_cropped <- terra::crop(ndvi_stack_scaled, crop_extent)
canopy_scaled_cropped <- terra::crop(canopy_scaled, crop_extent)
herby_scaled_cropped <- terra::crop(herby_scaled, crop_extent)

# to set the origin at (0,0)
ext(ndvi_stack_cropped) <- c(xmin - xmin, xmax - xmin, ymin - ymin, ymax - ymin)
ext(canopy_scaled_cropped) <- c(xmin - xmin, xmax - xmin, ymin - ymin, ymax - ymin)
ext(herby_scaled_cropped) <- c(xmin - xmin, xmax - xmin, ymin - ymin, ymax - ymin)

# Boyce index function

ecospat.boyce2 <- function (fit, obs, nclass = 0, window.w = "default", res = 100, 
                            PEplot = TRUE, rm.duplicate = TRUE, method = "spearman") 
{
  boycei <- function(interval, obs, fit) {
    pi <- sum(as.numeric(obs >= interval[1] & obs <= interval[2]))/length(obs)
    ei <- sum(as.numeric(fit >= interval[1] & fit <= interval[2]))/length(fit)
    return(round(pi/ei, 10))
  }
  if (inherits(fit, "RasterLayer")) {
    if (is.data.frame(obs) || is.matrix(obs)) {
      obs <- raster::extract(fit, obs)
    }
    fit <- getValues(fit)
    fit <- fit[!is.na(fit)]
    obs <- obs[!is.na(obs)]
  }
  mini <- min(fit, obs)
  maxi <- max(fit, obs)
  if (length(nclass) == 1) {
    if (nclass == 0) {
      if (window.w == "default") {
        window.w <- (max(fit) - min(fit))/10
      }
      vec.mov <- seq(from = mini, to = maxi - window.w, 
                     by = (maxi - mini - window.w)/res)
      vec.mov[res + 1] <- vec.mov[res + 1] + 1
      interval <- cbind(vec.mov, vec.mov + window.w)
    }
    else {
      vec.mov <- seq(from = mini, to = maxi, by = (maxi - 
                                                     mini)/nclass)
      interval <- cbind(vec.mov, c(vec.mov[-1], maxi))
    }
  }
  else {
    vec.mov <- c(mini, sort(nclass[!nclass > maxi | nclass < 
                                     mini]))
    interval <- cbind(vec.mov, c(vec.mov[-1], maxi))
  }
  f <- apply(interval, 1, boycei, obs, fit)
  to.keep <- which(f != "NaN")
  f <- f[to.keep]
  if (length(f) < 2) {
    b <- NA
  }
  else {
    r <- 1:length(f)
    if (rm.duplicate == TRUE) {
      r <- c(1:length(f))[f != c(f[-1], TRUE)]
    }
    b <- cor(f[r], vec.mov[to.keep][r], method = method)
  }
  HS <- apply(interval, 1, sum)/2
  if (length(nclass) == 1 & nclass == 0) {
    HS[length(HS)] <- HS[length(HS)] - 1
  }
  HS <- HS[to.keep]
  if (PEplot == TRUE) {
    plot(HS, f, xlab = "Habitat suitability", ylab = "Predicted/Expected ratio", 
         col = "grey", cex = 0.75)
    points(HS[r], f[r], pch = 19, cex = 0.75)
  }
  return(list(F.ratio = f, cor = round(b, 3), HS = HS))
}

# Formula with 1 pair of harmonics 

formula_twostep <- y ~ 
  
  ndvi_scaled +
  yday_s1:ndvi_scaled +
  yday_c1:ndvi_scaled +
  
  ndvi_2 +
  yday_s1:ndvi_2 +
  yday_c1:ndvi_2 +
  
  canopy_scaled +
  yday_s1:canopy_scaled +
  yday_c1:canopy_scaled +
  
  canopy_2 +
  yday_s1:canopy_2 +
  yday_c1:canopy_2 +
  
  herby_scaled +
  yday_s1:herby_scaled +
  yday_c1:herby_scaled +
  
  sl +
  yday_s1:sl +
  yday_c1:sl +
  
  log_sl +
  yday_s1:log_sl +
  yday_c1:log_sl +
  
  cos_ta +
  yday_s1:cos_ta +
  yday_c1:cos_ta +
  
  strata(step_id) +
  cluster(id)


# Jack knife validation

buffalo_ids <- unique(buffalo_CLR_year_harmonics$id)

data_minus1_list <- vector(mode = "list", length = 10)

# create a list of data frames that contain all individuals besides 1, sequentially for each id
for (i in 1:10) {
  data_minus1_list[[i]] <- buffalo_CLR_year_harmonics %>% filter(id != buffalo_ids[i])
}

# Using a loop instead of mapping over the function (mainly to assess progress)

indv_models <- vector(mode = "list", length = length(buffalo_ids))

for(i in 1:length(buffalo_ids)) {
  
  data_ssa <- data_minus1_list[[i]]
  
  tic()
  indv_models[[i]] <- Ts.estim(formula = formula_twostep,
                               data = data_ssa,
                               all.m.1 = TRUE,
                               D = "UN(1)",
                               itermax = 2000) 
  toc()
  
}

indv_models[[1]] # look at a single model

# Reconstructing coefficients with two harmonics with quadratics

yday <- seq(0,365,1)

yday_harmonics_df <- data.frame("yday_s1" = sin(2*pi*yday/365),
                                "yday_c1" = cos(2*pi*yday/365))

# get the names of the covariates
coef_names <- names(indv_models[[1]]$beta)[1:8]

# create list to store reconstructed harmonic cofficients from each model
indv_models_harm_list <- vector(mode = "list", length = length(coef_names))

for(j in 1:length(indv_models)) {
  
  # create temporary list to store coefficient values from each model
  harmonic_coefs_list <- vector(mode = "list", length = length(coef_names))
  
  for(i in 1:length(coef_names)) {
    
    harmonic_coefs_list[[i]] <- c(indv_models[[j]]$beta[which(names(indv_models[[j]]$beta) == paste0(coef_names[i]))],
                                  indv_models[[j]]$beta[which(names(indv_models[[j]]$beta) == paste0(coef_names[i], ":yday_s1"))],
                                  indv_models[[j]]$beta[which(names(indv_models[[j]]$beta) == paste0(coef_names[i], ":yday_c1"))])
    
  }
  
  indv_models_harm_list[[j]] <- data.frame("yday" = yday,
                                           "ndvi" = (harmonic_coefs_list[[1]][1] + as.matrix(yday_harmonics_df) %*% harmonic_coefs_list[[1]][2:3]),
                                           "ndvi_quad" = (harmonic_coefs_list[[2]][1] + as.matrix(yday_harmonics_df) %*% harmonic_coefs_list[[2]][2:3]),
                                           "canopy" = (harmonic_coefs_list[[3]][1] + as.matrix(yday_harmonics_df) %*% harmonic_coefs_list[[3]][2:3]),
                                           "canopy_quad" = (harmonic_coefs_list[[4]][1] + as.matrix(yday_harmonics_df) %*% harmonic_coefs_list[[4]][2:3]),
                                           "herby" = (harmonic_coefs_list[[5]][1] + as.matrix(yday_harmonics_df) %*% harmonic_coefs_list[[5]][2:3]),
                                           "sl" = (harmonic_coefs_list[[6]][1] + as.matrix(yday_harmonics_df) %*% harmonic_coefs_list[[6]][2:3]),
                                           "log_sl" = (harmonic_coefs_list[[7]][1] + as.matrix(yday_harmonics_df) %*% harmonic_coefs_list[[7]][2:3]),
                                           "cos_ta" = (harmonic_coefs_list[[8]][1] + as.matrix(yday_harmonics_df) %*% harmonic_coefs_list[[8]][2:3]))
  
}


# Reconstructing the Gamma and von Mises distributions from the tentative distributions (from Fieberg et al 2021: Appendix C)

# for individuals
# movement_parameters <- read_csv("outputs/buffalo_parametric_indv_movement_params_20230208.csv")
# movement_parameters_year <- movement_parameters %>% filter(id %in% buffalo_year_ids)

# for the population
# gamma$params$shape # 0.438167
# gamma$params$scale # 534.3507
# vonmises$params$kappa # 0.1848126

tentative_shape <- 0.438167
tentative_scale <- 534.3507
tentative_kappa <- 0.1848126

indv_models_harm_long_list <- vector(mode = "list", length = length(buffalo_ids))

for(i in 1:length(buffalo_ids)) {
  
  indv_models_harm_list[[i]] <- indv_models_harm_list[[i]] %>% mutate(shape = tentative_shape + log_sl,
                                                                      scale = 1/((1/tentative_scale) - sl),
                                                                      kappa = tentative_kappa + cos_ta)
  
  
  
  # write_csv(harmonic_coefs_df, "outputs/harmonic_coefs_df_with_movement_params_ndvi_canopy_quads_20230519.csv")
  
  # turning into a long data frame
  indv_models_harm_long_list[[i]] <- pivot_longer(indv_models_harm_list[[i]], cols = !1, names_to = "coef")
  
}

indv_models_harm_long <- cbind("id_excl" = factor(rep(buffalo_ids, each = nrow(indv_models_harm_long_list[[1]]))), bind_rows(indv_models_harm_long_list))

# write_csv(harmonic_coefs, "outputs/harmonic_coefs_df_long_with_movement_params_ndvi_canopy_quads_20230519.csv")

coefs <- unique(indv_models_harm_long$coef)
# coef_titles <- c("NDVI", "Slope", "Herbaceous vegetation", "Canopy cover", "Step length", "log(Step length)", "cos(Turning angle)")

ggplot() +
  geom_line(data = indv_models_harm_long %>% 
              # filter(id_excl == 2014) %>%
              filter(coef %in% unique(indv_models_harm_long$coef)[1:5]),
            aes(x = yday, y = value, colour = coef, group = interaction(id_excl, coef))) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(expression(Time-varying~parameter~values~beta)) +
  scale_x_continuous("Day of the year") +
  scale_color_discrete("Estimate") +
  # ggtitle(coef_titles[i]) +
  theme_classic() +
  theme(legend.position = "bottom")


# ggsave(paste("outputs/plots/clr_fitting/clr_1harmonics_all_minus1_20230608.png", sep = ""),
#   width=150, height=90, units="mm", dpi = 300)
# ggsave(paste("outputs/plots/clr_fitting/clr_1harmonic_all_pres_20230208.png", sep = ""),
#   width=300, height=180, units="mm", dpi = 300)

ggplot() +
  geom_path(data = indv_models_harm_long %>% 
              filter(coef %in% unique(indv_models_harm_long$coef)[c(9,11)]), 
            aes(x = yday, y = value, colour = coef, group = interaction(id_excl, coef))) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(expression(beta)) +
  scale_x_continuous("Day of the year") +
  scale_color_discrete("Estimate") +
  # ggtitle(coef_titles[i]) +
  theme_classic() +
  theme(legend.position = "bottom")

ggplot() +
  geom_path(data = indv_models_harm_long %>% 
              filter(coef %in% unique(indv_models_harm_long$coef)[10]), 
            aes(x = yday, y = value, colour = coef, group = interaction(id_excl, coef))) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(expression(beta)) +
  scale_x_continuous("Day of the year") +
  scale_color_discrete("Estimate") +
  # ggtitle(coef_titles[i]) +
  theme_classic() +
  theme(legend.position = "bottom")


# generating naive predictions - creating a data frame of monthly coefficients

years <- c(rep(2019, 7), rep(2018, 5))
months <- c(seq(1, 7, 1), seq(8, 12, 1))
day <- rep(1, 12)
dates <- make_datetime(year = years, month = months, day = day, tz = "Australia/Queensland")

# day of the middle of each month
mid_month <- c(16, 47, 75, 106, 136, 167, 197, 228, 259, 289, 320, 350)

monthly_coefs_list <- vector(mode = "list", length = length(buffalo_ids))
monthly_coefs_long_list <- vector(mode = "list", length = length(buffalo_ids))

for(i in 1:length(buffalo_ids)) {
  
  # get the coefficient at the middle of each month (15th day)
  monthly_coefs_list[[i]] <- data.frame("date" = dates,
                                        "ndvi" = indv_models_harm_list[[i]]$ndvi[mid_month],
                                        "ndvi_quad" = indv_models_harm_list[[i]]$ndvi_quad[mid_month],
                                        "canopy" = indv_models_harm_list[[i]]$canopy[mid_month],
                                        "canopy_quad" = indv_models_harm_list[[i]]$canopy_quad[mid_month],
                                        "herby" = indv_models_harm_list[[i]]$herby[mid_month],
                                        "shape" = indv_models_harm_list[[i]]$shape[mid_month],
                                        "scale" = indv_models_harm_list[[i]]$scale[mid_month],
                                        "kappa" = indv_models_harm_list[[i]]$kappa[mid_month])
  
  # write_csv(monthly_coefs, "outputs/monthly_coefs_wide_yday_with_movement_params_20230530.csv") 
  
  monthly_coefs_long_list[[i]] <- monthly_coefs_list[[i]] %>% pivot_longer(cols = !date)
  
  # sanity check to ensure that the coefficients line up to the correct months
  temp_plot <- monthly_coefs_long_list[[i]] %>% filter(!name %in% c("shape", "scale", "kappa")) %>% 
    ggplot(aes(x = month(date), y = value, colour = factor(name))) +
    geom_line() +
    scale_colour_viridis_d("Covariate") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_x_continuous(breaks = 1:12) +
    theme_classic()
  
  print(temp_plot)
  
}


# Predictions and cross-validation
## Naive prediction approach

boyce_naive_list <- vector(mode = "list", length = 12)

# outer loop for the exclusion of each individual
naive_pred_1harm_list <- vector(mode = "list", length = length(buffalo_ids))
naive_norm_1harm_list <- vector(mode = "list", length = length(buffalo_ids))
boyce_naive_list_list <- vector(mode = "list", length = length(buffalo_ids))

for(j in 1:length(buffalo_ids)) {
  
  # inner loop for each month
  naive_pred_stack <- c()
  naive_norm_stack <- c()
  
  for(i in 1:12) { 
    
    resources <- c(ndvi_stack_cropped[[which(time(ndvi_stack_cropped) == monthly_coefs_list[[j]]$date[i])]], 
                   canopy_scaled_cropped,
                   herby_scaled_cropped)
    
    # ndvi
    ndvi_lin <- resources[[1]]
    ndvi_lin_coef <- ndvi_lin * monthly_coefs_list[[j]]$ndvi[[i]]
    # plot(ndvi_lin)
    
    ndvi_quad <- resources[[1]]
    ndvi_quad_coef <- (ndvi_quad ^ 2) * monthly_coefs_list[[j]]$ndvi_quad[[i]] 
    # plot(ndvi_quad)
    
    ndvi_pred <- ndvi_lin_coef + ndvi_quad_coef
    # plot(ndvi_pred)
    
    # canopy cover 
    canopy_lin <- resources[[2]]
    canopy_lin_coef <- canopy_lin * monthly_coefs_list[[j]]$canopy[[i]]
    # plot(canopy_lin)
    
    canopy_quad <- resources[[2]]
    canopy_quad_coef <- (canopy_quad ^ 2) * monthly_coefs_list[[j]]$canopy_quad[[i]] 
    # plot(canopy_quad)
    
    canopy_pred <- canopy_lin_coef + canopy_quad_coef
    # plot(canopy_pred)
    
    # herby
    herby_lin <- resources[[3]]
    herby_pred <- herby_lin * monthly_coefs_list[[j]]$herby[[i]] 
    # plot(herby_pred)
    
    # combining
    naive_pred <- exp(ndvi_pred + canopy_pred + herby_pred)
    naive_norm <- naive_pred / global(naive_pred, fun = "sum", na.rm = TRUE)[[1]]
    
    plot(naive_norm)
    # points(as.numeric(buffalo_CLR_year_harmonics %>% filter(y == 1 & month == i & id == buffalo_ids[j]) %>% dplyr::select(x1) %>% unlist()), 
    #        as.numeric(buffalo_CLR_year_harmonics %>% filter(y == 1 & month == i & id == buffalo_ids[j]) %>% dplyr::select(y1) %>% unlist()))
    
    naive_pred_stack <- c(naive_pred_stack, naive_pred)
    naive_norm_stack <- c(naive_norm_stack, naive_norm)
    
    # setting up the buffalo observation data for the relevant month
    buffalo_obs <- buffalo_CLR_year_harmonics %>% filter(y == 1 & month == i & id == buffalo_ids[j]) %>%
      transmute(x = x1 - xmin, y = y1 - ymin)
    
    naive_raster <- raster(naive_norm)
    
    boyce_naive_list[[i]] <- ecospat.boyce2(naive_raster, buffalo_obs,
                                            method = "spearman"
                                            # method = "pearson",
                                            # method = "kendall"
    )
    
  }
  
  naive_pred_1harm_list[[j]] <- naive_pred_stack
  naive_norm_1harm_list[[j]] <- naive_norm_stack
  boyce_naive_list_list[[j]] <- boyce_naive_list
  
}

# for plotting the results

# for(j in 1:length(buffalo_ids)) {
#   for(i in 1:12) {
#       plot(naive_norm_1harm_list[[j]][[i]])
#   }
# }

# for(i in 1:12) {
#   # plot(naive_pred_stack[[i]])
#   plot(naive_pred_stack[[i]])
#   points(as.numeric(buffalo_CLR_year %>% filter(y == 1 & month == i) %>% dplyr::select(x1) %>% unlist()), 
#        as.numeric(buffalo_CLR_year %>% filter(y == 1 & month == i) %>% dplyr::select(y1) %>% unlist()))
# }


# assessing the spearman rank correlation coefficients

spearman_cor <- c()
spearman_cor_list <- vector(mode = "list", length = length(buffalo_ids))

for(j in 1:length(buffalo_ids)) {
  
  for(i in 1:12) {
    print(boyce_naive_list_list[[j]][[i]]$cor)
    spearman_cor[i] <- boyce_naive_list_list[[j]][[i]]$cor
  }
  
  spearman_cor_list[[j]] <- spearman_cor
  
}

for(j in 1:length(buffalo_ids)) {
  print(mean(spearman_cor_list[[j]]))
  print(sd(spearman_cor_list[[j]]))
}

# unlist this and find the overall mean - possibly also plot for each individual/month
# spearman_cor_list

# Using the Barnett-Moorcroft approximation prediction approach with Monte Carlo approximation

# Here we use the Barnett-Moorcroft approximation to estimate the UD, sampling from the movement kernel to approximate the integral in the numerator of:
  
 # BM equation

# inner loop
Beta_Z_z_list  <- vector(mode = "list", length = 12)
MCBM_preds_list <- vector(mode = "list", length = 12)
# outer loop
MCBM_preds_list_list <- vector(mode = "list", length = length(buffalo_ids))

# setting up the function
n_proposals <- 10
grid_res <- terra::res(naive_pred_stack[[1]])[1]
rows <- terra::nrow(naive_pred_stack[[1]])
cols <- terra::ncol(naive_pred_stack[[1]])
xmax <- ext(naive_pred_stack[[1]][[1]])[2] 
ymax <- ext(naive_pred_stack[[1]][[1]])[4]

for(j in 1:length(buffalo_ids)) {
  
  # create raster stack from list of naive predictions
  naive_preds <- rast(naive_pred_1harm_list[[j]])
  
  for(month_no in 1:12) {
    
    # for testing the loop
    # month_no <- 9
    
    tic()
    
    # Create grid of x and y points
    x_points <- rep(1:cols, each = rows)
    y_points <- rep(1:rows, times = cols)
    
    # Generate random angles and lengths
    ta <- runif(n_proposals * length(x_points), min = -pi, max = pi) # should be uniform
    sl <- rgamma(n_proposals * length(x_points), shape = monthly_coefs_list[[j]]$shape[month_no], scale = monthly_coefs_list[[j]]$scale[month_no])
    
    # Calculate proposal points
    x_proposal <- ((-(grid_res/2) + x_points * grid_res) + sl * sin(ta)) %% xmax
    y_proposal <- ((-(grid_res/2) + y_points * grid_res) + sl * cos(ta)) %% ymax
    
    # exp(beta * Z(x))
    # plot(naive_ud_cropped[[i]])
    
    # exp(beta * Z(z)) * psi(.) dz
    Beta_Z_z_proposed <- terra::extract(naive_preds[[month_no]], cbind(x_proposal, y_proposal))[,1]
    Beta_Z_z_array <- array(Beta_Z_z_proposed, dim = c(rows, cols , n_proposals))
    Beta_Z_z_matrix <- apply(Beta_Z_z_array, 1:2, mean, na.rm = TRUE)
    
    toc()
    
    Beta_Z_z_list[[month_no]] <- flip(terra::setValues(naive_preds[[month_no]], Beta_Z_z_matrix))
    # plot(naive_ud_cropped[[i]])
    # plot(Beta_Z_z_list[[month_no]])
    
    # exp(beta * Z(x)) * exp(beta * Z(z)) * psi(.) dz
    u_x_unnorm <- naive_preds[[month_no]] * Beta_Z_z_list[[month_no]]
    # plot(u_x_unnorm)
    u_x <- u_x_unnorm / as.numeric(terra::global(u_x_unnorm, fun = "sum", na.rm = TRUE))
    names(u_x) <- paste0("MCBM_month_", month_no)
    # plot(u_x)
    
    MCBM_preds_list[[month_no]] <- u_x
    
  }
  
  MCBM_preds_list_list[[j]] <- MCBM_preds_list
  
}

# save.image("MC-BM_approx_1harm_jacknife_20230609.RData")

# for(i in 1:12) plot(MCBM_preds_list_list[[1]][[i]])



# Plotting with points - naive and MC-BM

# for(month_no in 1:12) {
# 
# # exp(beta * Z(x))
# plot(naive_pred_stack[[month_no]])
# plot(naive_pred_stack[[month_no]])
# points(as.numeric(buffalo_CLR_year %>% 
#                     filter(y == 1 & month == month_no) %>% 
#                     dplyr::select(x1) %>% unlist()) - xmin, 
#        as.numeric(buffalo_CLR_year %>% 
#                     filter(y == 1 & month == month_no) %>% 
#                     dplyr::select(y1) %>% unlist()) - ymin)
# 
# # exp(beta * Z(z)) * psi(.) dz
# plot(Beta_Z_z_list[[month_no]])
# 
# # exp(beta * Z(x)) * exp(beta * Z(z)) * psi(.) dz
# plot(MCBM_preds_list[[month_no]])
# plot(MCBM_preds_list[[month_no]])
# points(as.numeric(buffalo_CLR_year %>%
#                     filter(y == 1 & month == month_no) %>%
#                     dplyr::select(x1) %>% unlist()) - xmin,
#        as.numeric(buffalo_CLR_year %>%
#                     filter(y == 1 & month == month_no) %>%
#                     dplyr::select(y1) %>% unlist()) - ymin)
# 
# }


# Boyce index 

# inner loop
boyce_naive_list <- vector(mode = "list", length = 12)
boyce_MCBM_list <- vector(mode = "list", length = 12)
# outer loop
boyce_naive_list_list <- vector(mode = "list", length = length(buffalo_ids))
boyce_MCBM_list_list <- vector(mode = "list", length = length(buffalo_ids))

for(j in 1:length(buffalo_ids)) {
  
  for(month_no in 1:12) {
    
    buffalo_obs <- buffalo_CLR_year_harmonics %>% 
      filter(y == 1 & month == month_no & id == buffalo_ids[[j]]) %>% 
      transmute(x = x1 - xmin, y = y1 - ymin)
    
    naive_raster <- raster(naive_norm_1harm_list[[j]][[month_no]])
    boyce_naive_list[[month_no]] <- ecospat.boyce2(naive_raster, buffalo_obs, method = "spearman")
    boyce_naive_list[[month_no]]
    
    MCBM_raster <- raster(MCBM_preds_list_list[[j]][[month_no]])
    boyce_MCBM_list[[month_no]] <- ecospat.boyce2(MCBM_raster, buffalo_obs, method = "spearman")
    boyce_MCBM_list[[month_no]]
    
  }
  
  boyce_naive_list_list[[j]] <- boyce_naive_list
  boyce_MCBM_list_list[[j]] <- boyce_MCBM_list
  
}

# extracting just the spearman rank correlation coefficient from each set of predictions

boyce_naive_cors <- c()
boyce_MCBM_cors <- c()

for(j in 1:length(buffalo_ids)) {
  for(month_no in 1:12) {
    boyce_naive_cors <- c(boyce_naive_cors, boyce_naive_list_list[[j]][[month_no]]$cor)
    boyce_MCBM_cors <- c(boyce_MCBM_cors, boyce_MCBM_list_list[[j]][[month_no]]$cor)
  }
}

boyce_correlation_df <- data.frame("id" = rep(buffalo_ids, each = 12),
                                   "month" = rep(1:12, length(buffalo_ids)),
                                   "Naive" = boyce_naive_cors,
                                   "MCBM" = boyce_MCBM_cors)

boyce_correlation_df_long <- boyce_correlation_df %>% pivot_longer(cols = !1:2)

ggplot(boyce_correlation_df_long) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  # hide outliers from boxplot as already included when jittering # , outlier.shape = NA
  geom_boxplot(aes(x = factor(month), y = value, colour = name, fill = name), alpha = 0.5) + 
  # geom_jitter(aes(x = factor(month), y = value, colour = name), 
  #             width = 0.1, alpha = 0.5, show.legend = FALSE) +
  # geom_path(aes(x = month, y = value, colour = factor(id), group = interaction(name, id))) +
  scale_y_continuous("Spearman rank correlation coefficient", limits = c(-1, 1)) +
  scale_x_discrete("Month") +
  scale_colour_viridis_d("Prediction approach", option = "D") +
  scale_fill_viridis_d("Prediction approach", option = "D") +
  theme_bw() +
  theme(legend.position = c(0.15, 0.15))

ggsave("outputs/plots/MCBM_naive_Boyce_jackknife_1harm_20230612.png",
       width=150, height=90, units="mm", dpi = 300)