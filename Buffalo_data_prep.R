# To load the buffalo GPS data used by many analysis scripts

library(tidyverse)
library(lubridate)

buffalo_data <- read_csv(file = "outputs/buffalo_parametric_popn_covs_20230208.csv")
buffalo_data <- buffalo_data %>% mutate(t1_ = lubridate::with_tz(buffalo_data$t1_, tzone = "Australia/Darwin"),
                                        t2_ = lubridate::with_tz(buffalo_data$t2_, tzone = "Australia/Darwin"))

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
         canopy_01 = canopy_cover/100,
         herby_scaled = scale(veg_herby),
         canopy_scaled = scale(canopy_01),
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


# Aligning step ids

# rounded_hours <- unique(round_date(buffalo_CLR_year$t1, "hour"))
# rounded_hours <- rounded_hours[order(rounded_hours)]
# unique_hours <- 1:length(unique(round_date(buffalo_CLR_year$t1, "hour")))
# hourly_steps <- data.frame(unique_hours, rounded_hours)
# 
# buffalo_CLR_year$step_aligned <- 0
# for(i in 1:length(unique_hours)) {
#   buffalo_CLR_year$step_aligned[which(buffalo_CLR_year$t1_rounded == hourly_steps$rounded_hours[i])] <- hourly_steps$unique_hours[i]
# }


# Add harmonic terms multiplied by covariates (essentially creating a design matrix)

buffalo_CLR_year_harmonics <- buffalo_CLR_year %>% 
  mutate(
    ndvi_2 = ndvi_scaled[,1]^2,
    canopy_2 = canopy_scaled[,1]^2
    
    # ndvi_s1 = ndvi_scaled[,1] * yday_s1,
    # ndvi_s2 = ndvi_scaled[,1] * yday_s2,
    # ndvi_s3 = ndvi_scaled[,1] * yday_s3,
    # ndvi_s4 = ndvi_scaled[,1] * yday_s4,
    # ndvi_c1 = ndvi_scaled[,1] * yday_c1,
    # ndvi_c2 = ndvi_scaled[,1] * yday_c2,
    # ndvi_c3 = ndvi_scaled[,1] * yday_c3,
    # ndvi_c4 = ndvi_scaled[,1] * yday_c4,
    
    
    # ndvi_2_s1 = ndvi_2 * yday_s1,
    # ndvi_2_s2 = ndvi_2 * yday_s2,
    # ndvi_2_s3 = ndvi_2 * yday_s3,
    # ndvi_2_s4 = ndvi_2 * yday_s4,
    # ndvi_2_c1 = ndvi_2 * yday_c1,
    # ndvi_2_c2 = ndvi_2 * yday_c2,
    # ndvi_2_c3 = ndvi_2 * yday_c3,
    # ndvi_2_c4 = ndvi_2 * yday_c4,
    
    # canopy_s1 = canopy_scaled[,1] * yday_s1,
    # canopy_s2 = canopy_scaled[,1] * yday_s2,
    # canopy_s3 = canopy_scaled[,1] * yday_s3,
    # canopy_s4 = canopy_scaled[,1] * yday_s4,
    # canopy_c1 = canopy_scaled[,1] * yday_c1,
    # canopy_c2 = canopy_scaled[,1] * yday_c2,
    # canopy_c3 = canopy_scaled[,1] * yday_c3,
    # canopy_c4 = canopy_scaled[,1] * yday_c4,
    
    # canopy_2_s1 = canopy_2 * yday_s1,
    # canopy_2_s2 = canopy_2 * yday_s2,
    # canopy_2_s3 = canopy_2 * yday_s3,
    # canopy_2_s4 = canopy_2 * yday_s4,
    # canopy_2_c1 = canopy_2 * yday_c1,
    # canopy_2_c2 = canopy_2 * yday_c2,
    # canopy_2_c3 = canopy_2 * yday_c3,
    # canopy_2_c4 = canopy_2 * yday_c4,
    # 
    # sl_s1 = sl * yday_s1,
    # sl_s2 = sl * yday_s2,
    # sl_s3 = sl * yday_s3,
    # sl_s4 = sl * yday_s4,
    # sl_c1 = sl * yday_c1,
    # sl_c2 = sl * yday_c2,
    # sl_c3 = sl * yday_c3,
    # sl_c4 = sl * yday_c4,
    # 
    # sl_scaled_s1 = sl_scaled * yday_s1,
    # sl_scaled_s2 = sl_scaled * yday_s2,
    # sl_scaled_s3 = sl_scaled * yday_s3,
    # sl_scaled_s4 = sl_scaled * yday_s4,
    # sl_scaled_c1 = sl_scaled * yday_c1,
    # sl_scaled_c2 = sl_scaled * yday_c2,
    # sl_scaled_c3 = sl_scaled * yday_c3,
    # sl_scaled_c4 = sl_scaled * yday_c4,
    # 
    # log_sl_s1 = log_sl * yday_s1,
    # log_sl_s2 = log_sl * yday_s2,
    # log_sl_s3 = log_sl * yday_s3,
    # log_sl_s4 = log_sl * yday_s4,
    # log_sl_c1 = log_sl * yday_c1,
    # log_sl_c2 = log_sl * yday_c2,
    # log_sl_c3 = log_sl * yday_c3,
    # log_sl_c4 = log_sl * yday_c4,
    # 
    # cos_ta_s1 = cos_ta * yday_s1,
    # cos_ta_s2 = cos_ta * yday_s2,
    # cos_ta_s3 = cos_ta * yday_s3,
    # cos_ta_s4 = cos_ta * yday_s4,
    # cos_ta_c1 = cos_ta * yday_c1,
    # cos_ta_c2 = cos_ta * yday_c2,
    # cos_ta_c3 = cos_ta * yday_c3,
    # cos_ta_c4 = cos_ta * yday_c4
    
  )

# Movement parameters from random sampling

# population
# gamma$params$shape # 0.438167
# gamma$params$scale # 534.3507
# vonmises$params$kappa # 0.1848126
# vonmises$params$mu # 0

movement_parameters <- read_csv("outputs/buffalo_parametric_indv_movement_params_20230208.csv")
movement_parameters_year <- movement_parameters %>% filter(id %in% buffalo_year_ids)