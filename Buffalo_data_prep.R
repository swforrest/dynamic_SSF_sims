# To load the buffalo GPS data used by many analysis scripts

library(tidyverse)
library(lubridate)

# buffalo_data_popn <- read_csv(file = "outputs/buffalo_parametric_popn_covs_20230208.csv")
# buffalo_data <- read_csv(file = "outputs/buffalo_parametric_indv_covs_20230208.csv")
buffalo_data <- read_csv(file = "outputs/buffalo_parametric_indv_hourly_stepgen_covs_2023-09-26.csv")

# buffalo_data %>% group_by(id) %>% filter(y == 1) %>% summarise(mean_sl = mean(sl_))
# buffalo_data %>% group_by(id) %>% filter(y == 0) %>% summarise(mean_sl = mean(sl_))

buffalo_data <- buffalo_data %>% mutate(t1_ = lubridate::with_tz(buffalo_data$t1_, tzone = "Australia/Darwin"),
                                        t2_ = lubridate::with_tz(buffalo_data$t2_, tzone = "Australia/Darwin"))

buffalo_CLR <- buffalo_data %>% 
  mutate(id_num = as.numeric(factor(id)), 
         step_id = step_id_, 
         x1 = x1_, x2 = x2_, 
         y1 = y1_, y2 = y2_, 
         t1 = t1_, 
         t1_rounded = round_date(buffalo_data$t1_, "hour"), 
         hour_t1 = hour(t1_rounded),
         t2 = t2_, 
         t2_rounded = round_date(buffalo_data$t2_, "hour"), 
         hour_t2 = hour(t2_rounded),
         hour = hour_t2,
         yday = yday(t1_),
         year = year(t1_), 
         month = month(t1_),
         sl = sl_, 
         # sl_scaled = scale(sl),
         log_sl = log(sl_), 
         ta = ta_, 
         cos_ta = cos(ta_),
         # ndvi_scaled = scale(ndvi_temporal),
         # ndvi_2 = ndvi_scaled[,1]^2,
         # canopy_01 = canopy_cover/100,
         # herby_scaled = scale(veg_herby),
         # canopy_scaled = scale(canopy_01),
         # canopy_2 = canopy_scaled[,1]^2,
         # elev_scaled = scale(DEM_H_end),
         # elev_2 = elev_scaled[,1]^2,
         # elev_delta_scaled = scale(elev_delta),
         # elev_0 = ifelse(DEM_H_end < 1, 1, DEM_H_end),
         # elev_log = log(elev_0),
         # slope_scaled = scale(slope_end),
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
         hour_c4 = cos(8*pi*hour/24)) #%>%
  # drop_na(c(cos_ta, sl)) %>% 
  # dplyr::select(!(x1_:case_))


# Filtering by animals that have more than a year of data

# buffalo_CLR %>% ggplot(aes(x = t1, y = factor(id), colour = factor(id))) +
#   geom_point(alpha = 0.1) +
#   scale_y_discrete("Buffalo ID") +
#   scale_x_datetime("Date") +
#   scale_colour_viridis_d() +
#   theme_bw() +
#   theme(legend.position = "none")
# 
# buffalo_CLR %>% dplyr::group_by(id) %>%  
#   summarise(min_time = min(t1), max_time = max(t1),
#             min_x = min(x2), max_x = max(x2),
#             min_y = min(y2), max_y = max(y2))

buffalo_year_ids <- c(2005, 2014, 2018, 2022, 2024, 2154, 2158, 2327, 2354, 2387)
buffalo_CLR_year <- buffalo_CLR %>% filter(id %in% buffalo_year_ids)
unique(buffalo_CLR_year$id)

min(buffalo_CLR_year$t1)
max(buffalo_CLR_year$t1)

buffalo_CLR_year_harmonics <- buffalo_CLR_year %>% filter(t1 < "2019-07-25 09:32:42 ACST")

# write_csv(buffalo_CLR_year_harmonics, paste0("outputs/buffalo_parametric_popn_covs_harmonics_", Sys.Date(), ".csv"))
write_csv(buffalo_CLR_year_harmonics, paste0("outputs/buffalo_parametric_indv_hourly_covs_harmonics_", Sys.Date(), ".csv"))

# buffalo_CLR_year_harmonics %>% ggplot(aes(x = t1, y = factor(id), colour = factor(id))) +
#   geom_point(alpha = 0.1) +
#   scale_y_discrete("Buffalo ID") +
#   scale_x_datetime("Date") +
#   scale_colour_viridis_d() +
#   theme_bw() +
#   theme(legend.position = "none")

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

# Movement parameters from random sampling

# population
# gamma$params$shape # 0.438167
# gamma$params$scale # 534.3507
# vonmises$params$kappa # 0.1848126
# vonmises$params$mu # 0

# movement_parameters <- read_csv("outputs/buffalo_parametric_indv_movement_params_20230208.csv")
# movement_parameters_year <- movement_parameters %>% filter(id %in% buffalo_year_ids)
