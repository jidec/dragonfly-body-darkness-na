# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(dplyr)
library(lubridate)
library(colorEvoHelpers)
library(moments)  # Add for skewness/kurtosis
library(ineq)     # Add for Gini coefficient
library(stringr)

# library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c("tibble", "moments", "ineq") # packages that your targets need to run
  # format = "qs", # Optionally set the default storage format. qs is fast.
  #
  # For distributed computing in tar_make(), supply a {crew} controller
  # as discussed at https://books.ropensci.org/targets/crew.html.
  # Choose a controller that suits your needs. For example, the following
  # sets a controller with 2 workers which will run as local R processes:
  #
  #   controller = crew::crew_controller_local(workers = 2)
  #
  # Alternatively, if you want workers to run on a high-performance computing
  # cluster, select a controller from the {crew.cluster} package. The following
  # example is a controller for Sun Grid Engine (SGE).
  # 
  #   controller = crew.cluster::crew_controller_sge(
  #     workers = 50,
  #     # Many clusters install R as an environment module, and you can load it
  #     # with the script_lines argument. To select a specific verison of R,
  #     # you may need to include a version string, e.g. "module load R/4.3.0".
  #     # Check with your system administrator if you are unsure.
  #     script_lines = "module load R"
  #   )
  #
  # Set other options as needed.
)

# tar_make_clustermq() is an older (pre-{crew}) way to do distributed computing
# in {targets}, and its configuration for your machine is below.
options(clustermq.scheduler = "multiprocess")

# tar_make_future() is an older (pre-{crew}) way to do distributed computing
# in {targets}, and its configuration for your machine is below.
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

# Run the R scripts in the R/ folder with your custom functions:
tar_source()
# source("other_functions.R") # Source other scripts as needed.

# note that flight_styles.csv is created using writePercherFlyerConsensus in /R beforehand
# wings.csv should not have rownames col
# Replace the target list below with your own:
list(
  # files
  tar_target(file1, "data/records_with_metrics.csv", format="file"),
  # filter for North America only
  tar_target(d_na, colorEvoHelpers::filterNorthAmerica(read.csv(file1))),
  # bind sex column using iNat annotations
  tar_target(d_sex, colorEvoHelpers::bindiNatSexAnnotation(d_na)),
  # parse iNat formatted datetimes (may be able to remove this later when datetime is consistent)
  tar_target(d_datetime, d_sex %>% mutate(datetime=lubridate::parse_date_time(date,orders = c(
    "Y-m-d H:M:S z", "Y-m-d H:M:S", "Y-m-d", "Y/m/d H:M p", "Y/m/d H:M:S p", "Y/m/d H:M:S",
    "a b d Y H:M:S z", "Y/m/d H:M p z", "Y/m/d H:M:S p z", "Y/m/d H:M p", "Y/m/d H:M:S p",
    "Y/m/d H:M:S", "a b d Y H:M:S", "Y/m/d H:M p", "Y/m/d H:M", "Y/m/d", "a b d Y H:M:S z"
  ))
  )),
  # bind yday, month, year, and season
  tar_target(d_temporal, bindDayMonthYearSeason(d_datetime)),
  # merge wings data into original (left join)
  tar_target(d_gridded, bindGridCells(d_temporal)),
  # bind prism data 
  #tar_target(d_prism, bindPrismData(d_gridded, prism_dir = "D:/GitProjects/inat-daily-activity-analysis/data/prism", variable="tmean",new_var_colname="prism_tmean")),
  
  tar_target(d_sp_cell_no_clim, d_gridded %>%
               dplyr::select(-geometry) %>%
               distinct(img_id, .keep_all = TRUE) %>%
               filter(yday > 0) %>% # remove erroneous ydays
               group_by(species, cell) %>%
               filter(n() >= 10) %>% # each species-cell must have at least 10 observations
               group_modify(~ {
                 # cap at 50 
                 if (nrow(.x) > 50) { # cap at 50 observations per species-cell - intended to reduce the impact of bias associated with very heavy sampling
                   .x <- slice_sample(.x, n = 50)
                 }
                 
                 # === ORIGINAL METRICS ===
                 lightness_vals <- .x$mean_cielab_lightness # lightness value is the mean cielab lightness
                 lightness_mean <- mean(lightness_vals)
                 lightness_median <- median(lightness_vals)
                 lightness_sd <- sd(lightness_vals)
                 
                 q <- quantile(lightness_vals, probs = c(0.05, 0.1, 0.25, 0.75, 0.9, 0.95))
                 
                 # === SENSITIVITY: IQR trim + use median ===
                 iqr <- IQR(lightness_vals)
                 q1 <- quantile(lightness_vals, 0.2) # 0.2 min
                 q3 <- quantile(lightness_vals, 0.8) # 0.8 max
                 
                 lower_fence <- q1 - 1.5 * iqr
                 upper_fence <- q3 + 1.5 * iqr
                 .x_trimmed <- .x %>% filter(mean_cielab_lightness >= lower_fence & 
                                               mean_cielab_lightness <= upper_fence)
                 
                 if (nrow(.x_trimmed) >= 5) {
                   black_20_sens <- median(.x_trimmed$hls_ch1_below0.2)
                   black_30_sens <- median(.x_trimmed$hls_ch1_below0.3)
                   mean_lightness_sens <- median(.x_trimmed$mean_cielab_lightness)
                   lightness_skew_sens <- moments::skewness(.x_trimmed$mean_cielab_lightness)
                   n_sens <- nrow(.x_trimmed)
                 } else {
                   black_20_sens <- NA_real_
                   black_30_sens <- NA_real_
                   mean_lightness_sens <- NA_real_
                   lightness_skew_sens <- NA_real_
                   n_sens <- NA_integer_
                 }
                 
                 tibble(
                   # === USED IN ANALYSIS ===
                   lightness_values = list(lightness_vals),
                   ids = list(.x$img_id),
                   urls = list(.x$img_url),
                   onset = quantile(.x$yday, 0.1),
                   cell_lat = Mode(.x$cell_lat),
                   cell_lon = Mode(.x$cell_lon),
                   year = Mode(.x$year),
                   black_20 = mean(.x$hls_ch1_below0.2),
                   black_30 = mean(.x$hls_ch1_below0.3),
                   mean_lightness = lightness_mean,
                   latitude = Mode(.x$cell_lat),
                   longitude = Mode(.x$cell_lon),
                   n = nrow(.x),
                   lightness_skew = moments::skewness(lightness_vals),
                   lightness_cv = lightness_sd / lightness_mean,
                   
                   # === SENSITIVITY METRICS ===
                   black_20_sens = black_20_sens,
                   black_30_sens = black_30_sens,
                   mean_lightness_sens = mean_lightness_sens,
                   lightness_skew_sens = lightness_skew_sens,
                   n_sens = n_sens,
                   n_removed = nrow(.x) - ifelse(is.na(n_sens), 0L, n_sens)
                   
                   # === COMMENTED OUT (bunch of other possible metrics not currently used) ===
                   # lightness_median = lightness_median,
                   # lightness_trimmed = mean(lightness_vals, trim = 0.1),
                   # lightness_var = var(lightness_vals),
                   # lightness_sd = lightness_sd,
                   # lightness_iqr = iqr,
                   # lightness_mad = mad(lightness_vals),
                   # brighter_var = if(length(brighter_half) > 1) var(brighter_half) else NA_real_,
                   # darker_var = if(length(darker_half) > 1) var(darker_half) else NA_real_,
                   # var_ratio_light = if(length(brighter_half) > 1 & length(darker_half) > 1) 
                   #   var(brighter_half) / var(darker_half) else NA_real_,
                   # lightness_kurtosis = moments::kurtosis(lightness_vals) - 3,
                   # tail_ratio_light = (q[["95%"]] - lightness_median) / (lightness_median - q[["5%"]]),
                   # upper_tail_weight = mean(lightness_vals > q[["75%"]]),
                   # lower_tail_weight = mean(lightness_vals < q[["25%"]]),
                   # prop_bright_outliers = mean(bright_outliers),
                   # prop_dark_outliers = mean(dark_outliers),
                   # outlier_asymmetry_light = sum(bright_outliers) - sum(dark_outliers),
                   # p90_p10_range = q[["90%"]] - q[["10%"]],
                   # p95_p05_range = q[["95%"]] - q[["5%"]],
                   # upper_decile_spread = q[["90%"]] - lightness_median,
                   # lower_decile_spread = lightness_median - q[["10%"]],
                   # decile_asymmetry = (q[["90%"]] - lightness_median) - (lightness_median - q[["10%"]]),
                   # bright_tail_excess = mean(lightness_vals > lightness_mean + lightness_sd),
                   # dark_tail_excess = mean(lightness_vals < lightness_mean - lightness_sd),
                   # tail_excess_ratio_light = if(sum(lightness_vals < lightness_mean - lightness_sd) > 0) 
                   #   mean(lightness_vals > lightness_mean + lightness_sd) / 
                   #   mean(lightness_vals < lightness_mean - lightness_sd) else NA_real_,
                   # lightness_gini = ineq::Gini(lightness_vals)
                 )
               }) %>%
               ungroup() %>%
               { 
                 m <- lm(onset ~ n, data = .)
                 mu <- mean(.$onset, na.rm = TRUE)
                 mutate(., 
                        bias = predict(m, newdata = .),
                        onset = onset - bias + mu) %>%
                   select(-bias)
               }),
  
  # Cell-level illumination metrics (streamlined)
  tar_target(cell_illumination_metrics, d_sp_cell_no_clim %>%
               group_by(cell) %>%
               summarise(
                 cell_skew_weighted = weighted.mean(lightness_skew, w = n, na.rm = TRUE),
                 cell_cv_weighted = weighted.mean(lightness_cv, w = n, na.rm = TRUE),
                 n_species = n(),
                 .groups = "drop"
               ) %>%
               mutate(
                 z_skew = scale(cell_skew_weighted)[,1],
                 z_cv = scale(cell_cv_weighted)[,1],
                 illumination_index = z_skew,  # simplified since var_ratio and tail_ratio removed
                 illumination_variance_index = (abs(z_skew) + z_cv) / 2
               )),
  
  # Merge illumination metrics back with climate data
  tar_target(d_sp_cell_with_illum_no_clim, d_sp_cell_no_clim %>%
               left_join(cell_illumination_metrics %>% 
                        select(cell, starts_with("cell_"), starts_with("illumination_"), 
                               starts_with("z_"), starts_with("mean_prop")), 
                        by = "cell")),
  
  # bind bioclim annual avg tmax and prec (now with illumination metrics)
  tar_target(d_sp_cell_clim, colorEvoHelpers::bindWorldClimData(d_sp_cell_with_illum_no_clim, var = "tavg", new_var_colname = "tmean", cellsize_km = 250) %>%
               colorEvoHelpers::bindWorldClimData(var = "srad", new_var_colname = "srad", cellsize_km = 250) %>% 
               colorEvoHelpers::bindWorldClimData(var = "vapr", new_var_colname = "vp", cellsize_km = 250) %>%
               colorEvoHelpers::bindWorldClimData(var = "bio", bioclim_var_num = 4, new_var_colname = "temp_seasonality") %>% 
               colorEvoHelpers::scaleDf(all_numerics=TRUE)),
  
  tar_target(d_sp_cell, d_sp_cell_clim %>% filter(!str_detect(species, "^(Dorocordulia|Cordulia|Somatochlora)"))),
  
  # read tree
  tar_target(tree, colorEvoHelpers::loadFixTree(tree_location="D:/GitProjects/dragonfly-body-darkness-na/data/dragonfly_tree.tre")),
  # trim data to tree and tree to data
  tar_target(d, colorEvoHelpers::trimDfToTree(d_sp_cell,tree)[[1]]),
  tar_target(tree_trimmed,colorEvoHelpers::trimDfToTree(d_sp_cell,tree)[[2]])
)