# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(dplyr)
library(lubridate)
library(colorEvoHelpers)

# library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c("tibble") # packages that your targets need to run
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
               filter(yday>0) %>%
               group_by(species, cell) %>%
               filter(n() >= 10) %>% # min of 10 per sp-cell
               group_modify(~ {
                 # cap at 50 
                 if (nrow(.x) > 50) {
                   .x <- slice_sample(.x, n = 50)
                 }
                 
                 tibble(
                   #onset = phenesse::quantile_ci(.x$yday, 0.1),
                   onset = quantile(.x$yday, 0.1),
                   cell_lat = Mode(.x$cell_lat),
                   cell_lon = Mode(.x$cell_lon),
                   year = Mode(.x$year),
                   black_20 = mean(.x$hls_ch1_below0.2),
                   black_30 = mean(.x$hls_ch1_below0.3),
                   mean_lightness = mean(.x$mean_cielab_lightness),
                   #wing_length = mean(.x$sp_winglength),
                   #flight_type = Mode(.x$flight_type),
                   latitude = cell_lat,
                   longitude = cell_lon,
                   n = nrow(.x)
                 )
               }) %>%
               ungroup() %>%
               { 
                 # fit across all groups & correct
                 m <- lm(onset ~ n, data = .)
                 mu <- mean(.$onset, na.rm = TRUE)
                 mutate(., 
                        bias = predict(m, newdata = .),
                        onset = onset - bias + mu) %>%
                   select(-bias)
               }),
  # bind bioclim annual avg tmax and prec
  tar_target(d_sp_cell,colorEvoHelpers::bindWorldClimData(d_sp_cell_no_clim, var = "tavg", new_var_colname = "tmean", cellsize_km = 250) %>%
               colorEvoHelpers::bindWorldClimData(var = "srad", new_var_colname = "srad", cellsize_km = 250) %>% 
               colorEvoHelpers::bindWorldClimData(var = "vapr", new_var_colname = "vp", cellsize_km = 250) %>%
               colorEvoHelpers::bindWorldClimData(var = "bio", bioclim_var_num = 4, new_var_colname = "temp_seasonality") %>% 
               colorEvoHelpers::scaleDf(all_numerics=TRUE)),
  # read tree
  tar_target(tree, colorEvoHelpers::loadFixTree(tree_location="D:/GitProjects/dragonfly-body-darkness-na/data/dragonfly_tree.tre")),
  # trim data to tree and tree to data
  tar_target(d, colorEvoHelpers::trimDfToTree(d_sp_cell,tree)[[1]]),
  tar_target(tree_trimmed,colorEvoHelpers::trimDfToTree(d_sp_cell,tree)[[2]])
)

