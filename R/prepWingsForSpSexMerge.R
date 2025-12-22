
# perform group, summarize, and pivot to yield a df with columns like
# "species", "sex", "fore_blackbrown", "fore_winglength"... etc.
prepWingsForSpSexMerge <- function(wings){
  library(dplyr)
  library(tidyr)
  df <- wings
  
  # Group by 'species', 'sex', 'wing_type' and calculate the means
  df_summary <- df %>%
    group_by(species, sex, wing_type) %>%
    summarize(
      sp_sex_blackbrown = mean(blackbrown, na.rm = TRUE),
      sp_sex_brownyellow = mean(brownyellow, na.rm = TRUE),
      sp_sex_winglength = mean(wing_length, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Calculate the mean values for each species, ignoring sex
  df_unknown <- df %>%
    group_by(species, wing_type) %>%
    summarize(
      sp_sex_blackbrown = mean(blackbrown, na.rm = TRUE),
      sp_sex_brownyellow = mean(brownyellow, na.rm = TRUE),
      sp_sex_winglength = mean(wing_length, na.rm = TRUE),
      sex = "unknown",
      .groups = 'drop'
    )
  
  # Combine the summary data with the unknown sex data
  df_combined <- bind_rows(df_summary, df_unknown)
  
  # Pivot the data to wide format
  df_wide <- df_combined %>%
    pivot_wider(
      names_from = wing_type,
      values_from = c(sp_sex_blackbrown, sp_sex_brownyellow),
      names_glue = "{.value}_{wing_type}"
    )
  
  # Remove any columns with "NA" in the column name
  df_wide <- df_wide %>% select(-contains("NA"))
  
  return(df_wide)
}
#wings_sp_sex <- prepWingsForSpSexMerge(wings)
