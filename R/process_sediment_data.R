# Extract and merge result tables
process_sediment_data <- function(dat) {

  # Execute the main function
  grad <- G2Sd::granstat(dat)

  # Separate textural descriptions
  grad$stat$fowa <- grad$stat$fowa |>
    tidyr::separate(
      col = sediment,
      into = c(
        "FWM_mean_description",
        "FWM_sorting_description",
        "FWM_skewness_description",
        "FWM_kurtosis_description"
      ),
      sep = "/"
    )

  # Group result datasets
  list_stats <- list(
    grad$stat$arith,
    grad$stat$geom,
    grad$stat$fowa,
    grad$index,
    grad$sedim$texture,
    grad$sedim$descript
  )

  # Full merge of tables by sample ID
  dat_final <- Reduce(function(x, y) {
    merge(x, y, by = "samples", all = TRUE)
  }, list_stats)

  # Organize statistical columns
  ordre_stats <- c(
    "arithmetic_mean"          = "mean.arith.um",
    "arithmetic_sorting"       = "sd.arith.um",
    "arithmetic_skewness"      = "skewness.arith.um",
    "arithmetic_kurtosis"      = "kurtosis.arith.um",

    "geometric_mean"           = "mean.geom.um",
    "geometric_sorting"        = "sd.geom.um",
    "geometric_skewness"       = "skewness.geom.um",
    "geometric_kurtosis"       = "kurtosis.geom.um",

    "logarithmic_mean"         = "mean.log.phi",
    "logarithmic_sorting"      = "sd.log.phi",
    "logarithmic_skewness"     = "skewness.log.phi",
    "logarithmic_kurtosis"     = "kurtosis.log.phi",

    "FWM_geometric_mean"       = "mean.fw.um",
    "FWM_geometric_sorting"    = "sd.fw.um",
    "FWM_geometric_skewness"   = "skewness.fw.um",
    "FWM_geometric_kurtosis"   = "kurtosis.fw.um",

    "FWM_logarithmic_mean"     = "mean.fw.phi",
    "FWM_logarithmic_sorting"  = "sd.fw.phi",
    "FWM_logarithmic_skewness" = "skewness.fw.phi",
    "FWM_logarithmic_kurtosis" = "kurtosis.fw.phi",

    "FWM_mean_description"     = "FWM_mean_description",
    "FWM_sorting_description"  = "FWM_sorting_description",
    "FWM_skewness_description" = "FWM_skewness_description",
    "FWM_kurtosis_description" = "FWM_kurtosis_description"
  )

  # Convert percentiles to the logarithmic Phi scale
  dat_final[,"D10.phi"] <- G2Sd:::.um2phi(dat_final$D10)
  dat_final[,"D25.phi"] <- G2Sd:::.um2phi(dat_final$D25)
  dat_final[,"D50.phi"] <- G2Sd:::.um2phi(dat_final$D50)
  dat_final[,"D75.phi"] <- G2Sd:::.um2phi(dat_final$D75)
  dat_final[,"D90.phi"] <- G2Sd:::.um2phi(dat_final$D90)
  dat_final[,"D90/D10.phi"] <- dat_final$D90.phi / dat_final$D10.phi
  dat_final[,"D90-D10.phi"] <- dat_final$D90.phi - dat_final$D10.phi
  dat_final[,"D75/D25.phi"] <- dat_final$D75.phi / dat_final$D25.phi
  dat_final[,"D75-D25.phi"] <- dat_final$D75.phi - dat_final$D25.phi

  # Organize percentiles columns
  ordre_percentile <- c(
    "percentile_10"                       = "D10", 
    "percentile_25"                       = "D25", 
    "percentile_50"                       = "D50", 
    "percentile_75"                       = "D75", 
    "percentile_90"                       = "D90", 
    "percentile_ratio_of_90_to_10"        = "D90/D10", 
    "interquartile_range_of_90_to_10"     = "D90-D10", 
    "percentile_ratio_of_75_to_25"        = "D75/D25", 
    "interquartile_range_of_75_to_25"     = "D75-D25",
    "percentile_10_phi"                   = "D10.phi", 
    "percentile_25_phi"                   = "D25.phi", 
    "percentile_50_phi"                   = "D50.phi", 
    "percentile_75_phi"                   = "D75.phi", 
    "percentile_90_phi"                   = "D90.phi", 
    "percentile_ratio_of_90_to_10_phi"    = "D90/D10.phi", 
    "interquartile_range_of_90_to_10_phi" = "D90-D10.phi", 
    "percentile_ratio_of_75_to_25_phi"    = "D75/D25.phi", 
    "interquartile_range_of_75_to_25_phi" = "D75-D25.phi"
  )

  # Organize texture columns
  ordre_texture <- c(
    "gravel" = "Gravel",
    "sand" = "Sand",
    "silt" = "silt",
    "clay" = "clay"
  )

  # Organize size classes columns
  ordre_classes <- c(
    "very_coarse_gravel" = "vcgravel", 
    "coarse_gravel"      = "cgravel", 
    "medium_gravel"      = "mgravel", 
    "fine_gravel"        = "fgravel", 
    "very_fine_gravel"   = "vfgravel", 
    "very_coarse_sand"   = "vcsand", 
    "coarse_sand"        = "csand", 
    "medium_sand"        = "msand", 
    "fine_sand"          = "fsand", 
    "very_fine_sand"     = "vfsand", 
    "very_coarse_silt"   = "vcsilt", 
    "coarse_silt"        = "csilt", 
    "medium_silt"        = "msilt", 
    "fine_silt"          = "fsilt", 
    "very_fine_silt"     = "vfsilt"
  )

  # Organize every columns
  ordre_cols <- c(
    "station_id" = "samples",
    "sediments_texture" = "texture",
    ordre_stats,
    ordre_percentile,
    ordre_texture,
    ordre_classes
  )

  # Final reorganization and summation of the global Silt class
  dat_final <- dat_final |>
    dplyr::mutate(
      silt = vcsilt + csilt + msilt + fsilt + vfsilt
    ) |>
    dplyr::select(tidyselect::any_of(ordre_cols))

  return(dat_final)
}
