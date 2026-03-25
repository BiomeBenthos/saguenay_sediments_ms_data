# Calculate statistical indices using the arithmetic Method of Moments
.moment.arith <- function(x) {
  x <- G2Sd:::.g2sd_tidy(x)

  all <- x |>
    dplyr::group_by(samples) |>
    dplyr::mutate(relative.value = value * 100 / sum(value)) |>
    dplyr::arrange(meshsize) |>
    # Calculate arithmetic midpoints between sieves for accurate representative particle sizes
    dplyr::mutate(meshsize.midpoint = (dplyr::lag(meshsize) + meshsize) / 2 / 1000)

  # Use the raw mesh size for the finest fraction where no lower boundary exists
  all[is.na(all$meshsize.midpoint), "meshsize.midpoint"] <- all[
    is.na(all$meshsize.midpoint),
    "meshsize"
  ] / 1000

  # Compute mean, sorting (SD), skewness, and kurtosis in metric units
  arith <- all |>
    dplyr::group_by(samples) |>
    dplyr::mutate(mean.arith.um = sum(relative.value * meshsize.midpoint) / 100) |>
    dplyr::mutate(
      sd.arith.um = sqrt(
        sum(relative.value * (meshsize.midpoint - mean.arith.um)^2) / 100
      )
    ) |>
    dplyr::mutate(
      skewness.arith.um = sum(
        relative.value * (meshsize.midpoint - mean.arith.um)^3
      ) / (100 * sd.arith.um^3)
    ) |>
    dplyr::mutate(
      kurtosis.arith.um = sum(
        relative.value * (meshsize.midpoint - mean.arith.um)^4
      ) / (100 * sd.arith.um^4)
    ) |>
    dplyr::mutate(
      mean.arith.um = mean.arith.um * 1000,
      sd.arith.um = sd.arith.um * 1000
    ) |>
    dplyr::select(-c("meshsize", "value", "relative.value", "meshsize.midpoint")) |>
    dplyr::distinct()

  return(arith)
}

# Calculate statistical indices using the geometric Method of Moments
.moment.geom <- function(x) {
  x <- G2Sd:::.g2sd_tidy(x)

  all <- x |>
    dplyr::group_by(samples) |>
    dplyr::mutate(relative.value = value * 100 / sum(value))
  
  # Prevent log(0) errors by setting a minimum floor for the finest mesh size
  all[all$meshsize == 0, "meshsize"] <- 0.1
  
  all <- all |>
    dplyr::group_by(samples) |>
    dplyr::arrange(meshsize) |>
    # Calculate geometric midpoints for logarithmic-based statistical calculations
    dplyr::mutate(meshsize.midpoint = sqrt(dplyr::lag(meshsize) * meshsize) / 1000)
  
  all[is.na(all$meshsize.midpoint), "meshsize.midpoint"] <- all[
    is.na(all$meshsize.midpoint),
    "meshsize"
  ] / 1000

  # Compute moments using natural logarithms (exp/log) for mathematical consistency
  geom <- all |>
    dplyr::group_by(samples) |>
    dplyr::mutate(
      mean.geom.um = exp(sum(relative.value * log(meshsize.midpoint)) / 100)
    ) |>
    dplyr::mutate(
      sd.geom.um = exp(sqrt(
        sum(relative.value * (log(meshsize.midpoint) - log(mean.geom.um))^2) / 100
      )),
      sd.log.phi = log2(sd.geom.um)
    ) |>
    dplyr::mutate(
      skewness.geom.um = sum(
        relative.value * (log(meshsize.midpoint) - log(mean.geom.um))^3
      ) / (100 * log(sd.geom.um)^3),
      skewness.log.phi = -1 * skewness.geom.um
    ) |>
    dplyr::mutate(
      kurtosis.geom.um = sum(
        relative.value * (log(meshsize.midpoint) - log(mean.geom.um))^4
      ) / (100 * log(sd.geom.um)^4),
      kurtosis.log.phi = kurtosis.geom.um
    ) |>
    dplyr::mutate(
      mean.geom.um = mean.geom.um * 1000,
      mean.log.phi = G2Sd:::.um2phi(mean.geom.um)
    ) |>
    dplyr::select(-c("meshsize", "value", "relative.value", "meshsize.midpoint")) |>
    dplyr::distinct()

  return(geom)
}

# Calculate graphical statistics according to Folk and Ward (1957)
.fowa.stat <- function(x, decreasing) {
  x = as.data.frame(x)
  mat.D <- G2Sd:::.percentile(x, decreasing)
  mat.D$percentile <- paste("D", mat.D$percentile, sep = "")
  
  # Prepare data in wide format for both Phi and Metric scales
  all.phi <- mat.D |>
    dplyr::select(-meshsize) |>
    tidyr::pivot_wider(names_from = percentile, values_from = phi)
  
  all.um <- mat.D |>
    dplyr::select(-phi) |>
    dplyr::mutate(meshsize = meshsize / 1000) |>
    tidyr::pivot_wider(names_from = percentile, values_from = meshsize)

  # Compute indices on the Phi scale (Standard Folk & Ward equations)
  all.phi <- all.phi |>
    dplyr::group_by(samples) |>
    dplyr::mutate(
      mean.fw.phi = (D16 + D50 + D84) / 3,
      sd.fw.phi = (D16 - D84) / 4 + (D5 - D95) / 6.6,
      skewness.fw.phi = (D16 + D84 - 2 * D50) / (2 * (D16 - D84)) + 
                        (D5 + D95 - 2 * D50) / (2 * (D5 - D95)),
      kurtosis.fw.phi = (D5 - D95) / (2.44 * (D25 - D75))
    ) |>
    dplyr::select(-c(tidyselect::starts_with("D")))

  # Compute indices on the metric scale using logarithmic transformations
  all.um <- all.um |>
    dplyr::group_by(samples) |>
    dplyr::mutate(
      mean.fw.um = exp((log(D16) + log(D50) + log(D84)) / 3) * 1000,
      sd.fw.um = exp((log(D84) - log(D16)) / 4 + (log(D95) - log(D5)) / 6.6),
      skewness.fw.um = (log(D16) + log(D84) - 2 * log(D50)) / (2 * (log(D84) - log(D16))) + 
                       (log(D5) + log(D95) - 2 * log(D50)) / (2 * (log(D95) - log(D5))),
      kurtosis.fw.um = (log(D95) - log(D5)) / (2.44 * (log(D75) - log(D25)))
    ) |>
    dplyr::select(-c(tidyselect::starts_with("D")))

  all <- dplyr::full_join(all.um, all.phi, by = "samples")

  # Apply qualitative classification labels based on index values
  all <- all |>
    dplyr::mutate(
      mean.descript = dplyr::case_when(
        mean.fw.phi <= -5 ~ "Very Coarse Gravel",
        mean.fw.phi > -5 & mean.fw.phi <= -4 ~ "Coarse Gravel",
        mean.fw.phi > -4 & mean.fw.phi <= -3 ~ "Medium Gravel",
        mean.fw.phi > -3 & mean.fw.phi <= -2 ~ "Fine Gravel",
        mean.fw.phi > -2 & mean.fw.phi <= -1 ~ "Very Fine Gravel",
        mean.fw.phi > -1 & mean.fw.phi <= 0 ~ "Very Coarse Sand",
        mean.fw.phi > 0 & mean.fw.phi <= 1 ~ "Coarse Sand",
        mean.fw.phi > 1 & mean.fw.phi <= 2 ~ "Medium Sand",
        mean.fw.phi > 2 & mean.fw.phi <= 3 ~ "Fine Sand",
        mean.fw.phi > 3 & mean.fw.phi <= 4 ~ "Very Fine Sand",
        mean.fw.phi > 4 & mean.fw.phi <= 5 ~ "Very Coarse Silt",
        mean.fw.phi > 5 & mean.fw.phi <= 6 ~ "Coarse Silt",
        mean.fw.phi > 6 & mean.fw.phi <= 7 ~ "Medium Silt",
        mean.fw.phi > 7 & mean.fw.phi <= 8 ~ "Fine Silt",
        mean.fw.phi > 8 & mean.fw.phi <= 9 ~ "Very Fine Silt",
        mean.fw.phi > 9 ~ "Clay",
        .default = "NA"
      ),
      sorting = dplyr::case_when(
        sd.fw.phi < 0.35 ~ "Very Well Sorted",
        sd.fw.phi >= 0.35 & sd.fw.phi < 0.5 ~ "Well Sorted",
        sd.fw.phi >= 0.5 & sd.fw.phi < 0.7 ~ "Moderately Well Sorted",
        sd.fw.phi >= 0.7 & sd.fw.phi < 1 ~ "Moderately Sorted",
        sd.fw.phi >= 1 & sd.fw.phi < 2 ~ "Poorly Sorted",
        sd.fw.phi >= 2 & sd.fw.phi < 4 ~ "Very Poorly Sorted",
        sd.fw.phi >= 4 ~ "Extremely Poorly Sorted",
        .default = "NA"
      ),
      skewness = dplyr::case_when(
        skewness.fw.phi >= 0.3 ~ "Very Fine Skewed",
        skewness.fw.phi >= 0.1 & skewness.fw.phi < 0.3 ~ "Fine Skewed",
        skewness.fw.phi > -0.1 & skewness.fw.phi < 0.1 ~ "Symmetrical",
        skewness.fw.phi > -0.3 & skewness.fw.phi <= -0.1 ~ "Coarse Skewed",
        skewness.fw.phi <= -0.3 ~ "Very Coarse Skewed",
        .default = "NA"
      ),
      kurtosis = dplyr::case_when(
        kurtosis.fw.phi < 0.67 ~ "Very Platykurtic",
        kurtosis.fw.phi >= 0.67 & kurtosis.fw.phi < 0.9 ~ "Platykurtic",
        kurtosis.fw.phi >= 0.9 & kurtosis.fw.phi <= 1.11 ~ "Mesokurtic",
        kurtosis.fw.phi > 1.11 & kurtosis.fw.phi <= 1.5 ~ "Leptokurtic",
        kurtosis.fw.phi > 1.5 & kurtosis.fw.phi <= 3 ~ "Very Leptokurtic",
        kurtosis.fw.phi > 3 ~ "Extremely Leptokurtic",
        .default = "NA"
      )
    )

  # Concatenate descriptions into a single textural string for overall classification
  folk.ward <- all |>
    tidyr::unite(
      mean.descript, sorting, skewness, kurtosis,
      col = "sediment", sep = "/"
    )

  return(folk.ward)
}