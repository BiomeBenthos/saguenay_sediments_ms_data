# Assign detailed grain-size classes
sed_classes <- function(x) {
  x <- G2Sd:::.g2sd_tidy(x)

  sedim <- c(63000, 31500, 2^c(4:-3) * 1000, 62.5, 31.25, 15.625, 7.8125, 3.90625, 1.953125)

  # Explicitly define all possible classes for the factor
  noms_classes <- c(
    "boulder", "vcgravel", "cgravel", "mgravel", "fgravel", "vfgravel",
    "vcsand", "csand", "msand", "fsand", "vfsand",
    "vcsilt", "csilt", "msilt", "fsilt", "vfsilt", "clay"
  )

  all <- x |>
    dplyr::group_by(samples) |>
    dplyr::mutate(relative.value = value * 100 / sum(value)) |>
    dplyr::mutate(
      class = dplyr::case_when(
        meshsize >= sedim[1] ~ "boulder",
        meshsize < sedim[1] & meshsize >= sedim[2] ~ "vcgravel",
        meshsize < sedim[2] & meshsize >= sedim[3] ~ "cgravel",
        meshsize < sedim[3] & meshsize >= sedim[4] ~ "mgravel",
        meshsize < sedim[4] & meshsize >= sedim[5] ~ "fgravel",
        meshsize < sedim[5] & meshsize >= sedim[6] ~ "vfgravel",
        meshsize < sedim[6] & meshsize >= sedim[7] ~ "vcsand",
        meshsize < sedim[7] & meshsize >= sedim[8] ~ "csand",
        meshsize < sedim[8] & meshsize >= sedim[9] ~ "msand",
        meshsize < sedim[9] & meshsize >= sedim[10] ~ "fsand",
        meshsize < sedim[10] & meshsize >= sedim[11] ~ "vfsand",
        meshsize < sedim[11] & meshsize >= sedim[12] ~ "vcsilt",
        meshsize < sedim[12] & meshsize >= sedim[13] ~ "csilt",
        meshsize < sedim[13] & meshsize >= sedim[14] ~ "msilt",
        meshsize < sedim[14] & meshsize >= sedim[15] ~ "fsilt",
        meshsize < sedim[15] & meshsize >= sedim[16] ~ "vfsilt",
        meshsize < sedim[16] ~ "clay",
        .default = "NA"
      )
    ) |>
    # Convert to factor to ensure existence of every category
    dplyr::mutate(class = factor(class, levels = noms_classes))

  sediment <- all |>
    # The .drop argument prevents dplyr from removing empty categories
    dplyr::group_by(samples, class, .drop = FALSE) |>
    dplyr::summarise(value = sum(relative.value), .groups = "drop") |>
    # The values_fill argument sets 0 instead of NA for missing classes
    tidyr::pivot_wider(names_from = class, values_from = value, values_fill = 0) |>
    dplyr::select(
      samples,
      tidyselect::all_of(noms_classes)
    )
    
  return(sediment)
}
