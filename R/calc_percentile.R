# Function to calculate specific percentiles from the cumulative grain-size distribution
.percentile <- function(x, decreasing) {
  x <- G2Sd:::.g2sd_tidy(x)
  
  # Replace zero values with a minimal threshold to prevent infinite values during Phi logarithmic conversion
  x[x$meshsize == 0, "meshsize"] <- 0.01
  x <- x |> dplyr::mutate(phi = G2Sd:::.um2phi(meshsize))

  # Compute relative and cumulative percentages, sorting from finest to coarsest particles
  all <- x |>
    dplyr::group_by(samples) |>
    dplyr::mutate(value.relative = value / sum(value)) |>
    dplyr::arrange(meshsize) |> 
    dplyr::mutate(cum.sum = cumsum(value.relative))

  # Define the standard target percentiles (from D5 to D95)
  D = c(0.05, 0.1, 0.16, 0.25, 0.5, 0.75, 0.84, 0.9, 0.95)
  mat.D <- NULL

  # Iterate through each target percentile to perform linear interpolation between known sieve sizes
  for (i in 1:length(D)) {
    # Identify the lower boundary (the largest sieve size smaller than the target percentile)
    Dmin <- all |> 
      dplyr::group_by(samples) |> 
      dplyr::filter(cum.sum <= D[i]) |>
      dplyr::filter(cum.sum == max(cum.sum)) |> 
      dplyr::filter(phi == max(phi)) |>
      dplyr::select(-c("value", "value.relative")) |>
      dplyr::rename(min = cum.sum, meshsize.min = meshsize, phi.min = phi)
    
    # Identify the upper boundary (the smallest sieve size larger than the target percentile)
    Dmax <- all |> 
      dplyr::group_by(samples) |> 
      dplyr::filter(cum.sum >= D[i]) |>
      dplyr::filter(cum.sum == min(cum.sum)) |> 
      dplyr::filter(phi == min(phi)) |>
      dplyr::select(-c("value", "value.relative")) |>
      dplyr::rename(max = cum.sum, meshsize.max = meshsize, phi.max = phi)
    
    # Calculate the exact position (ratio) of the percentile between the two boundaries
    Dall <- dplyr::full_join(Dmin, Dmax, by = "samples") |>
      dplyr::mutate(ratio = (D[i] - min) / (max - min)) |>
      dplyr::mutate(phi = ((phi.max - phi.min) * ratio) + phi.min)
    
    # Handle cases where the target percentile falls exactly on a sieve size
    Dall[(na.omit(Dall$phi.min) == na.omit(Dall$phi.max)), "phi"] <- Dall[(na.omit(Dall$phi.min) == na.omit(Dall$phi.max)), "phi.min"]
    
    # Convert interpolated Phi values back to metric units (micrometers)
    Dall <- Dall |>
      dplyr::mutate(meshsize = G2Sd:::.phi2um(phi)) |>
      dplyr::select(samples, phi, meshsize) |>
      dplyr::mutate(percentile = D[i] * 100)
      
    mat.D <- dplyr::bind_rows(mat.D, Dall)
  }
  return(mat.D)
}