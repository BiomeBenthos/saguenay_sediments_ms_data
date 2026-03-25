# Script to process sediments samples with G2Sd
library(G2Sd)

# Replace internal functions with custom versions
source("R/sed_classes.R")
source("R/calc_stats.R")
source("R/process_sediment_data.R")
source("R/calc_percentile.R")
assignInNamespace(
  x = ".sedim.descript", 
  value = sed_classes, 
  ns = "G2Sd"
)
utils::assignInNamespace(
  x = ".moment.arith", 
  value = .moment.arith, 
  ns = "G2Sd"
)
utils::assignInNamespace(
  x = ".moment.geom", 
  value = .moment.geom, 
  ns = "G2Sd"
)
utils::assignInNamespace(
  x = ".fowa.stat", 
  value = .fowa.stat, 
  ns = "G2Sd"
)
utils::assignInNamespace(
  x = ".percentile", 
  value = .percentile, 
  ns = "G2Sd"
)

#---------- 2022-2023 Data ----------#

dat2022_2023 <- readxl::read_xlsx(
  "data/raw/judith_granulo.xlsx", 
  sheet = "Merged 100pct"
) |>
  tibble::column_to_rownames(var = "...1") |>
  process_sediment_data()

dat2022_2023_events <- readxl::read_xlsx(
  "data/raw/enviro_Judith.xlsx"
) |>
  dplyr::select(
    Station, lat, long, date
  ) |>
  dplyr::mutate(
    date = as.Date(date)
  ) |>
  dplyr::rename(
    station_id = Station,
    latitude = lat,
    longitude = long,
    sampling_date = date
  ) |>
  dplyr::mutate(
    sample_id = sprintf("%s_vanveen_%s", station_id, sampling_date),
    measurement_id = sprintf("%s_grain_size", sample_id)
  )

dat2022_2023 <- dplyr::left_join(
  dat2022_2023,
  dat2022_2023_events,
  by = "station_id"
) |>
  dplyr::relocate(
    c(sample_id, measurement_id, sampling_date, latitude, longitude),
    .after = station_id
  )

#---------- 2024 Data ----------#

dat2024 <- readxl::read_xlsx(
  "data/raw/Rebecca Howman Sediment Data.xlsx",
  sheet = "Sheet1"
) |> 
  dplyr::arrange(dplyr::desc(`Size (μm)`)) |>
  tibble::column_to_rownames(var = "Size (μm)") |>
  process_sediment_data()

dat2024_events <- readxl::read_xlsx(
  "data/raw/Rebecca Howman Sediment Data.xlsx",
  sheet = "Sheet2"
) |>
  dplyr::rename(
    station_id = Station_ID,
    latitude = plan_latdd,
    longitude = plan_longdd,
    sampling_date = Date
  ) |>
  dplyr::mutate(
    sample_id = sprintf("%s_vanveen_%s", station_id, sampling_date),
    measurement_id = sprintf("%s_grain_size", sample_id)
  )

dat2024 <- dplyr::left_join(
  dat2024,
  dat2024_events,
  by = "station_id"
) |>
  dplyr::relocate(
    c(sample_id, measurement_id, sampling_date, latitude, longitude),
    .after = station_id
  )

#---------- MERGE ----------#

df_final <- dplyr::bind_rows(
  dat2022_2023,
  dat2024
) |>
  dplyr::mutate(
    sampling_protocol = "vanveen"
  ) |>
  dplyr::relocate(
    sampling_protocol, .after = sampling_date
  )

write.csv(
  df_final, 
  "data/processed/biome_saguenay_sediments_2022-2024.csv", 
  row.names = FALSE
)
