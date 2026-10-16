# PHASE 1 — Ghana renewable energy and spatial inequality research
# Data validation and geographic integration (no causal claims)
# Input folder: Desktop/PETER AND KONADU SSRN
# Outputs: same folder / Phase1_Outputs

required <- c('readr', 'dplyr', 'tibble', 'ggplot2', 'sf', 'terra')
missing_packages <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing_packages)) {
  message('Install these packages before running: ', paste(missing_packages, collapse=', '))
  install.packages(missing_packages, repos='https://cloud.r-project.org')
}
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(sf); library(terra)
})

home <- Sys.getenv('USERPROFILE')
folder <- file.path(home, 'Desktop', 'PETER AND KONADU SSRN')
if (!dir.exists(folder)) {
  candidates <- c(file.path(home,'OneDrive','Desktop','PETER AND KONADU SSRN'),
                  file.path(home,'OneDrive - Personal','Desktop','PETER AND KONADU SSRN'))
  found <- candidates[dir.exists(candidates)]
  if (length(found)) folder <- found[1]
}
if (!dir.exists(folder)) stop('Folder not found. Edit variable folder to your exact Desktop path: ', folder)
output <- file.path(folder, 'Phase1_Outputs')
dir.create(output, showWarnings=FALSE, recursive=TRUE)
message('Reading from: ', folder)
message('Saving results to: ', output)

csv_file <- file.path(folder, 'ghana_dre_atlas_settlements.csv')
if (!file.exists(csv_file)) stop('Missing required CSV: ', csv_file)
raw <- readr::read_csv(csv_file, show_col_types=FALSE, progress=FALSE)
required_vars <- c('geohash','lat','lon','population','mean_rwi', 'pv_value',
 'distance_to_existing_transmission_lines','has_nightlight', 'admin_cgaz_1','admin_cgaz_2')
missing_vars <- setdiff(required_vars, names(raw))
if (length(missing_vars)) stop('Missing fields: ', paste(missing_vars, collapse=', '))

# Variable-by-variable completeness report; NA plus blank strings
is_missing <- function(x) if (is.character(x)) is.na(x) | trimws(x)=='' else is.na(x)
variable_report <- tibble(variable=names(raw),
 missing_n=vapply(raw, function(x) sum(is_missing(x)), integer(1)),
 observed_n=vapply(raw, function(x) sum(!is_missing(x)), integer(1))) %>%
 mutate(missing_pct=round(100*missing_n/nrow(raw), 2)) %>%
 arrange(desc(missing_pct))
write_csv(variable_report, file.path(output,'01_variable_missingness.csv'))

# Flag rather than silently remove extreme/invalid records
frame <- raw %>% mutate(
 latitude=as.numeric(lat), longitude=as.numeric(lon),
 invalid_coordinates=is.na(latitude)|is.na(longitude)|
   latitude < 4 | latitude > 12 | longitude < -4 | longitude > 2,
 missing_rwi=is.na(mean_rwi), missing_solar=is.na(pv_value),
 invalid_population=is.na(population)|population < 0,
 missing_region=is.na(admin_cgaz_1)|trimws(as.character(admin_cgaz_1))=='',
 missing_district=is.na(admin_cgaz_2)|trimws(as.character(admin_cgaz_2))=='')

# Top population values are diagnostic only, not discarded.
write_csv(frame %>% arrange(desc(population)) %>%
 select(geohash,village_name,admin_cgaz_1,admin_cgaz_2,population,lat,lon) %>% head(30),
 file.path(output,'02_largest_population_records.csv'))
write_csv(frame %>% count(admin_cgaz_1, name='settlements') %>%
 arrange(desc(settlements)), file.path(output,'03_region_label_counts.csv'))
write_csv(frame %>% filter(invalid_coordinates|invalid_population|missing_region|missing_district) %>%
 select(geohash,village_name,lat,lon,admin_cgaz_1,admin_cgaz_2,
 invalid_coordinates,invalid_population,missing_region,missing_district),
 file.path(output,'04_flagged_settlements.csv'))

# Geohash duplicates warrant investigation (not automatic de-duplication)
duplicates <- frame %>% count(geohash, name='records') %>% filter(records > 1)
write_csv(duplicates, file.path(output,'05_duplicate_geohashes.csv'))

# Create points for valid coordinates; WGS84 EPSG:4326
valid <- frame %>% filter(!invalid_coordinates)
settlement_sf <- st_as_sf(valid, coords=c('longitude','latitude'), crs=4326, remove=FALSE)

# Extract administrative archives from the uploaded geofiles.zip.
# Inner shapefile ZIPs have different district vintages and must not be conflated.
geo_zip <- file.path(folder,'geofiles.zip')
if (file.exists(geo_zip)) {
  expanded <- file.path(output,'Geofiles_Unpacked')
  dir.create(expanded, showWarnings=FALSE, recursive=TRUE)
  unzip(geo_zip, exdir=expanded)
  inner_archives <- list.files(expanded, pattern='\\.zip$', recursive=TRUE, full.names=TRUE)
  for (archive in inner_archives) {
    destination <- file.path(expanded, tools::file_path_sans_ext(basename(archive)))
    dir.create(destination, showWarnings=FALSE, recursive=TRUE)
    unzip(archive, exdir=destination)
  }
  shp_files <- list.files(expanded, pattern='\\.shp$', recursive=TRUE, full.names=TRUE)
  join_reports <- list()
  if (length(shp_files)) for (path in shp_files) {
    nm <- tools::file_path_sans_ext(basename(path))
    message('Testing administrative boundary: ', nm)
    districts <- tryCatch(st_read(path, quiet=TRUE), error=function(e) NULL)
    if (is.null(districts)) next
    districts <- tryCatch(st_make_valid(st_transform(districts,4326)),error=function(e) NULL)
    if (is.null(districts)) next
    inside <- lengths(st_intersects(settlement_sf, districts))
    join_reports[[nm]] <- tibble(boundary_set=nm,
      polygon_count=nrow(districts),
      points_inside_any=sum(inside>0), points_outside_all=sum(inside==0),
      points_in_multiple=sum(inside>1))
    # Count-only: do not assign arbitrary district where polygons overlap.
  }
  if(length(join_reports)) write_csv(bind_rows(join_reports),
      file.path(output,'06_boundary_comparison.csv'))
} else message('Optional geofiles.zip not found; skipped polygon boundary validation.')

# Solar raster: year-round average-daily package preferred, avoiding doubled-up archives.
solar_zip <- file.path(folder,'Ghana_GISdata_LTAym_AvgDailyTotals_GlobalSolarAtlas-v2_GEOTIFF.zip')
if(file.exists(solar_zip)) {
  raster_dir <- file.path(output,'Solar_Unpacked')
  dir.create(raster_dir, showWarnings=FALSE, recursive=TRUE)
  zip_list <- unzip(solar_zip, list=TRUE)$Name
  # Only top-level rasters, not their monthly counterparts or XML/PDF metadata.
  chosen <- zip_list[grepl('/(PVOUT|GHI)\\.tif$',zip_list)]
  if(length(chosen)) unzip(solar_zip, files=chosen, exdir=raster_dir)
  tif_paths <- list.files(raster_dir, pattern='(PVOUT|GHI)\\.tif$',recursive=TRUE,full.names=TRUE)
  if(length(tif_paths)) {
    geo_points <- terra::vect(settlement_sf)
    for(tif in tif_paths) {
      label <- tools::file_path_sans_ext(basename(tif))
      rr <- tryCatch(terra::rast(tif),error=function(e) NULL)
      if(is.null(rr)) next
      message('Extracting ',label,'; check raster metadata to confirm units.')
      vals <- terra::extract(rr, terra::project(geo_points,terra::crs(rr)))
      valid[[paste0('raster_',tolower(label))]] <- vals[[2]]
    }
    # Note: raster values are not automatically interchangeable with pv_value.
    write_csv(valid %>% select(geohash,starts_with('raster_'),pv_value),
              file.path(output,'07_solar_raster_extraction.csv'))
  }
} else message('Solar GeoTIFF archive not found; skipped raster extraction.')

# Initial distribution plots only; never interpret correlations as causal estimates.
p <- ggplot(frame %>% filter(!is.na(mean_rwi),is.finite(mean_rwi)),
            aes(mean_rwi)) + geom_histogram(bins=40,fill='#156082',color='white') +
 labs(title='Distribution of settlement relative wealth index',
      subtitle='Observed settlements only; missing records excluded from this chart',
      x='Relative wealth index (source definition)',y='Number of settlements') + theme_minimal(base_size=12)
ggsave(file.path(output,'08_relative_wealth_distribution.png'),p,width=9,height=5,dpi=220)

summary_table <- tibble(
 metric=c('total_settlement_records','unique_nonmissing_geohashes','duplicate_geohash_groups',
          'invalid_or_outside_plausible_ghana_coordinates','missing_relative_wealth',
          'missing_atlas_pv_value','missing_region','missing_district','invalid_population'),
 value=c(nrow(frame),n_distinct(frame$geohash,na.rm=TRUE),nrow(duplicates),
 sum(frame$invalid_coordinates),sum(frame$missing_rwi),sum(frame$missing_solar),
 sum(frame$missing_region),sum(frame$missing_district),sum(frame$invalid_population)))
write_csv(summary_table,file.path(output,'00_validation_summary.csv'))
writeLines(c('PHASE 1 DATA VALIDATION — AUTOMATED NOTES',
 'All figures are calculated from local input files; no observations invented.',
 'Coordinates checked against broad plausibility bounds, not official national polygons.',
 'Inspect 06_boundary_comparison.csv before choosing a district boundary vintage.',
 'Do not interpret has_nightlight as confirmed electricity access.',
 'Do not interpret mean_rwi as directly measured household income.',
 'Do not interpret pv_value as actual electricity generation.',
 'PVOUT/GHI raster units must be confirmed from the included metadata.',
 'Population outliers have been flagged for review, not removed.',
 'The data are cross-sectional and cannot independently identify causal effects.'),
 file.path(output,'README_interpretation.txt'))
print(summary_table)
message('PHASE 1 COMPLETE. Open folder: ',output)