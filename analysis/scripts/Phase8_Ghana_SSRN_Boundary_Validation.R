# PHASE 8 — Ghana SSRN administrative boundary and district-label audit
# NO silent correction of geographic labels; both polygon vintages are tested independently.
# Cross-sectional associations are not causal effects or verified electricity access.
packages <- c('sf','readr','dplyr','tibble','ggplot2','stringr','tidyr')
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing_packages)) install.packages(missing_packages,repos='https://cloud.r-project.org')
suppressPackageStartupMessages({library(sf);library(readr);library(dplyr);library(tibble);library(ggplot2);library(stringr)})
home <- Sys.getenv('USERPROFILE')
candidates <- c(file.path(home,'OneDrive','Desktop','PETER AND KONADU SSRN'),
                file.path(home,'Desktop','PETER AND KONADU SSRN'),
                file.path(home,'OneDrive - Personal','Desktop','PETER AND KONADU SSRN'))
folder <- candidates[dir.exists(candidates)][1]
if(length(folder)==0 || is.na(folder)) stop('Research directory not found; edit folder manually.')
input_csv <- file.path(folder,'ghana_dre_atlas_settlements.csv')
geo_zip <- file.path(folder,'geofiles.zip')
if(!file.exists(input_csv) || !file.exists(geo_zip)) stop('Missing CSV or geofiles.zip in ',folder)
out <- file.path(folder,'Phase8_Outputs')
dir.create(out,showWarnings=FALSE,recursive=TRUE)
raw <- read_csv(input_csv,show_col_types=FALSE)
needed <- c('geohash','lat','lon','admin_cgaz_1','admin_cgaz_2')
if(!all(needed %in% names(raw))) stop('Missing required CSV columns: ',paste(setdiff(needed,names(raw)),collapse=', '))
if(anyDuplicated(raw$geohash)) stop('Duplicate geohash: resolve before polygon joins.')
raw <- raw %>% mutate(.source_row=row_number())
if(any(!is.finite(raw$lon)|!is.finite(raw$lat))) stop('Non-finite coordinates; audit before mapping.')
points <- st_as_sf(raw,coords=c('lon','lat'),crs=4326,remove=FALSE)
# Robust normalization for text agreement only (not a replacement for official harmonization).
normalize <- function(x) {
 x <- iconv(as.character(x),to='ASCII//TRANSLIT')
 x <- tolower(trimws(x))
 x <- gsub('\\s+region$','',x)
 x <- gsub('[^a-z0-9]+',' ',x)
 trimws(gsub('\\s+',' ',x))
}
# Decompress the outer ZIP and each district archive into distinct folders.
outer <- file.path(out,'archive_unpack')
dir.create(outer,showWarnings=FALSE,recursive=TRUE)
unzip(geo_zip,exdir=outer)
nested <- list.files(outer,pattern='\\.zip$',recursive=TRUE,full.names=TRUE)
selected <- nested[grepl('Districts_261|Districts_271',basename(nested),ignore.case=TRUE)]
if(length(selected)!=2) stop('Expected Districts_261.zip and Districts_271.zip; found: ',paste(basename(selected),collapse=', '))
results <- list(); detailed <- list()
for(archive in selected) {
 set_name <- if(grepl('261',basename(archive))) 'District_261' else 'District_272'
 unpack <- file.path(out,'archive_unpack',set_name)
 dir.create(unpack,recursive=TRUE,showWarnings=FALSE)
 unzip(archive,exdir=unpack)
 shp <- list.files(unpack,pattern='\\.shp$',recursive=TRUE,full.names=TRUE)
 if(length(shp)!=1) stop('Expected one shapefile in ',archive,'; got ',length(shp))
 poly <- st_read(shp,quiet=TRUE)
 if(!all(c('Region','District') %in% names(poly))) stop('Missing Region/District fields in ',shp)
 poly <- st_make_valid(st_transform(poly,4326))
 # Prefer explicit intersects; multiple hits can be shared-border effects.
 hits <- st_intersects(points,poly)
 nhit <- lengths(hits)
 ix <- vapply(hits,function(i) if(length(i)==1) i else NA_integer_,integer(1))
 region_poly <- rep(NA_character_,nrow(raw))
 district_poly <- rep(NA_character_,nrow(raw))
 match_one <- which(!is.na(ix))
 region_poly[match_one] <- as.character(poly$Region[ix[match_one]])
 district_poly[match_one] <- as.character(poly$District[ix[match_one]])
 comparison <- raw %>% transmute(geohash,lat,lon,source_region=admin_cgaz_1,
                 source_district=admin_cgaz_2) %>%
   mutate(boundary_set=set_name,polygon_matches=nhit,
          polygon_region=region_poly,polygon_district=district_poly,
          region_agrees=if_else(!is.na(source_region)&!is.na(polygon_region),
                                normalize(source_region)==normalize(polygon_region),NA),
          district_agrees=if_else(!is.na(source_district)&!is.na(polygon_district),
                                  normalize(source_district)==normalize(polygon_district),NA))
 # Regions may use alternate abbreviations; do not interpret every string mismatch as error.
 detailed[[set_name]] <- comparison
 results[[set_name]] <- tibble(boundary_set=set_name,polygon_count=nrow(poly),
  settlement_count=nrow(raw),one_polygon=sum(nhit==1),outside_all=sum(nhit==0),
  multiple_polygons=sum(nhit>1),
  n_region_comparable=sum(!is.na(comparison$region_agrees)),
  n_region_string_mismatch=sum(comparison$region_agrees==FALSE,na.rm=TRUE),
  n_district_comparable=sum(!is.na(comparison$district_agrees)),
  n_district_string_mismatch=sum(comparison$district_agrees==FALSE,na.rm=TRUE))
 write_csv(comparison %>% filter(polygon_matches!=1|region_agrees==FALSE|district_agrees==FALSE|is.na(source_region)|is.na(source_district)),
           file.path(out,paste0('02_flagged_',set_name,'.csv')))
 cat('Spatial polygon comparison complete:',set_name,' polygons:',nrow(poly),'\n')
}
summary <- bind_rows(results)
write_csv(summary,file.path(out,'01_boundary_coverage_comparison.csv'))
all_detailed <- bind_rows(detailed)
write_csv(all_detailed,file.path(out,'03_all_settlement_polygon_audit.csv'))
# Compare polygon labels only when uniquely assigned in both archives.
wide <- detailed[['District_261']] %>% select(geohash,
  hits261=polygon_matches,region261=polygon_region,district261=polygon_district) %>%
  inner_join(detailed[['District_272']] %>% select(geohash,
   hits272=polygon_matches,region272=polygon_region,district272=polygon_district),by='geohash') %>%
  mutate(comparable=hits261==1&hits272==1,
         district_names_disagree=if_else(comparable,
           normalize(district261)!=normalize(district272),NA),
         region_names_disagree=if_else(comparable,
           normalize(region261)!=normalize(region272),NA))
write_csv(wide %>% filter(!comparable|district_names_disagree==TRUE|region_names_disagree==TRUE),
          file.path(out,'04_boundary_vintage_disagreements.csv'))
# District-level mismatch rates: diagnosis only; names may change with administrative vintage.
by_region <- all_detailed %>% group_by(boundary_set,source_region) %>% summarise(
 n=n(),out_of_polygon=sum(polygon_matches==0),multi_polygon=sum(polygon_matches>1),
 region_string_mismatch=sum(region_agrees==FALSE,na.rm=TRUE),
 district_string_mismatch=sum(district_agrees==FALSE,na.rm=TRUE),.groups='drop')
write_csv(by_region,file.path(out,'05_region_level_audit.csv'))
# Check which of the 24 problematic labels fall inside mapped Ghana districts.
problem <- c('Zanzan','Plateaux Region')
flagged <- all_detailed %>% filter(is.na(source_region)|source_region %in% problem)
write_csv(flagged,file.path(out,'06_out_of_country_or_missing_region_records.csv'))
# Figures show coverage rather than claiming that mismatches are necessarily incorrect.
p <- summary %>% select(boundary_set,one_polygon,outside_all,multiple_polygons) %>%
 tidyr::pivot_longer(-boundary_set,names_to='membership',values_to='settlements') %>%
 ggplot(aes(boundary_set,settlements,fill=membership))+geom_col()+
 labs(title='Ghana settlement-to-district polygon membership audit',
      subtitle='Unique, absent, and ambiguous spatial matches; alternative boundary vintages',
      x='Boundary archive',y='Settlements',fill='Spatial matching',
      caption='Polygon membership is not proof that administrative labels or dates are harmonized.')+
 theme_minimal(base_size=11)
ggsave(file.path(out,'07_boundary_coverage.png'),p,width=9,height=5.5,dpi=300)
notes <- c('PHASE 8 — GEOGRAPHIC DATA INTEGRITY','',
 'Do not automatically choose a boundary vintage based solely on higher spatial coverage.',
 'Compare source metadata, reference year, official district definitions and geometry before selection.',
 'Polygon overlap is ambiguous: the script never assigns an arbitrary first polygon.',
 'Text mismatches can be naming differences; inspect source and polygon metadata.',
 'Border communities may have uncertain assignment due to coastline/boundary precision.',
 'The 24 unrecognized/missing region-label records need manual review in output 06.',
 'A district fixed-effects model using original labels is not equivalent to one using polygon-derived labels.',
 'If changing district definitions, rerun primary regression and spatial diagnostics before manuscript.',
 'The cross-sectional RWI/infrastructure findings are noncausal.')
writeLines(notes,file.path(out,'08_readme_interpretation.txt'))
cat('\nPHASE 8 COMPLETE:',out,'\n')
cat('\nBOUNDARY COVERAGE COMPARISON:\n');print(summary,width=Inf)
cat('\nBOUNDARY VINTAGE CONFLICT SUMMARY:\n');print(wide %>% summarise(
 unique_both=sum(comparable),different_district_names=sum(district_names_disagree,na.rm=TRUE),
 different_region_names=sum(region_names_disagree,na.rm=TRUE),unresolved_spatial_match=sum(!comparable)),width=Inf)
cat('\nFLAGGED SOURCE LABELS:\n');print(flagged %>% count(boundary_set,source_region,polygon_matches),n=Inf)