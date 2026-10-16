# Phase 2C: Audit non-Ghana region labels and repair HC3 inference.
# Original data are never overwritten; output is diagnostic, not a causal analysis.
packages <- c('readr','dplyr','tibble','sandwich')
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing)) install.packages(missing, repos='https://cloud.r-project.org')
suppressPackageStartupMessages({library(readr);library(dplyr);library(tibble)})
home <- Sys.getenv('USERPROFILE')
candidates <- c(file.path(home,'OneDrive','Desktop','PETER AND KONADU SSRN'),file.path(home,'Desktop','PETER AND KONADU SSRN'),file.path(home,'OneDrive - Personal','Desktop','PETER AND KONADU SSRN'))
folder <- candidates[dir.exists(candidates)][1]
if(is.na(folder)) stop('Folder not found; update candidates near the top of the script.')
out <- file.path(folder, 'Phase2C_Outputs');dir.create(out,showWarnings=FALSE,recursive=TRUE)
raw <- read_csv(file.path(folder,'ghana_dre_atlas_settlements.csv'),show_col_types=FALSE,progress=FALSE)
needed <- c('geohash','admin_cgaz_1','admin_cgaz_2','mean_rwi','pv_value','distance_to_existing_transmission_lines','dist_main_road_km','population')
stopifnot(all(needed %in% names(raw)))
# The officially recognized Ghana region names; not estimated from data.
ghana_regions <- c('Ahafo Region','Ashanti Region','Bono Region','Bono East Region','Central Region','Eastern Region','Greater Accra Region','North East Region','Northern Region','Oti Region','Savannah Region','Upper East Region','Upper West Region','Volta Region','Western Region','Western North Region')
d <- raw %>% transmute(geohash=as.character(geohash),region=na_if(trimws(as.character(admin_cgaz_1)),''),district=na_if(trimws(as.character(admin_cgaz_2)),''),rwi=as.numeric(mean_rwi),pv=as.numeric(pv_value),grid=as.numeric(distance_to_existing_transmission_lines),road=as.numeric(dist_main_road_km),pop=as.numeric(population)) %>% mutate(log_grid=log1p(grid),log_road=log1p(road),log_pop=log1p(pop),pv_100=pv/100,region_district=paste(region,district,sep=' | '),ghana_region=region %in% ghana_regions)
# Keep all unusual rows for independent geographic investigation.
flagged <- raw %>% mutate(region_check=trimws(as.character(admin_cgaz_1))) %>% filter(is.na(region_check) | !region_check %in% ghana_regions)
write_csv(flagged,file.path(out,'01_non_ghana_or_missing_region_FULL_RECORDS.csv'))
cat('\nUnrecognized/missing region records in original data:',nrow(flagged),'\n')
print(d %>% count(region,ghana_region,sort=TRUE))
x <- d %>% filter(is.finite(rwi),is.finite(pv_100),is.finite(log_grid),is.finite(log_road),is.finite(log_pop),!is.na(region),!is.na(district))
f <- rwi ~ log_grid + pv_100 + log_road + log_pop + region
summarize_fit <- function(df,label){
  df <- df %>% mutate(region=factor(region))
  fit <- lm(f,data=df)
  h <- hatvalues(fit)
  G <- dplyr::n_distinct(df$region_district)
  out1 <- tibble(sample=label,n=nobs(fit),clusters=G,rank=fit$rank,max_leverage=max(h),n_leverage_near_one=sum((1-h)<1e-8),n_leverage_above_099=sum(h>=.99))
  calc <- function(V,method,df_t=Inf){
    if(is.null(V)) return(tibble())
    b <- coef(fit); se <- sqrt(diag(V)); tst <- b/se
    p <- if(is.finite(df_t)) 2*pt(abs(tst),df=df_t,lower.tail=FALSE) else 2*pnorm(abs(tst),lower.tail=FALSE)
    tibble(sample=label,method=method,term=names(b),estimate=as.numeric(b),se=as.numeric(se),statistic=as.numeric(tst),p_value=as.numeric(p))
  }
  hc1 <- tryCatch(sandwich::vcovHC(fit,type='HC1'),error=function(e) NULL)
  hc3 <- tryCatch(sandwich::vcovHC(fit,type='HC3'),error=function(e) NULL)
  cl <- tryCatch(sandwich::vcovCL(fit,cluster=df$region_district,type='HC1'),error=function(e) NULL)
  list(diagnostic=out1,coef=bind_rows(calc(hc1,'HC1'),calc(hc3,'HC3'),calc(cl,'district_cluster_HC1_t_Gminus1',G-1)))
}
all_result <- summarize_fit(x,'original_complete_case')
# Prespecified geographical restriction: Ghana's 16 recognized regions.
ghana_x <- x %>% filter(ghana_region)
ghana_result <- summarize_fit(ghana_x,'ghana_16_regions_only')
checks <- bind_rows(all_result$diagnostic,ghana_result$diagnostic)
coeffs <- bind_rows(all_result$coef,ghana_result$coef)
write_csv(checks,file.path(out,'02_inference_and_leverage_comparison.csv'))
write_csv(coeffs,file.path(out,'03_regression_inference_comparison.csv'))
write_csv(coeffs %>% filter(term %in% c('log_grid','pv_100','log_road','log_pop')),file.path(out,'04_four_primary_coefficient_comparison.csv'))
write_csv(x %>% filter(!ghana_region),file.path(out,'05_excluded_from_geographic_sensitivity_complete_cases.csv'))
cat('\nINFERENCE AND LEVERAGE COMPARISON:\n');print(checks)
cat('\nFOUR MAIN COEFFICIENTS:\n');print(coeffs %>% filter(term %in% c('log_grid','pv_100','log_road','log_pop')),n=30,width=Inf)
cat('\nPhase 2C complete. Files saved at: ',out,'\n',sep='')
cat('IMPORTANT: Examine flagged coordinates and polygon membership before choosing final analysis sample.\n')
cat('HC3 becoming finite is not sufficient proof that all other model assumptions are satisfied.\n')