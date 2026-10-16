# PHASE 2B — Inference audit for Ghana SSRN analysis
# Goal: diagnose HC3 NaNs; compare HC1, HC2, HC3, district-cluster robust inference.
# No coefficients or findings are fabricated; run in RStudio against original CSV.
packages <- c('readr','dplyr','tibble','sandwich','lmtest','ggplot2')
missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) install.packages(missing,repos='https://cloud.r-project.org')
suppressPackageStartupMessages({library(readr);library(dplyr);library(tibble);library(ggplot2)})
home <- Sys.getenv('USERPROFILE')
candidates <- c(file.path(home,'OneDrive','Desktop','PETER AND KONADU SSRN'),file.path(home,'Desktop','PETER AND KONADU SSRN'),file.path(home,'OneDrive - Personal','Desktop','PETER AND KONADU SSRN'))
folder <- candidates[dir.exists(candidates)][1]
if(is.na(folder)) stop('Research folder not found. Edit candidates to match your PC.')
out <- file.path(folder,'Phase2B_Outputs');dir.create(out,showWarnings=FALSE,recursive=TRUE)
raw <- read_csv(file.path(folder,'ghana_dre_atlas_settlements.csv'),show_col_types=FALSE,progress=FALSE)
needed <- c('geohash','admin_cgaz_1','admin_cgaz_2','mean_rwi','pv_value','distance_to_existing_transmission_lines','dist_main_road_km','population')
stopifnot(all(needed %in% names(raw)))
d <- raw %>% transmute(geohash=as.character(geohash),region=na_if(trimws(as.character(admin_cgaz_1)),''),district=na_if(trimws(as.character(admin_cgaz_2)),''),rwi=as.numeric(mean_rwi),pv=as.numeric(pv_value),grid=as.numeric(distance_to_existing_transmission_lines),road=as.numeric(dist_main_road_km),pop=as.numeric(population)) %>% mutate(log_grid=log1p(grid),log_road=log1p(road),log_pop=log1p(pop),pv_100=pv/100,region_district=paste(region,district,sep=' | '))
x <- d %>% filter(is.finite(rwi),is.finite(pv_100),is.finite(log_grid),is.finite(log_road),is.finite(log_pop),!is.na(region),!is.na(district)) %>% mutate(region=factor(region))
f <- rwi ~ log_grid+pv_100+log_road+log_pop+region
fit <- lm(f,data=x)
lev <- hatvalues(fit); residual <- resid(fit)
k <- fit$rank; n <- nobs(fit)
cat('Model n:',n,'; matrix rank:',k,'; maximum leverage:',max(lev),'\n')
# A leverage value of exactly one makes the HC3 formula undefined.
lev_table <- tibble(geohash=x$geohash,region=as.character(x$region),district=x$district,population=x$pop,rwi=x$rwi,pv=x$pv,grid_distance_km=x$grid,road_distance_km=x$road,leverage=as.numeric(lev),residual=as.numeric(residual))
write_csv(lev_table %>% arrange(desc(leverage)) %>% head(100),file.path(out,'01_top_100_leverage_records.csv'))
write_csv(tibble(metric=c('n_complete_cases','design_rank','district_clusters','leverage_ge_0_99','leverage_ge_0_999','leverage_near_one','max_leverage','min_one_minus_leverage','zero_residual_count'),value=c(n,k,n_distinct(x$region_district),sum(lev>=.99),sum(lev>=.999),sum((1-lev)<1e-8),max(lev),min(1-lev),sum(abs(residual)<1e-10))),file.path(out,'02_leverage_diagnostics.csv'))
# Safely collect full inference matrix: NA indicates undefined/unavailable, never zero.
get_result <- function(V,label,df=Inf){
  b <- coef(fit); s <- sqrt(diag(V));stat <- b/s
  p <- if(is.finite(df)) 2*pt(abs(stat),df=df,lower.tail=FALSE) else 2*pnorm(abs(stat),lower.tail=FALSE)
  tibble(method=label,term=names(b),estimate=as.numeric(b),se=as.numeric(s),statistic=as.numeric(stat),p_value=as.numeric(p),df=as.numeric(df))
}
methods <- list()
for(tp in c('HC0','HC1','HC2','HC3')){
  z <- tryCatch(sandwich::vcovHC(fit,type=tp),error=function(e) {message(tp,': ',conditionMessage(e));NULL})
  if(!is.null(z)) methods[[tp]] <- get_result(z,tp)
}
G <- n_distinct(x$region_district)
vc <- sandwich::vcovCL(fit,cluster=x$region_district,type='HC1')
methods[['district_cluster_HC1_t_Gminus1']] <- get_result(vc,'district_cluster_HC1_t_Gminus1',df=G-1)
result <- bind_rows(methods)
write_csv(result,file.path(out,'03_all_inference_comparison.csv'))
# Evaluate whether singleton region levels explain exact-hat observations.
region_count <- x %>% count(region,name='settlements') %>% arrange(settlements)
write_csv(region_count,file.path(out,'04_region_sample_counts.csv'))
# Model sensitivity with no ultra-high leverage observations; diagnostic only, not preferred model.
if(any(lev>=.99)){
  y <- droplevels(x[lev < .99,,drop=FALSE])
  if(nrow(y)>100 && n_distinct(y$region)>1){
    fit_y <- lm(f,data=y)
    g_y <- n_distinct(y$region_district)
    v_y <- sandwich::vcovCL(fit_y,cluster=y$region_district,type='HC1')
    b <- coef(fit_y);s <- sqrt(diag(v_y));t <- b/s
    res_y <- tibble(term=names(b),estimate=as.numeric(b),cluster_se=as.numeric(s),cluster_p=as.numeric(2*pt(abs(t),df=g_y-1,lower.tail=FALSE)))
    write_csv(res_y,file.path(out,'05_high_leverage_exclusion_diagnostic_ONLY.csv'))
    write_csv(tibble(metric=c('excluded_records','retained_records','retained_clusters','max_retained_leverage'),value=c(sum(lev>=.99),nrow(y),g_y,max(hatvalues(fit_y)))),file.path(out,'06_high_leverage_exclusion_counts.csv'))
  }
}
# Distinctive results are always generated from user's data, never hardcoded.
cat('\nINFERENCE RESULTS — primary coefficients:\n')
print(result %>% filter(term %in% c('log_grid','pv_100','log_road','log_pop')))
cat('\nREGION COUNTS:\n'); print(region_count)
cat('\nPHASE 2B COMPLETE: ',out,'\n',sep='')
cat('If HC3 is not finite, do not report HC3 p-values. Investigate 01_top_100_leverage_records.csv and 02_leverage_diagnostics.csv.\n')
cat('Clustered inference assumes clusters are adequate; geographic residual dependence still needs spatial diagnostics.\n')