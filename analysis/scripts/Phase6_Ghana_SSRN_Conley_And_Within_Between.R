# PHASE 6 — Ghana SSRN: Spatial HAC and within-/between-district solar variation
# Cross-sectional ASSOCIATIONS only. Do not claim causal effects or electrification.
# Uses actual World Bank DRE Atlas settlement CSV (no imputed or fabricated values).
packages <- c('readr','dplyr','tibble','ggplot2','fixest')
new <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(new)) install.packages(new,repos='https://cloud.r-project.org')
suppressPackageStartupMessages({library(readr);library(dplyr);library(tibble);library(ggplot2);library(fixest)})
options(fixest_notes=FALSE)
home <- Sys.getenv('USERPROFILE')
possible <- c(file.path(home,'OneDrive','Desktop','PETER AND KONADU SSRN'),
              file.path(home,'Desktop','PETER AND KONADU SSRN'),
              file.path(home,'OneDrive - Personal','Desktop','PETER AND KONADU SSRN'))
folder <- possible[dir.exists(possible)][1]
if(is.na(folder)) stop('Desktop research folder not found. Edit folder path in script.')
output <- file.path(folder,'Phase6_Outputs')
dir.create(output,recursive=TRUE,showWarnings=FALSE)
raw <- read_csv(file.path(folder,'ghana_dre_atlas_settlements.csv'),show_col_types=FALSE,progress=FALSE)
need <- c('geohash','admin_cgaz_1','admin_cgaz_2','mean_rwi','pv_value',
          'distance_to_existing_transmission_lines','dist_main_road_km','population','lat','lon')
if(!all(need %in% names(raw))) stop('Missing required columns: ',paste(setdiff(need,names(raw)),collapse=', '))
regions <- c('Ahafo Region','Ashanti Region','Bono Region','Bono East Region',
             'Central Region','Eastern Region','Greater Accra Region','North East Region',
             'Northern Region','Oti Region','Savannah Region','Upper East Region',
             'Upper West Region','Volta Region','Western Region','Western North Region')
d <- raw %>% transmute(geohash=as.character(geohash),
    region=trimws(as.character(admin_cgaz_1)),district=trimws(as.character(admin_cgaz_2)),
    rwi=as.numeric(mean_rwi),pv=as.numeric(pv_value),
    grid=as.numeric(distance_to_existing_transmission_lines),
    road=as.numeric(dist_main_road_km),pop=as.numeric(population),
    lat=as.numeric(lat),lon=as.numeric(lon)) %>%
  mutate(region=na_if(region,''),district=na_if(district,''),
         cluster=paste(region,district,sep=' | '), pv_100=pv/100,
         log_grid=log1p(grid),log_road=log1p(road),log_pop=log1p(pop)) %>%
  filter(region %in% regions,!is.na(district),is.finite(rwi),is.finite(pv_100),
         is.finite(log_grid),is.finite(log_road),is.finite(log_pop),
         is.finite(lat),is.finite(lon),grid>=0,road>=0,pop>=0)
stopifnot(nrow(d)>1000,!anyDuplicated(d$geohash))
d <- d %>% mutate(lat_c=lat-mean(lat),lon_c=lon-mean(lon),
                  lat2=lat_c^2,lon2=lon_c^2,latlon=lat_c*lon_c) %>%
  group_by(cluster) %>% mutate(pv_district_mean_100=mean(pv_100),
                               pv_within_100=pv_100-pv_district_mean_100,
                               district_n=n(),district_pv_sd=sd(pv)) %>% ungroup()
# FE estimation: verify coefficients reproduce Phase 5 before extending analysis.
f_region <- feols(rwi ~ log_grid+pv_100+log_road+log_pop | region,data=d)
f_district <- feols(rwi ~ log_grid+pv_100+log_road+log_pop | cluster,data=d)
f_district_trend <- feols(rwi ~ log_grid+pv_100+log_road+log_pop +
                           lat_c+lon_c+lat2+lon2+latlon | cluster,data=d)
check <- tibble(model=c('region FE','district FE','district FE + spatial trend'),
                solar_estimate=c(coef(f_region)['pv_100'],coef(f_district)['pv_100'],
                                 coef(f_district_trend)['pv_100']))
write_csv(check,file.path(output,'00_reproduction_check.csv'))
print(check)
if(max(abs(check$solar_estimate-c(.0197,.0881,.111)))>.008) {
 warning('Coefficients differ materially from Phase 5; verify sample, formula, and file version.')
}
# Conley covariance: geographic-distance-based HAC, 25/50/100 km triangular kernel.
# Computational burden rises with cutoff; allow a large dataset sufficient time.
key <- c('log_grid','pv_100','log_road','log_pop')
models <- list('Region FE'=f_region,'District FE'=f_district,
               'District FE + trend'=f_district_trend)
extract_rows <- function(model_name,fit,method,V,df=Inf) {
  b <- coef(fit); keep <- intersect(key, names(b));
  s <- sqrt(pmax(0,diag(V)[keep])); est <- unname(b[keep]); se <- unname(s)
  crit <- if(is.infinite(df)) qnorm(.975) else qt(.975,df)
  p <- if(is.infinite(df)) 2*pnorm(-abs(est/se)) else 2*pt(-abs(est/se),df)
  tibble(model=model_name,method=method,term=keep,estimate=est,se=se,
         lower95=est-crit*se,upper95=est+crit*se,p_value=p)
}
inference <- list()
for(name in names(models)) {
  fit <- models[[name]]
  G <- n_distinct(d$cluster)
  Vcl <- vcov(fit,~cluster)
  inference[[length(inference)+1]] <- extract_rows(name,fit,'District clustered',Vcl,df=G-1)
  for(km in c(25,50,100)) {
    cat('Computing Conley',km,'km for',name,'...\n')
    V <- tryCatch(vcov(fit,fixest::vcov_conley(lat='lat',lon='lon',cutoff=km)),
                  error=function(e){warning(conditionMessage(e));NULL})
    if(!is.null(V)) inference[[length(inference)+1]] <-
      extract_rows(name,fit,paste0('Conley ',km,' km'),V)
  }
}
results <- bind_rows(inference)
write_csv(results,file.path(output,'01_spatial_HAC_comparison.csv'))
# Descriptive within/between decomposition. District-average solar terms are not causal.
# This decomposition helps show why district-FE and region-FE coefficients differ.
fb <- feols(rwi ~ log_grid+pv_within_100+pv_district_mean_100+log_road+log_pop | region,data=d)
G <- n_distinct(d$cluster)
Vb <- vcov(fb,~cluster)
b <- coef(fb); ix <- intersect(c('pv_within_100','pv_district_mean_100'),names(b))
se <- sqrt(diag(Vb)[ix]); dfr <- G-1
between <- tibble(term=ix,estimate=unname(b[ix]),se=unname(se),
  lower95=unname(b[ix]-qt(.975,dfr)*se),upper95=unname(b[ix]+qt(.975,dfr)*se),
  p_value=unname(2*pt(-abs(b[ix]/se),df=dfr)))
write_csv(between,file.path(output,'02_within_between_decomposition.csv'))
# Rank exposure to between/within contributions; weight only by observed settlements.
district_data <- d %>% group_by(region,cluster) %>% summarise(n=n(),
   sd_pv_kwh_kwp=sd(pv),mean_pv_kwh_kwp=mean(pv),
   mean_rwi=mean(rwi),.groups='drop')
write_csv(district_data,file.path(output,'03_district_solar_variation.csv'))
# Spatial HAC figure: no guarantee Conley removes geographic confounding.
p <- results %>% filter(term=='pv_100') %>%
 ggplot(aes(x=estimate,y=model,color=method)) +
 geom_vline(xintercept=0,linetype='dashed',color='grey50') +
 geom_errorbarh(aes(xmin=lower95,xmax=upper95),
                position=position_dodge(width=.6),height=.12) +
 geom_point(position=position_dodge(width=.6),size=2) +
 labs(title='Ghana solar–wealth association: spatial HAC sensitivity',
      subtitle='Modeled wealth index; 95% uncertainty intervals by geographic error structure',
      x='Coefficient per 100 kWh/kWp annual modeled PV potential',y=NULL,
      caption='Cross-sectional association only; Conley bandwidths are sensitivity choices') +
 theme_minimal(base_size=12)
ggsave(file.path(output,'04_solar_conley_inference.png'),p,width=12,height=6,dpi=300)
cat('\nPHASE 6 COMPLETE:',output,'\n')
cat('\nREPRODUCTION CHECK:\n');print(check,width=Inf)
cat('\nSPATIAL HAC SOLAR RESULTS:\n');print(results %>% filter(term=='pv_100'),n=Inf,width=Inf)
cat('\nWITHIN-BETWEEN DECOMPOSITION:\n');print(between,width=Inf)
cat('\nCAUTIONS: standard-error corrections do not resolve endogeneity or timing mismatches.\n')