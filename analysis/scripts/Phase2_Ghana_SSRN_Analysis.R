# PHASE 2 — Ghana SSRN spatial inequality working paper
# Descriptive patterns, audit and exploratory associations; NOT causal identification.
# Based on World Bank DRE Atlas settlement data and accompanying June 2025 data dictionary.
# Run the whole script in RStudio using Source.

packages <- c('readr','dplyr','tidyr','ggplot2','tibble','sandwich','lmtest')
needed <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
if(length(needed)) install.packages(needed, repos='https://cloud.r-project.org')
suppressPackageStartupMessages({library(readr); library(dplyr); library(tidyr); library(ggplot2); library(tibble)})

home <- Sys.getenv('USERPROFILE')
candidates <- c(file.path(home,'OneDrive','Desktop','PETER AND KONADU SSRN'),
                file.path(home,'Desktop','PETER AND KONADU SSRN'),
                file.path(home,'OneDrive - Personal','Desktop','PETER AND KONADU SSRN'))
folder <- candidates[dir.exists(candidates)][1]
if(is.na(folder)) stop('Folder missing. Set folder to the correct path manually.')
output <- file.path(folder,'Phase2_Outputs')
dir.create(output, recursive=TRUE,showWarnings=FALSE)
raw <- read_csv(file.path(folder,'ghana_dre_atlas_settlements.csv'),show_col_types=FALSE,progress=FALSE)
required <- c('geohash','admin_cgaz_1','admin_cgaz_2','mean_rwi','pv_value',
              'distance_to_existing_transmission_lines','dist_main_road_km',
              'population','has_nightlight','lat','lon')
stopifnot(all(required %in% names(raw)), nrow(raw)>0)

# Check supplied metadata, do not substitute guessed units.
cat('Source definitions: mean_rwi = modeled relative wealth index;\n',
    'pv_value = annual long-term PV potential (kWh/kWp, 1994–2015);\n',
    'distance_to_existing_transmission_lines = km;\n',
    'has_nightlight = VIIRS light-area overlap (Sep 2022–Aug 2023), NOT verified grid access.\n')

# Prepare observational analysis sample. Outliers remain until sensitivity analyses.
d <- raw %>% mutate(
 region=na_if(trimws(as.character(admin_cgaz_1)),''),
 district=na_if(trimws(as.character(admin_cgaz_2)),''),
 rwi=as.numeric(mean_rwi), pv=as.numeric(pv_value),
 grid_distance_km=as.numeric(distance_to_existing_transmission_lines),
 road_distance_km=as.numeric(dist_main_road_km),
 pop=as.numeric(population),
 nightlight=case_when(as.character(has_nightlight) %in% c('TRUE','True','true','1') ~ 1L,
                      as.character(has_nightlight) %in% c('FALSE','False','false','0') ~ 0L,
                      TRUE ~ NA_integer_),
 log_grid=log1p(grid_distance_km), log_road=log1p(road_distance_km),
 log_pop=log1p(pop),
 missing_rwi=is.na(rwi),
 region_known=!is.na(region),
 region_district=paste(coalesce(region,'Unknown'),coalesce(district,'Unknown'),sep=' | '))

# Check physical/semantic plausibility, do not silently exclude negative distances.
qc <- tibble(metric=c('total_records','missing_wealth','missing_pv','missing_region',
 'missing_grid_distance','negative_grid_distance','negative_road_distance',
 'invalid_population','missing_nightlight','population_above_one_million'),
 count=c(nrow(d),sum(is.na(d$rwi)),sum(is.na(d$pv)),sum(is.na(d$region)),
 sum(is.na(d$grid_distance_km)),sum(d$grid_distance_km<0,na.rm=TRUE),
 sum(d$road_distance_km<0,na.rm=TRUE),sum(is.na(d$pop)|d$pop<0),
 sum(is.na(d$nightlight)),sum(d$pop>1e6,na.rm=TRUE)))
write_csv(qc,file.path(output,'01_quality_and_sample_counts.csv'))
if(any(d$grid_distance_km<0,na.rm=TRUE) || any(d$road_distance_km<0,na.rm=TRUE))
  stop('Negative distance values: inspect 01_quality_and_sample_counts.csv first.')

# Region summary: unweighted settlement means/medians. Not population-level estimates.
region_summary <- d %>% group_by(region) %>% summarise(
 settlements=n(),missing_wealth=sum(is.na(rwi)),
 missing_wealth_pct=round(100*mean(is.na(rwi)),2),
 median_wealth=median(rwi,na.rm=TRUE),
 median_pv_kwh_kwp=median(pv,na.rm=TRUE),
 median_grid_distance_km=median(grid_distance_km,na.rm=TRUE),
 share_with_nightlight=mean(nightlight,na.rm=TRUE),.groups='drop') %>%
 mutate(across(where(is.numeric),~replace(.x,is.nan(.x),NA_real_))) %>% arrange(region)
write_csv(region_summary,file.path(output,'02_region_descriptives.csv'))

# Global quintiles of observed wealth, settlement-weighted and not representative of households.
qcut <- quantile(d$rwi,seq(0,1,0.2),na.rm=TRUE)
# RWI can contain ties: use ranks to prevent cut() from failing on non-unique quantiles.
observed <- d %>% filter(!is.na(rwi)) %>% arrange(rwi, geohash) %>%
 mutate(wealth_quintile=ntile(rwi,5),
        wealth_group=paste0('Q',wealth_quintile))
qsum <- observed %>% group_by(wealth_quintile,wealth_group) %>% summarise(
 settlements=n(), median_wealth=median(rwi),
 median_pv_kwh_kwp=median(pv,na.rm=TRUE),
 median_grid_distance_km=median(grid_distance_km,na.rm=TRUE),
 median_road_distance_km=median(road_distance_km,na.rm=TRUE),
 nightlight_share=mean(nightlight,na.rm=TRUE),.groups='drop')
write_csv(qsum,file.path(output,'03_wealth_quintile_patterns.csv'))
write_csv(tibble(quantile=c(0,.2,.4,.6,.8,1),value=as.numeric(qcut)),
          file.path(output,'04_wealth_quantile_cutpoints.csv'))

# Missingness is substantively important: compare records WITH and WITHOUT wealth.
missingness <- d %>% group_by(missing_rwi) %>% summarise(
 settlements=n(),median_population=median(pop,na.rm=TRUE),
 median_pv=median(pv,na.rm=TRUE),
 median_grid_distance_km=median(grid_distance_km,na.rm=TRUE),
 share_with_nightlight=mean(nightlight,na.rm=TRUE),.groups='drop')
write_csv(missingness,file.path(output,'05_wealth_missingness_comparison.csv'))

# Clear publication-readable graphics; titles do not claim effects or causal relationships.
p1 <- ggplot(qsum,aes(factor(wealth_quintile),median_grid_distance_km))+
 geom_col(fill='#22577A',width=.7)+
 geom_text(aes(label=sprintf('%.1f km',median_grid_distance_km)),vjust=-.35,size=3.6)+
 scale_y_continuous(expand=expansion(mult=c(0,.18)))+
 labs(title='Transmission-line proximity across relative-wealth groups',
 subtitle='Settlement-level medians; Q1 = lowest observed relative wealth',
 x='Relative wealth quintile',y='Median distance to existing transmission line (km)',
 caption='World Bank DRE Atlas; descriptive association, not electricity-access or causal effect.')+
 theme_minimal(base_size=12)
ggsave(file.path(output,'06_grid_distance_by_wealth_quintile.png'),p1,width=10,height=6,dpi=240)

p2 <- ggplot(qsum,aes(factor(wealth_quintile),median_pv_kwh_kwp))+
 geom_col(fill='#4D7C48',width=.7)+
 geom_text(aes(label=sprintf('%.0f',median_pv_kwh_kwp)),vjust=-.4,size=3.5)+
 scale_y_continuous(limits=c(0,NA),expand=expansion(mult=c(0,.14)))+
 labs(title='Modeled annual PV potential by relative-wealth group',
 subtitle='Long-term resource potential, not installed solar capacity or electricity generated',
 x='Relative wealth quintile',y='Median modeled PV potential (kWh/kWp/year)',
 caption='World Bank DRE Atlas; source PV averaging period: 1994–2015.')+
 theme_minimal(base_size=12)
ggsave(file.path(output,'07_solar_potential_by_wealth_quintile.png'),p2,width=10,height=6,dpi=240)

p3 <- ggplot(d,aes(x=lon,y=lat,color=grid_distance_km))+
 geom_point(size=.32,alpha=.72)+coord_equal()+
 scale_color_viridis_c(trans='sqrt',name='Grid distance (km)',na.value='grey70')+
 labs(title='Settlement geography and transmission-line proximity',
 subtitle='Point locations only; not a map of electrification or transmission ownership',
 x='Longitude',y='Latitude',caption='Source: World Bank DRE Atlas')+
 theme_minimal(base_size=11)
ggsave(file.path(output,'08_settlement_transmission_proximity_map.png'),p3,width=9,height=9,dpi=230)

# Regression: associations only, no spatial or causal interpretation.
# Region fixed effects account for broad between-region differences;
# district-clustered SE addresses dependence within named districts, but not all spatial autocorrelation.
model_data <- d %>% filter(is.finite(rwi),is.finite(pv),is.finite(log_grid),
 is.finite(log_road),is.finite(log_pop),!is.na(region),!is.na(district)) %>%
 mutate(region=factor(region))
write_csv(tibble(metric=c('model_complete_cases','excluded_from_model','distinct_district_labels',
 'n_geographic_clusters','population_p99_5'),
 value=c(nrow(model_data),nrow(d)-nrow(model_data),n_distinct(model_data$district),
 n_distinct(model_data$region_district),as.numeric(quantile(model_data$pop,.995,na.rm=TRUE)))),
 file.path(output,'09_model_sample.csv'))
if(nrow(model_data)<100 || n_distinct(model_data$region)<2) stop('Insufficient analysis sample.')
model_data <- model_data %>% mutate(pv_100=pv/100)
formula_main <- rwi ~ log_grid + pv_100 + log_road + log_pop + region
fit <- lm(formula_main,data=model_data)
coef_table <- function(model, label, clusters) {
  v_hc3 <- sandwich::vcovHC(model,type='HC3')
  hc <- lmtest::coeftest(model,vcov.=v_hc3)
  v_cl <- sandwich::vcovCL(model,cluster=clusters,type='HC1')
  cl <- lmtest::coeftest(model,vcov.=v_cl)
  tibble(model=label,term=rownames(hc),estimate=unname(hc[,1]),
         hc3_se=unname(hc[,2]),hc3_p=unname(hc[,4]),
         district_cluster_se=unname(cl[,2]),district_cluster_p=unname(cl[,4]))
}
base_results <- coef_table(fit,'main_complete_case',model_data$region_district)

# Sensitivity: remove extreme population tail for model stability (not as invalid data).
threshold <- quantile(model_data$pop,.995,na.rm=TRUE)
sens <- model_data %>% filter(pop<=threshold)
fit_sens <- lm(formula_main,data=sens)
sensitivity_results <- coef_table(fit_sens,'sensitivity_exclude_top_0.5_percent_population',
                                  sens$region_district)
write_csv(bind_rows(base_results,sensitivity_results),
          file.path(output,'10_exploratory_regression_coefficients.csv'))
write_csv(tibble(model=c('main','population_tail_sensitivity'),
 n=c(nobs(fit),nobs(fit_sens)),r_squared=c(summary(fit)$r.squared,summary(fit_sens)$r.squared),
 cluster_count=c(n_distinct(model_data$region_district),n_distinct(sens$region_district))),
 file.path(output,'11_regression_diagnostics.csv'))

writeLines(c(
 'PHASE 2 — INTERPRETATION GUIDE',
 'This is exploratory cross-sectional analysis; no causal effects are identified.',
 'Relative wealth index is modeled, NOT direct household income or expenditure.',
 'PV potential represents modeled long-run potential in kWh/kWp per year (1994–2015), not actual solar output.',
 'Grid distance is distance to existing transmission lines (unknown or mixed voltages); not proof of household grid connection.',
 'Night-light overlap is based on VIIRS Sep 2022–Aug 2023; not proof of access to electricity.',
 'Source variables have different observation/modeling periods; no longitudinal interpretation.',
 'Region and district labels must be verified against geographical boundaries before spatial inference.',
 'Regression standard errors use both HC3 and clusters built from region + district names.',
 'Clustering does not replace spatial residual dependence diagnostics (needed in later phase).',
 'Quintiles use rank-based grouping to accommodate RWI ties, and exclude missing wealth.',
 'Settlement summaries weight each settlement equally; they are not representative household estimates.',
 'RWI missingness must be assessed for selection bias; missingness is not treated as random.',
 'Review Phase1_Outputs/06_boundary_comparison.csv before picking District_261 versus District_272.',
 'Review regression results and outlier sensitivity before making substantive claims.'),
 file.path(output,'README_READ_BEFORE_INTERPRETATION.txt'))
cat('\nPHASE 2 COMPLETE. Outputs: ',output,'\n',sep='')
print(qc);print(qsum);print(missingness)
cat('\nModel sample n = ',nrow(model_data),'; number of district clusters = ',
 n_distinct(model_data$region_district),'\n',sep='')
cat('\nPrimary coefficients (associations only):\n')
print(base_results %>% filter(term %in% c('log_grid','pv_100','log_road','log_pop')))