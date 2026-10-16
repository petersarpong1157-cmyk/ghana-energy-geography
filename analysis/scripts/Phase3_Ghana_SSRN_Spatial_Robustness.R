# PHASE 3 — Ghana SSRN: spatial dependence and missing-wealth diagnostics
# Reproduces the geographically restricted Phase 2C model. No causal claims.
# Does not overwrite original data or prior phase outputs.
packages <- c('readr','dplyr','tibble','sandwich','FNN','ggplot2')
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing)) install.packages(missing,repos='https://cloud.r-project.org')
suppressPackageStartupMessages({library(readr);library(dplyr);library(tibble);library(ggplot2)})

home <- Sys.getenv('USERPROFILE')
candidates <- c(file.path(home,'OneDrive','Desktop','PETER AND KONADU SSRN'),
 file.path(home,'Desktop','PETER AND KONADU SSRN'),
 file.path(home,'OneDrive - Personal','Desktop','PETER AND KONADU SSRN'))
folder <- candidates[dir.exists(candidates)][1]
if(is.na(folder)) stop('Research folder not found: specify folder manually.')
out <- file.path(folder,'Phase3_Outputs')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
raw <- read_csv(file.path(folder,'ghana_dre_atlas_settlements.csv'),show_col_types=FALSE,progress=FALSE)
need <- c('geohash','admin_cgaz_1','admin_cgaz_2','mean_rwi','pv_value',
          'distance_to_existing_transmission_lines','dist_main_road_km','population','lat','lon')
stopifnot(all(need %in% names(raw)))
ghana_regions <- c('Ahafo Region','Ashanti Region','Bono Region','Bono East Region',
'Central Region','Eastern Region','Greater Accra Region','North East Region',
'Northern Region','Oti Region','Savannah Region','Upper East Region',
'Upper West Region','Volta Region','Western Region','Western North Region')
d <- raw %>% transmute(geohash=as.character(geohash),
 region=na_if(trimws(as.character(admin_cgaz_1)),''),
 district=na_if(trimws(as.character(admin_cgaz_2)),''),
 rwi=as.numeric(mean_rwi),pv=as.numeric(pv_value),
 grid=as.numeric(distance_to_existing_transmission_lines),
 road=as.numeric(dist_main_road_km),pop=as.numeric(population),
 lat=as.numeric(lat),lon=as.numeric(lon)) %>%
 mutate(log_grid=log1p(grid),log_road=log1p(road),
 log_pop=log1p(pop),pv_100=pv/100,
 cluster=paste(region,district,sep=' | '),
 missing_wealth=as.integer(is.na(rwi)),ghana_region=region %in% ghana_regions)
# Retain only valid covariates and recognized Ghana regions; no imputed wealth.
valid_x <- d %>% filter(ghana_region,!is.na(district),is.finite(pv_100),
 is.finite(log_grid),is.finite(log_road),is.finite(log_pop),
 is.finite(lat),is.finite(lon))
complete <- valid_x %>% filter(is.finite(rwi)) %>% mutate(region=factor(region))
f <- rwi ~ log_grid + pv_100 + log_road + log_pop + region
fit <- lm(f,data=complete)
lev <- hatvalues(fit)
if(max(lev) > 0.99999999) stop('Near-perfect leverage remains; investigate before inference.')
cl_v <- sandwich::vcovCL(fit,cluster=complete$cluster,type='HC1')
hc3_v <- sandwich::vcovHC(fit,type='HC3')
coef_tbl <- function(V,method){
 b <- coef(fit); se <- sqrt(diag(V))
 df <- if(method=='district_clustered') dplyr::n_distinct(complete$cluster)-1 else Inf
 stat <- b/se
 p <- if(is.finite(df)) 2*pt(abs(stat),df=df,lower.tail=FALSE) else 2*pnorm(abs(stat),lower.tail=FALSE)
 critical <- if(is.finite(df)) qt(.975,df=df) else qnorm(.975)
 tibble(method=method,term=names(b),estimate=unname(b),se=unname(se),
 lower_95=unname(b-critical*se),upper_95=unname(b+critical*se),p_value=unname(p))
}
coefs <- bind_rows(coef_tbl(hc3_v,'HC3'),coef_tbl(cl_v,'district_clustered'))
write_csv(coefs,file.path(out,'01_inference_with_95pct_CIs.csv'))

# Spatial Moran I: 8 nearest neighbors, row-standardized directed kNN weights.
# A degree-based flat mapping is approximated in km; note this is not survey design weighting.
# On observed point coordinates, kNN is invariant to varying settlement density.
set.seed(20261016)
coords <- cbind(complete$lon * 111.32*cos(mean(complete$lat)*pi/180),complete$lat * 110.57)
if(anyDuplicated(data.frame(x=coords[,1],y=coords[,2]))) {
  message('WARNING: Identical rounded point coordinates: check neighbor diagnostics.')
}
k <- 8L
nn <- FNN::get.knn(coords,k=k)$nn.index
res <- as.numeric(residuals(fit))
moran <- function(z){
 z <- z-mean(z)
 n <- length(z)
 # S0=n for row-standardized weights; Moran I=sum(z_i mean_neighbor_z)/sum(z_i^2)
 sum(z * rowMeans(matrix(z[as.vector(t(nn))],nrow=n,byrow=TRUE)))/sum(z*z)
}
# nn matrix is row-wise; as.vector(t(nn)) explicitly makes row-by-row sequence.
observed_I <- moran(res)
B <- 199L
permuted <- replicate(B,moran(sample(res,length(res),replace=FALSE)))
p_upper <- (1+sum(permuted>=observed_I))/(B+1)
p_two_sided <- (1+sum(abs(permuted)>=abs(observed_I)))/(B+1)
moran_summary <- tibble(n=nrow(complete),k_neighbors=k,permutations=B,
 moran_I=observed_I,mean_null_I=mean(permuted),
 p_upper_tailed=p_upper,p_two_sided=p_two_sided,
 method='Row-standardized directed 8-nearest-neighbor Moran I; Monte Carlo permutations')
write_csv(moran_summary,file.path(out,'02_residual_moran_test.csv'))

# Missingness: comparisons use all valid Ghana-region predictor rows, including missing RWI.
missing_summary <- valid_x %>% group_by(missing_wealth) %>% summarise(
 n=n(),median_population=median(pop),median_grid_km=median(grid),
 median_road_km=median(road),median_pv=median(pv),.groups='drop')
write_csv(missing_summary,file.path(out,'03_missing_vs_observed_wealth.csv'))
# Model probability of wealth being missing, rather than pretending missing is random.
# This assesses observed covariate dependence only, not missing-at-random validity.
miss_fit <- glm(missing_wealth ~ log_grid+pv_100+log_road+log_pop+factor(region),
 data=valid_x,family=binomial())
miss_coef <- summary(miss_fit)$coefficients
miss_tbl <- tibble(term=rownames(miss_coef),log_odds=miss_coef[,1],
 se=miss_coef[,2],p_value=miss_coef[,4],odds_ratio=exp(miss_coef[,1]))
write_csv(miss_tbl,file.path(out,'04_missing_wealth_logistic_regression.csv'))

# Population tail sensitivity, retaining previously chosen geographical restriction.
q995 <- quantile(complete$pop,.995,na.rm=TRUE)
tail_data <- complete %>% filter(pop<=q995)
tail_fit <- lm(f,data=tail_data)
get_cluster <- function(model,data,label){
 V <- sandwich::vcovCL(model,cluster=data$cluster,type='HC1')
 b <- coef(model);se <- sqrt(diag(V)); G <- n_distinct(data$cluster)
 tibble(model=label,n=nrow(data),clusters=G,term=names(b),
 estimate=unname(b),se=unname(se),p_value=2*pt(abs(b/se),df=G-1,lower.tail=FALSE))
}
tail_comparison <- bind_rows(get_cluster(fit,complete,'main'),
 get_cluster(tail_fit,tail_data,'remove_top_0.5pct_population'))
write_csv(tail_comparison,file.path(out,'05_population_tail_sensitivity.csv'))

# Visual for Moran permuted reference distribution: observational diagnostic only.
plot_m <- ggplot(tibble(I=permuted),aes(I))+geom_histogram(bins=28,fill='#52789C',color='white')+
 geom_vline(xintercept=observed_I,color='#C05A43',linewidth=1)+
 labs(title='Spatial dependence in Ghana settlement regression residuals',
 subtitle=paste0('Observed Moran I = ',round(observed_I,4),
 '; 8-nearest-neighbor permutations (B = ',B,')'),
 x='Moran I under randomly permuted residuals',y='Number of permutations',
 caption='Red line: observed residual Moran I. A test of residual spatial patterns, not causality.')+
 theme_minimal(base_size=12)
ggsave(file.path(out,'06_residual_spatial_autocorrelation.png'),plot_m,
 width=10,height=6,dpi=240)
writeLines(c(
 'PHASE 3 Interpretation safeguards',
 'The population and wealth measures are modeled atlas estimates, not direct household surveys.',
 'Spatial association cannot establish causality or individual access to the electricity grid.',
 'Moran I permutation test assesses residual spatial dependence under point exchangeability.',
 'A significant Moran I means district clustering may not fully capture spatial residual structure.',
 'kNN uses an approximate lat/lon-to-kilometer transformation; rerun with projected CRS in a refined analysis.',
 'Wealth missingness logistic regression only tests dependence on measured variables.',
 'If missingness is predictable, complete-case estimates may be selection-sensitive.',
 'The two population observations >1 million still deserve direct source verification.',
 'Boundary shapefile choice District_261 versus District_272 remains unresolved.',
 'Do not automatically use model-based spatial corrections without examining their assumptions.'),
 file.path(out,'READ_ME_INTERPRETATION.txt'))
cat('\nPHASE 3 COMPLETE:',out,'\n')
cat('Regression n:',nobs(fit),'clusters:',n_distinct(complete$cluster),
 'max leverage:',max(lev),'\n')
cat('\nPRIMARY COEFFICIENTS (district-clustered):\n')
print(coefs %>% filter(method=='district_clustered',term %in% c('log_grid','pv_100','log_road','log_pop')),
 width=Inf)
cat('\nRESIDUAL MORAN TEST:\n');print(moran_summary,width=Inf)
cat('\nWEALTH MISSINGNESS:\n');print(missing_summary,width=Inf)
cat('\nPOPULATION TAIL SENSITIVITY:\n')
print(tail_comparison %>% filter(term %in% c('log_grid','pv_100','log_road','log_pop')),width=Inf)