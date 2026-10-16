# PHASE 4: Ghana SSRN — spatial-block and missing-RWI sensitivity
# Exploratory cross-sectional associations; no causal interpretations.
# This script never overwrites original files or previous Phase outputs.
packages <- c('readr','dplyr','tibble','sandwich','ggplot2','FNN')
to_install <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
if(length(to_install)) install.packages(to_install, repos='https://cloud.r-project.org')
suppressPackageStartupMessages({library(readr);library(dplyr);library(tibble);library(ggplot2)})

home <- Sys.getenv('USERPROFILE')
options_dir <- c(file.path(home,'OneDrive','Desktop','PETER AND KONADU SSRN'),
                 file.path(home,'Desktop','PETER AND KONADU SSRN'),
                 file.path(home,'OneDrive - Personal','Desktop','PETER AND KONADU SSRN'))
folder <- options_dir[dir.exists(options_dir)][1]
if(is.na(folder)) stop('Folder not found. Edit the folder variable in this script.')
output <- file.path(folder,'Phase4_Outputs')
dir.create(output, recursive=TRUE, showWarnings=FALSE)
raw <- read_csv(file.path(folder,'ghana_dre_atlas_settlements.csv'),
                show_col_types=FALSE, progress=FALSE)
req <- c('geohash','admin_cgaz_1','admin_cgaz_2','mean_rwi','pv_value',
         'distance_to_existing_transmission_lines','dist_main_road_km','population','lat','lon')
if(!all(req %in% names(raw))) stop('Required column(s) missing: ',paste(setdiff(req,names(raw)),collapse=', '))
ghana_regions <- c('Ahafo Region','Ashanti Region','Bono Region','Bono East Region',
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
 log_grid=log1p(grid),log_road=log1p(road),log_pop=log1p(pop),
 pv_100=pv/100,cluster=paste(region,district,sep=' | '))
valid <- d %>% filter(region %in% ghana_regions,!is.na(district),
 is.finite(pv_100),is.finite(log_grid),is.finite(log_road),
 is.finite(log_pop),is.finite(lat),is.finite(lon))
complete <- valid %>% filter(is.finite(rwi))
# Approximate rectangular geographic blocks, deliberately independent of district labels.
# 0.5 and 1 degree are spatial robustness scales, not precise 50-/100-km circles.
block_id <- function(x,cell,shift=0) {
 paste(floor((x$lon+shift)/cell),floor((x$lat+shift)/cell),sep=':')
}
complete <- complete %>% mutate(block_05=block_id(.,0.5),block_05_shift=block_id(.,0.5,0.25),
                                block_1=block_id(.,1),block_1_shift=block_id(.,1,0.5))

base_formula <- rwi ~ log_grid + pv_100 + log_road + log_pop + factor(region)
district_formula <- rwi ~ log_grid + pv_100 + log_road + log_pop + factor(cluster)
main_fit <- lm(base_formula,data=complete)
# District fixed effects: compare within-district relationships rather than only region adjustments.
# R's lm may drop predictors if collinear with fixed effects; flag that in results.
district_fit <- lm(district_formula,data=complete)

extract <- function(fit,cl,method,model) {
  if(anyNA(coef(fit))) message('Aliased terms in ',model,'; inspect coefficients.')
  V <- sandwich::vcovCL(fit,cluster=cl,type='HC1',cadjust=TRUE)
  b <- coef(fit); b <- b[!is.na(b)]
  se <- sqrt(diag(V)[names(b)]); df <- length(unique(cl))-1
  if(df <= 1) stop('Insufficient clusters: ',method)
  tibble(model=model,method=method,n=nobs(fit),clusters=length(unique(cl)),
    term=names(b),estimate=as.numeric(b),se=as.numeric(se),
    lower_95=as.numeric(b-qt(.975,df)*se),upper_95=as.numeric(b+qt(.975,df)*se),
    p_value=as.numeric(2*pt(abs(b/se),df=df,lower.tail=FALSE)))
}
results <- dplyr::bind_rows(
 extract(main_fit,complete$cluster,'district clustered','region FE'),
 extract(main_fit,complete$block_05,'0.5-degree spatial blocks','region FE'),
 extract(main_fit,complete$block_05_shift,'shifted 0.5-degree blocks','region FE'),
 extract(main_fit,complete$block_1,'1-degree spatial blocks','region FE'),
 extract(main_fit,complete$block_1_shift,'shifted 1-degree blocks','region FE'),
 extract(district_fit,complete$cluster,'district clustered','district FE'),
 extract(district_fit,complete$block_05,'0.5-degree spatial blocks','district FE'),
 extract(district_fit,complete$block_1,'1-degree spatial blocks','district FE'))
write_csv(results,file.path(output,'01_spatial_block_inference.csv'))

# Moran I before and after adding district FE. Neighbor weights: directed row-standardized 8NN.
# Coordinates are projected approximately to km using sample mean latitude.
lat0 <- mean(complete$lat)
coords <- cbind(complete$lon *111.32*cos(lat0*pi/180),complete$lat*110.57)
if(anyDuplicated(data.frame(coords))) message('WARNING: coincident coordinates; inspect spatial-neighbor structure.')
nearest <- FNN::get.knn(coords,k=8)$nn.index
moran <- function(z) {
 z <- z-mean(z)
 sum(z*rowMeans(matrix(z[as.vector(t(nearest))],nrow=length(z),byrow=TRUE)))/sum(z^2)
}
set.seed(20261016)
B <- 999L
moran_test <- function(z,name){
 obs <- moran(z)
 perm <- replicate(B,moran(sample(z,length(z),replace=FALSE)))
 tibble(model=name,n=length(z),moran_I=obs,permutations=B,
 p_upper=(1+sum(perm>=obs))/(B+1),
 p_two_sided=(1+sum(abs(perm)>=abs(obs)))/(B+1))
}
moran_results <- bind_rows(moran_test(residuals(main_fit),'region FE'),
                            moran_test(residuals(district_fit),'district FE'))
write_csv(moran_results,file.path(output,'02_residual_moran_comparison.csv'))

# Missing-RWI sensitivity: inverse-probability weighting assumes missingness at random
# conditional on observed predictors. This assumption is untestable here.
# Weights are stabilized and capped to reduce influence, NOT a solution for MNAR.
valid <- valid %>% mutate(is_observed=as.integer(is.finite(rwi)))
obs_fit <- glm(is_observed ~ log_grid + pv_100 + log_road + log_pop + factor(region),
              data=valid,family=binomial())
valid$p_obs <- predict(obs_fit,type='response')
# Also save the observations with smallest fitted observation probability.
write_csv(valid %>% arrange(p_obs) %>%
 select(geohash,region,district,p_obs,is_observed,grid,road,pop,pv) %>% head(100),
 file.path(output,'03_low_observation_probability.csv'))
weighted <- valid %>% filter(is_observed==1) %>%
 mutate(raw_weight=mean(valid$is_observed)/pmax(p_obs,0.05),
        ipw=pmin(raw_weight,quantile(raw_weight,0.99,na.rm=TRUE)))
weighted_fit <- lm(base_formula,data=weighted,weights=ipw)
miss_results <- bind_rows(
 extract(main_fit,complete$cluster,'complete-case district cluster','unweighted'),
 extract(weighted_fit,weighted$cluster,'IPW district cluster','MAR sensitivity'))
write_csv(miss_results,file.path(output,'04_missingness_IPW_sensitivity.csv'))
write_csv(tibble(n_total=nrow(valid),n_observed=nrow(weighted),
 min_p_obs=min(valid$p_obs),p01_p_obs=as.numeric(quantile(valid$p_obs,.01)),
 max_raw_weight=max(weighted$raw_weight),max_capped_weight=max(weighted$ipw),
 effective_weighted_n=(sum(weighted$ipw)^2)/sum(weighted$ipw^2)),
 file.path(output,'05_IPW_diagnostics.csv'))

key <- c('log_grid','pv_100','log_road','log_pop')
p <- results %>% filter(term %in% key) %>%
 mutate(term=factor(term,levels=key)) %>%
 ggplot(aes(x=estimate,y=method,color=model)) +
 geom_vline(xintercept=0,linetype='dashed',color='grey55')+
 geom_errorbarh(aes(xmin=lower_95,xmax=upper_95),height=.15,
               position=position_dodge(width=.5))+
 geom_point(position=position_dodge(width=.5),size=1.8)+
 facet_wrap(~term,scales='free_x')+
 labs(title='Sensitivity of wealth associations to spatial inference',
      subtitle='Region and district fixed effects; intervals vary by clustering scheme',
      x='Regression coefficient and 95% confidence interval',y=NULL,color='Model')+
 theme_minimal(base_size=11)
ggsave(file.path(output,'06_spatial_inference_comparison.png'),p,width=13,height=8,dpi=300)

cat('\nPHASE 4 COMPLETE: ',output,'\n',sep='')
cat('\nSPATIAL INFERENCE — MAIN VARIABLES:\n')
print(results %>% filter(term %in% key) %>%
      select(model,method,clusters,term,estimate,se,p_value),n=Inf,width=Inf)
cat('\nMORAN COMPARISON:\n');print(moran_results,width=Inf)
cat('\nMISSINGNESS IPW SENSITIVITY:\n')
print(miss_results %>% filter(term %in% key) %>%
      select(model,method,term,estimate,se,p_value),n=Inf,width=Inf)
cat('\nIPW DIAGNOSTICS:\n')
print(read_csv(file.path(output,'05_IPW_diagnostics.csv'),show_col_types=FALSE),width=Inf)
cat('\nCAUTION: geographic blocks are an approximation, and IPW needs MAR/positivity.\n')
cat('Check exact geospatial polygon membership for flagged border records before final manuscript.\n')