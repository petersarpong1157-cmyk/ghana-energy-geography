# PHASE 5: Ghana SSRN — coefficient stability and spatial trend audit
# Inputs: original World Bank DRE settlements CSV. Nothing is overwritten.
# This phase tests sensitivity; it cannot establish causality or remove all spatial confounding.
packages <- c('readr','dplyr','tibble','ggplot2','sandwich')
new <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(new)) install.packages(new,repos='https://cloud.r-project.org')
suppressPackageStartupMessages({library(readr);library(dplyr);library(tibble);library(ggplot2)})

home <- Sys.getenv('USERPROFILE')
possible <- c(file.path(home,'OneDrive','Desktop','PETER AND KONADU SSRN'),
              file.path(home,'Desktop','PETER AND KONADU SSRN'),
              file.path(home,'OneDrive - Personal','Desktop','PETER AND KONADU SSRN'))
folder <- possible[dir.exists(possible)][1]
if(is.na(folder)) stop('Research folder not found. Set folder manually.')
output <- file.path(folder,'Phase5_Outputs')
dir.create(output,recursive=TRUE,showWarnings=FALSE)
raw <- read_csv(file.path(folder,'ghana_dre_atlas_settlements.csv'),show_col_types=FALSE,progress=FALSE)
need <- c('geohash','admin_cgaz_1','admin_cgaz_2','mean_rwi','pv_value',
          'distance_to_existing_transmission_lines','dist_main_road_km','population','lat','lon')
if(!all(need %in% names(raw))) stop('Missing columns: ',paste(setdiff(need,names(raw)),collapse=', '))
regions <- c('Ahafo Region','Ashanti Region','Bono Region','Bono East Region',
             'Central Region','Eastern Region','Greater Accra Region','North East Region',
             'Northern Region','Oti Region','Savannah Region','Upper East Region',
             'Upper West Region','Volta Region','Western Region','Western North Region')
d <- raw %>% transmute(geohash=as.character(geohash),
    region=trimws(as.character(admin_cgaz_1)),district=trimws(as.character(admin_cgaz_2)),
    rwi=suppressWarnings(as.numeric(mean_rwi)),pv=suppressWarnings(as.numeric(pv_value)),
    grid=suppressWarnings(as.numeric(distance_to_existing_transmission_lines)),
    road=suppressWarnings(as.numeric(dist_main_road_km)),
    pop=suppressWarnings(as.numeric(population)),
    lat=suppressWarnings(as.numeric(lat)),lon=suppressWarnings(as.numeric(lon))) %>%
  mutate(region=na_if(region,''),district=na_if(district,''),
         cluster=paste(region,district,sep=' | '),
         pv_100=pv/100,log_grid=log1p(grid),log_road=log1p(road),log_pop=log1p(pop)) %>%
  filter(region %in% regions,!is.na(district),is.finite(rwi),is.finite(pv_100),
         is.finite(log_grid),is.finite(log_road),is.finite(log_pop),
         is.finite(lat),is.finite(lon),grid>=0,road>=0,pop>=0)
stopifnot(nrow(d)>1000,all(!duplicated(d$geohash)))
d <- d %>% mutate(block05=paste(floor(lon/0.5),floor(lat/0.5),sep=':'),
                  block1=paste(floor(lon),floor(lat),sep=':'))
# Geographic trend coordinates centered to ease polynomial interpretation.
d <- d %>% mutate(lat_c=lat-mean(lat),lon_c=lon-mean(lon),
                  lat2=lat_c^2,lon2=lon_c^2,latlon=lat_c*lon_c)

f_region <- rwi ~ log_grid + pv_100 + log_road + log_pop + factor(region)
f_region_geo <- update(f_region,. ~ . + lat_c + lon_c + lat2 + lon2 + latlon)
f_district <- rwi ~ log_grid + pv_100 + log_road + log_pop + factor(cluster)
f_district_geo <- update(f_district,. ~ . + lat_c + lon_c + lat2 + lon2 + latlon)

get_stats <- function(fit,dat,model,scheme='district') {
  cl <- switch(scheme,district=dat$cluster,block05=dat$block05,block1=dat$block1)
  G <- length(unique(cl))
  if(G < 10) stop('Too few clusters: ',G)
  V <- sandwich::vcovCL(fit,cluster=cl,type='HC1',cadjust=TRUE)
  b <- coef(fit); b <- b[!is.na(b)]
  keep <- intersect(c('log_grid','pv_100','log_road','log_pop'),names(b))
  s <- sqrt(diag(V)[keep]); b <- b[keep]; df <- G-1
  tibble(model=model,inference=scheme,n=nobs(fit),clusters=G,term=keep,
         estimate=unname(b),se=unname(s),
         lower95=unname(b-qt(.975,df)*s),upper95=unname(b+qt(.975,df)*s),
         p_value=unname(2*pt(abs(b/s),df=df,lower.tail=FALSE)))
}
models <- list('Region FE'=f_region,'Region FE + spatial trend'=f_region_geo,
               'District FE'=f_district,'District FE + spatial trend'=f_district_geo)
fits <- lapply(models,function(f) lm(f,data=d))
results <- dplyr::bind_rows(lapply(names(fits),function(name){
  fit <- fits[[name]]
  bind_rows(lapply(c('district','block05','block1'),function(s) get_stats(fit,d,name,s)))
}))
write_csv(results,file.path(output,'01_specification_comparison.csv'))

# Does the district-FE coefficient rely on districts with little local PV variation?
within <- d %>% group_by(cluster,region) %>% summarise(n=n(),
  pv_sd=sd(pv,na.rm=TRUE),pv_min=min(pv),pv_max=max(pv),
  pv_range=pv_max-pv_min,.groups='drop') %>% arrange(pv_sd)
write_csv(within,file.path(output,'02_district_within_pv_variation.csv'))
# Sensitivity to dropping one district with unusually low PV variation (descriptive only).
threshold <- 10 # kWh/kWp: prespecified diagnostic cutoff, not a scientifically validated threshold
eligible <- within %>% filter(n>=10,is.finite(pv_sd),pv_sd>=threshold) %>% pull(cluster)
restricted <- d %>% filter(cluster %in% eligible)
restricted_results <- if(length(unique(restricted$cluster))>=20 && nrow(restricted)>1000){
  f <- lm(f_district,data=restricted)
  bind_rows(get_stats(f,restricted,'District FE; districts n>=10 & PV SD>=10','district'),
            get_stats(f,restricted,'District FE; districts n>=10 & PV SD>=10','block05'))
} else tibble()
write_csv(restricted_results,file.path(output,'03_within_district_variation_sensitivity.csv'))

# Leave-one-region-out (LOO) refits. Track how each excluded region moves coefficients.
loo <- bind_rows(lapply(regions,function(r){
  dat <- d %>% filter(region!=r)
  fit <- lm(f_district,data=dat)
  get_stats(fit,dat,paste0('Omit ',r),'district') %>% mutate(omitted_region=r)
}))
write_csv(loo,file.path(output,'04_leave_one_region_out.csv'))
loo_summary <- loo %>% group_by(term) %>% summarise(min_estimate=min(estimate),
   max_estimate=max(estimate),min_p=min(p_value),max_p=max(p_value),
   n_p_under_005=sum(p_value<.05),n_models=n(),.groups='drop')
write_csv(loo_summary,file.path(output,'05_LOO_summary.csv'))

# Direct visual audit: solar PV estimate across FE and inference choices.
plot1 <- results %>% filter(term=='pv_100') %>%
 ggplot(aes(x=estimate,y=model,color=inference))+
 geom_vline(xintercept=0,linetype='dashed',color='grey50')+
 geom_errorbarh(aes(xmin=lower95,xmax=upper95),height=.14,position=position_dodge(width=.55))+
 geom_point(position=position_dodge(width=.55),size=2.2)+
 labs(title='Solar-potential coefficient changes with geographic controls',
      subtitle='Association with modeled relative wealth; 95% cluster-based intervals',
      x='RWI coefficient per 100 kWh/kWp annual modeled PV potential',y=NULL,color='Cluster scheme',
      caption='No causal interpretation; block clustering is approximate')+
 theme_minimal(base_size=12)
ggsave(file.path(output,'06_solar_specification_sensitivity.png'),plot1,width=12,height=6.5,dpi=300)
plot2 <- loo %>% filter(term=='pv_100') %>%
 ggplot(aes(x=reorder(omitted_region,estimate),y=estimate))+
 geom_hline(yintercept=0,linetype='dashed',color='grey55')+
 geom_pointrange(aes(ymin=lower95,ymax=upper95))+
 coord_flip()+labs(title='Leave-one-region-out stability: district fixed effects',
   subtitle='Each point excludes all settlements in the named region',
   x='Omitted region',y='Solar PV coefficient and 95% district-clustered interval')+
 theme_minimal(base_size=11)
ggsave(file.path(output,'07_solar_leave_one_region_out.png'),plot2,width=10,height=7,dpi=300)

cat('\nPHASE 5 COMPLETE:',output,'\n')
cat('\nSPECIFICATION COMPARISON (solar potential):\n')
print(results %>% filter(term=='pv_100') %>%
  select(model,inference,n,clusters,estimate,se,lower95,upper95,p_value),n=Inf,width=Inf)
cat('\nDISTRICT SOLAR VARIATION SUMMARY:\n')
print(tibble(districts=nrow(within),districts_n_below_10=sum(within$n<10),
  districts_pv_sd_below_10=sum(is.finite(within$pv_sd)&within$pv_sd<10),
  median_within_district_pv_sd=median(within$pv_sd,na.rm=TRUE),
  districts_used_restricted=length(eligible)),width=Inf)
cat('\nRESTRICTED DISTRICT RESULTS:\n');print(restricted_results,width=Inf)
cat('\nLEAVE-ONE-REGION-OUT SUMMARY:\n');print(loo_summary,width=Inf)
cat('\nNOTE: Remaining residual spatial autocorrelation, timing mismatch and modeled indicators limit interpretation.\n')