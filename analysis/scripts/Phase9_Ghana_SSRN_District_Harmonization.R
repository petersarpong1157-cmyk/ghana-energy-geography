# PHASE 9 | Ghana SSRN: geographic label harmonization and model sensitivity
# Requires original settlement CSV and Phase8_Outputs/03_all_settlement_polygon_audit.csv.
# Never silently overwrite labels. Results are observational, not causal.
packages <- c('readr','dplyr','tidyr','stringr','sandwich','ggplot2')
new <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(new)) install.packages(new,repos='https://cloud.r-project.org')
suppressPackageStartupMessages({library(readr);library(dplyr);library(tidyr);library(stringr);library(ggplot2)})
userhome <- Sys.getenv('USERPROFILE')
paths <- c(file.path(userhome,'OneDrive','Desktop','PETER AND KONADU SSRN'),file.path(userhome,'Desktop','PETER AND KONADU SSRN'),file.path(userhome,'OneDrive - Personal','Desktop','PETER AND KONADU SSRN'))
folder <- paths[dir.exists(paths)][1]
if(length(folder)==0 || is.na(folder)) stop('Cannot locate Desktop/PETER AND KONADU SSRN. Set `folder` manually.')
source_csv <- file.path(folder,'ghana_dre_atlas_settlements.csv')
audit_csv <- file.path(folder,'Phase8_Outputs','03_all_settlement_polygon_audit.csv')
if(!all(file.exists(c(source_csv,audit_csv)))) stop('Missing input CSV or Phase 8 audit; check folder paths.')
out <- file.path(folder,'Phase9_Outputs');dir.create(out,recursive=TRUE,showWarnings=FALSE)
raw <- read_csv(source_csv,show_col_types=FALSE,progress=FALSE)
audit <- read_csv(audit_csv,show_col_types=FALSE,progress=FALSE)
required <- c('geohash','admin_cgaz_1','admin_cgaz_2','mean_rwi','pv_value','distance_to_existing_transmission_lines','dist_main_road_km','population','lat','lon')
if(!all(required %in% names(raw))) stop('Missing columns in source: ',paste(setdiff(required,names(raw)),collapse=', '))
require_audit <- c('geohash','boundary_set','polygon_matches','polygon_region','polygon_district')
if(!all(require_audit %in% names(audit))) stop('Phase 8 audit missing: ',paste(setdiff(require_audit,names(audit)),collapse=', '))
if(anyDuplicated(raw$geohash)) stop('Duplicate source geohash values')
if(anyDuplicated(audit[c('boundary_set','geohash')])) stop('Duplicate boundary_set/geohash audit values')
if(!setequal(unique(audit$boundary_set),c('District_261','District_272'))) stop('Unexpected boundary sets in audit.')
norm <- function(x){x<-iconv(as.character(x),to='ASCII//TRANSLIT');x<-tolower(trimws(x));x<-gsub('\\b(municipal|metropolitan|district|assembly|region)\\b',' ',x);trimws(gsub(' +',' ',gsub('[^a-z0-9]+',' ',x)))}
# Geography source - polygon joins use unique match only.
wide <- audit %>% select(geohash,boundary_set,polygon_matches,polygon_region,polygon_district) %>%
 pivot_wider(names_from=boundary_set,values_from=c(polygon_matches,polygon_region,polygon_district),names_sep='__')
if(nrow(wide)!=nrow(raw)) stop('Phase 8 audit and settlement data have different ID coverage.')
d <- raw %>% left_join(wide,by='geohash') %>%
 transmute(geohash, source_region=trimws(as.character(admin_cgaz_1)), source_district=trimws(as.character(admin_cgaz_2)),
   n261=polygon_matches__District_261,n272=polygon_matches__District_272,
   region261=polygon_region__District_261,region272=polygon_region__District_272,
   district261=polygon_district__District_261,district272=polygon_district__District_272,
   rwi=as.numeric(mean_rwi),pv=as.numeric(pv_value),
   grid=as.numeric(distance_to_existing_transmission_lines),road=as.numeric(dist_main_road_km),
   pop=as.numeric(population),lat=as.numeric(lat),lon=as.numeric(lon)) %>%
 mutate(across(c(source_region,source_district,region261,region272,district261,district272),~na_if(.x,'')),
   matches_both=n261==1 & n272==1,
   source261_equal=norm(source_district)==norm(district261),
   source272_equal=norm(source_district)==norm(district272),
   names261272_equal=norm(district261)==norm(district272),
   pv_100=pv/100,log_grid=log1p(grid),log_road=log1p(road),log_pop=log1p(pop))
# Diagnostic examples, including top mismatch combinations (do not automatically reassign source districts).
comparable <- d %>% filter(matches_both,!is.na(source_district),!is.na(district261),!is.na(district272))
by_pair <- comparable %>% count(source_region,source_district,district261,district272,source261_equal,source272_equal,sort=TRUE)
write_csv(by_pair,file.path(out,'01_source_vs_polygon_name_pairs.csv'))
write_csv(comparable %>% filter(!source261_equal | !source272_equal | !names261272_equal) %>%
 select(geohash,lat,lon,source_region,source_district,region261,district261,region272,district272,source261_equal,source272_equal,names261272_equal),
 file.path(out,'02_geographic_mismatch_locations.csv'))
write_csv(d %>% filter(!matches_both | is.na(source_region) | !(source_region %in% c('Ahafo Region','Ashanti Region','Bono Region','Bono East Region','Central Region','Eastern Region','Greater Accra Region','North East Region','Northern Region','Oti Region','Savannah Region','Upper East Region','Upper West Region','Volta Region','Western Region','Western North Region'))) %>%
 select(geohash,lat,lon,source_region,source_district,n261,n272,region261,district261,region272,district272,rwi,pv),
 file.path(out,'03_problematic_geographies_for_manual_review.csv'))
# Common sample deliberately restricted to Ghana recognized region names, both polygon assignments unique, and complete model fields.
regions <- c('Ahafo Region','Ashanti Region','Bono Region','Bono East Region','Central Region','Eastern Region','Greater Accra Region','North East Region','Northern Region','Oti Region','Savannah Region','Upper East Region','Upper West Region','Volta Region','Western Region','Western North Region')
model <- d %>% filter(matches_both,source_region %in% regions,!is.na(source_district),!is.na(region261),!is.na(region272),!is.na(district261),!is.na(district272),
 is.finite(rwi),is.finite(pv_100),is.finite(log_grid),is.finite(log_road),is.finite(log_pop),is.finite(lat),is.finite(lon),grid>=0,road>=0,pop>=0) %>%
 mutate(source_cluster=paste(source_region,source_district,sep=' | '),
  polygon261_cluster=paste(region261,district261,sep=' | '),polygon272_cluster=paste(region272,district272,sep=' | '))
if(nrow(model)<1000) stop('Too few complete observations. Investigate geographic join.')
summary <- tibble(metric=c('original_rows','unique_spatial_matches_both','comparable_source_labels','source_261_name_mismatches','source_272_name_mismatches','261_272_name_mismatches','matched_complete_model_n','source_clusters_common_n','polygon261_clusters_common_n','polygon272_clusters_common_n'),
 value=c(nrow(d),sum(d$matches_both,na.rm=TRUE),nrow(comparable),sum(!comparable$source261_equal,na.rm=TRUE),sum(!comparable$source272_equal,na.rm=TRUE),sum(!comparable$names261272_equal,na.rm=TRUE),nrow(model),dplyr::n_distinct(model$source_cluster),dplyr::n_distinct(model$polygon261_cluster),dplyr::n_distinct(model$polygon272_cluster)))
write_csv(summary,file.path(out,'04_harmonization_summary.csv'))
# Report stability on SAME OBSERVATIONS; cluster choices independently varied.
main_terms <- c('log_grid','pv_100','log_road','log_pop')
fe_specs <- list(source='source_cluster',polygon261='polygon261_cluster',polygon272='polygon272_cluster')
fits <- lapply(names(fe_specs),function(nm){var<-fe_specs[[nm]];lm(as.formula(paste('rwi ~ log_grid + pv_100 + log_road + log_pop + factor(',var,')',sep='')),data=model)})
names(fits) <- names(fe_specs)
get_coef <- function(fit,cluster,modelname,clustername){
 G <- length(unique(cluster)); V <- sandwich::vcovCL(fit,cluster=cluster,type='HC1',cadjust=TRUE)
 b <- coef(fit)[main_terms]; s <- sqrt(diag(V)[main_terms]);df <- G-1
 tibble(model_FE=modelname,inference_clusters=clustername,n=nobs(fit),clusters=G,term=main_terms,
 estimate=unname(b),se=unname(s),lower95=unname(b-qt(.975,df)*s),upper95=unname(b+qt(.975,df)*s),p_value=unname(2*pt(abs(b/s),df=df,lower.tail=FALSE)))
}
results <- bind_rows(lapply(names(fits),function(nm){
 fit <- fits[[nm]]
 bind_rows(lapply(names(fe_specs),function(clnm) get_coef(fit,model[[fe_specs[[clnm]]]],nm,clnm)))
}))
write_csv(results,file.path(out,'05_FE_and_cluster_sensitivity_common_sample.csv'))
# Which geographic labels the original study clusters merged/split?
assignment <- model %>% distinct(source_cluster,polygon261_cluster,polygon272_cluster)
source_splits <- assignment %>% group_by(source_cluster) %>% summarise(n_polygon261=n_distinct(polygon261_cluster),n_polygon272=n_distinct(polygon272_cluster),.groups='drop') %>% arrange(desc(n_polygon261))
poly_merges <- assignment %>% group_by(polygon261_cluster) %>% summarise(n_source=n_distinct(source_cluster),.groups='drop') %>% arrange(desc(n_source))
write_csv(source_splits,file.path(out,'06_source_labels_mapping_multiple_polygon_districts.csv'))
write_csv(poly_merges,file.path(out,'07_polygon_districts_mapping_multiple_source_labels.csv'))
# Plot direct comparison for *fixed* inference choice: each FE uses polygon261 clustering.
plotdata <- results %>% filter(inference_clusters=='polygon261',term %in% c('log_grid','pv_100'))
fig <- ggplot(plotdata,aes(x=estimate,y=model_FE,color=model_FE))+
 geom_vline(xintercept=0,linetype='dashed',color='grey55')+
 geom_errorbarh(aes(xmin=lower95,xmax=upper95),height=.15)+geom_point(size=2.4)+
 facet_wrap(~term,scales='free_x')+theme_minimal(base_size=12)+
 labs(title='Ghana wealth regression: district-label sensitivity',subtitle='Same settlements; different fixed effects; inference clustered by 261-district polygon labels',x='Coefficient and 95% confidence interval',y='Fixed effects',color=NULL,caption='Labels and vintage require provenance review; associations are noncausal.')
ggsave(file.path(out,'08_boundary_label_model_sensitivity.png'),fig,width=10,height=5.5,dpi=300)
cat('\nPHASE 9 COMPLETE:',out,'\n')
cat('\nGEOGRAPHIC HARMONIZATION SUMMARY:\n');print(summary,n=Inf)
cat('\nMODEL SENSITIVITY (POLYGON261 CLUSTERED):\n');print(results %>% filter(inference_clusters=='polygon261'),n=Inf,width=Inf)
cat('\nMOST COMMON SOURCE-POLYGON NAME PAIRS:\n');print(head(by_pair,15),n=15,width=Inf)
cat('\nCAUTION: Both shapefiles require authenticated boundary dates, provenance, and district-code harmonization.\n')
cat('This phase does not automatically declare polygon labels more correct than source labels.\n')