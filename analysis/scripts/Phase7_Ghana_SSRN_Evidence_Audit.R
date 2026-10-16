# PHASE 7 — Ghana SSRN: evidence strength audit, ALL covariates
# Reads completed Phase 6 outputs; does not fabricate/re-estimate missing results.
# Interpretations remain cross-sectional associations, not causal effects or verified electrification.
packages <- c('readr','dplyr','tidyr','ggplot2','tibble')
new <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(new)) install.packages(new,repos='https://cloud.r-project.org')
suppressPackageStartupMessages({library(readr);library(dplyr);library(tidyr);library(ggplot2);library(tibble)})
home <- Sys.getenv('USERPROFILE')
possible <- c(file.path(home,'OneDrive','Desktop','PETER AND KONADU SSRN'),
              file.path(home,'Desktop','PETER AND KONADU SSRN'),
              file.path(home,'OneDrive - Personal','Desktop','PETER AND KONADU SSRN'))
folder <- possible[dir.exists(possible)][1]
if(is.na(folder)) stop('Research folder not located; set folder manually.')
infile <- file.path(folder,'Phase6_Outputs','01_spatial_HAC_comparison.csv')
if(!file.exists(infile)) stop('Missing Phase 6 results: ',infile)
out <- file.path(folder,'Phase7_Outputs')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
d <- read_csv(infile,show_col_types=FALSE)
required <- c('model','method','term','estimate','se','lower95','upper95','p_value')
if(!all(required %in% names(d))) stop('Expected Phase 6 columns absent: ',paste(setdiff(required,names(d)),collapse=', '))
keys <- c('log_grid','pv_100','log_road','log_pop')
mods <- c('Region FE','District FE','District FE + trend')
methods <- c('District clustered','Conley 25 km','Conley 50 km','Conley 100 km')
d <- d %>% filter(term %in% keys,model %in% mods,method %in% methods)
if(anyDuplicated(d[c('model','method','term')])) stop('Duplicate model/method/term combinations.')
if(nrow(d)!=length(keys)*length(mods)*length(methods)) warning('Not all 48 expected comparisons are present; results will show available combinations only.')
if(any(!is.finite(d$estimate)|!is.finite(d$se)|d$se<=0)) stop('Non-finite or nonpositive standard errors detected.')
# 95% interval excludes zero; do not conflate p-value cutoffs across distinct inference methods.
d <- d %>% mutate(direction=case_when(lower95>0~'positive',upper95<0~'negative',TRUE~'uncertain'),
                  interval_excludes_zero=(lower95>0|upper95<0),
                  ci_width=upper95-lower95)
write_csv(d,file.path(out,'01_all_predictors_HAC_intervals.csv'))
summary <- d %>% group_by(model,term) %>% summarise(
  n_methods=n(),n_intervals_exclude_zero=sum(interval_excludes_zero),
  n_positive=sum(direction=='positive'),n_negative=sum(direction=='negative'),
  n_uncertain=sum(direction=='uncertain'),
  min_se=min(se),max_se=max(se),.groups='drop')
write_csv(summary,file.path(out,'02_interval_consistency_summary.csv'))
# Focus 100km comparisons: strongest sensitivity in selected distance grid, NOT universally best.
long <- d %>% filter(method=='Conley 100 km') %>%
 select(model,term,estimate,se,lower95,upper95,p_value,direction)
write_csv(long,file.path(out,'03_conley_100km_main_results.csv'))
# For audit, verify phase 6 reproduction table matches visible coefficients within rounding.
check_path <- file.path(folder,'Phase6_Outputs','00_reproduction_check.csv')
if(file.exists(check_path)) {
 check <- read_csv(check_path,show_col_types=FALSE)
 write_csv(check,file.path(out,'04_reproduction_check_copy.csv'))
}
# Small multiple plot, every predictor allowed its own axis; avoids mixing coefficient magnitudes.
labels <- c(log_grid='log(1 + transmission distance, km)',
            pv_100='PV potential per 100 kWh/kWp/year',
            log_road='log(1 + road distance, km)',
            log_pop='log(1 + population)')
p <- d %>% mutate(term=factor(term,levels=keys,labels=labels[keys]),
                  model=factor(model,levels=rev(mods)),method=factor(method,levels=methods)) %>%
 ggplot(aes(x=estimate,y=model,color=method))+
 geom_vline(xintercept=0,linetype='dashed',color='grey45')+
 geom_errorbarh(aes(xmin=lower95,xmax=upper95),height=.12,
                position=position_dodge(width=.68))+
 geom_point(position=position_dodge(width=.68),size=1.7)+
 facet_wrap(~term,scales='free_x',ncol=2)+
 labs(title='Ghana settlement wealth associations: inference sensitivity',
      subtitle='95% intervals; coefficient estimates change with FE model, not SE method',
      x='Coefficient (RWI units)',y=NULL,color='Inference',
      caption='Observational associations only. Conley cutoffs are sensitivity choices, not causal correction.')+
 theme_minimal(base_size=11)+theme(legend.position='bottom')
ggsave(file.path(out,'05_all_predictors_inference_sensitivity.png'),p,width=13,height=8,dpi=300)
# Publication-quality prose must be tied to exact result table; write an audit memo, not invented claims.
lines <- c('PHASE 7 EVIDENCE AUDIT — Ghana SSRN',
 'Source: Phase6_Outputs/01_spatial_HAC_comparison.csv',
 'Interpretation: modeled Relative Wealth Index is the outcome.',
 'A distance coefficient is NOT a causal effect of grid or road construction.',
 'Transmission-line proximity is NOT verified household electrification.',
 'PV value is modeled long-run output potential, NOT actual installed generation.',
 'Spatial HAC protects inference under its assumptions but DOES NOT remove confounding.',
 'The largest tested cutoff (100km) is a sensitivity choice, not a validated true spatial range.',
 'Check polygon district membership and temporal mismatch before manuscript submission.',
 'Compare all four predictors under region and district FE, rather than selecting only significant results.',
 'Reported signs and 95% interval conclusions are in CSV files 02 and 03.',
 'Do not claim causal impacts, investment effects, or SSRN acceptance from these results.')
writeLines(lines,file.path(out,'06_interpretation_and_limits.txt'))
cat('\nPHASE 7 COMPLETE:',out,'\n')
cat('\n100 KM CONLEY RESULTS — ALL MAIN VARIABLES:\n');print(long,n=Inf,width=Inf)
cat('\nCONFIDENCE-INTERVAL CONSISTENCY:\n');print(summary,n=Inf,width=Inf)