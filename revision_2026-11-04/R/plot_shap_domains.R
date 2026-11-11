Sys.setlocale("LC_CTYPE","Chinese (Simplified)_China.utf8")
.libPaths(unique(c(Sys.getenv("STROKE_R_LIBRARY", ""),.libPaths())))
suppressPackageStartupMessages(library(ggplot2))
root<-file.path(getwd(),"revision_2026-11-04");out<-file.path(root,"shap_domain_results")
dir.create(file.path(out,"figures"),showWarnings=FALSE)
theme<-theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),panel.grid.major.y=element_blank(),
  legend.position="bottom",plot.caption=element_text(hjust=0,size=9),panel.spacing=grid::unit(1.2,"lines"))
colors<-c("Demographic"="#7B3294","Metabolic/renal"="#008837","Inflammation-related hematologic"="#C47C22","Vital signs"="#2166AC")
label<-function(x){
 n<-c(age_years="Age",sex="Sex",glucose="Glucose",sodium="Sodium",potassium="Potassium",chloride="Chloride",bicarbonate="Bicarbonate",
 bun="Blood urea nitrogen",creatinine="Creatinine",calcium="Calcium",wbc="White blood cell count",platelets="Platelet count",
 hemoglobin="Hemoglobin",hematocrit="Hematocrit",heart_rate="Heart rate",respiratory_rate="Respiratory rate",spo2="SpO2",sbp="Systolic BP",dbp="Diastolic BP",map="Mean arterial pressure")
 unname(n[x])
}
save_fig<-function(p,name,w,h){
 ggsave(file.path(out,"figures",paste0(name,".png")),p,width=w,height=h,dpi=300)
 ggsave(file.path(out,"figures",paste0(name,".svg")),p,width=w,height=h,device=grDevices::svg)
}
d<-read.csv(file.path(out,"conceptual_variable_ranking.csv"));d<-d[d$rank<=15,]
d$panel<-paste0(d$window_hours,"h / ",ifelse(d$dataset=="MIMIC_internal_test","MIMIC internal test","eICU external validation"))
d$panel<-factor(d$panel,levels=c("24h / MIMIC internal test","24h / eICU external validation","6h / MIMIC internal test","6h / eICU external validation"))
d$key<-paste(d$panel,d$name,sep="::");d$key<-factor(d$key,levels=d$key[order(d$estimate)])
labsmap<-setNames(label(d$name),as.character(d$key))
p<-ggplot(d,aes(estimate,key,fill=domain))+geom_col(width=.72)+geom_errorbar(aes(xmin=lower,xmax=upper),orientation="y",width=.2,linewidth=.35)+
 facet_wrap(~panel,ncol=2,scales="free_y")+scale_y_discrete(labels=labsmap)+scale_fill_manual(values=colors)+
 labs(title="Frozen XGBoost: validation-set SHAP contributions",subtitle="Top 15 conceptual variables; measurements and corresponding missingness flags grouped by signed sum",
 x="Mean absolute grouped SHAP (log-odds)",y=NULL,fill=NULL,
 caption="Bars: mean |sum of feature SHAP within variable|. Whiskers: 95% conditional cluster-bootstrap intervals (1,000 replicates).\nNative exact TreeSHAP uses training tree path covers; aggregation does not recompute coalition Shapley values.\nThe two landmarks concern different eligible populations. Contributions are predictive allocations, not causal effects.")+theme+
 guides(fill=guide_legend(nrow=2,byrow=TRUE))
save_fig(p,"SHAP_validation_top_variables",12.5,11.5)
d<-read.csv(file.path(out,"external_domain_contributions.csv"))
metrics<-c("sum_absolute","absolute_signed_sum","sum_absolute_per_encoded_feature")
d<-d[d$metric %in% metrics,];d$window<-factor(paste0(d$window_hours,"-hour landmark"),levels=c("24-hour landmark","6-hour landmark"))
d$definition<-factor(d$metric,levels=metrics,labels=c("Sum of mean absolute feature SHAP","Mean absolute signed domain sum","Sum divided by encoded feature count"))
d$domain<-factor(d$domain,levels=rev(names(colors)))
p<-ggplot(d,aes(estimate,domain,fill=domain))+geom_col(width=.65)+geom_errorbar(aes(xmin=lower,xmax=upper),orientation="y",width=.15,linewidth=.4)+
 facet_grid(window~definition,scales="free_x")+scale_fill_manual(values=colors)+scale_y_discrete(labels=function(x)gsub("Inflammation-related hematologic","Inflammation-related\nhematologic",x))+
 labs(title="eICU: domain rankings depend on aggregation",x="Contribution magnitude (log-odds)",y=NULL,
 caption="Three descriptive summaries of the same native feature SHAP values. Summed absolute magnitudes depend on domain size.\nSigned contributions can cancel within a domain. Normalization answers a different question; none establishes intrinsic domain importance.\nWhiskers: 95% conditional connected-hospital/patient cluster-bootstrap intervals; frozen models, no recalibration.")+theme+theme(legend.position="none")
save_fig(p,"SHAP_external_domain_definitions",14,6.8)
d<-read.csv(file.path(out,"domain_permutation_summary.csv"));d<-d[d$dataset=="eICU_external_validation" & d$metric=="auroc",]
d$window<-factor(paste0(d$window_hours,"-hour landmark"),levels=c("24-hour landmark","6-hour landmark"))
d$domain<-factor(d$domain,levels=rev(names(colors)))
p<-ggplot(d,aes(mean,domain,fill=domain))+geom_vline(xintercept=0,color="grey55",linetype="dashed")+
 geom_col(width=.65)+geom_errorbar(aes(xmin=q025,xmax=q975),orientation="y",width=.15,linewidth=.45)+
 facet_grid(window~scope)+scale_fill_manual(values=colors)+scale_y_discrete(labels=function(x)gsub("Inflammation-related hematologic","Inflammation-related\nhematologic",x))+
 labs(title="eICU: joint domain permutation sensitivity",x="AUROC decrease after joint shuffling",y=NULL,
 caption="Mean of 50 shuffles; whiskers show 2.5th-97.5th percentiles of shuffle variation, NOT patient-sampling confidence intervals.\nMeasurements and missingness flags move together within a domain. Shuffling can break cross-domain associations.\nWithin-hospital shuffles preserve hospital membership. These diagnostics do not establish causal or independent domain effects.")+theme+theme(legend.position="none")
save_fig(p,"Domain_joint_permutation_external",11.8,6.8)
