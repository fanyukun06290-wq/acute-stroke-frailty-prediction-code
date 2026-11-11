Sys.setlocale("LC_CTYPE","Chinese (Simplified)_China.utf8")
.libPaths(unique(c(Sys.getenv("STROKE_R_LIBRARY", ""),.libPaths())))
suppressPackageStartupMessages(library(ggplot2))
root<-file.path(getwd(),"revision_2026-11-04");out<-file.path(root,"subgroup_results")
dir.create(file.path(out,"figures"),showWarnings=FALSE)
theme<-theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),panel.grid.major.y=element_blank(),plot.caption=element_text(hjust=0,size=9),legend.position="bottom")
save_fig<-function(p,name,w,h){
 ggsave(file.path(out,"figures",paste0(name,".png")),p,width=w,height=h,dpi=300)
 ggsave(file.path(out,"figures",paste0(name,".svg")),p,width=w,height=h,device=grDevices::svg)
}
pop<-read.csv(file.path(out,"subgroup_populations.csv"));m<-read.csv(file.path(out,"subgroup_metric_intervals.csv"))
dd<-merge(m[m$metric=="auroc",],pop,by=c("window_hours","dataset","family","subgroup"),suffixes=c("","_pop"))
dd<-dd[dd$status=="Reported",]
order<-c("All records","<65",">=65","Female","Male","AIS_STRICT","ICH_STRICT","SAH","Mixed/ambiguous","Non-mixed AIS/ICH/SAH","ICH/SAH combined")
labels<-c("All records"="Overall","<65"="Age <65",">=65"="Age >=65","Female"="Female","Male"="Male","AIS_STRICT"="Ischemic (AIS_STRICT)","ICH_STRICT"="Intracerebral hemorrhage","SAH"="Subarachnoid hemorrhage","Mixed/ambiguous"="Mixed/ambiguous","Non-mixed AIS/ICH/SAH"="Non-mixed AIS/ICH/SAH","ICH/SAH combined"="ICH/SAH combined")
dd$y<-factor(dd$subgroup,levels=rev(order))
dd$panel<-factor(paste0(dd$window_hours,"h / ",ifelse(dd$dataset=="MIMIC_internal_test","MIMIC internal test","eICU external validation")),
 levels=c("24h / MIMIC internal test","24h / eICU external validation","6h / MIMIC internal test","6h / eICU external validation"))
dd$count<-sprintf("N=%s; E=%s%s",dd$records,dd$events,ifelse(dd$events<20|dd$non_events<20|dd$clusters<20|dd$records<100," *",""))
p<-ggplot(dd,aes(estimate,y,color=family))+geom_vline(xintercept=.5,linetype="dotted",color="grey70")+
 geom_errorbar(aes(xmin=lower,xmax=upper),orientation="y",width=.18,linewidth=.55)+geom_point(size=2.1)+
 geom_vline(xintercept=1.015,color="grey85",linewidth=.3)+geom_text(aes(x=1.19,label=count),hjust=1,color="#444444",size=2.7)+facet_wrap(~panel,ncol=2)+
 scale_y_discrete(labels=labels)+scale_x_continuous(limits=c(.40,1.20),breaks=c(.4,.5,.6,.7,.8,.9,1))+
 scale_color_manual(values=c("Overall"="#222222","Age"="#7B3294","Sex"="#008837","Phenotype"="#2166AC","Restricted evaluation"="#C47C22"))+
 labs(title="Frozen XGBoost: subgroup discrimination and precision",x="AUROC (95% conditional cluster-bootstrap interval)",y=NULL,color=NULL,
 caption="N: ICU records; E: observed endpoint events. * Small-sample flag. Empty/small-cell sex-unknown groups are omitted; counts appear in a separate right column.\nPatients are resampled internally; connected hospital/patient components externally. Same whole-cohort weights used across subgroups (1,000 replicates).\nRestricted evaluations use the unchanged pooled model. AUROC differences are descriptive, not tests of effect modification or fairness.")+theme
save_fig(p,"Subgroup_AUROC_forest",14.5,9)
d<-read.csv(file.path(out,"phenotype_case_mix.csv"));d<-d[d$partition %in% c("Training","All_analysis_records") & !(d$database=="mimic" & d$partition=="All_analysis_records"),]
d$cohort<-factor(ifelse(d$database=="mimic","MIMIC training","eICU external validation"),levels=c("MIMIC training","eICU external validation"))
d$window<-factor(paste0(d$window_hours,"-hour landmark"),levels=c("24-hour landmark","6-hour landmark"))
p<-ggplot(d,aes(cohort,record_fraction,fill=phenotype))+geom_col(width=.65)+
 geom_text(aes(label=sprintf("%.1f%%",100*record_fraction)),position=position_stack(vjust=.5),size=3.2)+
 facet_wrap(~window,ncol=2)+scale_y_continuous(labels=function(x)paste0(round(100*x),"%"),expand=expansion(mult=c(0,.01)))+
 scale_fill_manual(values=c("AIS_STRICT"="#2166AC","ICH_STRICT"="#B2182B","SAH"="#C47C22","Mixed/ambiguous"="#999999"))+
 labs(title="Stroke case mix differs between development and external validation",x=NULL,y="Proportion of ICU records",fill=NULL,
 caption="Mutually exclusive retrospective coding categories; record-based denominators, not unique-patient proportions.\nMixed/ambiguous includes multiple acute groups, ambiguous secondary codes and mixed CVST. Phenotype was not a model input.\nThe two landmarks have different eligible populations. Differences may affect discrimination and calibration during transport.")+theme
save_fig(p,"Phenotype_case_mix",10.5,6)
