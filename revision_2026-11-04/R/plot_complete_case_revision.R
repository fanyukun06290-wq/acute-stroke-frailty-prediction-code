Sys.setlocale("LC_CTYPE","Chinese (Simplified)_China.utf8")
.libPaths(unique(c(Sys.getenv("STROKE_R_LIBRARY", ""),.libPaths())))
suppressPackageStartupMessages(library(ggplot2))
root<-file.path(getwd(),"revision_2026-11-04");out<-file.path(root,"complete_case_results")
dir.create(file.path(out,"figures"),showWarnings=FALSE)
save_fig<-function(p,name,w,h){
 ggsave(file.path(out,"figures",paste0(name,".png")),p,width=w,height=h,dpi=300)
 ggsave(file.path(out,"figures",paste0(name,".svg")),p,width=w,height=h,device=grDevices::svg)
}
labels<-c(age_years="Age",sex="Sex",glucose="Glucose",sodium="Sodium",potassium="Potassium",chloride="Chloride",bicarbonate="Bicarbonate",
 bun="Blood urea nitrogen",creatinine="Creatinine",calcium="Calcium",wbc="White blood cell count",platelets="Platelet count",
 hemoglobin="Hemoglobin",hematocrit="Hematocrit",heart_rate="Heart rate",respiratory_rate="Respiratory rate",spo2="SpO2",sbp="Systolic BP",dbp="Diastolic BP",map="Mean arterial pressure")
window<-function(x)factor(paste0(x,"-hour landmark"),levels=c("24-hour landmark","6-hour landmark"))
d<-read.csv(file.path(out,"raw_input_missingness.csv"))
d<-d[d$population!="MIMIC_all_analysis",]
d$variable<-sub("_(mean_24h|first_6h)$","",d$raw_input)
d$variable<-factor(d$variable,levels=rev(names(labels)))
d$cohort<-factor(d$population,levels=c("MIMIC_training","MIMIC_internal_test","eICU_external_validation"),labels=c("MIMIC\ntraining","MIMIC\ninternal test","eICU\nexternal validation"))
d$window<-window(d$window_hours)
d$text<-ifelse(d$missing>0 & d$missing_fraction<.01,sprintf("%.2f%%",100*d$missing_fraction),sprintf("%.1f%%",100*d$missing_fraction))
p<-ggplot(d,aes(cohort,variable,fill=missing_fraction))+geom_tile(color="white",linewidth=.4)+
 geom_text(aes(label=text,color=missing_fraction>.6),size=2.8)+scale_color_manual(values=c(`FALSE`="#222222",`TRUE`="white"),guide="none")+
 scale_fill_gradient(low="#F7FBFF",high="#08519C",limits=c(0,1),labels=function(x)paste0(round(100*x),"%"))+
 scale_y_discrete(labels=labels)+facet_wrap(~window,ncol=2)+labs(title="Missing raw inputs before training-derived imputation",x=NULL,y=NULL,fill="Unavailable",
 caption="Missingness follows the original time-availability, physiological-range and extraction rules; it does not necessarily mean no measurement was taken.\nComplete cases require all 20 raw inputs, including known F/M sex; derived missingness flags are not additional completeness criteria.\nRates are record-based. 6-hour values are first eligible measurements; 24-hour values are means, so the windows are not interchangeable.")+
 theme_minimal(base_size=11)+theme(panel.grid=element_blank(),legend.position="bottom",plot.caption=element_text(hjust=0,size=9),axis.text.x=element_text(size=9))
save_fig(p,"Raw_input_missingness_heatmap",11.8,8.7)
d<-read.csv(file.path(out,"complete_case_metric_intervals.csv"));d<-d[d$dataset=="eICU_external_validation" & d$metric %in% c("auroc","calibration_intercept"),]
d$window<-window(d$window_hours);d$metric<-factor(d$metric,levels=c("auroc","calibration_intercept"),labels=c("AUROC","Calibration-in-the-large intercept"))
d$key<-paste(d$window,d$subset,sep="::")
keyorder<-unlist(lapply(levels(d$window),function(w)paste(w,c("Incomplete cases","Complete cases","All records"),sep="::")))
d$key<-factor(d$key,levels=keyorder)
lut<-setNames(paste0(d$subset,"\nN=",d$records,"; E=",d$events),as.character(d$key))
vline<-data.frame(metric=factor("Calibration-in-the-large intercept",levels=levels(d$metric)),x=0)
p<-ggplot(d,aes(estimate,key,color=subset))+geom_vline(data=vline,aes(xintercept=x),inherit.aes=FALSE,linetype="dashed",color="grey60")+
 geom_errorbar(aes(xmin=lower,xmax=upper),orientation="y",width=.15,linewidth=.6)+geom_point(size=2.2)+
 facet_grid(window~metric,scales="free")+scale_y_discrete(labels=lut)+
 scale_color_manual(values=c("All records"="#222222","Complete cases"="#2166AC","Incomplete cases"="#C47C22"))+
 labs(title="eICU: frozen-model performance after complete-case restriction",x="Estimate (95% conditional cluster-bootstrap interval)",y=NULL,
 caption="N: ICU records; E: endpoint events. Same frozen XGBoost predictions, with no refitting or recalibration.\nIntervals use whole-cohort connected-hospital/patient resampling, then subset restriction (1,000 replicates).\nDifferences concern selected populations, not the effect of imputation or model improvement. A CI crossing zero does not establish adequate calibration.")+
 theme_bw(base_size=11)+theme(legend.position="none",panel.grid.minor=element_blank(),panel.grid.major.y=element_blank(),plot.caption=element_text(hjust=0,size=9))
save_fig(p,"Complete_case_external_performance",12.3,7.4)
