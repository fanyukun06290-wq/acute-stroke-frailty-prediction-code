Sys.setlocale("LC_CTYPE","Chinese (Simplified)_China.utf8")
.libPaths(unique(c(Sys.getenv("STROKE_R_LIBRARY", ""),.libPaths())))
suppressPackageStartupMessages(library(ggplot2))
root<-file.path(getwd(),"revision_2026-11-04");out<-file.path(root,"endpoint_results")
dir.create(file.path(out,"figures"),showWarnings=FALSE)
theme<-theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),legend.position="bottom",plot.caption=element_text(hjust=0,size=9))
save_fig<-function(p,name,w,h){
 ggsave(file.path(out,"figures",paste0(name,".png")),p,width=w,height=h,dpi=300)
 ggsave(file.path(out,"figures",paste0(name,".svg")),p,width=w,height=h,device=grDevices::svg)
}
d<-read.csv(file.path(out,"competing_event_curves.csv"))
d$panel<-factor(paste0(d$window_hours,"h / ",ifelse(d$dataset=="MIMIC_internal_test","MIMIC internal test","eICU external validation")),levels=c("24h / MIMIC internal test","24h / eICU external validation","6h / MIMIC internal test","6h / eICU external validation"))
a<-d[d$metric!="one_minus_km_censor_discharge",]
labels<-c(death_cif="In-hospital death",alive_discharge_cif="Alive hospital discharge",remaining_hospitalized="Alive and still hospitalized")
a$label<-unname(labels[a$metric]);cols<-c("In-hospital death"="#B2182B","Alive hospital discharge"="#2166AC","Alive and still hospitalized"="#666666")
p<-ggplot(a,aes(icu_day,estimate,color=label,fill=label))+geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.08,color=NA)+geom_line(linewidth=.7)+facet_wrap(~panel,ncol=2)+
 scale_color_manual(values=cols)+scale_fill_manual(values=cols)+scale_y_continuous(labels=function(x)paste0(round(100*x),"%"),limits=c(0,1))+
 labs(title="Competing events after the ICU landmark",x="Days since ICU admission",y="State probability / cumulative incidence",color=NULL,fill=NULL,
 caption="Death and alive hospital discharge are absorbing events; remaining observations are administratively censored at ICU day 30.\nRibbons: 95% pointwise conditional cluster-bootstrap intervals (1,000 replicates), not simultaneous bands. Curves start at each landmark.\nAt day 30, event probabilities reflect the original strict <30-day boundary. No post-discharge mortality is observed.")+theme
save_fig(p,"Competing_event_state_probabilities",11.8,8.5)
a<-d[d$metric %in% c("death_cif","one_minus_km_censor_discharge"),]
a$label<-ifelse(a$metric=="death_cif","AJ death CIF (discharge competes)","1-KM (discharge censored; diagnostic only)")
p<-ggplot(a,aes(icu_day,estimate,color=label,linetype=label))+geom_line(linewidth=.8)+facet_wrap(~panel,ncol=2)+
 scale_color_manual(values=c("AJ death CIF (discharge competes)"="#B2182B","1-KM (discharge censored; diagnostic only)"="#666666"))+
 scale_linetype_manual(values=c("AJ death CIF (discharge competes)"="solid","1-KM (discharge censored; diagnostic only)"="dashed"))+
 scale_y_continuous(labels=function(x)paste0(round(100*x),"%"),limits=c(0,1))+
 labs(title="Censoring alive discharge changes the estimated quantity",x="Days since ICU admission",y="Estimated probability",color=NULL,linetype=NULL,
 caption="AJ: Aalen-Johansen cumulative incidence of death before discharge. Final AJ CIF equals the binary endpoint event fraction in these retained records.\nDischarge-censored 1-KM is a different, hypothetical quantity; it is not an estimate of observed in-hospital absolute risk or complete all-cause mortality.\nNeither method recovers post-discharge deaths or removes differences in discharge practices. Frozen predictions are unchanged.")+theme
save_fig(p,"AJ_vs_discharge_censored_KM",11.8,8.2)
r<-read.csv(file.path(out,"fixed_risk_bin_cif.csv"));r<-r[r$dataset=="eICU_external_validation",]
key<-c("window_hours","dataset","risk_bin","records","events")
obs<-r[r$metric=="death_cif",c(key,"estimate","lower","upper")]
pr<-r[r$metric=="mean_predicted",c(key,"estimate")];names(pr)[ncol(pr)]<-"mean_predicted"
a<-merge(obs,pr,by=key);a$window<-factor(paste0(a$window_hours,"-hour landmark"),levels=c("24-hour landmark","6-hour landmark"))
a$risk_bin<-factor(a$risk_bin,levels=c("<0.05","0.05-<0.10","0.10-<0.20",">=0.20"))
p<-ggplot(a,aes(mean_predicted,estimate,color=risk_bin))+geom_abline(slope=1,intercept=0,linetype="dashed",color="grey60")+
 geom_errorbar(aes(ymin=lower,ymax=upper),width=.008)+geom_point(size=2.5)+facet_wrap(~window,ncol=2)+
 scale_x_continuous(breaks=seq(0,.5,.1))+scale_y_continuous(breaks=seq(0,.5,.1))+coord_equal(xlim=c(0,.5),ylim=c(0,.5))+
 labs(title="eICU: observed death CIF in fixed predicted-risk strata",x="Mean frozen predicted probability",y="Observed death CIF before ICU day 30",color="Fixed risk bin",
 caption="Bins: <0.05, 0.05-<0.10, 0.10-<0.20, >=0.20; no external optimization. Vertical bars: 95% pointwise cluster-bootstrap intervals.\nAxes zoomed to 0-0.50; all bins retained. The final death CIF equals the event fraction here, and does not add independent calibration evidence.\nThis plot does not establish adequate local calibration or clinical readiness; post-discharge mortality remains unavailable.")+theme
save_fig(p,"External_risk_bin_CIF_calibration",10.8,5.8)
