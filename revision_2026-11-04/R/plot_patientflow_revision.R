Sys.setlocale("LC_CTYPE","Chinese (Simplified)_China.utf8")
.libPaths(unique(c(Sys.getenv("STROKE_R_LIBRARY", ""),.libPaths())))
suppressPackageStartupMessages(library(ggplot2))
root<-file.path(getwd(),"revision_2026-11-04");out<-file.path(root,"patientflow_results")
dir.create(file.path(out,"figures"),showWarnings=FALSE)
save_fig<-function(p,name,w,h){
 ggsave(file.path(out,"figures",paste0(name,".png")),p,width=w,height=h,dpi=300)
 ggsave(file.path(out,"figures",paste0(name,".svg")),p,width=w,height=h,device=grDevices::svg)
}
flow<-read.csv(file.path(out,"cohort_flow_counts.csv"));nodes<-edges<-list();nn<-ne<-0L
panel_levels<-c("24h / MIMIC","24h / eICU","6h / MIMIC","6h / eICU")
fmt<-function(x)format(x,big.mark=",",scientific=FALSE,trim=TRUE)
for(h in c(24L,6L))for(db in c("mimic","eicu")){
 q<-flow[flow$window_hours==h & flow$database==db,];panel<-paste0(h,"h / ",if(db=="mimic")"MIMIC" else "eICU")
 box<-function(stage,label,x,y,width=1.9,role="Retained",events=FALSE){
  z<-q[q$stage==stage,];stopifnot(nrow(z)==1)
  text<-paste0(label,"\nN=",fmt(z$records),"; P=",fmt(z$unique_patients),if(events)paste0("\nE=",fmt(z$events)) else "")
  nn<<-nn+1L;nodes[[nn]]<<-data.frame(panel=panel,xmin=x-width/2,xmax=x+width/2,ymin=y-.60,ymax=y+.60,x=x,y=y,label=text,role=role)
 }
 edge<-function(x,y,xe,ye){ne<<-ne+1L;edges[[ne]]<<-data.frame(panel=panel,x=x,y=y,xend=xe,yend=ye)}
 box("Source coded adult index stays","Coded adult index stays",1,8)
 box("Valid chronology","Pass core time checks",1,6)
 box("Landmark eligible","Alive in ICU and hospital\nat landmark",1,4)
 box("Analysis cohort","Known-outcome analysis cohort",1,2,events=TRUE)
 box("Excluded invalid chronology","Excluded: failed core\ntime checks",3.25,7,role="Excluded")
 box("Excluded before/not reaching landmark","Excluded: not reaching\nlandmark in ICU/hospital",3.25,5,role="Excluded")
 box("Excluded unknown/conflicting outcome","Excluded: outcome\nunknown/conflicting",3.25,3,role="Excluded")
 for(y in c(8,6,4))edge(1,y-.60,1,y-1.40)
 for(y in c(7,5,3))edge(1.1,y+.40,2.30,y)
 if(db=="mimic"){
  box("Training partition","Training",.40,0,width=1.65,role="Partition",events=TRUE)
  box("Internal test partition","Internal test",2.20,0,width=1.65,role="Partition",events=TRUE)
  edge(1,1.4,.40,.60);edge(1,1.4,2.20,.60)
 }else{
  box("External validation","External validation",1,0,role="Partition",events=TRUE);edge(1,1.4,1,.60)
 }
}
n<-do.call(rbind,nodes);e<-do.call(rbind,edges);n$panel<-factor(n$panel,levels=panel_levels);e$panel<-factor(e$panel,levels=panel_levels)
p<-ggplot()+geom_segment(data=e,aes(x,y,xend=xend,yend=yend),arrow=grid::arrow(length=grid::unit(.10,"inches")),linewidth=.4,color="#555555")+
 geom_rect(data=n,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=role),color="#555555",linewidth=.35)+
 geom_text(data=n,aes(x,y,label=label),size=3.0,lineheight=.95)+facet_wrap(~panel,ncol=2)+
 scale_fill_manual(values=c(Retained="#EAF2F8",Excluded="#F8EFEA",Partition="#EDF4E8"))+coord_cartesian(xlim=c(-.5,4.4),ylim=c(-.9,8.9),clip="off")+
 labs(title="Source-to-analysis cohort accounting by ICU landmark",
  caption="N: ICU/index-hospitalization records; P: unique longitudinal patients; E: endpoint events. Unique-patient counts across branches are not additive.\nSource boxes start with extracted adult stroke-coded index stays; earlier all-database screening counts are not available in these local files.\nMIMIC partitions are disjoint by patient within each window. eICU uses uniquepid; both windows require survival and continued ICU/hospital stay at the landmark.")+theme_void(base_size=11)+
 theme(legend.position="none",strip.text=element_text(size=12,face="bold"),panel.spacing=grid::unit(1.3,"lines"),plot.caption=element_text(hjust=0,size=9),plot.title=element_text(size=16),plot.margin=margin(12,12,12,12))
save_fig(p,"Reconciled_cohort_flow",14.5,11)
m<-read.csv(file.path(out,"repeat_admission_metric_intervals.csv"));pop<-read.csv(file.path(out,"repeat_admission_populations.csv"))
d<-merge(m[m$metric=="auroc",],pop,by=c("window_hours","dataset","analysis"))
d$label<-factor(d$analysis,levels=c("equal_patient_weight","one_record_per_patient","primary"),labels=c("Equal patient weight","One record per patient","All eligible records"))
d$panel<-factor(paste0(d$window_hours,"h / ",ifelse(d$dataset=="MIMIC_internal_test","MIMIC internal test","eICU external validation")),
 levels=c("24h / MIMIC internal test","24h / eICU external validation","6h / MIMIC internal test","6h / eICU external validation"))
d$count<-paste0("N=",d$records,"; P=",d$unique_patients)
p<-ggplot(d,aes(estimate,label,color=analysis))+geom_errorbar(aes(xmin=lower,xmax=upper),orientation="y",width=.15,linewidth=.6)+geom_point(size=2.2)+
 geom_vline(xintercept=1.015,color="grey85",linewidth=.3)+geom_text(aes(x=1.19,label=count),hjust=1,color="#444444",size=2.7)+
 facet_wrap(~panel,ncol=2)+scale_x_continuous(limits=c(.5,1.2),breaks=c(.5,.6,.7,.8,.9,1))+
 scale_color_manual(values=c(primary="#222222",one_record_per_patient="#2166AC",equal_patient_weight="#008837"))+
 labs(title="Repeated-admission sensitivity of frozen XGBoost validation",x="AUROC (95% conditional cluster-bootstrap interval)",y=NULL,
 caption="N: evaluated records; P: unique patients. One-record selection uses a fixed outcome-blind hash priority, not admission chronology.\nEqual weighting gives each patient total weight 1 across their eligible records; raw record and event counts remain unchanged.\nModels and patient partitions are frozen. Selection is not repeated in bootstrap, so intervals omit selection and development uncertainty.")+
 theme_bw(base_size=11)+theme(legend.position="none",panel.grid.minor=element_blank(),panel.grid.major.y=element_blank(),plot.caption=element_text(hjust=0,size=9))
save_fig(p,"Repeat_admission_AUROC",12.8,6.5)
