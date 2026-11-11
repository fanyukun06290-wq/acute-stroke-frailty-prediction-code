Sys.setlocale("LC_CTYPE","Chinese (Simplified)_China.utf8")
.libPaths(unique(c(Sys.getenv("STROKE_R_LIBRARY", ""),.libPaths())))
suppressPackageStartupMessages(library(ggplot2))
root<-file.path(getwd(),"revision_2026-11-04");out<-file.path(root,"dca_threshold_results")
d<-read.csv(file.path(out,"decision_curve_intervals.csv"));a<-read.csv(file.path(out,"default_strategy_intervals.csv"))
d<-d[d$dataset=="eICU_external_validation",];a<-a[a$dataset=="eICU_external_validation",]
window<-function(h) factor(paste0(h,"-hour landmark"),levels=c("24-hour landmark","6-hour landmark"))
d$window<-window(d$window_hours);a$window<-window(a$window_hours)
raw<-d[d$metric=="net_benefit" & d$model %in% c("xgboost","random_forest","elastic_net"),]
labels<-c(xgboost="XGBoost",random_forest="Random forest",elastic_net="Elastic net")
raw$label<-unname(labels[raw$model])
all<-data.frame(window=a$window,threshold=a$threshold,estimate=a$all_flagged,label="Flag all")
none<-data.frame(window=a$window,threshold=a$threshold,estimate=0,label="Flag none")
lines<-rbind(raw[,c("window","threshold","estimate","label")],all,none)
lines$label<-factor(lines$label,levels=c("XGBoost","Random forest","Elastic net","Flag all","Flag none"))
theme<-theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),legend.position="bottom",
    plot.caption=element_text(hjust=0,size=9),panel.spacing=grid::unit(1,"lines"))
fig<-ggplot()+geom_ribbon(data=raw[raw$model=="xgboost",],aes(threshold,ymin=lower,ymax=upper),fill="#2166AC",alpha=.14)+
    geom_line(data=lines,aes(threshold,estimate,color=label,linetype=label),linewidth=.65)+
    geom_vline(xintercept=.10,color="#999999",linetype="dotted")+facet_wrap(~window,nrow=1)+
    scale_color_manual(values=c("XGBoost"="#2166AC","Random forest"="#1B7837","Elastic net"="#762A83","Flag all"="#333333","Flag none"="#777777"))+
    scale_linetype_manual(values=c("XGBoost"="solid","Random forest"="solid","Elastic net"="solid","Flag all"="dashed","Flag none"="dotted"))+
    coord_cartesian(xlim=c(.01,.50),ylim=c(-.04,.20),expand=FALSE)+
    labs(title="Exploratory decision curves in eICU",x="Risk threshold probability",y="Net benefit per record",color=NULL,linetype=NULL,
        caption=paste("XGBoost ribbon: 95% pointwise cluster-bootstrap intervals (not a simultaneous band). Vertical marker: inherited 0.10 operating point.",
            "\nY-axis zoom: -0.04 to 0.20; full numerical values are retained. The two landmarks concern different conditional populations.",
            "\nFlagging represents hypothetical further risk review, not a validated intervention, treatment benefit or life-support policy."))+theme
dir.create(file.path(out,"figures"),showWarnings=FALSE)
ggsave(file.path(out,"figures/Decision_curve_primary_models.png"),fig,width=11.5,height=5.7,dpi=300)
ggsave(file.path(out,"figures/Decision_curve_primary_models.svg"),fig,width=11.5,height=5.7,device=grDevices::svg)
gain<-d[d$model=="xgboost" & d$metric=="difference_vs_best_default",]
fig<-ggplot(gain,aes(threshold,estimate))+geom_ribbon(aes(ymin=lower,ymax=upper),fill="#2166AC",alpha=.16)+
    geom_line(color="#2166AC",linewidth=.7)+geom_hline(yintercept=0,linetype="dashed",color="#777777")+
    geom_vline(xintercept=.10,color="#999999",linetype="dotted")+facet_wrap(~window,nrow=1)+
    labs(title="XGBoost net-benefit difference from the better default strategy",x="Risk threshold probability",y="Net benefit difference per record",
        caption=paste("The better of flag-all and flag-none is recomputed within every paired bootstrap replicate.",
            "\nShading: 95% pointwise conditional intervals; grid comparisons are exploratory, not a validated threshold-selection rule.",
            "\nPositive theoretical net benefit is not proof of improved clinical decisions or patient outcomes."))+theme
ggsave(file.path(out,"figures/Decision_curve_best_default_difference.png"),fig,width=11.5,height=5.4,dpi=300)
ggsave(file.path(out,"figures/Decision_curve_best_default_difference.svg"),fig,width=11.5,height=5.4,device=grDevices::svg)
cat("DCA_THRESHOLD_PLOTS_WRITTEN\n")
