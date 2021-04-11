## to be sourced from main analysis script

ggp_theme_default <- theme(

    ## panel parameters
    panel.border=element_blank(),
    panel.background=element_rect(fill="white"),
    panel.grid.major=element_line(color="white"),
    panel.grid.minor=element_line(color="white"),

    ## plot parameters
    plot.background=element_rect(fill="white"),
    plot.margin=margin(t=0.5, r=0.5, b=0.5, l=0.5, unit="lines"),
    plot.title=element_text(size=12, face="plain", hjust=0.5),
    plot.subtitle=element_text(size=10, face="plain", hjust=0.5),

    ## axis parameters
    axis.line.x.bottom=element_line(),
    axis.line.y.left=element_line(),
    axis.title.x=element_text(size=10, face="plain"),
    axis.title.y=element_text(size=10, face="plain"),
    axis.text.x=element_text(size=9, face="plain"),
    axis.text.y=element_text(size=9, face="plain"),
    # axis.ticks.x=element_blank(),
    # axis.ticks.y=element_blank(),

    ## legend parameters
    legend.key=element_blank(),
    legend.key.size=unit(0.3, "cm"),
    legend.key.width=unit(0.3, "cm"),
    legend.position="none",
    legend.title=element_blank(),
    legend.text=element_text(size=10, margin=margin(t=0.2)),
    legend.spacing.x=unit(0.2, "cm"),
    legend.spacing.y=unit(0.2, "cm"),

    ## facet parameters
    strip.background=element_rect(fill="black"),
    strip.text=element_text(colour="white")

)
