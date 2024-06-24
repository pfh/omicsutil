
# Pipeline utilities

library(tidyverse)

# Arbitrary code in a pipeline
f <- \(a,b) eval.parent(substitute(\(.)b))(a)

# Special assignments in a pipeline
set <- \(a,b) eval.parent(substitute(\(.){b;.}))(a)



# Reporting utilities

give_figure <- function(f, name, width, height, res=200, out_dir="output/qmd") {
    #ggsave(paste0("output/qmd/",name,".png"), p, width=width,height=height, dpi=200, limitsize=FALSE)
    #ggsave(paste0("output/qmd/",name,".svg"), p, width=width,height=height, limitsize=FALSE)
    
    png(paste0(out_dir,"/",name,".png"), width=width,height=height,units="in",res=res)
    f()
    dev.off()
    
    svg(paste0(out_dir,"/",name,".svg"), width=width,height=height)
    f()
    dev.off()
    
    cat(paste0("Save-able images: <a href=\"",name,".png\">[PNG]</a>\n"))
    cat(paste0("<a href=\"",name,".svg\">[SVG]</a>\n"))
    cat(paste0("<br><img src=\"",name,".png\" width=",width*100,">\n\n"))
}
