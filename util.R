
# Pipeline utilities

library(tidyverse)

# Arbitrary code in a pipeline
f <- \(a,b) eval.parent(substitute(\(.)b))(a)

# Special assignments in a pipeline
set <- \(a,b) eval.parent(substitute(\(.){b;.}))(a)



# Reporting utilities

nice_name <- function(...) {
    paste0(...) |> str_replace_all("[ :/\\\\]","_")
}

give_figure <- function(f, name, width, height, res=200) {
    filename <- paste0(out_prefix,nice_name(name))
    
    png(paste0(out_dir,"/",filename,".png"), width=width,height=height,units="in",res=res)
    f()
    dev.off()
    
    svg(paste0(out_dir,"/",filename,".svg"), width=width,height=height)
    f()
    dev.off()
    
    cat(paste0("<b>",name,"</b> "))
    cat(paste0("<a href=\"",filename,".png\">[PNG]</a>"))
    cat(paste0("<a href=\"",filename,".svg\">[SVG]</a>\n"))
    cat(paste0("<br><img src=\"",name,".png\" width=",width*100,">\n\n"))
}

give_table <- function(df, name, what="rows") {
    filename <- paste0(out_prefix,nice_name(name))
    write_csv(df, paste0(out_dir,"/",filename,".csv"))
    cat(paste0("* <a href=\"",filename,".csv\">", name, "</a> (csv, ", nrow(df), " ", what, ")\n\n"))
}
