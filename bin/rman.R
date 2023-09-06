#!/usr/bin/env Rscript

#args = commandArgs(trailingOnly=TRUE)

extractRd <- function (args) {
    for (file in args) {
        if (!file.exists(file)) {
            stop(paste("Error: File",file,"does not exists"))
        }
        fin  = file(file, "r")
        flag = FALSE
        while(length((line = readLines(fin,n=1)))>0) {
            if (grepl("^#' *\\\\name\\{.+\\}",line)) {
                name = gsub("\\$","_",gsub("^#' *\\\\name\\{(.+)\\}.*","\\1",line))
                print(paste("name = '",name,"'",sep=''))
                flag = TRUE
                mandir = sub("R/.+","man/",file)
                fout = file(paste(mandir,name,".Rd",sep=""),"w")
                
            } else if (flag & (grepl('^[a-zA-Z" ].*',line) | grepl("^$",line))) {
                flag = FALSE
                close(fout) 
                next
            } 
            if (flag) {
                line = sub("^#' ?","",line)
                cat(line,"\n",file=fout)
            }
        }
        close(fin)
    }
}
main <- function (argv) {
    extractRd(argv)
}
if (sys.nframe() == 0L && !interactive()) {
    main(commandArgs(trailingOnly=TRUE))
}

# TODO: Add a parser for commonly written field names
# Idea: Define a look-up table or map, to map the field 
# name to the LaTeX field name.
# Define a look-up table to map the common field names
# tbl <- c(
#     "Name", "Alias", "Title", "Usage", "Description",
#     "Arguments", "Return", "Examples", "Value", "See also"
# )
# Convert from #'[[Nn]ame:] to '\name{', i.e. exchange first
# letter with `\` + first letter lowercase, and ':'  with '{'
# and before 