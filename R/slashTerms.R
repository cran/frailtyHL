slashTerms <-
structure(function (x) 
{
    if (!("/" %in% all.names(x))) 
        return(x)
    if (x[[1]] != as.name("/")) 
        stop("unparseable formula for grouping factor")
    list(slashTerms(x[[2]]), slashTerms(x[[3]]))
}, source = c("function (x) ", "{", "    if (!(\"/\" %in% all.names(x))) ", 
"        return(x)", "    if (x[[1]] != as.name(\"/\")) ", "        stop(\"unparseable formula for grouping factor\")", 
"    list(slashTerms(x[[2]]), slashTerms(x[[3]]))", "}"))
