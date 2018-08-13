makeInteraction <- function(x)
### from a list of length 2 return recursive interaction terms
{
    if (length(x) < 2) return(x)
    trm1 <- makeInteraction(x[[1]])
    trm11 <- if(is.list(trm1)) trm1[[1]] else trm1
    list(substitute(foo:bar, list(foo=x[[2]], bar = trm11)), trm1)
}
