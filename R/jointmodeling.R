jointmodeling <-
function(Model="mean",RespDist="gaussian",Link=NULL,LinPred="constant",RandDist=NULL,
Offset=NULL) {
    if (Model=="mean" && is.null(Link)) Link="identity"
    if (Model=="dispersion" && is.null(Link)) Link="log"
    res<-list(Model,Link,LinPred,RandDist,Offset)
    return(res)
}
