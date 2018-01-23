trace <-
function(x) {stopifnot(is.matrix(x), nrow(x)==ncol(x)); sum(diag(x))}
