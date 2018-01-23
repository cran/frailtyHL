p.sum <-
function(...) {
x <-list(...)
eval(parse(text=paste("x[[", 1:length(x), "]]", sep="", collapse="+")))
}
