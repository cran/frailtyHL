print.hlike <-
function(object) {
print(object$beta); cat("\n"); print(object$theta);
cat("\n",ifelse(object$converge==1, "Successfully Converged",
"Failed to Converge"),"\n")
}
