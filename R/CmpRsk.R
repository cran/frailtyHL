CmpRsk <-
function(time, index) {
# Verify Arugments
if(any(time<0)) stop("Invalid time")
if(any(index<0, diff(sort(unique(index)))!=1)) stop("Invalid index")
# Create object
temp <- cbind(time, index)
class(temp) = "CmpRsk"
return(temp)
}
