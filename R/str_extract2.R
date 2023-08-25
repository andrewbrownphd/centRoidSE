#Helper function that prevents needing an entire package dependency
#for a single function
str_extract2 <- function(x) {
  hee <- (strsplit(as.character(x), split=".", fixed=T))[[1]][2]
  return(ifelse(is.na(hee), NA, paste(".",hee,sep="")))
}