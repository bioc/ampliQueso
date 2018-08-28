.loadLibs <-function(packageName){
  
  if( do.call("require",list(packageName) ) ){
    print( paste(packageName, "is loaded correctly",sep=" ") )	
  } 
  else {
    print( paste("Trying to install",packageName, sep=" " ) )
    do.call("install.packages",list(packageName) )
    if( do.call("require",list(packageName) ) )
    	{}
    else	
	{
	   print("Package not found in CRAN, trying Biocondutor repository")
	   if (!requireNamespace("BiocManager", quietly=TRUE))
    	   install.packages("BiocManager")
	 #  .biocLite(packageName )
	}
    if( do.call("require",list(packageName) ) ){
        print( paste(packageName,"installed and loaded",sep=" ") )
    } 
    else {
        stop( paste("Could not install",packageName, sep=" ") )
    }
  }  


}


#escape special characters for LaTeX output
.escSpecialChars <-function(str){
  if(is.null(str))
    return("")
  
  str <-gsub("_", "\\_", str, fixed = TRUE)
  str <-gsub("&", "\\&", str, fixed = TRUE)
  str <-gsub("$", "\\$", str, fixed = TRUE)
  
  return(str)
}

.setdiff.data.frame <-
  function(A,B) A[ !duplicated( rbind(B,A) )[ -seq_len(nrow(B))] , ]