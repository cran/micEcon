## -------- system formulas -------------------
snqProfitSystem <- function(nNetput,nFix) {
   system <- list()
   for(i in 1:nNetput) {
      system[[i]] <- paste("q",as.character(i)," ~",sep="")
      for(j in 1:nNetput) {
         if(j==1) {
            system[[i]] <- paste(system[[i]]," pr", as.character(j),sep="")
         } else {
            system[[i]] <- paste(system[[i]]," + pr", as.character(j),sep="")
         }
      }
      for(j in 1:nNetput) {
         for(k in 1:nNetput) {
            system[[i]] <- paste(system[[i]], " + pq", as.character(i), ".",
                                 as.character(j), ".", as.character(k), sep="")
         }
      }
      if( nFix > 0 ) {
         for(j in 1:nFix) {
            system[[i]] <- paste(system[[i]]," + f", as.character(j),sep="")
         }
         for(j in 1:nFix) {
            for(k in 1:nFix) {
               system[[i]] <- paste(system[[i]], " + fq", as.character(i), ".",
                                    as.character(j), ".", as.character(k), sep="")
            }
         }
      }
      system[[i]] <- as.formula(system[[i]])
   }
   return( system )
}
