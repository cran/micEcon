translogCheckMono <- function( xNames, data, coef, increasing = TRUE,
   strict = FALSE, quadHalf = TRUE, dataLogged = FALSE,
   tol = 2 * .Machine$double.eps ) {

   result <- list()

   deriv <- translogDeriv( xNames = xNames, data = data, coef = coef,
      quadHalf = quadHalf, dataLogged = dataLogged )$deriv

   if( increasing ) {
      if( strict ) {
         result$exog <- deriv > 0
      } else {
         result$exog <- deriv >= - tol
      }
   } else {
      if( strict ) {
         result$exog <- deriv < 0
      } else {
         result$exog <- deriv <= tol
      }
   }

   result$obs <- rowSums( !result$exog ) == 0
   result$increasing <- increasing
   result$strict     <- strict

   class( result ) <- "translogCheckMono"
   return( result )
}

