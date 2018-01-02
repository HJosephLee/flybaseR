#' Open FlyBase gene reports using a web browser
#'
#' The function takes FlyBase IDs (e.g. FBgn#######) or Gene symbols as an input, and opens FlyBase gene report from the default web browser (e.g. Jump2Gene).
#' @param x a character or vector. FlyBase IDs or symbols. 
#' @param n The maximum number of genes to be opened.  Default = 10.
#' @keywords flybase
#' @export
#' @examples
#' j2g("e2f")
#' j2g(c("e2f", "e2f2", "dsx"))
#' j2g(c("FBgn0000504"))


j2g <- function(x, n){
      
      if ( !requireNamespace("rvest", quietly = T) ) {
            stop("'rvest' package should be installed.", call. = F)
      }
      
      if ( missing(x) ){
            stop("No FlyBase IDs or Gene Symbols entered.", call. = F)
      }
      
      if ( missing(n) ){
            n <- 10
      }
      
      if (length(x) > n){
            stop("You are trying to open too many tabs on your browser.  You can set the 'n' parameter if you really want to open +10 tabs." , call. = F)
      }
      
      require(rvest)
      
      session <- html_session("http://flybase.org")
      form <- html_form(session)[[1]]
      
      for (i in 1:as.character(x)){
            result <- suppressMessages(submit_(session, set_values(form, context = i)))
            browseURL(url = result$url)
      }
      
}
