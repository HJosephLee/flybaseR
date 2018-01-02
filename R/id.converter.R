#' Converting or Updating FlyBase IDs or Symbols
#'
#' The function takes FlyBase IDs (e.g. FBgn#######) or Gene symbols as an input, and converts it into updated IDs or symbols.
#' The function accesses FlyBase, so requires internet-connection.  FlyBase ID inputs are bundled as 1,000.  Symbols are 100.
#' FlyBase IDs for genes that are split into multiple genes will be concatenated with two colons (::).  Genes that does not have matching IDs will be shown as "unknown".
#' Note: Certain gene symbols would appear as "unknown" even if the gene exists, and have FlyBase IDs. This is because the ID converter in FlyBase website cannot convert the gene.  For example, CG31976 cannot be converted by FlyBase, although you can find the gene from Jump2Genes. 

#' @param x a vector. FlyBase IDs or names to be converted.
#' @param output Output-types in either FlyBase IDs or symbols. "Name", "Symbol", "n", or "s" produces gene symbols as an output.  FlyBase IDs for any other characters. Default="fbid" (FlyBase IDs).  
#' @param DmelOnly if True, non-melanogaster gene IDs will be ignored.  Default = T.
#' @keywords flybase
#' @export
#' @examples
#' id.converter(x, output="symbol", DmelOnly=T)

id.converter <- function(x, output, DmelOnly){
      
      if ( !requireNamespace("rvest", quietly = T) ) {
            stop("'rvest' package should be installed.", call. = F)
      }
      
      require(rvest)
      
      #######################################################
      # Handling arguments, or assigning default arguments.
      #######################################################
      
      if ( missing(x) ){
            stop("No FlyBase IDs or Gene Symbols to be converted!!", call. = F)
      }
      
      if ( missing(output) ){
            output.type <- "fbid"
            message("FlyBase IDs or Gene Symbols will be updated to the most recent version of FlyBase IDs..")
      } else {
            if ( toupper(output) == "N" | toupper(output) == "NAME" | toupper(output) == "S" | toupper(output) == "SYMBOL" ) {
                  output.type <- "symbol"
                  message("FlyBase IDs or Gene Symbols will be updated to the most recent version of Gene Symbols..")
            } else {
                  output.type <- "fbid"
                  message("FlyBase IDs or Gene Symbols will be updated to the most recent version of FlyBase IDs..")
            }
      }
      
      
      if ( missing(DmelOnly) ){
            DmelOnly <- T
            message("non-melanogaster IDs and Symbols will be ignored.")
      }
      
      
      #################################################
      # ID conversion to the current FBIDs or symbols
      #################################################
      
      session <- html_session("http://flybase.org/convert/id")
      form <- html_form(session)[[2]]
      
      if ( length(x[ !grepl("FBgn", x)]) >= 100 ){ bundle <- 100 } else { bundle <- 1000 }
       
      for (i in 1:ceiling(length(x)/bundle)){
            
            if (i != ceiling(length(x)/bundle)){ temp.x <- as.character(x[ (bundle*(i-1)+1):(bundle*i) ])
            } else {
                  temp.x <- as.character(x[ (bundle*(i-1)+1):length(x) ])
            }
            
            message(paste("Processing ", prettyNum(min((bundle*i), length(x)), big.mark=",", big.interval=3), " genes out of ", prettyNum(length(x), big.mark=",", big.interval=3), sep=""))
            
            form <- set_values(form, ids = paste(as.character(temp.x), collapse = "\n"))
            conversion.table <- html_table(suppressMessages(submit_form(session, form)))[[1]]
            conversion.table <- conversion.table[-1, ]
            colnames(conversion.table) <- c("submitted", "current", "converted", "symbols")
            conversion.table$submitted <- gsub(" - unknown ID", "", conversion.table$submitted, fixed=T)
            conversion.table$converted <- gsub("^.+ - unknown ID", "unknown", conversion.table$converted, perl=T)
            conversion.table$symbols <- gsub("^.+ - unknown ID", "unknown", conversion.table$symbols, perl=T)
            
            
            if ( DmelOnly == T ){ 
                  conversion.table <- conversion.table[ grepl("^D.+\\\\", conversion.table$symbols, perl=T) == F, ] # e.g. "Dyak\"
            }
            
            if ( output.type == "fbid" ){
                  temp.result <- sapply(as.character(temp.x), function(y){ paste(conversion.table[ y == conversion.table$submitted, 3], collapse = "::") }, simplify = T, USE.NAMES = F)
            } else {
                  temp.result <- sapply(as.character(temp.x), function(y){ paste(conversion.table[ y == conversion.table$submitted, 4], collapse = "::") }, simplify = T, USE.NAMES = F)
            }
            
            temp.result[ temp.result == "" ] <- "unknown"
            temp.result <- gsub("unknown::", "", gsub("::unknown", "", temp.result, fixed = T), fixed=T) # when more than 2 genes are concatenaterd with an unknown, remove the unknown tag. 
            
            if ( i == 1 ){ result <- temp.result } else { result <- c(result, temp.result) }
            
      } # for i
      
      
      ############################
      # Reporting
      ############################
      
      message(paste("Total # of ID mismatches : ", length(setdiff(x, result)), " (will show no match if ID -> symbol conversion)" , sep=""))
      message(paste("Total # unknowns : ", length(result[ result == "unknown" ]), sep=""))
      message(paste("Total # split IDs : ", length(result[ grepl("::", result) ]), sep=""))
      message(paste("Total # duplicated or merged IDs : ", length(result[ duplicated(result) & result != "unknown" ]), sep=""))
      message("Unknown IDs were updated into 'unknown'.  If a gene split in the updated gene model, split IDs were concatenated with colons (::).")
      message("If a gene has been merged in the updated gene model, they were duplicated in the final result")
      message("Done.")
      
      return(result)
      
}
