#' Converting or Updating FlyBase IDs or Symbols
#'
#' The function takes FlyBase IDs (e.g. FBgn0000003) or Gene symbols as an input, and converts it into updated IDs or symbols using the FlyBase ID converter (web).
#' The function accesses FlyBase, so requires internet-connection.  FlyBase ID inputs are bundled as 1,000.  100 for symbols.
#' FlyBase IDs for genes that are split into multiple genes will be concatenated with two colons (::).  Genes that does not have matching IDs will be shown as "unknown".
#' Certain gene symbols would appear as "unknown" even if the gene exists, and have FlyBase IDs. This is because the ID converter in FlyBase website cannot convert the gene.  For example, CG31976 cannot be converted by FlyBase, although you can find the gene from the gene report.  Setting diehard.symbols = T will look for gene report pages of such unconvertible genes one by one.  The process is essentially slow because it accesses FlyBase for each gene. 

#' @param x a vector. FlyBase IDs or names to be converted.
#' @param symbols Logical.  If TRUE, the output will be gene symbols, rather than FlyBase IDs.  Default = F
#' @param bundle.size Numeric.  The number of FlyBase IDs or symbols to be submitted to FlyBase at once. Default is 1,000 if there are less than 100 symbols; 100 if more than 1,000 symbols.  Reduce the number down if Timeout error occurs.
#' @param DmelOnly Logical.  If TRUE, non-melanogaster gene IDs will be ignored.  Default = T.
#' @param polite.access Numeric.  Intervals between FlyBase access for each bundle as seconds.  Default = 0.
#' @param diehard.symbols Logical.  If TRUE, ntervals between FlyBase access for each bundle as seconds.  Default = 0.
#' @param convert.into "genes", "transcripts", or "polypeptides". "g", "t", or "p" is also possible. If missing, the IDs will be updated to the most recent IDs only.
#' @keywords flybase
#' @export
#' @examples
#' id.converter(x, symbols = T)
#' id.converter(x, bundle.size = 50, be.polite = 10, convert.into = "transcripts")
#' id.converter(x, symbols = T, bundle.size = 50, diehard.symbols = T)


id.converter <- function(x, symbols, bundle.size, DmelOnly, polite.access, diehard.symbols, convert.into){
      
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
      
      if ( missing(symbols) ){
            symbols <- F
      } 
      
      if (symbols == F){
            message("The result will have the most recent version of FlyBase IDs.")
      } else {
            message("The result will have the most recent version of Gene Symbols.") 
      }

      if ( missing(bundle.size) ){
            if ( length(x[ !grepl("FBgn", x)]) >= 100 ){ bundle.size <- 100 } else { bundle.size <- 1000 }
      }
      
      if ( missing(DmelOnly) ){
            DmelOnly <- T
            message("non-melanogaster IDs and Symbols will be ignored.")
      }
      
      if ( missing(polite.access) ){
            polite.access <- 0
      }
      
      if ( missing(diehard.symbols) ){
            diehard.symbols <- F
      }
      
      
      #################################################
      # ID conversion to the current FBIDs or symbols
      #################################################
      
      session <- html_session("http://flybase.org/convert/id")
      form.original <- html_form(session)[[2]]
      
      for (i in 1:ceiling(length(x)/bundle.size)){
            
            if (i != ceiling(length(x)/bundle.size)){ temp.x <- as.character(x[ (bundle.size*(i-1)+1):(bundle.size*i) ])
            } else {
                  temp.x <- as.character(x[ (bundle.size*(i-1)+1):length(x) ])
            }
            
            message(paste("Processing ", prettyNum(min((bundle.size*i), length(x)), big.mark=",", big.interval=3), " genes out of ", prettyNum(length(x), big.mark=",", big.interval=3), sep=""))
            
            if ( missing(convert.into) | toupper(convert.into) %in% c( "G", "GENE", "GENES", "T", "TRANSCRIPT", "TRANSCRIPTS", "RNA", "P", "POLYPEPTIDE", "POLYPEPTIDES", "PROTEIN", "PROTEINS" ) == F ){
                  form <- set_values(form.original, ids = paste(as.character(temp.x), collapse = "\n"))
            } else {
                  if ( toupper(convert.into) %in% c( "G", "GENE", "GENES" ) ){
                        form <- set_values(form.original, ids = paste(as.character(temp.x), collapse = "\n"), mode = "convert", convert = "fbgn")
                  }
                  if ( toupper(convert.into) %in% c( "T", "TRANSCRIPT", "TRANSCRIPTS", "RNA" ) ){
                        form <- set_values(form.original, ids = paste(as.character(temp.x), collapse = "\n"), mode = "convert", convert = "fbtr")
                  }
                  if ( toupper(convert.into) %in% c( "P", "POLYPEPTIDE", "POLYPEPTIDES", "PROTEIN", "PROTEINS" )){
                        form <- set_values(form.original, ids = paste(as.character(temp.x), collapse = "\n"), mode = "convert", convert = "fbpp")
                  }
            }
            
            
            conversion.table <- html_table(suppressMessages(submit_form(session, form)))[[1]]
            conversion.table <- conversion.table[-1, ]
            colnames(conversion.table) <- c("submitted", "current", "converted", "symbols")
            conversion.table$submitted <- gsub(" - unknown ID", "", conversion.table$submitted, fixed=T)
            conversion.table$converted <- gsub("^.+ - unknown ID", "unknown", conversion.table$converted, perl=T)
            conversion.table$symbols <- gsub("^.+ - unknown ID", "unknown", conversion.table$symbols, perl=T)
            
            
            if ( DmelOnly == T ){ 
                  conversion.table <- conversion.table[ grepl("^D.+\\\\", conversion.table$symbols, perl=T) == F, ] # e.g. "Dyak\"
            }
            
            if ( symbols == F ){
                  temp.result <- sapply(as.character(temp.x), function(y){ paste(conversion.table[ y == conversion.table$submitted, 3], collapse = "::") }, simplify = T, USE.NAMES = F)
            } else {
                  temp.result <- sapply(as.character(temp.x), function(y){ paste(conversion.table[ y == conversion.table$submitted, 4], collapse = "::") }, simplify = T, USE.NAMES = F)
            }
            
            temp.result[ temp.result == "" ] <- "unknown"
            temp.result <- gsub("unknown::", "", gsub("::unknown", "", temp.result, fixed = T), fixed=T) # when more than 2 genes are concatenaterd with an unknown, remove the unknown tag. 
            
            if ( i == 1 ){ result <- temp.result } else { result <- c(result, temp.result) }
            
            # being polite.
            Sys.sleep(time=as.numeric(polite.access))
            
      } # for i
      
      
      #################################
      # Converting the diehard symbols
      #################################
      
      if ( diehard.symbols == T & length(result[ result == "unknown"]) > 0 ){
            
            message("Checking FlyBase gene reports for genes that are not converted by the FlyBase ID converter.  This process is very slow and takes a long time!")
            temp.df <- data.frame( symbols = x[ which(result == "unknown")], id = NA, newsymbols = NA, stringsAsFactors = F)
            
            session2 <- html_session("http://flybase.org")
            form2 <- html_form(session2)[[1]]
            
            for (j in 1:nrow(temp.df)){
                  temp.df[j, 2] <- suppressMessages(submit_form(session2, set_values(form2, context = temp.df[j,1])))$url
                  
                  if ( j %% 10 == 0 | j == nrow(temp.df)){
                        message(paste("Processing ", prettyNum(j, big.mark=",", big.interval=3), " genes out of ", prettyNum(nrow(temp.df), big.mark=",", big.interval=3), sep=""))
                  }
                  
            } # j
            
            temp.df$id <- gsub("^.+FBgn", "FBgn", temp.df$id, perl = T)
            temp.df$id[ !grepl("FBgn", temp.df$id, fixed = T) ] <- "unknown"
            
            # dealing with the four most notorious symbols that FlyBase cannot handle easily. 
            
            temp.df[ temp.df$symbols == "E2f", 2 ] <- "FBgn0011766" # e2f indicates E2f1, but FlyBase cannot decide if it is E2f1 or E2f2 when only "E2f" is entered with a capital E. 
            temp.df[ temp.df$symbols == "dm", 2 ] <- "FBgn0262656" # Myc
            temp.df[ temp.df$symbols == "th", 2 ] <- "FBgn0260635" # Diap1 = thread 
            temp.df[ temp.df$symbols == "W", 2 ] <- "FBgn0003997" # hid.  This gene used be called as Wrinkled. 
            
            if (symbols == F){ result[ which(result == "unknown") ] <- temp.df$id 
            
            } else {
                  
                  form3 <- html_form(session)[[2]]
                  form3 <- set_values(form3, ids = paste(temp.df$id, collapse = "\n"))
                  conversion.table2 <- html_table(suppressMessages(submit_form(session, form3)))[[1]]
                  conversion.table2 <- conversion.table2[-1, ]
                  colnames(conversion.table2) <- c("submitted", "current", "converted", "symbols")
                  conversion.table2$submitted <- gsub(" - unknown ID", "", conversion.table2$submitted, fixed=T)
                  conversion.table2$converted <- gsub("^.+ - unknown ID", "unknown", conversion.table2$converted, perl=T)
                  conversion.table2$symbols <- gsub("^.+ - unknown ID", "unknown", conversion.table2$symbols, perl=T)
                  
                  if ( DmelOnly == T ){ 
                        conversion.table2 <- conversion.table2[ grepl("^D.+\\\\", conversion.table2$symbols, perl=T) == F, ] # e.g. "Dyak\"
                  }
                  
                  temp.df$newsymbols <- sapply( temp.df$id, function(y){ paste(conversion.table2[ y == conversion.table2$submitted, 4], collapse = "::") }, simplify = T, USE.NAMES = F)
                  result[ which(result == "unknown") ] <- temp.df$newsymbols
                  
            }
      }
      
      
      
      ############################
      # Reporting
      ############################
      
      message(paste("Total # of ID mismatches : ", length(setdiff(x, result)), " (will show no match if ID -> symbol or gene -> protein conversion)" , sep=""))
      message(paste("Total # unknowns : ", length(result[ result == "unknown" ]), sep=""))
      message(paste("Total # split IDs : ", length(result[ grepl("::", result) ]), sep=""))
      message(paste("Total # duplicated or merged IDs : ", length(result[ duplicated(result) & result != "unknown" ]), sep=""))
      message("Unknown IDs were updated into 'unknown'.  If a gene split in the updated gene model, split IDs were concatenated with colons (::).")
      message("If a gene has been merged in the updated gene model, they were duplicated in the final result")
      message("Done.")
      
      return(result)
      
}
