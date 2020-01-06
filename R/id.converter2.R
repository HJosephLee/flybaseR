#'  Updating FlyBase IDs to a certain version.
#'
#' The function takes FlyBase IDs (e.g. FBgn#######) as an input, and converts it into certain versions of IDs. This function is not able to handle gene symbols.
#' The function accesses the FlyBase FTP site, so requires internet-connection.
#' FlyBase IDs for genes that are split into multiple genes will be concatenated with two colons (::).  Genes that does not have matching IDs will be shown as "unknown".

#' @param x a vector. FlyBase IDs.
#' @param version FlyBase ID version for the updated result. It should be either version numbers (e.g. 6.12) or FlyBase Release Numbers (e.g. FB2017_04).  The default is "current".  
#' @param thread The function itself is slow, so in order to speed up you can use multiple CPU threads for parallelization. Default : thread = 1.
#' @keywords flybase
#' @export
#' @examples
#' id.converter2(x, version=6.12)
#' id.converter2(x, version="FB2016_04", thread=4)

id.converter2 <- function(x, version, thread){
      
      
      if ( !requireNamespace("RCurl", quietly = T) ) {
            stop("RCurl should be installed.", call. = F)
      }
      
      require(RCurl)
      
      #######################################################
      # Handling arguments, or assigning default arguments.
      #######################################################
      
      if ( missing(x) ){
            stop("No Flybase IDs to be converted!!", call. = F)
      }
      
      if ( missing(version) ){
            version <- "current"
            message("IDs will be updated to the most recent Flybase version..")
      }
      
      if ( missing(thread) ){
            thread <- 1
            message("Running on a single core..")
      }
      
      ##########################
      # Handling versions
      ##########################
      
      version.table <- unlist( strsplit( getURL("ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/", dirlistonly=T), split = "\n" ) )
      
      version <- gsub("^.+FB20", "FB20", version.table[ grepl(as.character(version), version.table) ])
      if (length(version) > 1) { 
            stop(paste("Multiple versions for the query :", version, collapse = " ")) 
      }
      
      if (grepl("FB", version) == F & version != "current"){
            stop("fbid.updater2 does not support Pre-2006 genomes")
      }
      
      if (version != "current"){
            version <- substring(version, nchar(version)-8, nchar(version))
      }
      
      #########################################
      # Get the gene ID list file from FlyBase
      #########################################
      
      message("Downloading fbgn_annotation_ID file from Flybase..")
      
      temp <- tempfile()
      download.file(paste("ftp://ftp.flybase.net/releases/",version, "/precomputed_files/genes/fbgn_annotation_ID.tsv.gz", sep=""), temp, method="curl", quiet = T)
      conversion.table <- read.delim( gzfile(temp), blank.lines.skip = T, comment.char = "#", header=F, stringsAsFactors = F)
      conversion.table <- conversion.table[ !grepl("\\", conversion.table$V1, fixed=T), ] # removing non-melanogastor genes
      
      unlink(temp, force=T)
      
      #########################
      # ID conversion
      #########################
      
      message("Updating IDs..")
      
      if (thread == 1){
            temp.x <- as.character(x)
            temp.x[ which(temp.x %in% conversion.table$V3 == F) ] <-  sapply(temp.x[ temp.x %in% conversion.table$V3 == F ], 
                                                                             function(y){ paste(conversion.table[ grepl(y, conversion.table$V4, fixed = T), 3], collapse = "::") }, simplify = T, USE.NAMES = F) # So, if a gene ID became splitted, then conmessageenate the names.
            
      } else {
            
            if ( !requireNamespace("parallel", quietly = TRUE) ) {
                  stop("'parallel' package is required for multi-threading.", call. = FALSE)
            }
            require(parallel)
            
            clusters <- makeCluster(thread)
            temp.x <- as.character(x)
            temp.x[ which(temp.x %in% conversion.table$V3 == F) ] <-  parSapply(clusters, temp.x[ temp.x %in% conversion.table$V3 == F ], 
                                                                                function(y){ paste(conversion.table[ grepl(y, conversion.table$V4, fixed = T), 3], collapse = "::") }, simplify = T, USE.NAMES = F) # So, if a gene ID became splitted, then conmessageenate the names.

            stopCluster(clusters)
            
      }
      
      temp.x[ temp.x == "" ] <- "unknown"
      
      ############################
      # Reporting
      ############################
      
      message(paste("Total # of ID mismatches : ", length(setdiff(x, temp.x)), sep=""))
      message(paste("Total # unknowns : ", length(temp.x[ temp.x == "unknown" ]), sep=""))
      message(paste("Total # split IDs : ", length(temp.x[ grepl("::", temp.x) ]), sep=""))
      message(paste("Total # duplicated or merged IDs : ", length(temp.x[ duplicated(temp.x) & temp.x != "unknown" ]), sep=""))
      message("Unknown IDs were updated into 'unknown'.  If a gene split in the updated gene model, split IDs were concatenated with colons (::).")
      message("If a gene has been merged in the updated gene model, they were duplicated in the final result")
      message("Done.")

      return(temp.x)
      
} # function
