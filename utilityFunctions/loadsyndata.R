

#' SynID <- 'syn23583548'
#' A function
#' 
#' @param SynID A vector
#' @param names A numeric
#' 
#' @return The dataframe pulled from synapse
#' @importFrom synapser synGet
#' @importFrom data.table fread
#' @export
TMT_Express_Load <- function( SynID, names ){
  #'@param  SynID a synapse ID of a proteomics csv matrix eg. syn21266454 or syn21266454
  #'@param  names number corresponding to the row that are the rownames (Use 0 if none apply )
  if( !is.character(SynID) ){
    return("Invalid Synapse ID: Input SynID is not a character string")
  }
  if( !grepl('syn', SynID) | !((nchar(SynID) == 11) | (nchar(SynID) == 10 ))  ){
    return("Invalid Synapse ID: Input SynID is not a valid synapse ID")
  }
  #if(  grepl('Benefactor not found for', as.character( try( synapser::synGetPermissions(SynID), silent = TRUE ) )) ){
  #  return("Syn ID does not exist")
  #}
  if(  'You do not have READ permission for the requested entity' %in% as.character( try( syn_temp$getPermissions(SynID), silent = TRUE ) ) ){
    return("User does not have READ access")
  }
  if( !grepl('.csv$', syn_temp$get(SynID)$path ) ){
    return("File to import must be a CSV File")
  }
  if( !is.numeric(names) ){
    return("names is not a number")
  }
  
  if( names == 0 ){
    import <- data.frame( data.table::fread( syn_temp$get(SynID)$path, header=T, sep=',' ))
  }else{
    import <- data.frame( data.table::fread( syn_temp$get(SynID)$path, header=T, sep=',' ), row.names = 1)
  }
  return( import )
}
