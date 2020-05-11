#' @title Generate table of momentum component combinations
#' @param psqmax Integer, maximum p^2 = px^2 + py^2 + pz^2 to be included in momentum list
#'
#' @return
#' Returns a \link{data.frame} with all possible momentum combinations.
#' 
#' @export
mom_combinations <- function(psqmax){
  pmax <- ceiling(sqrt(psqmax))+1
  pseq <- (-pmax):pmax
  moms <- NULL
  for( px in pseq ){
    for( py in pseq ){
      for( pz in pseq ){
        if( px^2 + py^2 + pz^2 <= psqmax )
          moms <- rbind(moms,
                        data.frame(px=px, py=py, pz=pz))
      }
    }
  }
  srt_idcs <- order(moms$px^2 + moms$py^2 + moms$pz^2,
                    abs(moms$px), abs(moms$py), abs(moms$pz))
  moms <- moms[srt_idcs,]
  rownames(moms) <- NULL

  return( moms )
}
