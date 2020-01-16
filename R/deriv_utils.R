#' @title create list of chains of displacements
#' Multilpe covariant displacements, when applied in order, form
#' a list of displacments. Each consists of a direction and a dimension.
#' @param max_depth Positive integer, number of displacement combinations
#'                  to construct.
#' @param dims Integer vector, which lattice dimensions to consider. Default 0:3
#' @param dirs Integer vector, which displacement directions to consider. Default
#'             forward and backward <-> c(0,1)
#' @return List of data frames, each with columns 'dim' and 'dir' of 'max_depth' rows.
create_displ_chains <- function(max_depth, dims=c(0:3), dirs=c(0,1) ){
  stopifnot(max_depth > 0)
  # there are 4 dimensions, 2 directions
  # -> factor of 8 per level
  # -> in general, length(dims)*length(dirs) per level
  n_chains <- (length(dims)*length(dirs))^max_depth

  displ <- list()
  for( i in 1:max_depth ){
    displ[[length(displ)+1]] <- dims
  }
  for( i in 1:max_depth ){
    displ[[length(displ)+1]] <- dirs
  }

  displ_table <- as.matrix(expand.grid(displ))
  colnames(displ_table) <- NULL
  rownames(displ_table) <- NULL
  dimnames(displ_table) <- NULL
 
  displ_chains <- list()
  for( i in 1:n_chains ){
    displ_chains[[i]] <- data.frame(dim=displ_table[i, 1:max_depth],
                                    dir=displ_table[i, (max_depth+1):(2*max_depth)])
  } 
 
  return(displ_chains) 

}

