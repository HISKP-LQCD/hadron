# functions for simple error propagation

compute_square <- function(x,name=NA,debug=FALSE) {
  rval <- list( val=x$val^2, 
                 dval=2*x$val*x$dval,
                 name=name )
  if(debug) {
    print(sprintf("compute_square: %s",as.character(name)))
    print(rval)
  }
  return(rval)
}

compute_sqrt <- function(x,name=NA,debug=FALSE){
  nidx <- which( x$val < 0)
  if(length(nidx)>0){
    x["val",nidx] <- abs(x["val",nidx])
    warning(sprintf("compute_sqrt: Warning, negative value replaced by absolute value for %s!\n",name))
  }
  rval <- list( val=sqrt(x$val),
                dval=0.5*x$dval/sqrt(x$val),
                name=name )
  if(debug) {
    print(sprintf("compute_sqrt: %s",as.character(name)))
    print(rval)
  }
  return(rval)
}

compute_ratio <- function(dividend,divisor,name=NA,debug=FALSE) {
  rval <- list( val=dividend$val / divisor$val, 
                dval=sqrt( (dividend$dval/divisor$val)^2 + (divisor$dval*dividend$val/divisor$val^2)^2 ), 
                name=name )
  if(debug) {
    print(sprintf("compute_ratio: %s",as.character(name)))
    print(rval)
  }
  return(rval)
}

compute_product <- function(a,b,name=NA,debug=FALSE) {
  rval <- list( val=a$val * b$val, 
                      dval=sqrt( (a$dval*b$val)^2 + (b$dval*a$val)^2 ), 
                      name=name )
  if(debug) {
    print(sprintf("compute_product: %s",as.character(name)))
    print(rval)
  }
  return(rval)
}

compute_sum <- function(a,b,name=NA,debug=FALSE) {
  rval <- list( val=a$val + b$val,
                      dval=sqrt( a$dval^2 + b$dval^2 ),
                      name=name )
  if(debug) {
    message(sprintf("compute_sum: %s\n",as.character(name)))
    print(rval)
  }
  return(rval)
}

compute_difference <- function(pos,neg,name=NA,debug=FALSE) {
  neg$val <- -neg$val
  rval <- compute_sum(a=pos,b=neg,name=name,debug=FALSE)
  if(debug) {
    message(sprintf("compute_difference: %s\n",as.character(name)))
    print(rval)
  }
  return(rval)
}

