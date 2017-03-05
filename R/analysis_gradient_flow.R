analysis_gradient_flow <- function(path,basename,read.data=TRUE,plot=FALSE,skip=0,start=0,scale=1,dbg=FALSE) {
  if(missing(basename)){
    path.strings <- strsplit(x=path,split='/')[[1]]
    basename <- path.strings[length(path.strings)]
  }
  if(read.data) {
    raw.gradflow <- readgradflow(path=path,skip=skip)
    save(raw.gradflow,file=sprintf("%s.raw.gradflow.Rdata",basename),compress=FALSE)
  }else{
    cat(sprintf("Warning, reading data from %s.raw.gradflow.Rdata, if the number of samples changed, set read.data=TRUE to reread all output files\n",basename))
    load(sprintf("%s.raw.gradflow.Rdata",basename))
  }
  if(dbg==TRUE) print(raw.gradflow)

  t_vec <- unique(raw.gradflow$t)
  Ncol <- ncol(raw.gradflow)
  Nrow <- length(t_vec)
  # allocate some memory, for each observable, have space for the value, the error and the autocorrelation time
  # create a list with NULL rownams and sensible column names for this purpose
  subnames <- c("value","dvalue","ddvalue","tauint","dtauint")
  cnames <- colnames(raw.gradflow[,3:Ncol])
  outer.cnames <- t(outer(cnames,subnames,FUN=function(x,y) { paste(x,y,sep=".") }))
  grad.dimnames <- list()
  grad.dimnames[[1]] <- NULL
  grad.dimnames[[2]] <- c("t",as.vector(outer.cnames) )
  
  gradflow <- as.data.frame(matrix(data=NA,nrow=Nrow,ncol=(length(subnames)*(Ncol-2)+1),dimnames=grad.dimnames))
  for(i_row in 1:length(t_vec)){
    uwerr.gradflow <- apply(X=raw.gradflow[which(raw.gradflow$t==t_vec[i_row]),3:Ncol],MARGIN=2,FUN=uwerrprimary)
    summaryvec <- c(t_vec[i_row])
    for(i_col in 1:length(cnames)) {
      obs <- uwerr.gradflow[[cnames[i_col]]]
      summaryvec <- c(summaryvec, obs$value, obs$dvalue, obs$ddvalue, obs$tauint*scale, obs$dtauint*scale)
    }
    gradflow[i_row,] <- summaryvec
  }
  
  gf.scales <- list()
  
  gf.scales[["sqrt_t0_Eplaq"]] <- list(val=rep(0,3), ref.val=0.1416, ref.dval=0.0008, approx.idx=0,
                                       label="\\sqrt{ t_0^{\\mathrm{plaq}} }", obs="tsqEplaq")

  gf.scales[["sqrt_t0_Esym"]] <- list(val=rep(0,3), ref.val=0.1416, ref.dval=0.0008, approx.idx=0,
                                      label="\\sqrt{ t_0^{\\mathrm{sym}} }", obs="tsqEsym")

  gf.scales[["w0_Wsym"]] <- list(val=rep(0,3), ref.val=0.1755, ref.dval=0.0019, approx.idx=0,
                                 label="\\w_0^{\\mathrm{sym}}", obs="Wsym")

  for( i in 1:length(gf.scales) ){
    val <- gradflow[,sprintf("%s.value", gf.scales[[i]]$obs)]
    dval <- gradflow[,sprintf("%s.dvalue", gf.scales[[i]]$obs)] 
    tref <- c( approx( x=val+dval,y=gradflow$t,xout=0.3)$y, 
               approx( x=val, y=gradflow$t, xout=0.3 )$y, 
               approx( x=val-dval, y=gradflow$t, xout=0.3)$y )

    sqrt_tref <- sqrt(tref)
  }

  save(gradflow,file=sprintf("%s.result.gradflow.Rdata",basename),compress=FALSE)
   
  # find w_0 and its lower and upper values, note how x and y are reversed for this purpose
  w0sq <- c( approx(x=gradflow$Wsym.value+gradflow$Wsym.dvalue,y=gradflow$t,xout=0.3)$y, 
             approx( x=gradflow$Wsym.value, y=gradflow$t, xout=0.3 )$y, 
             approx( x=gradflow$Wsym.value-gradflow$Wsym.dvalue, y=gradflow$t, xout=0.3)$y )
 
  # find t value closest to w0sq
  w0sq_approx <- w0sq[2]
  for( i in 1:length(t_vec) ){
    if( t_vec[i] >= w0sq_approx ) { w0sq_approx <- t_vec[i]; break }
  }

  cat(sprintf("w0^2/a^2: %f (+%f -%f)\n",w0sq[2],w0sq[3]-w0sq[2],w0sq[2]-w0sq[1]))
   
  w0 <- sqrt(w0sq)
  cat(sprintf("w0/a: %f (+%f -%f)\n",w0[2],w0[3]-w0[2],w0[2]-w0[1]))

  # TODO: add database of inputs to hadron, such that one can choose from multiple values of w0
  #       to compute the lattice spacing from
  w0fm <- list(val=0.1755,dval=0.0019,name="w_0(fm)")
  a <- matrix(nrow=3,ncol=2)
  for(i in 1:length(w0)){
    a_fm <- compute_ratio(dividend=w0fm,divisor=list(val=w0[i],dval=0,name="w_0"),name="a(fm)")
    a[i,] <- c(a_fm$val,a_fm$dval)
  }
  cat(sprintf("Using w0 = %f (%f)\n",w0fm$val,w0fm$dval))
  cat("a(fm):", a[2,1], " (", a[3,1]-a[2,1], ",", a[1,1]-a[2,1]," ) ( ", a[2,2], ")\n")
  cat("a(fm) ~ ", a[2,1], "(", sqrt( 0.5*(abs(a[3,1]-a[2,1])+abs(a[1,1]-a[2,1]))^2 + a[2,2]^2 ), ")\n");

  
  if(plot) {
    tikzfiles <- tikz.init(basename=sprintf("%s.gradflow",basename),width=4,height=4)
    # set up plot
    plot(x=gradflow$t, y=gradflow$Wsym.value,
         type='n',xlim=c(0,1.25*w0sq[2]),ylim=c(0.0,0.4),
         xlab="$t/a^2$",ylab="$W(t)$",las=1)
    # draw errorband
    poly.col <- rgb(red=1.0,green=0.0,blue=0.0,alpha=0.6)
    poly.x <- c(gradflow$t,rev(gradflow$t))
    poly.y <- c(gradflow$Wsym.value+gradflow$Wsym.dvalue,rev(gradflow$Wsym.value-gradflow$Wsym.dvalue))
    polygon(x=poly.x,y=poly.y,col=poly.col)
    abline(h=0.3)
    abline(v=w0sq)
    
    # plot MD history of Wsym at w0
    plot(y=raw.gradflow[which(raw.gradflow$t==w0sq_approx),"Wsym"],
         x=start+c( 0:( length(raw.gradflow[which(raw.gradflow$t==w0sq_approx),"Wsym"]) - 1 ) )*scale,
         type='l',lwd=3,
         main="",xlab="$N_\\mathrm{conf}$",ylab=sprintf("$W\\left( t/a^2 = %s \\right)$",w0sq_approx),las=1)
    
    tikz.finalize(tikzfiles)
  }
  
}

