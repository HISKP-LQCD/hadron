bla <- read.table("scatteringlenght.dat")
L <- bla$V2
mps <- bla$V4
fps <- bla$V5
rpi <-  mps^2/(4.0*pi*fps)^2*g1( L*mps )

data.frame((1-0.5*rpi), (1+2*rpi))
