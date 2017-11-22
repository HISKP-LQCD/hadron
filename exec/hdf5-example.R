library(rhdf5)

## assume you have an hdf5 file "h5file"

## find out the hdf5 datasets in this file

file <- "h5file"

h5ls(file)

## decide which one you need, say

dsname <- "correlator"

## then you may

x <- readbinarycf(files=file, T=96, hdf5format=TRUE, hdf5name=dsname,
                  hdf5index=c(1,2)) 

## the data types used in our contraction code are complex. with
## hdf5index we access real and imaginary part (1,2), respectively.
## per default

x$cf

## contains the real part and

x$icf

## the imaginary part
## You may change this with

x <- readbinarycf(files=file, T=96, hdf5format=TRUE, hdf5name=dsname,
                  hdf5index=c(2,1))
