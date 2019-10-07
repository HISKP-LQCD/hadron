#args <- list(
#    ens = c("A40.24", "A40.32"),
#    all.dirs = list(
##            A40.20 = c("p0/T1u", "p1/A1g"),
#            A40.24 = c("p0/T1u", "p1/A1g", "p1/Ep1g"),
#            A40.32 = c("p0/T1u", "p1/A1g", "p1/Ep1g")
#            ),
#    pcs = list(
##            A40.20 = c(1,3),
#            A40.24 = c(2,3,2),
#            A40.32 = c(2,3,2)
#            ),
#    data.paths = c("/hiskp2/werner/pipi_I1/data/A40.24/5_fit-data/", "/hiskp2/werner/pipi_I1/data/A40.32/5_fit-data/"),
##    data.paths = c("/hiskp2/werner/pipi_I1/data/A40.20/5_fit-data/", "/hiskp2/werner/pipi_I1/data/A40.24/5_fit-data/", "/hiskp2/werner/pipi_I1/data/A40.32/5_fit-data/"),
#    output.path = "/hiskp2/werner/pipi_I1/data/A40"
##    boot.R = 150,
#    )

args <- list(
    ens = c("A30.32"),
    all.dirs = list(
#            A30.32 = c("p0/T1u", "p1/A1", "p1/E", "p2/A1", "p2/B1", "p2/B2", "p3/A1", "p3/E", "p4/A1", "p4/E")
            A30.32 = c("p1/A1", "p1/E", "p2/A1", "p2/B1", "p2/B2", "p3/A1", "p3/E", "p4/A1", "p4/E")
            ),
    pcs = list(
#            A30.32 = c(3,3,2,3,2,2,2,1,1,1)
            A30.32 = c(3,2,3,2,2,2,1,1,1)
            ),
    data.paths = c("/hiskp4/werner/pipi_I1/data/A30.32/5_fit-data/"),
    output.path = "/hiskp4/werner/pipi_I1/data/A30.32"
#    boot.R = 150,
    )
    
