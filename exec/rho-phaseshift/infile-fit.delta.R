args <- list(
    ens = c("A40.24", "A40.32"),
    all.dirs = list(
#            A40.20 = c("p0/T1u", "p1/A1g"),
            A40.24 = c("p0/T1u", "p1/A1g", "p1/Ep1g"),
            A40.32 = c("p0/T1u", "p1/A1g", "p1/Ep1g")
            ),
    pcs = list(
#            A40.20 = c(1,3),
            A40.24 = c(2,3,2),
            A40.32 = c(2,3,2)
            ),
    data.paths = c("/hiskp2/werner/pipi_I1/data/A40.24/5_fit-data/", "/hiskp2/werner/pipi_I1/data/A40.32/5_fit-data/"),
#    data.paths = c("/hiskp2/werner/pipi_I1/data/A40.20/5_fit-data/", "/hiskp2/werner/pipi_I1/data/A40.24/5_fit-data/", "/hiskp2/werner/pipi_I1/data/A40.32/5_fit-data/"),
    output.path = "/hiskp2/werner/pipi_I1/data/A40"
#    boot.R = 150,
    )
    
