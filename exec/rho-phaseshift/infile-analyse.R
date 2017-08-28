## Input file with default parameters
args <- list(
    ens = "A40.24",
    path.to.data = c("/hiskp2/werner/pipi_I1/data/A40.24/3_gevp-data/"),
    output.path = "/hiskp2/werner/pipi_I1/data/A40.24/",
    disp = c("lat"),
    maxpcs = 3,
    L = 24,
    T = 48,
    path.to.hadron = "/hiskp2/werner/code/hadron/",
    t0 = 2,
    t10 = c(7, 5, 5),
    t11 = c(11, 11, 11),
    t21 = c(19, 17, 15),
    t.step = 2,
#    disp = "lat",
#    type = "subtracted",
    reread = FALSE,
    dirs = c("p1/A1g"),
#    dirs = c("p0/T1u", "p1/A1g", "p1/Ep1g", "p2/A1g", "p2/A2g", "p2/A2u", "p3/A1g", "p3/Ep1g"),
    boot.R = 1500,
    boot.l = 4,
    seed = 1234
    )
    
