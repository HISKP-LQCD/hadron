## Input file with default parameters
args <- list(
    ens = "A30.32",
    path.to.data = c("/home/maow/Build/sLapH-projection/A30.32-for-testing/3_gevp-data/"),
    output.path = "/home/maow/Build/sLapH-projection/A30.32-for-testing/",
    disp = c("lat"),
    maxpcs = 3,
    L = 32,
    T = 64,
    t0 = 2,
    t10 = c(7, 5, 5),
    t11 = c(11, 11, 11),
    t21 = c(19, 17, 15),
    t.step = 2,
#    disp = "lat",
#    type = "subtracted",
    reread = FALSE,
    dirs = c("p0/T1u"),
#    dirs = c("p0/T1u", "p1/A1g", "p1/Ep1g", "p2/A1g", "p2/A2g", "p2/A2u", "p3/A1g", "p3/Ep1g"),
    boot.R = 100,
    boot.l = 4,
    seed = 1234
    )
    
