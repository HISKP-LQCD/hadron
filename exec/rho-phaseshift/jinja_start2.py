#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright © 2013, 2017 Martin Ueding <dev@martin-ueding.de>

import jinja2
import math
import os
import os.path
import subprocess
from subprocess import call

def setup_runs():

    outpath = "/hiskp2/werner/pipi_I1/data/"
    lattice = "A40.24"
    L = 24
    T = 48
#    irreps = {    0 : ["T1u"]}
    irreps = {    0 : ["T1u"], 
                  1 : ["A1g", "Ep1g"], 
                  2 : ["A1g", "A2g", "A2u"], 
                  3 : ["A1g", "Ep1g"]
                }

    # Set up jinja environment
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(".templates")
    )

    # Rendering parameters.R with values and creating folder structure
    template = env.get_template("template-parameters.j2")
    for p, irreps_p in sorted(irreps.items()):
        for irrep in irreps_p:
            folder = "%s/%s/5_fit-data/p%1d/%s"%(outpath,lattice,p,irrep)
            if not os.path.isdir(folder):
                os.makedirs(folder)
            with open(folder+"/"+"parameters.R", "w") as f:
                f.write(template.render())

    # Rendering run.sh with values.
    template = env.get_template("template-run.j2")
    with open("analyse.sh", "w") as f:
        f.write(template.render(
            irrepmoms = irreps,
            outpath = outpath,
            lattice = lattice,
          ))

    os.chmod("analyse.sh", 0o770)


def start_runs():

    subprocess.call("./analyse.sh", shell=True)

#    cd /hiskp2/werner/pipi_I1/data/A40.32/5_fit-data
#    for d in p0/T1u p1/A1g p1/Ep1g p2/A2g p2/A2u p3/A1g p3/Ep1g; do
#      echo doing $d
#      cd $d
#      Rscript /hiskp2/werner/code/hadron/exec/rho-phaseshift/analyse.R >& analyse.log
#      Rscript /hiskp2/werner/code/hadron/exec/rho-phaseshift/average.data.R >& average.log
#      cd ../../
#    done

def main():
    setup_runs()
    start_runs()

if __name__ == "__main__":
    main()


