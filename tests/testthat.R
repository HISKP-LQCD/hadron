library(testthat)
library(hadron)

# The developers of this package only use Linux or the Linux subsystem on
# Windows. Therefore we do not want to invest time in getting all the unit
# tests working on Windows. If the operating system is not linux, we just don't
# run our unit tests.
if (Sys.info()['sysname'] == 'Linux') {
  test_check("hadron")
}
