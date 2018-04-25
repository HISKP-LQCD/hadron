context("test-hdf5_datasetname.R")

test_that("to df and back", {
  datasetnames <- c('C40B_uuuu_p110.d000.g5_p1-1-1.d000.g5_p-100.d000.g5_p-101.d000.g5',
                    'C40B_uuuu_p111.d000.g5_p-1-1-1.d000.g5_p000.d000.g5_p000.d000.g5')
  
  name_df <- datasetname_to_df(datasetnames)
  back <- df_to_datasetname(name_df)
  
  expect_equal(back, datasetnames)
})
