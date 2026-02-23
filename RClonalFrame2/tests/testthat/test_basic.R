#Testing basic package functions
context("Testing basic package functions")

test_that("Testing basic package functions.", {
  set.seed(0)
  expect_silent(tre<-rtree(10))
  expect_is(tre,'phylo')
})


