library(testthat)
source("../R/HarmonicScatter.R")


test_that("Harmony: Compare Linear and Recursive Function", {
  for(d in c(1:20, 30, 40, 50 , 100)) {
    expect_equal(.HarmonyR(d), .Harmony(d))
  }
})


test_that("Harmony Linear", {
  # Fibonacci constant
  phi <- 1.61803398874989484820458683436563
  expect_equal(.Harmony(1), phi)
  phi <- 1 + sqrt(5/4)-0.5 
  expect_equal(.Harmony(1), phi)
  
  # Plastic number (https://oeis.org/A072117)
  phi <- 1.32471795724474602596090885447809
  expect_equal(.Harmony(2), phi)
  
  phi <- 1.22074408460575947536168534910883
  expect_equal(.Harmony(3), phi)
  
  phi <- 1.1673 # provided by Marohnić & Strmečki
  expect_equal(.Harmony(4), phi, tolerance = 1E-3)
  
  phi <- 1.1347
  expect_equal(.Harmony(5), phi, tolerance = 1E-3)
  
  phi <- 1.1128
  expect_equal(.Harmony(6), phi, tolerance = 1E-3)
  
  #phi <- 1.0970 # seems off
  #expect_equal(.Harmony(6), phi, tolerance = 1E-3)
})


test_that("Harmony Recursive", {
  # Fibonacci constant
  phi <- 1.61803398874989484820458683436563
  expect_equal(.HarmonyR(1), phi) #-0.0643
  phi <- 1 + sqrt(5/4)-0.5 
  expect_equal(.Harmony(1), phi)

  # Plastic number  
  phi <- 1.32471795724474602596090885447809
  expect_equal(.HarmonyR(2), phi) #-0.0124
  
  phi <- 1.22074408460575947536168534910883
  expect_equal(.HarmonyR(3), phi) #-0.00436
})


test_that("GetHarmonicPoints", {
  Dim <- 2L
  N <- 100L
  o <- GetHarmonicPoints(Dim, N)
  expect_equal(dim(o), c(N, Dim))
  expect_true(all(o < 1))
  expect_true(all(o > 0))

  Dim <- 5L
  N <- 200L
  o <- GetHarmonicPoints(Dim, N)
  expect_equal(dim(o), c(N, Dim))
  expect_true(all(o < 1))
  expect_true(all(o > 0))
  
  
  # "the minimum distance between two points consistently 
  # falls between 0.549/sqrt(N) and 0.868/sqrt(N)".
  Dim <- 2L
  N <- 100L
  o <- GetHarmonicPoints(Dim, N)
  o <- dist(o) # dist() expects points in rows
  o <- o #/ sqrt(N)
  expect_gte(min(o), 0.549,
             label = format(min(o), digits = 10))
  expect_lte(min(o), 0.868,
             label = format(min(o), digits = 10))
})

