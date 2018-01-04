
context("Check functions")

test_that("eta.bound.par forces parameters to be in admissible range", {
  
  lower <- c(0, 5)
  upper <- c(1, 6)
  
  f <- gpt:::eta.bound.par(c(-1,Inf), lower, upper)
  expect(all(f>lower & f<upper & abs(f) != Inf),
         "eta not forced to range")
  
  expect_warning(gpt:::eta.bound.par(c(-1,Inf), lower, upper, scale=.4))
  expect_error(gpt:::eta.bound.par(c(-1,Inf), lower, upper-5))
})


test_that("check.input.par returns parameters that match labels", {
  
  l <- letters[1:10]
  u <- runif(10)
  n <- gpt:::check.input.par(u, l)
  
  expect_silent(gpt:::check.input.par(u, l))
  expect_error(gpt:::check.input.par(u, "a"))
  
  expect_named(n, l)
  
  names(u) <- letters[sample(10)]
  m <- gpt:::check.input.par(u, l)
  expect_named(m, l)
  expect_equal(u[l],m)
})  


test_that("data.check returns error for invalid data", {
  file <- paste0(path.package("gpt"), "/models/2htm_exgauss.txt")
  mpt <- new("mpt", file=file)
  
  ##### if everything is ok:
  
  x <- sample(mpt@cat.names, 10, replace = TRUE)
  y <- exp(rnorm(10, 1,1))
  res <- gpt:::data.check(mpt, x, y)
  for(case in 1:3){
    if (case == 2) {
      res <- gpt:::data.check(mpt, "xx", "yy", data.frame(xx=x,yy=y))
      expect_null(res$group)
    }else if(case == 3){
      res <- gpt:::data.check(mpt, "xx", "yy", data.frame(xx=x,yy=y), group=rep(letters[1:2],5))
      expect_type(res$group, "character")
    }
    expect_named(res, c("x","y","group"))
    expect_type(res$x, "integer")
    expect_type(res$y, "double")
    expect_length(dim(res$y), 2)
    
  }
  
  #### errors:
  expect_error(gpt:::data.check(mpt, sample(letters, 10), y))
  expect_error(gpt:::data.check(mpt, x[1:5], y), all=FALSE)
  expect_error(gpt:::data.check(mpt, "x","yy", data.frame(xx=x,yy=y)), all=FALSE)
  expect_error(gpt:::data.check(mpt, 1:10, y))
  
})