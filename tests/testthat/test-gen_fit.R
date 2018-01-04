context("Data generation and fitting for GPT models")


test_that("gpt_gen generates valid GPT data", {
  
  ############# 1D exGaussian 
  n <- c(targets=300, lures=300)
  theta <- c(do=.7,dn=.5, g=.4)
  eta <- c(mu=400, sig=50, lambda_do=300, 
           lambda_go=500, lambda_gn=500, 
           lambda_dn=300)   
  
  file <- paste0(path.package("gpt"), "/models/2htm_exgauss.txt")
  gen <- gpt_gen(n=n, theta=theta, eta=eta, latent="exgauss", file=file)
  
  expect_s3_class(gen, "data.frame")
  expect_identical(names(gen), c("tree","x","y","state"))
  expect_s3_class(gen$tree, "factor")
  expect_s3_class(gen$x, "factor")
  expect_type(gen$y, "double")
  expect_type(gen$state, "integer")
  
  means <- c(by(gen$y, gen$state %in% c(1,6), mean))
  expect_equal(means, eta["mu"]+eta[c("lambda_go","lambda_do")], 
               tolerance=.1, check.attributes=FALSE)
    
  ############# 2D normal 
  eta <- c(m1_d=100, m1_g=120, s1=20, 
           m2_d=.4, m2_g=.6, s2=.1)   
  file <- paste0(path.package("gpt"), "/models/2htm_2normal.txt")
  gen <- gpt_gen(n=n, theta=theta, eta=eta, 
                 latent=c("normal","normal"), file=file)
  
  expect_s3_class(gen, "data.frame")
  expect_identical(names(gen), c("tree","x","y.1","y.2","state"))
  expect_s3_class(gen$tree, "factor")
  expect_s3_class(gen$x, "factor")
  expect_type(gen$y.1, "double")
  expect_type(gen$state, "integer")
  
  means <- c(by(gen$y.1, gen$state %in% c(1,6), mean))
  expect_equal(means, eta[c("m1_g","m1_d")], 
               tolerance=.1, check.attributes=FALSE)
})


test_that("gpt_fit recovers true exGaussian parameters in 2HTM", {
  
  
  n <- c(targets=130, lures=130)
  theta <- c(do=.7,dn=.5, g=.4)
  eta <- c(mu=400, sig=50, lambda_do=300, 
           lambda_go=500, lambda_gn=500, 
           lambda_dn=300)   
  
  file <- paste0(path.package("gpt"), "/models/2htm_exgauss.txt")
  gen <- gpt_gen(n=n, theta=theta, eta=eta, latent="exgauss", file=file)
  fit <- gpt_fit(x=gen$x, y=gen$y, latent="exgauss", file=file, n.fit = c(1,1),
                 restrictions=list("do=dn", "lambda_do=lambda_dn", 
                                   "lambda_go=lambda_gn"))
  
  
  expect_equal(fit$fit.EM$par[fit$gpt@theta], 
               fit$fit.grad$par[fit$gpt@theta], tolerance=.1)
  expect_equal(theta[fit$gpt@theta], 
               fit$fit.grad$par[fit$gpt@theta], tolerance=.25)
  expect_equal(eta[fit$gpt@eta], 
               fit$fit.grad$par[fit$gpt@eta], tolerance=.25)
  
  
})