#########################
#### Reparameterized eta Parameters
#########################

htm <- "
# tree ; cat  ; equation      ;  normal
old    ; hit  ; do            ;  m_d,sig    
old    ; hit  ; (1-do)*g      ;  m_g,sig    
old    ; miss ; (1-do)*(1-g)  ;  m_g,sig    
       
new    ; fa   ; (1-dn)*g      ;  m_g,sig    
new    ; cr   ; (1-dn)*(1-g)  ;  m_g,sig    
new    ; cr   ; dn            ;  m_d,sig
"
htm.r <- "
# tree ; cat  ; equation      ;  normal
old    ; hit  ; do            ;  m_d,5*sig    
old    ; hit  ; (1-do)*g      ;  m_d + add,5*sig    
old    ; miss ; (1-do)*(1-g)  ;  m_d + add,5*sig    
       
new    ; fa   ; (1-dn)*g      ;  m_d + add,5*sig    
new    ; cr   ; (1-dn)*(1-g)  ;  m_d + add,5*sig    
new    ; cr   ; dn            ;  m_d,5*sig
"

n <- c(targets=.5, lures=.5)  * 1000 
theta <- c(do=.6, g=.5)
eta <-   c(m_d=5, m_g=10, sig=1)
eta.r <- c(m_d=5, add=5, sig=1/5)

context("Model definition")

test_that("eta parameters are correctly reparameterized", {
  
  set.seed(123)
  
  expect_silent(hr <- read_gpt(htm.r, "normal"))
  expect_identical(hr@eta, c("add", "m_d", "sig"))
  
  expect_silent(hr2 <- read_gpt(htm.r, "normal", restrictions = list("m_d = add", "sig = 1")))
  expect_identical(hr2@eta, c("add"))
  
  # replace only variables (ok: d->D but not: lambda_d => lambda_D)
  expect_silent(hr3 <- read_gpt(htm.r, "normal", restrictions = list("d = xx", "m = 1")))
  expect_identical(hr3@eta, c("add", "m_d", "sig"))
  
  
  gen <- gpt_gen(n=n, theta=theta, eta=eta, latent=c("normal"),  
                 restrictions = list("do=dn"), file= htm)
  gen2 <- gpt_gen(n=n, theta=theta, eta=eta.r, latent=c("normal"),  
                  restrictions = list("do=dn"), file= htm.r)
  expect_true(ks.test(gen$y, gen2$y)$p> .001)
  
  #### GPT
  fit <- gpt_fit(x="x", y="y", data=gen, latent=c("normal"), 
                 n.fit = c(2,1), maxit = c(1000,1000), 
                 starting.values = c(eta, theta),
                 restrictions = list("do=dn"), file=htm)
  expect_warning(fit2 <- gpt_fit(x="x", y="y", data=gen, latent=c("normal"), 
                                 eta.lower =c("add" = 0, "sig" = 0),
                                 # eta.upper = c("m_d" = 1000), 
                                 n.fit = c(2,1), maxit = c(1000,1000), 
                                 starting.values = c(eta.r, theta),
                                 restrictions = list("do=dn"), file=htm.r))
  
  expect_equal(fit$fit.grad$loglik, fit2$fit.grad$loglik, tolerance = .1)
  expect_equal(fit$fit.grad$par[-c(4,5)], fit2$fit.grad$par[-c(3,5)], tolerance = .01)
  expect_equal(unname(fit$fit.grad$par[4]), sum(fit2$fit.grad$par[3:4]), tolerance = .01)
  expect_equal(fit$fit.grad$par[5], 5*fit2$fit.grad$par[5], tolerance = .01)
  
})