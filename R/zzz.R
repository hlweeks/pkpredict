# Define default pipercillin prior when loading this package
.onLoad <- function(libname, pkgname) {
  
  options(pkpredict.pip.default.prior = list("log_pk_mean" = c(lv_1 = 2.7760840,
                                                               lk_10 = -0.8777951,
                                                               lk_12 = -4.7726861,
                                                               lk_21 = -3.9570281),
                                             "log_pk_vcov" =  matrix(c(0.203362856, -0.009494482, -0.0025150202,  0.0011462905,
                                                                     -0.009494482,  0.153998510, -0.0027744820, -0.0064576674,
                                                                     -0.0025150202, -0.0027744820,  0.0013094598,  0.0003188164,
                                                                     0.0011462905, -0.0064576674,  0.0003188164,  0.0012682863), 
                                                                     nrow=4,ncol=4,byrow=TRUE),
                                             "log_err_mean" = -0.2741756,
                                             "log_err_sdev" = 0.02111675
                                             ))
  
  invisible()
}

# Remove default pipercillin prior when unloading this package
.onUnload <- function(libname, pkgname) {
  
  options(pip.default.prior = NULL)
  
  invisible()
}